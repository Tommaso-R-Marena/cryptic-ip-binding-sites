"""Batch infrastructure for large AlphaFold proteome analyses."""

from __future__ import annotations

import csv
import hashlib
import json
import multiprocessing as mp
import os
import sqlite3
import tempfile
import time
import logging
import threading
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

from ..utils.profiling import timed
from ..utils.resources import set_memory_limit
from ..utils.input_validation import (
    BatchDownloadConfig,
    ExportConfig,
    ParallelProcessorConfig,
    parse_or_raise,
)
from ..utils.logging_utils import configure_logging, log_with_context
from ..errors import (
    BatchItemProcessingError,
    CacheOperationError,
    NetworkRetryError,
    OperationTimeoutError,
    RecoveryStateError,
)
from .alphafold_client import AlphaFoldClient


class AlphaFoldBatchDownloader:
    """Download AlphaFold structures for entire UniProt proteomes with resume support."""

    UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

    def __init__(
        self,
        output_dir: Path | str,
        state_path: Optional[Path | str] = None,
        requests_per_second: float = 4.0,
        max_retries: int = 5,
        backoff_seconds: float = 1.5,
        session: Optional[requests.Session] = None,
        timeout_seconds: float = 30.0,
        log_dir: Path | str = "logs",
    ) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.state_path = Path(state_path) if state_path else self.output_dir / "download_state.json"
        config = parse_or_raise(
            BatchDownloadConfig,
            {
                "requests_per_second": requests_per_second,
                "max_retries": max_retries,
                "backoff_seconds": backoff_seconds,
                "timeout_seconds": timeout_seconds,
            },
            "Invalid batch downloader configuration",
        )
        self.requests_per_second = config.requests_per_second
        self.max_retries = config.max_retries
        self.backoff_seconds = config.backoff_seconds
        self.timeout_seconds = config.timeout_seconds
        self.session = session or requests.Session()
        self.af_client = AlphaFoldClient(cache_dir=self.output_dir)
        self.logger = configure_logging(log_dir, "batch_downloader")

    def _load_state(self) -> Dict[str, Any]:
        if self.state_path.exists():
            try:
                return json.loads(self.state_path.read_text())
            except json.JSONDecodeError as exc:
                raise RecoveryStateError(
                    f"Checkpoint file '{self.state_path}' is corrupted and cannot be resumed."
                ) from exc
        return {"completed": {}, "failed": {}}

    def _save_state(self, state: Dict[str, Any]) -> None:
        self.state_path.write_text(json.dumps(state, indent=2, sort_keys=True))

    def _rate_limit(self, started_at: float) -> None:
        minimum_interval = 1.0 / self.requests_per_second
        elapsed = time.time() - started_at
        if elapsed < minimum_interval:
            time.sleep(minimum_interval - elapsed)

    def _request_with_retry(self, method: str, url: str, **kwargs: Any) -> requests.Response:
        for attempt in range(1, self.max_retries + 1):
            try:
                response = self.session.request(method, url, timeout=self.timeout_seconds, **kwargs)
                response.raise_for_status()
                return response
            except requests.Timeout as exc:
                if attempt == self.max_retries:
                    raise OperationTimeoutError(
                        f"Request to {url} timed out after {self.max_retries} attempts."
                    ) from exc
                time.sleep(self.backoff_seconds * (2 ** (attempt - 1)))
            except requests.RequestException as exc:
                if attempt == self.max_retries:
                    raise NetworkRetryError(
                        f"Request to {url} failed after {self.max_retries} attempts."
                    ) from exc
                time.sleep(self.backoff_seconds * (2 ** (attempt - 1)))
        raise RuntimeError("Unreachable retry branch")

    @timed()
    def fetch_proteome_uniprot_ids(self, proteome_id: str) -> List[str]:
        """List all UniProt accessions for a proteome using UniProt pagination."""
        cursor: Optional[str] = None
        uniprot_ids: List[str] = []

        while True:
            params = {
                "query": f"proteome:{proteome_id}",
                "fields": "accession",
                "format": "json",
                "size": "500",
            }
            if cursor:
                params["cursor"] = cursor

            started = time.time()
            response = self._request_with_retry("GET", self.UNIPROT_SEARCH_URL, params=params)
            payload = response.json()
            uniprot_ids.extend(entry["primaryAccession"] for entry in payload.get("results", []))
            cursor = payload.get("nextCursor")
            self._rate_limit(started)

            if not cursor:
                break

        return uniprot_ids

    @timed()
    def download_proteomes(self, proteome_ids: Sequence[str], resume: bool = True) -> Dict[str, int]:
        """Download AlphaFold structures for every protein in each proteome."""
        state = self._load_state() if resume else {"completed": {}, "failed": {}}
        summary = {"downloaded": 0, "skipped": 0, "failed": 0}

        for proteome_id in proteome_ids:
            all_ids = self.fetch_proteome_uniprot_ids(proteome_id)
            completed = set(state["completed"].get(proteome_id, []))
            failed = set(state["failed"].get(proteome_id, []))

            progress = tqdm(all_ids, desc=f"{proteome_id}", unit="protein")
            for uniprot_id in progress:
                if uniprot_id in completed:
                    summary["skipped"] += 1
                    continue

                started = time.time()
                try:
                    self.af_client.fetch_structure(uniprot_id)
                    completed.add(uniprot_id)
                    failed.discard(uniprot_id)
                    summary["downloaded"] += 1
                except Exception as exc:
                    log_with_context(
                        self.logger,
                        logging.WARNING,
                        "Skipping failed structure download",
                        proteome_id=proteome_id,
                        uniprot_id=uniprot_id,
                        error=str(exc),
                    )
                    failed.add(uniprot_id)
                    summary["failed"] += 1
                finally:
                    state["completed"][proteome_id] = sorted(completed)
                    state["failed"][proteome_id] = sorted(failed)
                    self._save_state(state)
                    self._rate_limit(started)
                    progress.set_postfix(done=len(completed), failed=len(failed))

        return summary


class AnalysisCache:
    """SQLite-backed cache for expensive structure analysis results."""

    def __init__(
        self,
        db_path: Path | str,
        pipeline_version: str,
        pipeline_params: Optional[Dict[str, Any]] = None,
        query_cache_size: int = 2048,
        max_db_retries: int = 5,
        log_dir: Path | str = "logs",
    ):
        self.db_path = Path(db_path)
        self.pipeline_version = pipeline_version
        self.pipeline_params = pipeline_params or {}
        self.params_hash = self._hash_pipeline_params(self.pipeline_params)
        self.connection = sqlite3.connect(self.db_path, check_same_thread=False)
        self.connection.row_factory = sqlite3.Row
        self.query_cache_size = query_cache_size
        self.max_db_retries = max_db_retries
        self.logger = configure_logging(log_dir, "analysis_cache")
        self._db_lock = threading.Lock()
        self._query_cache: OrderedDict[Tuple[str, str, str], Dict[str, Any]] = OrderedDict()
        self._initialize_schema()

    def _with_db_retry(self, operation: Callable[[], Any]) -> Any:
        for attempt in range(1, self.max_db_retries + 1):
            try:
                with self._db_lock:
                    return operation()
            except sqlite3.OperationalError as exc:
                if "locked" not in str(exc).lower() or attempt == self.max_db_retries:
                    raise CacheOperationError(
                        f"Database operation failed after {attempt} attempts: {exc}"
                    ) from exc
                time.sleep(0.1 * (2 ** (attempt - 1)))

    @staticmethod
    def _hash_pipeline_params(pipeline_params: Dict[str, Any]) -> str:
        serialized = json.dumps(pipeline_params, sort_keys=True, default=str)
        return hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    def _initialize_schema(self) -> None:
        self.connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS analysis_cache (
                uniprot_id TEXT NOT NULL,
                organism TEXT NOT NULL,
                analysis_results TEXT NOT NULL,
                timestamp TEXT NOT NULL,
                pipeline_version TEXT NOT NULL,
                params_hash TEXT NOT NULL,
                PRIMARY KEY (uniprot_id, pipeline_version, params_hash)
            );

            CREATE INDEX IF NOT EXISTS idx_cache_lookup
            ON analysis_cache (uniprot_id, organism, pipeline_version, params_hash);
            """
        )
        self.connection.commit()

    def _cache_key(self, uniprot_id: str, organism: str) -> Tuple[str, str, str]:
        return (uniprot_id, organism, self.params_hash)

    def _update_query_cache(self, key: Tuple[str, str, str], value: Dict[str, Any]) -> None:
        self._query_cache[key] = value
        self._query_cache.move_to_end(key)
        while len(self._query_cache) > self.query_cache_size:
            self._query_cache.popitem(last=False)

    @timed()
    def get_cached_result(self, uniprot_id: str, organism: str) -> Optional[Dict[str, Any]]:
        key = self._cache_key(uniprot_id, organism)
        if key in self._query_cache:
            return self._query_cache[key]

        row = self._with_db_retry(
            lambda: self.connection.execute(
                """
                SELECT analysis_results
                FROM analysis_cache
                WHERE uniprot_id = ?
                  AND organism = ?
                  AND pipeline_version = ?
                  AND params_hash = ?
                """,
                (uniprot_id, organism, self.pipeline_version, self.params_hash),
            ).fetchone()
        )

        if row is None:
            return None

        result = json.loads(row["analysis_results"])
        self._update_query_cache(key, result)
        return result

    @timed()
    def set_cached_result(self, uniprot_id: str, organism: str, analysis_results: Dict[str, Any]) -> None:
        timestamp = datetime.now(timezone.utc).isoformat()
        self._with_db_retry(
            lambda: self.connection.execute(
                """
                INSERT INTO analysis_cache (
                    uniprot_id, organism, analysis_results, timestamp, pipeline_version, params_hash
                ) VALUES (?, ?, ?, ?, ?, ?)
                ON CONFLICT(uniprot_id, pipeline_version, params_hash)
                DO UPDATE SET
                    organism = excluded.organism,
                    analysis_results = excluded.analysis_results,
                    timestamp = excluded.timestamp
                """,
                (
                    uniprot_id,
                    organism,
                    json.dumps(analysis_results),
                    timestamp,
                    self.pipeline_version,
                    self.params_hash,
                ),
            )
        )
        self._with_db_retry(self.connection.commit)
        self._update_query_cache(self._cache_key(uniprot_id, organism), analysis_results)


    @timed()
    def set_cached_results_batch(self, rows: Sequence[Tuple[str, str, Dict[str, Any]]]) -> int:
        """Batch upsert cache rows in one transaction for throughput."""
        timestamp = datetime.now(timezone.utc).isoformat()
        payload = [
            (
                uniprot_id,
                organism,
                json.dumps(analysis_results),
                timestamp,
                self.pipeline_version,
                self.params_hash,
            )
            for uniprot_id, organism, analysis_results in rows
        ]
        def _batch_upsert() -> None:
            with self.connection:
                self.connection.executemany(
                    """
                    INSERT INTO analysis_cache (
                        uniprot_id, organism, analysis_results, timestamp, pipeline_version, params_hash
                    ) VALUES (?, ?, ?, ?, ?, ?)
                    ON CONFLICT(uniprot_id, pipeline_version, params_hash)
                    DO UPDATE SET
                        organism = excluded.organism,
                        analysis_results = excluded.analysis_results,
                        timestamp = excluded.timestamp
                    """,
                    payload,
                )

        self._with_db_retry(_batch_upsert)
        for uniprot_id, organism, analysis_results in rows:
            self._update_query_cache(self._cache_key(uniprot_id, organism), analysis_results)
        return len(rows)

    def invalidate_outdated_cache(self) -> int:
        """Remove cached entries that do not match current pipeline settings."""
        cursor = self.connection.execute(
            """
            DELETE FROM analysis_cache
            WHERE pipeline_version != ? OR params_hash != ?
            """,
            (self.pipeline_version, self.params_hash),
        )
        self.connection.commit()
        return cursor.rowcount

    @timed()
    def export_results(self, output_path: Path | str, export_format: str) -> Path:
        export = parse_or_raise(
            ExportConfig,
            {"export_format": export_format},
            "Invalid export configuration",
        )
        export_path = Path(output_path)
        data = pd.read_sql_query("SELECT * FROM analysis_cache", self.connection)

        if export.export_format == "csv":
            data.to_csv(export_path, index=False)
        elif export.export_format == "json":
            data.to_json(export_path, orient="records", indent=2)
        elif export.export_format == "hdf5":
            data.to_hdf(export_path, key="analysis_cache", mode="w")

        return export_path

    @timed()
    def vacuum(self) -> None:
        self.connection.execute("VACUUM")

    def close(self) -> None:
        self.connection.close()


_WORKER_FUNCTION: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None


def _init_worker(analyze_function: Callable[[Dict[str, Any]], Dict[str, Any]]) -> None:
    global _WORKER_FUNCTION
    _WORKER_FUNCTION = analyze_function


def _run_worker(item: Dict[str, Any]) -> Dict[str, Any]:
    if _WORKER_FUNCTION is None:
        raise RuntimeError("Worker function is not initialized")
    return _WORKER_FUNCTION(item)


class ParallelProcessor:
    """Process proteins in parallel with chunking + checkpoint resume."""

    def __init__(
        self,
        analyze_function: Callable[[Dict[str, Any]], Dict[str, Any]],
        workers: Optional[int] = None,
        chunk_size: int = 100,
        checkpoint_path: Path | str = "processing_checkpoint.json",
        memory_limit_gb: Optional[float] = None,
        use_memmap: bool = False,
        timeout_per_item_seconds: Optional[float] = None,
        checkpoint_every: int = 1,
        log_dir: Path | str = "logs",
    ) -> None:
        config = parse_or_raise(
            ParallelProcessorConfig,
            {
                "workers": workers or max(1, mp.cpu_count() - 1),
                "chunk_size": chunk_size,
                "checkpoint_every": checkpoint_every,
            },
            "Invalid parallel processor configuration",
        )
        self.analyze_function = analyze_function
        self.workers = config.workers
        self.chunk_size = config.chunk_size
        self.checkpoint_every = config.checkpoint_every
        self.checkpoint_path = Path(checkpoint_path)
        self.memory_limit_gb = memory_limit_gb
        self.use_memmap = use_memmap
        self.timeout_per_item_seconds = timeout_per_item_seconds
        self.logger = configure_logging(log_dir, "parallel_processor")
        self._memmap_paths: List[Path] = []

    def _load_checkpoint(self) -> Dict[str, Any]:
        if self.checkpoint_path.exists():
            try:
                return json.loads(self.checkpoint_path.read_text())
            except json.JSONDecodeError as exc:
                raise RecoveryStateError(
                    f"Checkpoint file '{self.checkpoint_path}' is invalid JSON."
                ) from exc
        return {"processed": [], "started_at": time.time()}

    def _save_checkpoint(self, checkpoint: Dict[str, Any]) -> None:
        self.checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True))

    def _chunked(self, items: Sequence[Dict[str, Any]]) -> Iterator[Sequence[Dict[str, Any]]]:
        for index in range(0, len(items), self.chunk_size):
            yield items[index : index + self.chunk_size]

    @timed()
    def run(self, items: Sequence[Dict[str, Any]], resume: bool = True) -> List[Dict[str, Any]]:
        if self.memory_limit_gb:
            set_memory_limit(self.memory_limit_gb)

        checkpoint = self._load_checkpoint() if resume else {"processed": [], "started_at": time.time()}
        processed = set(checkpoint["processed"])
        remaining = [item for item in items if item["uniprot_id"] not in processed]

        if not remaining:
            return []

        results: List[Dict[str, Any]] = []
        total = len(remaining)
        progress = tqdm(total=total, desc="Analyzing", unit="protein")
        start = time.time()

        with mp.Pool(processes=self.workers, initializer=_init_worker, initargs=(self.analyze_function,)) as pool:
            for chunk in self._chunked(remaining):
                iterator = pool.imap_unordered(_run_worker, chunk)
                while True:
                    try:
                        result = iterator.next(timeout=self.timeout_per_item_seconds)
                    except mp.TimeoutError as exc:
                        raise OperationTimeoutError(
                            f"Parallel item processing exceeded timeout ({self.timeout_per_item_seconds}s)."
                        ) from exc
                    except StopIteration:
                        break
                    except Exception as exc:
                        log_with_context(
                            self.logger,
                            logging.WARNING,
                            "Skipping failed item during parallel processing",
                            error=str(exc),
                        )
                        continue

                    if "uniprot_id" not in result:
                        raise BatchItemProcessingError("unknown", "result missing uniprot_id")

                    results.append(result)
                    processed.add(result["uniprot_id"])
                    if len(processed) % self.checkpoint_every == 0:
                        checkpoint["processed"] = sorted(processed)
                        self._save_checkpoint(checkpoint)

                    progress.update(1)
                    elapsed = time.time() - start
                    rate = progress.n / elapsed if elapsed else 0
                    eta = (total - progress.n) / rate if rate else 0
                    progress.set_postfix(eta_seconds=f"{eta:.1f}")

        checkpoint["processed"] = sorted(processed)
        self._save_checkpoint(checkpoint)
        progress.close()
        return results

    @timed()
    def run_streaming(self, items: Sequence[Dict[str, Any]], resume: bool = True) -> Iterator[Dict[str, Any]]:
        """Yield results incrementally to minimize peak memory usage."""
        if self.memory_limit_gb:
            set_memory_limit(self.memory_limit_gb)

        checkpoint = self._load_checkpoint() if resume else {"processed": [], "started_at": time.time()}
        processed = set(checkpoint["processed"])
        remaining = [item for item in items if item["uniprot_id"] not in processed]
        if not remaining:
            return

        with mp.Pool(processes=self.workers, initializer=_init_worker, initargs=(self.analyze_function,)) as pool:
            for chunk in self._chunked(remaining):
                for result in pool.imap_unordered(_run_worker, chunk):
                    processed.add(result["uniprot_id"])
                    checkpoint["processed"] = sorted(processed)
                    self._save_checkpoint(checkpoint)
                    yield result

    @timed()
    def create_result_memmap(self, size: int, columns: int = 8) -> np.memmap:
        """Create a temporary memory-mapped array for large numeric outputs."""
        memmap_path = Path(tempfile.gettempdir()) / f"cryptic_ip_{os.getpid()}_{int(time.time())}.mmap"
        self._memmap_paths.append(memmap_path)
        return np.memmap(memmap_path, dtype="float32", mode="w+", shape=(size, columns))

    def cleanup_temp_files(self) -> None:
        for path in self._memmap_paths:
            if path.exists():
                path.unlink()
        self._memmap_paths.clear()


class ResourceMonitor:
    """Background monitor that captures process resource usage for long jobs."""

    def __init__(self, logger_name: str = "resource_monitor", interval_seconds: float = 5.0):
        self.interval_seconds = interval_seconds
        self._stop = threading.Event()
        self.logger = configure_logging("logs", logger_name)
        self._thread: Optional[threading.Thread] = None

    def start(self) -> None:
        def _loop() -> None:
            while not self._stop.is_set():
                rss_bytes = 0
                try:
                    import resource

                    rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                except Exception:
                    rss_bytes = 0
                log_with_context(self.logger, logging.INFO, "resource_snapshot", rss_kb=rss_bytes)
                self._stop.wait(self.interval_seconds)

        self._thread = threading.Thread(target=_loop, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2)


@timed()
def append_results_to_file(results: Iterable[Dict[str, Any]], output_path: Path | str) -> Path:
    """Append incremental run outputs to CSV in a checkpoint-safe way."""
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    rows = list(results)
    if not rows:
        return output_file

    fieldnames = sorted({key for row in rows for key in row.keys()})
    write_header = not output_file.exists()

    with output_file.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)

    return output_file
