"""Batch infrastructure for large AlphaFold proteome analyses."""

from __future__ import annotations

import csv
import hashlib
import json
import multiprocessing as mp
import sqlite3
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import pandas as pd
import requests
from tqdm import tqdm

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
    ) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.state_path = Path(state_path) if state_path else self.output_dir / "download_state.json"
        self.requests_per_second = max(requests_per_second, 0.1)
        self.max_retries = max_retries
        self.backoff_seconds = backoff_seconds
        self.session = session or requests.Session()
        self.af_client = AlphaFoldClient(cache_dir=self.output_dir)

    def _load_state(self) -> Dict[str, Any]:
        if self.state_path.exists():
            return json.loads(self.state_path.read_text())
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
                response = self.session.request(method, url, timeout=30, **kwargs)
                response.raise_for_status()
                return response
            except requests.RequestException:
                if attempt == self.max_retries:
                    raise
                time.sleep(self.backoff_seconds * (2 ** (attempt - 1)))
        raise RuntimeError("Unreachable retry branch")

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
                except Exception:
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

    def __init__(self, db_path: Path | str, pipeline_version: str, pipeline_params: Optional[Dict[str, Any]] = None):
        self.db_path = Path(db_path)
        self.pipeline_version = pipeline_version
        self.pipeline_params = pipeline_params or {}
        self.params_hash = self._hash_pipeline_params(self.pipeline_params)
        self.connection = sqlite3.connect(self.db_path)
        self.connection.row_factory = sqlite3.Row
        self._initialize_schema()

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

    def get_cached_result(self, uniprot_id: str, organism: str) -> Optional[Dict[str, Any]]:
        row = self.connection.execute(
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

        if row is None:
            return None
        return json.loads(row["analysis_results"])

    def set_cached_result(self, uniprot_id: str, organism: str, analysis_results: Dict[str, Any]) -> None:
        timestamp = datetime.now(timezone.utc).isoformat()
        self.connection.execute(
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
        self.connection.commit()

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

    def export_results(self, output_path: Path | str, export_format: str) -> Path:
        export_path = Path(output_path)
        data = pd.read_sql_query("SELECT * FROM analysis_cache", self.connection)

        if export_format == "csv":
            data.to_csv(export_path, index=False)
        elif export_format == "json":
            data.to_json(export_path, orient="records", indent=2)
        elif export_format == "hdf5":
            data.to_hdf(export_path, key="analysis_cache", mode="w")
        else:
            raise ValueError(f"Unsupported export format: {export_format}")

        return export_path

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
    ) -> None:
        self.analyze_function = analyze_function
        self.workers = workers or max(1, mp.cpu_count() - 1)
        self.chunk_size = chunk_size
        self.checkpoint_path = Path(checkpoint_path)

    def _load_checkpoint(self) -> Dict[str, Any]:
        if self.checkpoint_path.exists():
            return json.loads(self.checkpoint_path.read_text())
        return {"processed": [], "started_at": time.time()}

    def _save_checkpoint(self, checkpoint: Dict[str, Any]) -> None:
        self.checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True))

    def _chunked(self, items: Sequence[Dict[str, Any]]) -> Iterator[Sequence[Dict[str, Any]]]:
        for index in range(0, len(items), self.chunk_size):
            yield items[index : index + self.chunk_size]

    def run(self, items: Sequence[Dict[str, Any]], resume: bool = True) -> List[Dict[str, Any]]:
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
                for result in pool.imap_unordered(_run_worker, chunk):
                    results.append(result)
                    processed.add(result["uniprot_id"])
                    checkpoint["processed"] = sorted(processed)
                    self._save_checkpoint(checkpoint)

                    progress.update(1)
                    elapsed = time.time() - start
                    rate = progress.n / elapsed if elapsed else 0
                    eta = (total - progress.n) / rate if rate else 0
                    progress.set_postfix(eta_seconds=f"{eta:.1f}")

        progress.close()
        return results


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
