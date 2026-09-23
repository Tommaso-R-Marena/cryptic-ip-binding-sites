"""Concurrent, verified downloads for every remote data source.

Every structure and annotation the pipeline uses is fetched through here, so
the properties that decide whether a dataset can be trusted are implemented
once:

* **Nothing unverified is kept.** A response is checked (gzip integrity,
  expected format markers, minimum size) before it is written, and it is
  written to a temporary file and renamed into place, so a download cut off
  half-way can never be mistaken for a cached file on the next run.
* **Transient failures are retried; permanent ones are not.** Connection
  errors, timeouts, 429 and 5xx are retried with exponential backoff and
  jitter, honouring ``Retry-After``. A 404 is recorded as absent at once.
* **Resumable.** A destination that already exists and passes validation is
  not downloaded again.
* **Concurrent but polite.** Requests run on one event loop with a global
  concurrency bound and a per-host bound, and an optional per-host rate limit.
* **Accountable.** Every job yields a :class:`FetchResult` - status, bytes,
  SHA-256, attempts, error - which callers write to a manifest.

``asyncio`` suits this workload: downloads are latency-bound, so dozens can be
in flight on one thread without the overhead or the thread-safety questions of
a thread pool sharing a ``requests.Session``.
"""

from __future__ import annotations

import asyncio
import gzip
import hashlib
import json
import logging
import os
import random
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence
from urllib.parse import urlsplit

LOGGER = logging.getLogger(__name__)

DEFAULT_USER_AGENT = (
    "cryptic-ip-binding-sites/1.0 "
    "(+https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites)"
)
RETRY_STATUSES = frozenset({408, 425, 429, 500, 502, 503, 504})


class ValidationError(ValueError):
    """Downloaded content is not what the job expected."""


Validator = Callable[[bytes], None]


# ---------------------------------------------------------------- validators
def validate_pdb(payload: bytes) -> None:
    """A PDB-format coordinate file: has coordinates and an END record.

    The whole file is searched: a large entry's header (REMARK, SEQRES, ...)
    can run to several hundred kilobytes before the first coordinate record.
    """
    if b"\nATOM  " not in payload and b"\nHETATM" not in payload and not payload.startswith((b"ATOM  ", b"HETATM")):
        raise ValidationError("no ATOM/HETATM records")
    if b"\nEND" not in payload[-4096:] and not payload.rstrip().endswith(b"END"):
        raise ValidationError("no END record: truncated PDB file")


def validate_mmcif(payload: bytes) -> None:
    """An mmCIF coordinate file: a data block with an atom_site loop."""
    if not payload.lstrip().startswith(b"data_"):
        raise ValidationError("does not start with a data_ block")
    if b"_atom_site." not in payload:
        raise ValidationError("no _atom_site category")


def validate_json(payload: bytes) -> None:
    try:
        json.loads(payload)
    except ValueError as exc:
        raise ValidationError(f"invalid JSON: {exc}") from exc


def validate_nonempty(payload: bytes) -> None:
    if not payload.strip():
        raise ValidationError("empty response")


def for_path(path: Path) -> Optional[Validator]:
    """The validator implied by a destination's extension."""
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".pdb") or name.endswith(".ent"):
        return validate_pdb
    if name.endswith(".cif"):
        return validate_mmcif
    if name.endswith(".json"):
        return validate_json
    return validate_nonempty


# ---------------------------------------------------------------- jobs
@dataclass
class FetchJob:
    """One download.

    Attributes:
        key: Caller's identifier (e.g. accession), echoed in the result.
        url: Source URL.
        dest: Destination path. ``None`` returns the payload in memory.
        decompress: Gunzip the payload before validating and writing
            (``"auto"`` decompresses when the payload is gzip).
        validator: Content check; defaults to one chosen by ``dest``'s suffix.
    """

    key: str
    url: str
    dest: Optional[Path] = None
    decompress: object = "auto"
    validator: Optional[Validator] = None


@dataclass
class FetchResult:
    key: str
    url: str
    ok: bool
    status: Optional[int] = None
    path: Optional[str] = None
    n_bytes: int = 0
    sha256: str = ""
    attempts: int = 0
    seconds: float = 0.0
    from_cache: bool = False
    not_found: bool = False
    error: str = ""
    payload: Optional[bytes] = field(default=None, repr=False)

    def to_record(self) -> Dict[str, object]:
        record = asdict(self)
        record.pop("payload", None)
        return record


def _decode(payload: bytes, decompress: object) -> bytes:
    gz = payload[:2] == b"\x1f\x8b"
    if decompress is True or (decompress == "auto" and gz):
        try:
            return gzip.decompress(payload)
        except (OSError, EOFError) as exc:
            # A truncated or corrupt gzip stream is a failed download, never a
            # file to be stored as it is.
            raise ValidationError(f"corrupt gzip stream: {exc}") from exc
    return payload


def _check(payload: bytes, job: FetchJob) -> None:
    validator = job.validator or (for_path(job.dest) if job.dest is not None else validate_nonempty)
    if validator is not None:
        validator(payload)


def _cached(job: FetchJob) -> Optional[FetchResult]:
    """A destination that exists and validates is a completed download."""
    if job.dest is None or not job.dest.exists() or job.dest.stat().st_size == 0:
        return None
    data = job.dest.read_bytes()
    try:
        _check(data, job)
    except ValidationError:
        LOGGER.warning("Cached %s fails validation; downloading again", job.dest)
        return None
    return FetchResult(
        key=job.key, url=job.url, ok=True, path=str(job.dest), n_bytes=len(data),
        sha256=hashlib.sha256(data).hexdigest(), from_cache=True,
    )


def _write_atomic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.part")
    tmp.write_bytes(data)
    os.replace(tmp, path)


# ---------------------------------------------------------------- engine
class _HostLimiter:
    """Per-host concurrency bound and minimum interval between requests."""

    def __init__(self, per_host: int, min_interval: float) -> None:
        self.per_host = per_host
        self.min_interval = min_interval
        self._sems: Dict[str, asyncio.Semaphore] = {}
        self._locks: Dict[str, asyncio.Lock] = {}
        self._last: Dict[str, float] = {}

    def semaphore(self, host: str) -> asyncio.Semaphore:
        if host not in self._sems:
            self._sems[host] = asyncio.Semaphore(self.per_host)
            self._locks[host] = asyncio.Lock()
        return self._sems[host]

    async def pace(self, host: str) -> None:
        if self.min_interval <= 0:
            return
        async with self._locks[host]:
            wait = self._last.get(host, 0.0) + self.min_interval - time.monotonic()
            if wait > 0:
                await asyncio.sleep(wait)
            self._last[host] = time.monotonic()


async def _fetch_one(session, job: FetchJob, limiter: _HostLimiter, sem: asyncio.Semaphore,
                     *, max_retries: int, backoff: float, max_backoff: float) -> FetchResult:
    import aiohttp

    cached = _cached(job)
    if cached is not None:
        return cached

    started = time.monotonic()
    host = urlsplit(job.url).netloc
    result = FetchResult(key=job.key, url=job.url, ok=False)
    for attempt in range(1, max_retries + 1):
        result.attempts = attempt
        retry_after: Optional[float] = None
        try:
            async with sem, limiter.semaphore(host):
                await limiter.pace(host)
                async with session.get(job.url) as response:
                    result.status = response.status
                    if response.status == 404 or response.status == 410:
                        result.not_found = True
                        result.error = f"HTTP {response.status}"
                        break
                    if response.status in RETRY_STATUSES:
                        header = response.headers.get("Retry-After")
                        retry_after = float(header) if header and header.isdigit() else None
                        raise aiohttp.ClientResponseError(
                            response.request_info, response.history,
                            status=response.status, message="retryable status",
                        )
                    if response.status >= 400:
                        result.error = f"HTTP {response.status}"
                        break
                    raw = await response.read()
            data = _decode(raw, job.decompress)
            _check(data, job)
            if job.dest is not None:
                _write_atomic(job.dest, data)
                result.path = str(job.dest)
            else:
                result.payload = data
            result.ok = True
            result.error = ""
            result.n_bytes = len(data)
            result.sha256 = hashlib.sha256(data).hexdigest()
            break
        except ValidationError as exc:
            # Wrong content from a 200 response: retried, since a proxy or
            # server hiccup can return an error page or a cut stream.
            result.error = f"validation: {exc}"
        except (aiohttp.ClientError, asyncio.TimeoutError, OSError) as exc:
            result.error = f"{type(exc).__name__}: {exc}"[:300]
        if attempt < max_retries:
            delay = retry_after if retry_after is not None else min(max_backoff, backoff * 2 ** (attempt - 1))
            await asyncio.sleep(delay * (0.75 + 0.5 * random.random()))
    result.seconds = round(time.monotonic() - started, 3)
    return result


async def fetch_all_async(
    jobs: Sequence[FetchJob],
    *,
    concurrency: int = 16,
    per_host: int = 8,
    min_interval: float = 0.0,
    max_retries: int = 6,
    backoff: float = 1.0,
    max_backoff: float = 60.0,
    timeout: float = 180.0,
    user_agent: str = DEFAULT_USER_AGENT,
    headers: Optional[Mapping[str, str]] = None,
    progress: Optional[Callable[[FetchResult], None]] = None,
) -> List[FetchResult]:
    """Run ``jobs`` concurrently; results are returned in job order."""
    import aiohttp

    limiter = _HostLimiter(per_host=max(1, per_host), min_interval=min_interval)
    sem = asyncio.Semaphore(max(1, concurrency))
    client_timeout = aiohttp.ClientTimeout(total=timeout, sock_connect=30)
    connector = aiohttp.TCPConnector(limit=max(1, concurrency), limit_per_host=max(1, per_host))
    hdrs = {"User-Agent": user_agent, **dict(headers or {})}
    async with aiohttp.ClientSession(
        timeout=client_timeout, connector=connector, headers=hdrs, trust_env=True
    ) as session:

        async def run(job: FetchJob) -> FetchResult:
            result = await _fetch_one(
                session, job, limiter, sem,
                max_retries=max_retries, backoff=backoff, max_backoff=max_backoff,
            )
            if progress is not None:
                progress(result)
            return result

        return list(await asyncio.gather(*(run(job) for job in jobs)))


def fetch_all(jobs: Sequence[FetchJob], **kwargs) -> List[FetchResult]:
    """Synchronous entry point for :func:`fetch_all_async`."""
    if not jobs:
        return []
    return asyncio.run(fetch_all_async(list(jobs), **kwargs))


def write_manifest(results: Iterable[FetchResult], path: Path) -> Dict[str, int]:
    """Record every download's outcome; return counts by outcome."""
    records = [result.to_record() for result in results]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(records, indent=2), encoding="utf-8")
    summary = {
        "total": len(records),
        "downloaded": sum(1 for r in records if r["ok"] and not r["from_cache"]),
        "cached": sum(1 for r in records if r["from_cache"]),
        "not_found": sum(1 for r in records if r["not_found"]),
        "failed": sum(1 for r in records if not r["ok"] and not r["not_found"]),
    }
    return summary


# ---------------------------------------------------------------- AlphaFold
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


def _pick_alphafold_entry(entries: object, accession: str) -> Optional[Mapping[str, object]]:
    """The canonical full-length (F1) model for ``accession``.

    The prediction endpoint can list several entries - isoforms, fragments of
    very long proteins - so the entry is chosen by accession and fragment
    rather than taken first.
    """
    if not isinstance(entries, list):
        entries = [entries] if isinstance(entries, dict) else []
    candidates = [e for e in entries if isinstance(e, dict)]
    exact = [e for e in candidates if str(e.get("uniprotAccession", "")).upper() == accession.upper()]
    pool = exact or candidates
    for entry in pool:
        if str(entry.get("entryId", "")).endswith("-F1"):
            return entry
    return pool[0] if pool else None


def fetch_alphafold_models(
    accessions: Sequence[str],
    out_dir: Path,
    *,
    fmt: str = "pdb",
    concurrency: int = 16,
    manifest: Optional[Path] = None,
    **kwargs,
) -> Dict[str, FetchResult]:
    """Download the current AlphaFold model for each accession.

    The model URL comes from the prediction API (``pdbUrl``/``cifUrl``), so the
    file fetched is the current release whatever its version suffix, and a
    stale cached version is not mistaken for it: the destination carries the
    version the API reports.

    Returns:
        Mapping accession -> result of the model download (or of the API
        lookup, when that failed).
    """
    out_dir = Path(out_dir)
    unique = list(dict.fromkeys(a.strip().upper() for a in accessions if a and a.strip()))
    lookups = fetch_all(
        [FetchJob(key=a, url=ALPHAFOLD_API.format(accession=a), validator=validate_json) for a in unique],
        concurrency=concurrency, **kwargs,
    )
    url_key = "pdbUrl" if fmt == "pdb" else "cifUrl"
    jobs: List[FetchJob] = []
    results: Dict[str, FetchResult] = {}
    for lookup in lookups:
        if not lookup.ok:
            results[lookup.key] = lookup
            continue
        entry = _pick_alphafold_entry(json.loads(lookup.payload or b"null"), lookup.key)
        url = entry.get(url_key) if entry else None
        if not url:
            lookup.ok = False
            lookup.not_found = True
            lookup.error = f"no {url_key} in AlphaFold API response"
            lookup.payload = None
            results[lookup.key] = lookup
            continue
        name = str(url).rsplit("/", 1)[-1]
        jobs.append(FetchJob(key=lookup.key, url=str(url), dest=out_dir / name))
    for result in fetch_all(jobs, concurrency=concurrency, **kwargs):
        results[result.key] = result
    if manifest is not None:
        write_manifest(results.values(), manifest)
    return results


# ---------------------------------------------------------------- RCSB
RCSB_FILE = "https://files.rcsb.org/download/{pdb_id}.{fmt}.gz"


def fetch_rcsb_structures(
    pdb_ids: Sequence[str],
    out_dir: Path,
    *,
    prefer: Sequence[str] = ("cif", "pdb"),
    concurrency: int = 8,
    manifest: Optional[Path] = None,
    **kwargs,
) -> Dict[str, FetchResult]:
    """Download RCSB entries, falling back through ``prefer`` formats.

    mmCIF first: entries with more than 62 chains or 99,999 atoms have no
    legacy PDB file. An entry is retried in the next format only when the
    previous one was absent (404), not when it failed transiently.
    """
    out_dir = Path(out_dir)
    pending = list(dict.fromkeys(p.strip().upper() for p in pdb_ids if p and p.strip()))
    results: Dict[str, FetchResult] = {}
    for fmt in prefer:
        if not pending:
            break
        jobs = [
            FetchJob(key=p, url=RCSB_FILE.format(pdb_id=p, fmt=fmt), dest=out_dir / f"{p}.{fmt}", decompress=True)
            for p in pending
        ]
        next_pending = []
        for result in fetch_all(jobs, concurrency=concurrency, **kwargs):
            results[result.key] = result
            if result.not_found:
                next_pending.append(result.key)
        pending = next_pending
    if manifest is not None:
        write_manifest(results.values(), manifest)
    return results


# ---------------------------------------------------------------- large files
async def _download_segmented(
    url: str,
    dest: Path,
    *,
    segments: int,
    min_segment: int,
    max_retries: int,
    timeout: float,
    user_agent: str,
    progress: Optional[Callable[[int, int], None]],
) -> FetchResult:
    import aiohttp

    started = time.monotonic()
    result = FetchResult(key=dest.name, url=url, ok=False)
    tmp = dest.with_name(f".{dest.name}.part")
    hdrs = {"User-Agent": user_agent}
    client_timeout = aiohttp.ClientTimeout(total=None, sock_connect=30, sock_read=timeout)
    async with aiohttp.ClientSession(timeout=client_timeout, headers=hdrs, trust_env=True) as session:
        async with session.head(url, allow_redirects=True) as head:
            result.status = head.status
            if head.status in (404, 410):
                result.not_found = True
                result.error = f"HTTP {head.status}"
                return result
            head.raise_for_status()
            size = int(head.headers.get("Content-Length", "0") or 0)
            ranged = head.headers.get("Accept-Ranges", "").lower() == "bytes" and size > 0

        n = max(1, min(segments, size // max(1, min_segment))) if ranged else 1
        bounds = [(i * size // n, (i + 1) * size // n - 1) for i in range(n)] if ranged else [(0, -1)]
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(tmp, "wb") as handle:
            if size:
                handle.truncate(size)
        done = 0

        async def segment(lo: int, hi: int) -> None:
            nonlocal done
            offset = lo
            for attempt in range(1, max_retries + 1):
                headers = {"Range": f"bytes={offset}-{hi}"} if hi >= 0 else {}
                try:
                    async with session.get(url, headers=headers) as response:
                        if hi >= 0 and response.status != 206:
                            raise aiohttp.ClientResponseError(
                                response.request_info, response.history,
                                status=response.status, message="range not honoured",
                            )
                        response.raise_for_status()
                        with open(tmp, "r+b") as handle:
                            handle.seek(offset)
                            async for chunk in response.content.iter_chunked(1 << 20):
                                handle.write(chunk)
                                offset += len(chunk)
                                done += len(chunk)
                                if progress is not None:
                                    progress(done, size)
                    if hi < 0 or offset > hi:
                        return
                except (aiohttp.ClientError, asyncio.TimeoutError, OSError) as exc:
                    if attempt == max_retries:
                        raise
                    LOGGER.warning("segment %d-%d: %s; resuming at %d", lo, hi, exc, offset)
                    await asyncio.sleep(min(60.0, 2 ** attempt) * (0.75 + 0.5 * random.random()))
            raise RuntimeError(f"segment {lo}-{hi} incomplete")

        try:
            await asyncio.gather(*(segment(lo, hi) for lo, hi in bounds))
        except Exception as exc:
            result.error = f"{type(exc).__name__}: {exc}"[:300]
            tmp.unlink(missing_ok=True)
            return result

    actual = tmp.stat().st_size
    if size and actual != size:
        result.error = f"size mismatch: {actual} of {size} bytes"
        tmp.unlink(missing_ok=True)
        return result
    os.replace(tmp, dest)
    result.ok = True
    result.path = str(dest)
    result.n_bytes = actual
    result.seconds = round(time.monotonic() - started, 3)
    return result


def download_large_file(
    url: str,
    dest: Path,
    *,
    segments: int = 16,
    min_segment: int = 8 << 20,
    max_retries: int = 8,
    timeout: float = 300.0,
    user_agent: str = DEFAULT_USER_AGENT,
    progress: Optional[Callable[[int, int], None]] = None,
) -> FetchResult:
    """Download one large file in parallel byte ranges.

    Many servers - the EBI FTP mirror among them - throttle each connection,
    so a multi-gigabyte archive over one stream can take hours; several
    ranges in parallel take minutes. A segment that drops resumes from the byte
    it reached. The file is assembled under a temporary name and renamed only
    once its size matches ``Content-Length``. Servers that do not honour ranges
    get a single stream.
    """
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return FetchResult(key=dest.name, url=url, ok=True, path=str(dest),
                           n_bytes=dest.stat().st_size, from_cache=True)
    return asyncio.run(
        _download_segmented(
            url, dest, segments=segments, min_segment=min_segment, max_retries=max_retries,
            timeout=timeout, user_agent=user_agent, progress=progress,
        )
    )
