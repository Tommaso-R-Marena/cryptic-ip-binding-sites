"""Rigorous RCSB PDB client for exhaustive ground-truth data collection.

Design goals
------------
The validation dataset is the foundation of every benchmark in this project, so
the client that builds it has to be exhaustive, restartable and auditable:

**Exhaustive pagination.** The RCSB search API caps a single response at 10 000
rows. A one-shot request therefore silently truncates large result sets. This
client reads ``total_count`` and pages until the result set is exhausted.

**Batched metadata via GraphQL.** Per-entry REST calls cost one request per
entry plus one per polymer entity, which is thousands of requests for a
proteome-scale query. The GraphQL endpoint returns entry, polymer-entity,
non-polymer-entity and instance-level annotations for up to a few hundred
identifiers per request, cutting request counts by two orders of magnitude and
making a full collection feasible.

**Content-addressed cache.** Every response is cached on disk keyed by a hash of
the request. Re-running a collection is free, interrupted runs resume, and the
cache doubles as the provenance record of exactly what the API returned.

**Auditable downloads.** Structure files are written atomically, hashed with
SHA-256, and recorded in a manifest together with the retrieval timestamp so a
dataset can be re-verified byte-for-byte later.

Network access is required. Every public method raises
:class:`RcsbUnavailableError` with an actionable message when the API cannot be
reached, so callers can degrade gracefully instead of producing a silently empty
dataset.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, Optional, Sequence

import requests

LOGGER = logging.getLogger(__name__)

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
RCSB_REST_URL = "https://data.rcsb.org/rest/v1/core"
RCSB_FILES_URL = "https://files.rcsb.org/download"

#: The search API rejects ``rows`` above this value.
MAX_ROWS_PER_PAGE = 10_000

#: GraphQL identifier batch size. Large batches reduce request counts but risk
#: gateway timeouts; 100 is a reliable compromise for entry-level queries.
GRAPHQL_BATCH_SIZE = 100


class RcsbUnavailableError(RuntimeError):
    """Raised when the RCSB API cannot be reached or returns a fatal error."""


@dataclass
class RequestStats:
    """Counters describing the network traffic of a collection run."""

    requests: int = 0
    cache_hits: int = 0
    retries: int = 0
    bytes_downloaded: int = 0
    errors: int = 0

    def to_dict(self) -> Dict[str, int]:
        """Return a JSON-serialisable snapshot of the counters."""
        return {
            "requests": self.requests,
            "cache_hits": self.cache_hits,
            "retries": self.retries,
            "bytes_downloaded": self.bytes_downloaded,
            "errors": self.errors,
        }


@dataclass
class DownloadRecord:
    """Provenance record for one retrieved structure file."""

    identifier: str
    path: str
    url: str
    sha256: str
    n_bytes: int
    retrieved_at_utc: str
    from_cache: bool

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serialisable representation."""
        return {
            "identifier": self.identifier,
            "path": self.path,
            "url": self.url,
            "sha256": self.sha256,
            "n_bytes": self.n_bytes,
            "retrieved_at_utc": self.retrieved_at_utc,
            "from_cache": self.from_cache,
        }


def sha256_bytes(payload: bytes) -> str:
    """Return the hex SHA-256 digest of a byte string."""
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    """Return the hex SHA-256 digest of a file, read in chunks."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


class ResponseCache:
    """Content-addressed on-disk cache for API responses.

    Keys are SHA-256 digests of ``(kind, payload)`` so a cache hit is only
    possible for a byte-identical request. Values are JSON documents.
    """

    def __init__(self, root: Optional[Path]) -> None:
        """Initialise the cache.

        Args:
            root: Directory for cache entries. ``None`` disables caching.
        """
        self.root = Path(root) if root is not None else None
        if self.root is not None:
            self.root.mkdir(parents=True, exist_ok=True)

    @property
    def enabled(self) -> bool:
        """Whether the cache is active."""
        return self.root is not None

    def key(self, kind: str, payload: Any) -> str:
        """Compute the cache key for a request."""
        blob = json.dumps({"kind": kind, "payload": payload}, sort_keys=True, default=str)
        return sha256_bytes(blob.encode("utf-8"))

    def _path(self, key: str) -> Path:
        assert self.root is not None
        return self.root / key[:2] / f"{key}.json"

    def get(self, key: str) -> Optional[Any]:
        """Return the cached value for ``key``, or ``None`` on a miss."""
        if not self.enabled:
            return None
        path = self._path(key)
        if not path.exists():
            return None
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            LOGGER.warning("Discarding corrupt cache entry %s: %s", path, exc)
            path.unlink(missing_ok=True)
            return None

    def put(self, key: str, value: Any) -> None:
        """Store ``value`` under ``key`` using an atomic replace."""
        if not self.enabled:
            return
        path = self._path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(value, sort_keys=True), encoding="utf-8")
        os.replace(tmp, path)


class RcsbClient:
    """HTTP client for the RCSB search, data and file services.

    Args:
        cache_dir: Directory for the response cache. ``None`` disables caching.
        max_retries: Attempts per request before raising.
        timeout_s: Per-request timeout in seconds.
        min_interval_s: Minimum delay between requests (politeness throttle).
        user_agent: ``User-Agent`` header; identifying the client is requested
            by the RCSB terms of use.
    """

    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        *,
        max_retries: int = 5,
        timeout_s: float = 60.0,
        min_interval_s: float = 0.05,
        user_agent: str = "cryptic-ip-binding-sites/1.0 (structural bioinformatics pipeline)",
    ) -> None:
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": user_agent, "Accept": "application/json"})
        self.cache = ResponseCache(cache_dir)
        self.max_retries = max(1, int(max_retries))
        self.timeout_s = float(timeout_s)
        self.min_interval_s = float(min_interval_s)
        self.stats = RequestStats()
        self._last_request_at = 0.0

    # ------------------------------------------------------------------ HTTP

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_request_at
        if elapsed < self.min_interval_s:
            time.sleep(self.min_interval_s - elapsed)
        self._last_request_at = time.monotonic()

    def _request(
        self,
        method: str,
        url: str,
        *,
        allow_404: bool = False,
        **kwargs: Any,
    ) -> Optional[requests.Response]:
        """Issue a request with retry, backoff and ``Retry-After`` handling.

        Args:
            method: HTTP method.
            url: Absolute URL.
            allow_404: Return ``None`` instead of raising when the server
                answers 404 (used for optional resources).
            **kwargs: Passed to :meth:`requests.Session.request`.

        Returns:
            The response, or ``None`` for a tolerated 404.

        Raises:
            RcsbUnavailableError: When all attempts fail.
        """
        last_error: Optional[BaseException] = None
        for attempt in range(1, self.max_retries + 1):
            self._throttle()
            try:
                response = self.session.request(method, url, timeout=self.timeout_s, **kwargs)
                self.stats.requests += 1
            except requests.RequestException as exc:
                last_error = exc
                self.stats.errors += 1
                self.stats.retries += 1
                backoff = min(2.0**attempt, 30.0)
                LOGGER.warning(
                    "%s %s failed (attempt %d/%d): %s; retrying in %.1fs",
                    method,
                    url,
                    attempt,
                    self.max_retries,
                    exc,
                    backoff,
                )
                time.sleep(backoff)
                continue

            if response.status_code == 404 and allow_404:
                return None
            if response.status_code == 429:
                retry_after = response.headers.get("Retry-After")
                wait = float(retry_after) if retry_after else min(2.0**attempt, 30.0)
                self.stats.retries += 1
                LOGGER.warning("Rate limited by %s; sleeping %.1fs", url, wait)
                time.sleep(wait)
                continue
            if 500 <= response.status_code < 600:
                self.stats.retries += 1
                backoff = min(2.0**attempt, 30.0)
                LOGGER.warning(
                    "Server error %d from %s (attempt %d/%d); retrying in %.1fs",
                    response.status_code,
                    url,
                    attempt,
                    self.max_retries,
                    backoff,
                )
                time.sleep(backoff)
                continue

            try:
                response.raise_for_status()
            except requests.HTTPError as exc:
                self.stats.errors += 1
                raise RcsbUnavailableError(
                    f"{method} {url} returned {response.status_code}: {response.text[:200]}"
                ) from exc
            return response

        raise RcsbUnavailableError(
            f"{method} {url} failed after {self.max_retries} attempts. "
            "Check network connectivity and any outbound proxy policy."
        ) from last_error

    def _post_json(self, url: str, payload: Mapping[str, Any], *, cache_kind: str) -> Any:
        key = self.cache.key(cache_kind, {"url": url, "payload": payload})
        cached = self.cache.get(key)
        if cached is not None:
            self.stats.cache_hits += 1
            return cached
        response = self._request("POST", url, json=payload)
        assert response is not None
        data = response.json()
        self.cache.put(key, data)
        return data

    def _get_json(self, url: str, *, cache_kind: str, allow_404: bool = True) -> Optional[Any]:
        key = self.cache.key(cache_kind, {"url": url})
        cached = self.cache.get(key)
        if cached is not None:
            self.stats.cache_hits += 1
            return cached if cached != {"__missing__": True} else None
        response = self._request("GET", url, allow_404=allow_404)
        if response is None:
            self.cache.put(key, {"__missing__": True})
            return None
        data = response.json()
        self.cache.put(key, data)
        return data

    # ---------------------------------------------------------------- search

    def search(
        self,
        query: Mapping[str, Any],
        *,
        return_type: str = "entry",
        rows_per_page: int = 5_000,
        max_results: Optional[int] = None,
        request_options: Optional[Mapping[str, Any]] = None,
    ) -> List[str]:
        """Run a search query and page through the *entire* result set.

        Args:
            query: A RCSB search query node.
            return_type: ``"entry"``, ``"polymer_entity"``,
                ``"non_polymer_entity"``, ``"assembly"`` or
                ``"polymer_instance"``.
            rows_per_page: Rows per request; clamped to the API maximum.
            max_results: Optional cap on total identifiers returned.
            request_options: Extra request options merged into each page.

        Returns:
            Sorted, de-duplicated identifiers.

        Raises:
            RcsbUnavailableError: When the search service is unreachable.
        """
        rows = max(1, min(int(rows_per_page), MAX_ROWS_PER_PAGE))
        identifiers: List[str] = []
        seen: set = set()
        start = 0
        total: Optional[int] = None

        while True:
            options: Dict[str, Any] = {
                "paginate": {"start": start, "rows": rows},
                "results_verbosity": "compact",
                # Identifier sort makes pagination deterministic; relevance
                # scores are not stable across pages.
                "sort": [{"sort_by": "rcsb_entry_info.deposited_atom_count", "direction": "desc"}],
            }
            if request_options:
                options.update(dict(request_options))
            payload = {
                "query": dict(query),
                "return_type": return_type,
                "request_options": options,
            }
            data = self._post_json(RCSB_SEARCH_URL, payload, cache_kind="search")
            if not data:
                break

            if total is None:
                total = int(data.get("total_count", 0))
                LOGGER.info("Search reports %d %s hits", total, return_type)

            page = _extract_identifiers(data)
            if not page:
                break
            for identifier in page:
                upper = identifier.upper()
                if upper not in seen:
                    seen.add(upper)
                    identifiers.append(upper)
            start += rows
            if max_results is not None and len(identifiers) >= max_results:
                identifiers = identifiers[:max_results]
                break
            if total is not None and start >= total:
                break

        LOGGER.info("Collected %d unique %s identifiers", len(identifiers), return_type)
        return sorted(identifiers)

    def search_chemcomp_full_text(self, term: str, *, max_results: int = 500) -> List[str]:
        """Full-text search of the chemical component dictionary.

        Args:
            term: Search phrase, e.g. ``"inositol phosphate"``.
            max_results: Maximum component identifiers to return.

        Returns:
            Chemical component identifiers.
        """
        query = {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": term},
        }
        return self.search(
            query,
            return_type="mol_definition",
            rows_per_page=min(max_results, MAX_ROWS_PER_PAGE),
            max_results=max_results,
            request_options={"sort": [{"sort_by": "score", "direction": "desc"}]},
        )

    def search_entries_with_components(
        self,
        comp_ids: Sequence[str],
        *,
        experimental_methods: Optional[Sequence[str]] = None,
        max_resolution: Optional[float] = None,
        max_results: Optional[int] = None,
    ) -> List[str]:
        """Find every entry containing at least one of the given components.

        Args:
            comp_ids: Chemical component identifiers to match.
            experimental_methods: Restrict to these ``exptl.method`` values.
                ``None`` (the default) accepts **all** methods, which is
                deliberate: restricting to X-ray discards cryo-EM structures of
                large IP-dependent assemblies.
            max_resolution: Optional upper bound on reported resolution (Å).
            max_results: Optional cap on identifiers returned.

        Returns:
            Entry identifiers.
        """
        if not comp_ids:
            return []

        component_node = {
            "type": "group",
            "logical_operator": "or",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": (
                            "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"
                        ),
                        "operator": "exact_match",
                        "value": comp_id.upper(),
                    },
                }
                for comp_id in sorted({c.upper() for c in comp_ids})
            ],
        }

        nodes: List[Dict[str, Any]] = [component_node]
        if experimental_methods:
            nodes.append(
                {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "exptl.method",
                                "operator": "exact_match",
                                "value": method,
                            },
                        }
                        for method in experimental_methods
                    ],
                }
            )
        if max_resolution is not None:
            nodes.append(
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": float(max_resolution),
                    },
                }
            )

        query: Dict[str, Any] = (
            nodes[0]
            if len(nodes) == 1
            else {"type": "group", "logical_operator": "and", "nodes": nodes}
        )
        return self.search(query, return_type="entry", max_results=max_results)

    def search_decoy_entries(
        self,
        *,
        exclude_comp_ids: Sequence[str],
        max_resolution: float = 2.5,
        min_residues: int = 80,
        max_results: int = 2_000,
    ) -> List[str]:
        """Find high-resolution entries that contain **no** inositol phosphate.

        These entries supply protein-level negatives: pockets drawn from proteins
        with no known IP ligand. Without them, every negative pocket comes from an
        IP-binding protein, which makes the benchmark easier than a proteome
        screen and inflates apparent performance.

        Args:
            exclude_comp_ids: Components whose presence disqualifies an entry.
            max_resolution: Resolution ceiling in Å.
            min_residues: Minimum deposited polymer residue count.
            max_results: Maximum identifiers to return.

        Returns:
            Entry identifiers with no matching component.
        """
        nodes: List[Dict[str, Any]] = [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": float(max_resolution),
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_polymer_monomer_count",
                    "operator": "greater_or_equal",
                    "value": int(min_residues),
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                    "operator": "greater_or_equal",
                    "value": 1,
                },
            },
        ]
        for comp_id in sorted({c.upper() for c in exclude_comp_ids}):
            nodes.append(
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": (
                            "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"
                        ),
                        "operator": "exact_match",
                        "value": comp_id,
                        "negation": True,
                    },
                }
            )
        query = {"type": "group", "logical_operator": "and", "nodes": nodes}
        return self.search(query, return_type="entry", max_results=max_results)

    # ------------------------------------------------------------ metadata

    def fetch_chemcomp(self, comp_id: str) -> Optional[Dict[str, Any]]:
        """Fetch a chemical component definition.

        Args:
            comp_id: Component identifier.

        Returns:
            The component document, or ``None`` when unknown.
        """
        return self._get_json(
            f"{RCSB_REST_URL}/chemcomp/{comp_id.upper()}", cache_kind="chemcomp"
        )

    def fetch_entry_metadata(
        self, pdb_ids: Sequence[str], *, batch_size: int = GRAPHQL_BATCH_SIZE
    ) -> Dict[str, Dict[str, Any]]:
        """Fetch rich entry metadata for many entries using batched GraphQL.

        One request covers up to ``batch_size`` entries and returns entry-level
        experimental details, polymer entity annotations (UniProt cross
        references, source organism, EC numbers) and non-polymer entity
        instance annotations - everything the dataset builder needs.

        Args:
            pdb_ids: Entry identifiers.
            batch_size: Identifiers per GraphQL request.

        Returns:
            Mapping from upper-case entry identifier to its metadata document.
        """
        unique = sorted({pdb_id.upper() for pdb_id in pdb_ids})
        out: Dict[str, Dict[str, Any]] = {}
        for chunk in _chunked(unique, max(1, int(batch_size))):
            payload = {"query": _ENTRY_GRAPHQL_QUERY, "variables": {"ids": list(chunk)}}
            try:
                data = self._post_json(RCSB_GRAPHQL_URL, payload, cache_kind="graphql_entries")
            except RcsbUnavailableError as exc:
                LOGGER.warning("GraphQL batch failed (%d ids): %s", len(chunk), exc)
                continue
            if data.get("errors"):
                LOGGER.warning("GraphQL reported errors: %s", str(data["errors"])[:300])
            for entry in (data.get("data") or {}).get("entries") or []:
                if not entry:
                    continue
                identifier = str(entry.get("rcsb_id", "")).upper()
                if identifier:
                    out[identifier] = entry
        LOGGER.info("Fetched metadata for %d/%d entries", len(out), len(unique))
        return out

    # ------------------------------------------------------------- download

    def download_structure(
        self,
        pdb_id: str,
        out_dir: Path,
        *,
        prefer_format: str = "cif",
        decompress: bool = True,
    ) -> Optional[DownloadRecord]:
        """Download a structure file, preferring mmCIF over legacy PDB.

        mmCIF is preferred because entries with more than 62 chains or 99 999
        atoms - common for the large assemblies where inositol phosphates act as
        structural cofactors - have no legacy PDB format file at all. Falling
        back to PDB only when mmCIF is unavailable keeps those entries in the
        dataset instead of silently dropping them.

        Args:
            pdb_id: Entry identifier.
            out_dir: Destination directory.
            prefer_format: ``"cif"`` or ``"pdb"``.
            decompress: Store the decompressed file (gzip is used in transit).

        Returns:
            A :class:`DownloadRecord`, or ``None`` when no format could be
            retrieved.
        """
        identifier = pdb_id.upper()
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        formats = ("cif", "pdb") if prefer_format == "cif" else ("pdb", "cif")
        for fmt in formats:
            suffix = f".{fmt}"
            target = out_dir / f"{identifier}{suffix}"
            if target.exists() and target.stat().st_size > 0:
                return DownloadRecord(
                    identifier=identifier,
                    path=str(target),
                    url=f"{RCSB_FILES_URL}/{identifier}{suffix}.gz",
                    sha256=sha256_file(target),
                    n_bytes=target.stat().st_size,
                    retrieved_at_utc=datetime.fromtimestamp(
                        target.stat().st_mtime, tz=timezone.utc
                    ).isoformat(),
                    from_cache=True,
                )

            url = f"{RCSB_FILES_URL}/{identifier}{suffix}.gz"
            try:
                response = self._request("GET", url, allow_404=True, stream=True)
            except RcsbUnavailableError as exc:
                LOGGER.warning("Download failed for %s (%s): %s", identifier, fmt, exc)
                continue
            if response is None:
                LOGGER.debug("No %s file for %s", fmt, identifier)
                continue

            raw = response.content
            self.stats.bytes_downloaded += len(raw)
            try:
                content = gzip.decompress(raw) if decompress else raw
            except (OSError, EOFError):
                content = raw

            tmp = target.with_suffix(target.suffix + ".part")
            tmp.write_bytes(content)
            os.replace(tmp, target)
            return DownloadRecord(
                identifier=identifier,
                path=str(target),
                url=url,
                sha256=sha256_bytes(content),
                n_bytes=len(content),
                retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
                from_cache=False,
            )

        LOGGER.warning("No structure file available for %s", identifier)
        return None

    def download_structures(
        self,
        pdb_ids: Sequence[str],
        out_dir: Path,
        *,
        max_workers: int = 4,
        prefer_format: str = "cif",
    ) -> List[DownloadRecord]:
        """Download many structures concurrently.

        Args:
            pdb_ids: Entry identifiers.
            out_dir: Destination directory.
            max_workers: Concurrent download threads. Kept modest to respect
                the RCSB fair-use policy.
            prefer_format: Preferred file format.

        Returns:
            One record per successfully retrieved entry.
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed

        unique = sorted({pdb_id.upper() for pdb_id in pdb_ids})
        records: List[DownloadRecord] = []
        workers = max(1, min(int(max_workers), 8))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(
                    self.download_structure, pdb_id, out_dir, prefer_format=prefer_format
                ): pdb_id
                for pdb_id in unique
            }
            for future in as_completed(futures):
                pdb_id = futures[future]
                try:
                    record = future.result()
                except Exception as exc:  # noqa: BLE001 - one bad entry must not abort
                    LOGGER.warning("Download of %s raised: %s", pdb_id, exc)
                    continue
                if record is not None:
                    records.append(record)
        records.sort(key=lambda rec: rec.identifier)
        LOGGER.info("Retrieved %d/%d structure files", len(records), len(unique))
        return records


def _extract_identifiers(data: Mapping[str, Any]) -> List[str]:
    """Extract identifiers from a search response in either verbosity format."""
    result_set = data.get("result_set") or []
    identifiers: List[str] = []
    for item in result_set:
        if isinstance(item, str):
            identifiers.append(item)
        elif isinstance(item, Mapping):
            value = item.get("identifier")
            if value:
                identifiers.append(str(value))
    return identifiers


def _chunked(items: Sequence[str], size: int) -> Iterator[Sequence[str]]:
    """Yield consecutive chunks of at most ``size`` items."""
    for start in range(0, len(items), size):
        yield items[start : start + size]


#: GraphQL document retrieving the entry metadata the dataset builder needs.
#:
#: Only long-stable schema fields are requested. Metadata is treated as
#: best-effort enrichment throughout the pipeline: ligand geometry, burial and
#: labels are computed from the coordinate file, which is authoritative, so a
#: GraphQL schema change degrades annotation richness without invalidating any
#: measurement.
_ENTRY_GRAPHQL_QUERY = """
query EntryMetadata($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    exptl { method }
    refine { ls_R_factor_R_free ls_R_factor_R_work }
    rcsb_accession_info { initial_release_date deposit_date }
    rcsb_entry_info {
      resolution_combined
      experimental_method
      deposited_polymer_monomer_count
      deposited_atom_count
      polymer_entity_count_protein
      nonpolymer_entity_count
    }
    polymer_entities {
      rcsb_id
      entity_poly { rcsb_sample_sequence_length }
      rcsb_polymer_entity { pdbx_description }
      rcsb_polymer_entity_container_identifiers {
        auth_asym_ids
        reference_sequence_identifiers { database_name database_accession }
      }
      rcsb_entity_source_organism { ncbi_scientific_name ncbi_taxonomy_id }
      rcsb_ec_lineage { id name }
    }
    nonpolymer_entities {
      rcsb_id
      nonpolymer_comp { chem_comp { id name formula formula_weight } }
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}
"""
