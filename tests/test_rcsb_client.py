"""Tests for the RCSB client's pagination, caching and failure handling.

Network calls are stubbed, so these run anywhere. The behaviours under test are
the ones that silently corrupt a dataset when wrong: truncated pagination, cache
collisions, and errors swallowed into an empty result.
"""

import json
from pathlib import Path

import pytest
import requests

from cryptic_ip.database.rcsb_client import (
    MAX_ROWS_PER_PAGE,
    RcsbClient,
    RcsbUnavailableError,
    ResponseCache,
    sha256_bytes,
    sha256_file,
)


class FakeResponse:
    """Stand-in for a requests.Response."""

    def __init__(self, payload=None, status_code=200, content=b"", headers=None):
        self._payload = payload
        self.status_code = status_code
        self.content = content
        self.headers = headers or {}
        self.text = json.dumps(payload) if payload is not None else ""

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")


def test_response_cache_round_trip(tmp_path: Path):
    cache = ResponseCache(tmp_path)
    key = cache.key("search", {"a": 1})
    assert cache.get(key) is None
    cache.put(key, {"value": 42})
    assert cache.get(key) == {"value": 42}


def test_response_cache_keys_differ_for_different_requests(tmp_path: Path):
    cache = ResponseCache(tmp_path)
    assert cache.key("search", {"a": 1}) != cache.key("search", {"a": 2})
    assert cache.key("search", {"a": 1}) != cache.key("entry", {"a": 1})


def test_response_cache_discards_corrupt_entries(tmp_path: Path):
    cache = ResponseCache(tmp_path)
    key = cache.key("search", {"a": 1})
    cache.put(key, {"value": 1})
    corrupt = tmp_path / key[:2] / f"{key}.json"
    corrupt.write_text("not json", encoding="utf-8")
    assert cache.get(key) is None
    assert not corrupt.exists()


def test_disabled_cache_is_a_no_op():
    cache = ResponseCache(None)
    assert not cache.enabled
    cache.put(cache.key("x", {}), {"a": 1})
    assert cache.get(cache.key("x", {})) is None


def test_search_pages_through_the_entire_result_set(monkeypatch):
    """A result set larger than one page must not be silently truncated."""
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    total = 12_000
    identifiers = [f"E{i:05d}" for i in range(total)]

    def fake_request(method, url, **kwargs):
        start = kwargs["json"]["request_options"]["paginate"]["start"]
        rows = kwargs["json"]["request_options"]["paginate"]["rows"]
        return FakeResponse(
            {"total_count": total, "result_set": identifiers[start : start + rows]}
        )

    monkeypatch.setattr(client, "_request", fake_request)
    found = client.search({"type": "terminal"}, rows_per_page=5000)
    assert len(found) == total
    assert found[0] == "E00000"


def test_search_clamps_rows_to_the_api_maximum(monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    seen = {}

    def fake_request(method, url, **kwargs):
        seen["rows"] = kwargs["json"]["request_options"]["paginate"]["rows"]
        return FakeResponse({"total_count": 1, "result_set": ["AAAA"]})

    monkeypatch.setattr(client, "_request", fake_request)
    client.search({"type": "terminal"}, rows_per_page=999_999)
    assert seen["rows"] == MAX_ROWS_PER_PAGE


def test_search_deduplicates_and_uppercases(monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)

    def fake_request(method, url, **kwargs):
        start = kwargs["json"]["request_options"]["paginate"]["start"]
        if start > 0:
            return FakeResponse({"total_count": 2, "result_set": []})
        return FakeResponse({"total_count": 2, "result_set": ["abcd", "ABCD", "efgh"]})

    monkeypatch.setattr(client, "_request", fake_request)
    assert client.search({"type": "terminal"}) == ["ABCD", "EFGH"]


def test_search_handles_verbose_result_objects(monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)

    def fake_request(method, url, **kwargs):
        start = kwargs["json"]["request_options"]["paginate"]["start"]
        if start > 0:
            return FakeResponse({"total_count": 1, "result_set": []})
        return FakeResponse({"total_count": 1, "result_set": [{"identifier": "1ZY7"}]})

    monkeypatch.setattr(client, "_request", fake_request)
    assert client.search({"type": "terminal"}) == ["1ZY7"]


def test_component_query_builds_one_node_per_identifier(monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    captured = {}

    def fake_request(method, url, **kwargs):
        captured["query"] = kwargs["json"]["query"]
        return FakeResponse({"total_count": 0, "result_set": []})

    monkeypatch.setattr(client, "_request", fake_request)
    client.search_entries_with_components(["ihp", "I3P"], experimental_methods=["X-RAY DIFFRACTION"])

    query = captured["query"]
    assert query["logical_operator"] == "and"
    component_node = query["nodes"][0]
    values = {node["parameters"]["value"] for node in component_node["nodes"]}
    assert values == {"IHP", "I3P"}


def test_component_query_without_identifiers_returns_nothing(monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    assert client.search_entries_with_components([]) == []


def test_decoy_query_negates_every_ligand(monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    captured = {}

    def fake_request(method, url, **kwargs):
        captured["query"] = kwargs["json"]["query"]
        return FakeResponse({"total_count": 0, "result_set": []})

    monkeypatch.setattr(client, "_request", fake_request)
    client.search_decoy_entries(exclude_comp_ids=["IHP"], max_results=5)

    negated = [
        node
        for node in captured["query"]["nodes"]
        if node["parameters"].get("negation") is True
    ]
    assert len(negated) == 1
    assert negated[0]["parameters"]["value"] == "IHP"


def test_request_failure_raises_an_actionable_error(monkeypatch):
    client = RcsbClient(cache_dir=None, max_retries=2, min_interval_s=0.0)

    def always_fail(method, url, timeout=None, **kwargs):
        raise requests.ConnectionError("no route to host")

    monkeypatch.setattr(client.session, "request", always_fail)
    monkeypatch.setattr("time.sleep", lambda *_: None)

    with pytest.raises(RcsbUnavailableError, match="failed after 2 attempts"):
        client._request("GET", "https://example.invalid")
    assert client.stats.errors >= 2


def test_download_prefers_mmcif_and_records_a_checksum(tmp_path: Path, monkeypatch):
    """mmCIF first, because large assemblies have no legacy PDB file at all."""
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    payload = b"data_1ZY7\n#\n"
    requested = []

    def fake_request(method, url, allow_404=False, **kwargs):
        requested.append(url)
        return FakeResponse(content=__import__("gzip").compress(payload))

    monkeypatch.setattr(client, "_request", fake_request)
    record = client.download_structure("1zy7", tmp_path)

    assert record is not None
    assert requested[0].endswith("1ZY7.cif.gz")
    assert Path(record.path).read_bytes() == payload
    assert record.sha256 == sha256_bytes(payload)
    assert record.from_cache is False


def test_download_falls_back_to_pdb_when_mmcif_is_absent(tmp_path: Path, monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)

    def fake_request(method, url, allow_404=False, **kwargs):
        if url.endswith(".cif.gz"):
            return None
        return FakeResponse(content=__import__("gzip").compress(b"ATOM\n"))

    monkeypatch.setattr(client, "_request", fake_request)
    record = client.download_structure("1ABC", tmp_path)
    assert record is not None
    assert record.path.endswith(".pdb")


def test_download_returns_none_when_no_format_exists(tmp_path: Path, monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    monkeypatch.setattr(client, "_request", lambda *a, **k: None)
    assert client.download_structure("ZZZZ", tmp_path) is None


def test_existing_file_is_reused_and_flagged_as_cached(tmp_path: Path, monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    target = tmp_path / "1ZY7.cif"
    target.write_bytes(b"data_1ZY7\n")

    def fail(*args, **kwargs):
        raise AssertionError("cached file must not be re-downloaded")

    monkeypatch.setattr(client, "_request", fail)
    record = client.download_structure("1ZY7", tmp_path)
    assert record is not None
    assert record.from_cache is True
    assert record.sha256 == sha256_file(target)


def test_entry_metadata_batches_and_merges(monkeypatch):
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    calls = []

    def fake_request(method, url, **kwargs):
        ids = kwargs["json"]["variables"]["ids"]
        calls.append(len(ids))
        return FakeResponse(
            {"data": {"entries": [{"rcsb_id": identifier} for identifier in ids]}}
        )

    monkeypatch.setattr(client, "_request", fake_request)
    metadata = client.fetch_entry_metadata([f"E{i:04d}" for i in range(250)], batch_size=100)
    assert len(metadata) == 250
    assert calls == [100, 100, 50]


def test_entry_metadata_survives_a_failed_batch(monkeypatch):
    """One bad batch must not lose the whole collection's annotations."""
    client = RcsbClient(cache_dir=None, min_interval_s=0.0)
    state = {"calls": 0}

    def fake_request(method, url, **kwargs):
        state["calls"] += 1
        if state["calls"] == 1:
            raise RcsbUnavailableError("gateway timeout")
        ids = kwargs["json"]["variables"]["ids"]
        return FakeResponse(
            {"data": {"entries": [{"rcsb_id": identifier} for identifier in ids]}}
        )

    monkeypatch.setattr(client, "_request", fake_request)
    metadata = client.fetch_entry_metadata(["A001", "B002"], batch_size=1)
    assert len(metadata) == 1
