"""Tests for the async download layer, against a local HTTP server."""

from __future__ import annotations

import asyncio
import gzip
import json
import threading
from collections import Counter
from pathlib import Path

import pytest

aiohttp = pytest.importorskip("aiohttp")
from aiohttp import web  # noqa: E402

from cryptic_ip.database import async_fetch  # noqa: E402
from cryptic_ip.database.async_fetch import FetchJob, fetch_all, write_manifest  # noqa: E402

PDB = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n"
    "ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 90.00           C\n"
    "END\n"
).encode()
CIF = b"data_1ABC\nloop_\n_atom_site.id\n1\n"


class Server:
    """A local server whose routes misbehave on purpose."""

    def __init__(self):
        self.hits = Counter()
        self.flaky_left = 2
        app = web.Application()
        app.router.add_get("/{name}", self.handle)
        app.router.add_get("/download/{name}", self.handle)
        self.runner = web.AppRunner(app)
        self.loop = asyncio.new_event_loop()
        self.thread = threading.Thread(target=self.loop.run_forever, daemon=True)

    async def handle(self, request):
        name = request.match_info["name"]
        self.hits[name] += 1
        if name == "ok.pdb":
            return web.Response(body=PDB)
        if name == "ok.pdb.gz":
            return web.Response(body=gzip.compress(PDB))
        if name == "GOOD.pdb":
            return web.Response(body=PDB)
        if name == "ONLYPDB.pdb.gz":
            return web.Response(body=gzip.compress(PDB))
        if name == "ok.cif.gz":
            return web.Response(body=gzip.compress(CIF))
        if name == "flaky.pdb":
            if self.flaky_left > 0:
                self.flaky_left -= 1
                return web.Response(status=503)
            return web.Response(body=PDB)
        if name == "ratelimited.pdb":
            if self.hits[name] == 1:
                return web.Response(status=429, headers={"Retry-After": "0"})
            return web.Response(body=PDB)
        if name == "truncated.pdb.gz":
            return web.Response(body=gzip.compress(PDB)[:-10])
        if name == "noend.pdb":
            return web.Response(body=PDB.replace(b"END\n", b""))
        if name == "html.pdb":
            return web.Response(body=b"<html>Service unavailable</html>")
        if name == "forbidden.pdb":
            return web.Response(status=403)
        if name.startswith("api-"):
            accession = name[4:]
            base = str(request.url.origin())
            return web.json_response([
                {"uniprotAccession": accession + "-2", "entryId": f"AF-{accession}-2-F1",
                 "pdbUrl": f"{base}/wrong.pdb"},
                {"uniprotAccession": accession, "entryId": f"AF-{accession}-F1",
                 "pdbUrl": f"{base}/AF-{accession}-F1-model_v6.pdb"},
            ])
        if name.startswith("AF-") and name.endswith(".pdb"):
            return web.Response(body=PDB)
        return web.Response(status=404)

    def __enter__(self):
        self.thread.start()
        asyncio.run_coroutine_threadsafe(self.runner.setup(), self.loop).result()
        site = web.TCPSite(self.runner, "127.0.0.1", 0)
        asyncio.run_coroutine_threadsafe(site.start(), self.loop).result()
        port = self.runner.addresses[0][1]
        self.base = f"http://127.0.0.1:{port}"
        return self

    def __exit__(self, *exc):
        asyncio.run_coroutine_threadsafe(self.runner.cleanup(), self.loop).result()
        self.loop.call_soon_threadsafe(self.loop.stop)


@pytest.fixture
def server():
    with Server() as s:
        yield s


FAST = dict(backoff=0.01, max_backoff=0.02, max_retries=3)


def test_download_decompress_validate_and_resume(server, tmp_path):
    jobs = [
        FetchJob("a", f"{server.base}/ok.pdb", tmp_path / "a.pdb"),
        FetchJob("b", f"{server.base}/ok.pdb.gz", tmp_path / "b.pdb"),
        FetchJob("c", f"{server.base}/ok.cif.gz", tmp_path / "c.cif"),
    ]
    results = fetch_all(jobs, **FAST)
    assert all(r.ok and not r.from_cache for r in results)
    assert (tmp_path / "b.pdb").read_bytes() == PDB  # stored decompressed
    assert [r.key for r in results] == ["a", "b", "c"]  # job order

    again = fetch_all(jobs, **FAST)
    assert all(r.from_cache for r in again)
    assert server.hits["ok.pdb"] == 1  # resumed, not re-downloaded


def test_transient_failures_are_retried(server, tmp_path):
    results = fetch_all(
        [
            FetchJob("f", f"{server.base}/flaky.pdb", tmp_path / "f.pdb"),
            FetchJob("r", f"{server.base}/ratelimited.pdb", tmp_path / "r.pdb"),
        ],
        **FAST,
    )
    assert all(r.ok for r in results)
    assert results[0].attempts == 3
    assert results[1].attempts == 2


def test_absent_is_not_retried_and_forbidden_is_not_retried(server, tmp_path):
    missing, forbidden = fetch_all(
        [
            FetchJob("m", f"{server.base}/missing.pdb", tmp_path / "m.pdb"),
            FetchJob("x", f"{server.base}/forbidden.pdb", tmp_path / "x.pdb"),
        ],
        **FAST,
    )
    assert missing.not_found and not missing.ok and missing.attempts == 1
    assert not forbidden.ok and not forbidden.not_found and forbidden.attempts == 1
    assert not (tmp_path / "m.pdb").exists()


@pytest.mark.parametrize(
    "route, reason",
    [("truncated.pdb.gz", "corrupt gzip"), ("noend.pdb", "no END"), ("html.pdb", "no ATOM")],
)
def test_bad_content_is_never_written(server, tmp_path, route, reason):
    (result,) = fetch_all([FetchJob("t", f"{server.base}/{route}", tmp_path / "t.pdb")], **FAST)
    assert not result.ok
    assert reason in result.error
    assert result.attempts == FAST["max_retries"]
    assert not (tmp_path / "t.pdb").exists()
    assert not list(tmp_path.glob(".*.part"))


def test_corrupt_cached_file_is_replaced(server, tmp_path):
    dest = tmp_path / "a.pdb"
    dest.write_bytes(PDB[:40])  # a download that was cut off
    (result,) = fetch_all([FetchJob("a", f"{server.base}/ok.pdb", dest)], **FAST)
    assert result.ok and not result.from_cache
    assert dest.read_bytes() == PDB


def test_manifest(server, tmp_path):
    results = fetch_all(
        [
            FetchJob("a", f"{server.base}/ok.pdb", tmp_path / "a.pdb"),
            FetchJob("m", f"{server.base}/missing.pdb", tmp_path / "m.pdb"),
        ],
        **FAST,
    )
    summary = write_manifest(results, tmp_path / "manifest.json")
    assert summary == {"total": 2, "downloaded": 1, "cached": 0, "not_found": 1, "failed": 0}
    records = json.loads((tmp_path / "manifest.json").read_text())
    assert records[0]["sha256"] and "payload" not in records[0]


def test_alphafold_picks_the_canonical_f1_entry(server, tmp_path, monkeypatch):
    monkeypatch.setattr(async_fetch, "ALPHAFOLD_API", server.base + "/api-{accession}")
    results = async_fetch.fetch_alphafold_models(
        ["p78563", "P78563"], tmp_path, manifest=tmp_path / "m.json", **FAST
    )
    assert set(results) == {"P78563"}  # normalised and de-duplicated
    result = results["P78563"]
    assert result.ok
    assert Path(result.path).name == "AF-P78563-F1-model_v6.pdb"  # isoform entry not taken
    assert server.hits["wrong.pdb"] == 0


def test_rcsb_falls_back_only_when_absent(server, tmp_path, monkeypatch):
    monkeypatch.setattr(async_fetch, "RCSB_FILE", server.base + "/{pdb_id}.{fmt}.gz")
    results = async_fetch.fetch_rcsb_structures(["onlypdb", "nothere"], tmp_path, **FAST)
    # mmCIF absent -> legacy PDB fetched and stored decompressed.
    assert results["ONLYPDB"].ok
    assert Path(results["ONLYPDB"].path).name == "ONLYPDB.pdb"
    assert (tmp_path / "ONLYPDB.pdb").read_bytes() == PDB
    assert server.hits["ONLYPDB.cif.gz"] == 1
    # Absent in every format: reported as not found, nothing written.
    assert results["NOTHERE"].not_found
    assert not list(tmp_path.glob("NOTHERE*"))


def test_alphafold_client_fetches_current_release(server, tmp_path, monkeypatch):
    from cryptic_ip.database.alphafold_client import AlphaFoldClient

    monkeypatch.setattr(async_fetch, "ALPHAFOLD_API", server.base + "/api-{accession}")
    client = AlphaFoldClient(cache_dir=tmp_path)
    (tmp_path / "AF-P78563-F1-model_v4.pdb").write_bytes(PDB)  # a stale release
    path = client.fetch_structure("P78563")
    assert path.name == "AF-P78563-F1-model_v6.pdb"


def test_alphafold_client_offline_uses_newest_cached_by_number(tmp_path, monkeypatch):
    from cryptic_ip.database.alphafold_client import AlphaFoldClient

    # Nothing listens on port 9: the API is unreachable.
    monkeypatch.setattr(async_fetch, "ALPHAFOLD_API", "http://127.0.0.1:9/{accession}")
    monkeypatch.setattr(async_fetch, "fetch_all", _fast(async_fetch.fetch_all))
    client = AlphaFoldClient(cache_dir=tmp_path)
    for version in (4, 10):
        (tmp_path / f"AF-P78563-F1-model_v{version}.pdb").write_bytes(PDB)
    # v10 sorts before v4 as a string; the numeric version must decide.
    assert client.fetch_structure("P78563").name == "AF-P78563-F1-model_v10.pdb"


def _fast(fn):
    def wrapped(jobs, **kwargs):
        kwargs.update(FAST)
        return fn(jobs, **kwargs)

    return wrapped


def test_large_file_downloads_in_ranges_and_matches(tmp_path):
    """Served with Range support: assembled from segments, byte-identical."""
    payload = bytes(range(256)) * 40_000  # ~10 MB
    source = tmp_path / "src.bin"
    source.write_bytes(payload)
    seen_ranges = []

    async def handle(request):
        seen_ranges.append(request.headers.get("Range"))
        return web.FileResponse(source)

    app = web.Application()
    app.router.add_route("*", "/archive.tar", handle)
    runner = web.AppRunner(app)
    loop = asyncio.new_event_loop()
    thread = threading.Thread(target=loop.run_forever, daemon=True)
    thread.start()
    asyncio.run_coroutine_threadsafe(runner.setup(), loop).result()
    site = web.TCPSite(runner, "127.0.0.1", 0)
    asyncio.run_coroutine_threadsafe(site.start(), loop).result()
    port = runner.addresses[0][1]
    try:
        result = async_fetch.download_large_file(
            f"http://127.0.0.1:{port}/archive.tar", tmp_path / "out.tar",
            segments=4, min_segment=1 << 20,
        )
    finally:
        asyncio.run_coroutine_threadsafe(runner.cleanup(), loop).result()
        loop.call_soon_threadsafe(loop.stop)
    assert result.ok, result.error
    assert (tmp_path / "out.tar").read_bytes() == payload
    assert len([r for r in seen_ranges if r]) == 4
    assert not list(tmp_path.glob(".*.part"))


def test_fetch_cli_exit_status(server, tmp_path, monkeypatch, capsys):
    """CI relies on this: failures fail the step; absent entries only with --require-all."""
    from scripts import fetch_structures

    monkeypatch.setattr(async_fetch, "RCSB_FILE", server.base + "/{pdb_id}.{fmt}.gz")
    monkeypatch.setattr(async_fetch, "fetch_all", _fast(async_fetch.fetch_all))

    assert fetch_structures.main(["rcsb", "--out", str(tmp_path), "onlypdb"]) == 0
    assert (tmp_path / "ONLYPDB.pdb").exists()
    manifest = json.loads((tmp_path / "rcsb_manifest.json").read_text())
    assert manifest[0]["ok"] and manifest[0]["sha256"]

    # Absent: reported, not fatal - unless every identifier is required.
    assert fetch_structures.main(["rcsb", "--out", str(tmp_path), "nothere"]) == 0
    assert fetch_structures.main(["rcsb", "--out", str(tmp_path), "--require-all", "nothere"]) == 1
    assert "absent NOTHERE" in capsys.readouterr().out

    # A server that keeps failing: fatal.
    monkeypatch.setattr(async_fetch, "RCSB_FILE", server.base + "/flaky.pdb?{pdb_id}{fmt}")
    server.flaky_left = 100
    assert fetch_structures.main(["rcsb", "--out", str(tmp_path / "x"), "--prefer", "pdb", "zzzz"]) == 1


def test_fetch_cli_reads_identifier_lists(tmp_path):
    from scripts import fetch_structures

    (tmp_path / "ids.csv").write_text("pdb_id,other\n1abc,x\n1ABC,y\n2def,z\n")
    (tmp_path / "ck.json").write_text(json.dumps({"processed": ["P1", "p2"]}))
    import argparse

    args = argparse.Namespace(
        ids=["3ghi"], ids_file=None, ids_csv=str(tmp_path / "ids.csv"), column="pdb_id",
        ids_json=str(tmp_path / "ck.json"), json_key="processed",
    )
    assert fetch_structures._identifiers(args) == ["3GHI", "1ABC", "2DEF", "P1", "P2"]


def test_pdb_client_validates_and_refetches_a_cut_file(server, tmp_path, monkeypatch):
    from cryptic_ip.database.pdb_client import PDBClient

    monkeypatch.setattr(async_fetch, "fetch_all", _fast(async_fetch.fetch_all))
    monkeypatch.setattr(PDBClient, "BASE_URL", server.base)
    client = PDBClient(cache_dir=tmp_path)
    (tmp_path / "GOOD.pdb").write_bytes(PDB[:30])  # a cut-off earlier download
    assert client.fetch_structure("good") == tmp_path / "GOOD.pdb"
    assert (tmp_path / "GOOD.pdb").read_bytes() == PDB  # re-fetched, not reused
    assert server.hits["GOOD.pdb"] == 1
    client.fetch_structure("GOOD")
    assert server.hits["GOOD.pdb"] == 1  # a valid cached file is reused
    with pytest.raises(ValueError):
        client.fetch_structure("MISSING")



def test_pdb_with_a_long_header_is_valid():
    """7SNQ's header runs past 200 KB before its first ATOM record."""
    from cryptic_ip.database.async_fetch import validate_pdb

    header = b"".join(b"REMARK 999 %-60d\n" % i for i in range(5000))  # ~370 KB
    assert len(header) > 300_000
    validate_pdb(header + PDB)  # must not raise
    with pytest.raises(async_fetch.ValidationError):
        validate_pdb(header + b"END\n")
