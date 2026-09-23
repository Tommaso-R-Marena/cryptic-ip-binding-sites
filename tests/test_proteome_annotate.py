"""The proteome screen's UniProt annotation download, against a local server."""

from __future__ import annotations

import asyncio
import gzip
import threading

import pytest

aiohttp = pytest.importorskip("aiohttp")
from aiohttp import web  # noqa: E402


def test_uniprot_annotation_download(tmp_path, monkeypatch):
    """Gzipped TSV from the stream endpoint is decompressed, checked and renamed."""
    from scripts import proteome_screen

    tsv = (
        "Entry\tGene Names (primary)\tProtein names\tKeywords\tSubcellular location [CC]\tBinding site\tFunction [CC]\n"
        "P78563\tADARB1\tDouble-stranded RNA-specific editase 1\tRNA-binding\tNucleus\t"
        'BINDING 376; /ligand="1D-myo-inositol hexakisphosphate"\t\n'
    ).encode()
    seen = {}

    async def handle(request):
        seen.update(request.query)
        return web.Response(body=gzip.compress(tsv))

    app = web.Application()
    app.router.add_get("/stream", handle)
    runner = web.AppRunner(app)
    loop = asyncio.new_event_loop()
    threading.Thread(target=loop.run_forever, daemon=True).start()
    asyncio.run_coroutine_threadsafe(runner.setup(), loop).result()
    site = web.TCPSite(runner, "127.0.0.1", 0)
    asyncio.run_coroutine_threadsafe(site.start(), loop).result()
    monkeypatch.setattr(proteome_screen, "UNIPROT_STREAM", f"http://127.0.0.1:{runner.addresses[0][1]}/stream")
    try:
        assert proteome_screen.main(["annotate", "--organism", "human", "--output-dir", str(tmp_path)]) == 0
    finally:
        asyncio.run_coroutine_threadsafe(runner.cleanup(), loop).result()
        loop.call_soon_threadsafe(loop.stop)
    assert seen["query"] == "proteome:UP000005640" and seen["compressed"] == "true"
    import pandas as pd

    frame = pd.read_csv(tmp_path / "human_uniprot.tsv", sep="\t", dtype=str)
    assert list(frame.columns)[:2] == ["uniprot_id", "gene"]
    assert frame.loc[0, "uniprot_id"] == "P78563"
