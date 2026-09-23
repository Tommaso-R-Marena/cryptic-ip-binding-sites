"""Tests for proteome archive resolution and extraction."""

from __future__ import annotations

import io
import tarfile

import pytest

from cryptic_ip.database.downloader import extract_models, newest_archive

LISTING = """
<a href="UP000002195_44689_DICDI_v4.tar">UP000002195_44689_DICDI_v4.tar</a>
<a href="UP000002195_44689_DICDI_v6.tar">UP000002195_44689_DICDI_v6.tar</a>
<a href="UP000002311_559292_YEAST_v6.tar">UP000002311_559292_YEAST_v6.tar</a>
<a href="UP000002311_559292_YEAST_v10.tar">UP000002311_559292_YEAST_v10.tar</a>
"""


def test_newest_archive_uses_the_listing_not_the_organism_key():
    # Dictyostelium's mnemonic is DICDI, not DICTYOSTELIUM.
    assert newest_archive(LISTING, "UP000002195") == "UP000002195_44689_DICDI_v6.tar"
    # Versions compare numerically: v10 is newer than v6.
    assert newest_archive(LISTING, "UP000002311") == "UP000002311_559292_YEAST_v10.tar"
    with pytest.raises(RuntimeError):
        newest_archive(LISTING, "UP000005640")


def test_extract_models_keeps_pdb_and_cannot_escape(tmp_path):
    archive = tmp_path / "a.tar"
    with tarfile.open(archive, "w") as tar:
        for name, data in (
            ("AF-P1-F1-model_v6.pdb.gz", b"pdb"),
            ("AF-P1-F1-model_v6.cif.gz", b"cif"),
            ("../../evil-F1-model_v6.pdb.gz", b"x"),
        ):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    out = tmp_path / "out"
    assert extract_models(archive, out) == 2
    assert sorted(p.name for p in out.iterdir()) == ["AF-P1-F1-model_v6.pdb.gz", "evil-F1-model_v6.pdb.gz"]
    assert not (tmp_path.parent / "evil-F1-model_v6.pdb.gz").exists()
