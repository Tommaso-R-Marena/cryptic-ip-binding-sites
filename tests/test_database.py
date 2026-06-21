"""Tests for proteome database manager."""

import pytest
from pathlib import Path

from cryptic_ip.database.manager import ProteomeManager


@pytest.fixture
def proteome_dir(tmp_path):
    """Minimal proteome directory with one AlphaFold-like PDB."""
    pdb = tmp_path / "AF-P78563-F1-model_v6.pdb"
    pdb.write_text(
        "ATOM      1  N   MET A   1      11.104  13.207  14.503  1.00 90.00           N\n"
        "ATOM      2  CA  MET A   1      11.551  11.823  14.240  1.00 90.00           C\n"
        "TER\nEND\n",
        encoding="utf-8",
    )
    return tmp_path


class TestProteomeManager:
    def test_build_catalog(self, proteome_dir):
        manager = ProteomeManager(str(proteome_dir))
        catalog = manager.build_catalog(force=True)
        assert len(catalog) == 1
        assert catalog.iloc[0]["uniprot_id"] == "P78563"

    def test_get_structure_path(self, proteome_dir):
        manager = ProteomeManager(str(proteome_dir))
        manager.build_catalog(force=True)
        path = manager.get_structure_path("P78563")
        assert path is not None
        assert Path(path).exists()
