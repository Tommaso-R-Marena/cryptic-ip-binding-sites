"""Tests for structure validation."""

import pytest
from pathlib import Path

from cryptic_ip.database.alphafold_client import AlphaFoldClient
from cryptic_ip.database.pdb_client import PDBClient

pytestmark = pytest.mark.requires_network


class TestStructureDownload:
    """Test downloading real structures via project clients."""

    def test_download_adar2_alphafold(self, tmp_path):
        client = AlphaFoldClient(cache_dir=tmp_path)
        path = client.fetch_structure("P78563")
        assert path.exists()
        assert path.stat().st_size > 10000

    def test_download_adar2_crystal(self, tmp_path):
        client = PDBClient(cache_dir=tmp_path)
        path = client.fetch_structure("1ZY7")
        assert path.exists()
        assert path.stat().st_size > 10000


class TestStructureParsing:
    """Test parsing cached crystal structures with BioPython."""

    @pytest.fixture
    def adar2_structure(self):
        path = Path("data/validation/1ZY7.pdb")
        if not path.exists():
            pytest.skip("ADAR2 crystal structure not cached")
        from Bio.PDB import PDBParser

        return PDBParser(QUIET=True).get_structure("adar2", str(path))

    def test_structure_has_residues(self, adar2_structure):
        residues = list(adar2_structure.get_residues())
        assert len(residues) > 50

    def test_structure_has_atoms(self, adar2_structure):
        atoms = list(adar2_structure.get_atoms())
        assert len(atoms) > 500

    def test_find_basic_residues(self, adar2_structure):
        basic = {"ARG", "LYS", "HIS"}
        count = sum(1 for r in adar2_structure.get_residues() if r.get_resname() in basic)
        assert count > 10
