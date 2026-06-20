"""Unit tests for structure parsing helpers."""

from pathlib import Path

import pytest

from cryptic_ip.validation.structure_context import (
    ligand_context,
    ligand_total_sasa,
    load_structure,
    parse_ligand_resnames,
)


@pytest.fixture
def adar2_pdb() -> Path:
    path = Path("data/validation/1ZY7.pdb")
    if not path.exists():
        pytest.skip("ADAR2 crystal structure not available")
    return path


class TestLoadStructure:
    def test_load_pdb(self, adar2_pdb):
        structure = load_structure(adar2_pdb)
        assert len(list(structure.get_atoms())) > 0


class TestLigandContext:
    def test_adar2_has_buried_ligand(self, adar2_pdb):
        residues, ligand_sasa, centroid = ligand_context(adar2_pdb)
        assert len(residues) >= 4
        assert ligand_sasa is not None
        assert ligand_sasa < 10.0
        assert centroid is not None

    def test_parse_ligand_resnames(self, adar2_pdb):
        found = parse_ligand_resnames(adar2_pdb)
        assert "IHP" in found or "IP6" in found


class TestLigandSasa:
    def test_ligand_total_sasa(self, adar2_pdb):
        resnames = parse_ligand_resnames(adar2_pdb)
        ligand_id = next(iter(resnames))
        total = ligand_total_sasa(adar2_pdb, ligand_id)
        assert total >= 0.0
