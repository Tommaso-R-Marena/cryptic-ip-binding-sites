"""Tests for the chain-aware structure array view.

The regressions covered here are the two that silently corrupted earlier
results: collapsing residues from different chains onto one key, and counting
every alternate conformation of an atom.
"""

from pathlib import Path

import numpy as np
import pytest

from cryptic_ip.analysis.structure_arrays import (
    STANDARD_AMINO_ACIDS,
    load_structure_arrays,
    phosphate_group_indices,
)

TWO_CHAIN_PDB = """HEADER    TWO CHAINS SHARING RESIDUE NUMBERS
ATOM      1  N   ARG A  42       0.000   0.000   0.000  1.00 50.00           N
ATOM      2  CA  ARG A  42       1.500   0.000   0.000  1.00 50.00           C
ATOM      3  NH1 ARG A  42       3.000   0.000   0.000  1.00 50.00           N
ATOM      4  N   ARG B  42      20.000   0.000   0.000  1.00 60.00           N
ATOM      5  CA  ARG B  42      21.500   0.000   0.000  1.00 60.00           C
ATOM      6  NH1 ARG B  42      23.000   0.000   0.000  1.00 60.00           N
HETATM    7  O   HOH A 201      10.000  10.000  10.000  1.00 30.00           O
TER
END
"""

ALTLOC_PDB = """HEADER    ALTERNATE CONFORMATIONS
ATOM      1  N   SER A   1       0.000   0.000   0.000  1.00 50.00           N
ATOM      2  CA ASER A   1       1.500   0.000   0.000  0.70 50.00           C
ATOM      3  CA BSER A   1       1.600   0.500   0.000  0.30 50.00           C
TER
END
"""

PHOSPHATE_PDB = """HEADER    PHOSPHATE GEOMETRY
HETATM    1  C1  IHP A 101       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  O1  IHP A 101       1.430   0.000   0.000  1.00 20.00           O
HETATM    3  P1  IHP A 101       3.030   0.000   0.000  1.00 20.00           P
HETATM    4  O1A IHP A 101       4.530   0.000   0.000  1.00 20.00           O
HETATM    5  O1B IHP A 101       3.030   1.500   0.000  1.00 20.00           O
HETATM    6  O9Z IHP A 101      20.000   0.000   0.000  1.00 20.00           O
TER
END
"""


def _write(tmp_path: Path, name: str, content: str) -> Path:
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    return path


def test_residues_in_different_chains_get_distinct_keys(tmp_path: Path):
    """Residue 42 of chain A and chain B must not collapse into one entry."""
    arrays = load_structure_arrays(_write(tmp_path, "two_chain.pdb", TWO_CHAIN_PDB))
    keys = [key for key in arrays.residue_keys]
    assert len(keys) == 2
    assert {key[1] for key in keys} == {"A", "B"}
    assert all(key[2] == 42 for key in keys)


def test_masking_one_residue_selects_only_that_chain(tmp_path: Path):
    arrays = load_structure_arrays(_write(tmp_path, "two_chain.pdb", TWO_CHAIN_PDB))
    mask = arrays.mask_residue((0, "A", 42, ""))
    assert mask.sum() == 3
    assert set(arrays.chain_ids[mask]) == {"A"}


def test_water_is_excluded_by_default(tmp_path: Path):
    arrays = load_structure_arrays(_write(tmp_path, "two_chain.pdb", TWO_CHAIN_PDB))
    assert not np.any(arrays.resnames == "HOH")
    arrays_with_water = load_structure_arrays(
        _write(tmp_path, "two_chain2.pdb", TWO_CHAIN_PDB), keep_solvent=True
    )
    assert np.any(arrays_with_water.resnames == "HOH")
    assert arrays_with_water.is_solvent.sum() == 1


def test_alternate_conformations_are_deduplicated_by_occupancy(tmp_path: Path):
    """Only the highest-occupancy altloc survives, so atoms are not double counted."""
    arrays = load_structure_arrays(_write(tmp_path, "altloc.pdb", ALTLOC_PDB))
    ca_mask = arrays.atom_names == "CA"
    assert ca_mask.sum() == 1
    assert arrays.occupancies[ca_mask][0] == pytest.approx(0.70)


def test_subset_recomputes_residue_keys(tmp_path: Path):
    arrays = load_structure_arrays(_write(tmp_path, "two_chain.pdb", TWO_CHAIN_PDB))
    subset = arrays.subset(arrays.chain_ids == "B")
    assert subset.n_atoms == 3
    assert subset.residue_keys == [(0, "B", 42, "")]
    assert set(subset.residue_index.tolist()) == {0}


def test_phosphate_groups_are_found_by_bonding_not_atom_names(tmp_path: Path):
    """A distant oxygen named like a phosphate oxygen must not be included."""
    arrays = load_structure_arrays(_write(tmp_path, "phos.pdb", PHOSPHATE_PDB))
    mask = arrays.mask_resnames(["IHP"])
    indices = phosphate_group_indices(arrays, mask)
    names = sorted(str(arrays.atom_names[i]) for i in indices)
    # P1 plus the three oxygens within bonding distance; O9Z is 17 A away.
    assert "P1" in names
    assert "O9Z" not in names
    assert "C1" not in names


def test_phosphate_indices_empty_without_phosphorus(tmp_path: Path):
    arrays = load_structure_arrays(_write(tmp_path, "two_chain.pdb", TWO_CHAIN_PDB))
    assert phosphate_group_indices(arrays, np.ones(arrays.n_atoms, dtype=bool)).size == 0


def test_missing_file_raises(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        load_structure_arrays(tmp_path / "absent.pdb")


def test_polymer_flag_tracks_standard_residues(tmp_path: Path):
    arrays = load_structure_arrays(_write(tmp_path, "phos.pdb", PHOSPHATE_PDB))
    assert not arrays.is_polymer.any()
    assert arrays.is_hetero.all()
    assert "ARG" in STANDARD_AMINO_ACIDS
