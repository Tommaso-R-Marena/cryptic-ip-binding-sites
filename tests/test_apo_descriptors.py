"""Tests for computing pocket descriptors on apo structures.

The classifier is deployed on AlphaFold models, which never carry a ligand, but
was trained on deposited holo structures. When descriptors were computed with
the ligand present, the ligand covered the residues lining its own pocket, so
those residues read as buried because the ligand was there - a signature of an
*occupied* pocket that an apo target cannot show. These tests pin the fix and
the evidence for it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from cryptic_ip.analysis.features import PocketFeatureExtractor
from cryptic_ip.analysis.inositol_detection import detect_inositol_residues
from cryptic_ip.analysis.structure_arrays import (
    is_polymer_residue,
    load_structure_arrays,
    write_apo_structure,
)
from cryptic_ip.testing.synthetic import default_benchmark_specs, write_synthetic_structure
from scripts.extract_pocket_features import cache_fingerprint

#: Three residues: an ordinary alanine, a selenomethionine recorded as HETATM
#: (as the PDB does), and a phosphate ion standing in for a bound ligand.
MSE_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
HETATM    4  N   MSE A   2       3.300   1.600   0.000  1.00 20.00           N
HETATM    5  CA  MSE A   2       3.900   2.900   0.000  1.00 20.00           C
HETATM    6 SE   MSE A   2       5.600   3.200   1.000  1.00 20.00          SE
HETATM    7  P   PO4 A 101      10.000  10.000  10.000  1.00 20.00           P
HETATM    8  O1  PO4 A 101      11.500  10.000  10.000  1.00 20.00           O
END
"""


class TestPolymerRule:
    @pytest.mark.parametrize(
        "het_flag, resname, expected",
        [
            (" ", "ALA", True),
            ("H_MSE", "MSE", True),  # modified residue recorded as HETATM
            ("H_SEP", "SEP", True),
            ("H_ARG", "ARG", False),  # a free amino acid bound as a ligand
            ("H_IHP", "IHP", False),
            ("W", "HOH", False),
        ],
    )
    def test_polymer_membership(self, het_flag, resname, expected):
        assert is_polymer_residue(het_flag, resname) is expected

    def test_selenomethionine_is_loaded_as_protein(self, tmp_path):
        """MSE used to be classified as ligand, leaving a hole in the chain."""
        path = tmp_path / "mse.pdb"
        path.write_text(MSE_PDB, encoding="utf-8")
        arrays = load_structure_arrays(path)

        mse = arrays.resnames == "MSE"
        assert mse.any()
        assert arrays.is_polymer[mse].all()
        assert not arrays.is_polymer[arrays.resnames == "PO4"].any()


class TestApoStructure:
    def test_ligands_removed_and_modified_residues_kept(self, tmp_path):
        holo = tmp_path / "holo.pdb"
        holo.write_text(MSE_PDB, encoding="utf-8")
        apo = load_structure_arrays(write_apo_structure(holo, tmp_path / "apo.pdb"))

        assert set(apo.resnames) == {"ALA", "MSE"}
        assert apo.is_hetero.sum() == 0

    def test_protein_is_preserved_in_the_same_frame(self, tmp_path):
        """Holo ligand coordinates must still label the apo pockets."""
        spec = default_benchmark_specs(n_buried=1, n_surface=0, n_decoy=0, seed=0)[0]
        holo_path = write_synthetic_structure(spec, tmp_path)
        holo = load_structure_arrays(holo_path)
        apo = load_structure_arrays(write_apo_structure(holo_path, tmp_path / "apo.pdb"))

        assert detect_inositol_residues(apo) == []
        assert apo.n_atoms == int(holo.is_polymer.sum())
        assert np.allclose(
            np.sort(holo.coords[holo.is_polymer], axis=0), np.sort(apo.coords, axis=0), atol=1e-3
        )


class TestTheLeakTheApoModeRemoves:
    """Recorded evidence for computing descriptors without the ligand."""

    @pytest.fixture(scope="class")
    def descriptors(self, tmp_path_factory):
        root = tmp_path_factory.mktemp("leak")
        spec = default_benchmark_specs(n_buried=1, n_surface=0, n_decoy=0, seed=0)[0]
        holo_path = write_synthetic_structure(spec, root)
        holo = load_structure_arrays(holo_path)
        apo = load_structure_arrays(write_apo_structure(holo_path, root / "apo.pdb"))
        site = holo.coords[detect_inositol_residues(holo)[0].atom_indices].mean(axis=0)
        return (
            PocketFeatureExtractor(holo, n_points=128).extract(1, site).features,
            PocketFeatureExtractor(apo, n_points=128).extract(1, site).features,
        )

    def test_the_ligand_makes_its_own_pocket_look_buried(self, descriptors):
        holo, apo = descriptors
        # With the ligand present the lining residues are almost fully covered.
        assert holo["mean_relative_sasa"] < 0.05
        assert apo["mean_relative_sasa"] > 5 * holo["mean_relative_sasa"]
        # Depth inherits the same bias through its definition of "exposed".
        assert holo["burial_depth"] > 2 * apo["burial_depth"]

    def test_enclosure_is_unaffected_because_it_ignores_the_ligand(self, descriptors):
        holo, apo = descriptors
        assert holo["enclosure"] == pytest.approx(apo["enclosure"])


class TestCacheFingerprint:
    BASE = {
        "apo": True, "require_cryptic": False, "sasa_points": 256,
        "min_alpha_spheres": 3, "skip_electrostatics": True,
    }

    def test_identical_settings_share_a_fingerprint(self):
        assert cache_fingerprint(dict(self.BASE), None) == cache_fingerprint(dict(self.BASE), None)

    @pytest.mark.parametrize(
        "change",
        [{"apo": False}, {"require_cryptic": True}, {"sasa_points": 128}],
    )
    def test_settings_that_change_the_rows_change_the_fingerprint(self, change):
        """A cached row computed under other settings must not be reused."""
        assert cache_fingerprint({**self.BASE, **change}, None) != cache_fingerprint(
            dict(self.BASE), None
        )

    def test_explicit_ligand_identifiers_change_the_fingerprint(self):
        assert cache_fingerprint(dict(self.BASE), ["IHP"]) != cache_fingerprint(
            dict(self.BASE), None
        )
