"""Tests for coordinate-based inositol phosphate identification.

The point of :mod:`cryptic_ip.analysis.inositol_detection` is that it decides
from atoms rather than from a list of identifiers, so these tests build ligands
geometrically and check that the classification follows the chemistry: a
cyclohexane ring that is not hydroxylated is rejected, a hydroxylated one with no
phosphate is InsP0, and adding phosphates walks the series up. None of the
assertions depend on a residue name, which is the property under test.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.analysis.inositol_detection import (
    detect_inositol_residues,
    detected_comp_ids,
    summarise_detection,
)
from cryptic_ip.analysis.structure_arrays import StructureArrays

#: Cyclohexane C-C bond length; the ring is built as a planar hexagon, which is
#: geometrically sufficient for bond-graph detection.
CC_BOND = 1.54


def _ring_coords(n: int = 6, bond: float = CC_BOND) -> np.ndarray:
    """Return coordinates of a planar n-membered carbon ring."""
    radius = bond / (2.0 * np.sin(np.pi / n))
    angles = np.arange(n) * 2.0 * np.pi / n
    return np.column_stack(
        [radius * np.cos(angles), radius * np.sin(angles), np.zeros(n)]
    )


def _build_arrays(
    coords: np.ndarray,
    elements: list,
    resname: str = "LIG",
    *,
    is_polymer: bool = False,
) -> StructureArrays:
    """Wrap a single hetero residue's atoms in a StructureArrays view."""
    n = len(elements)
    return StructureArrays(
        coords=np.asarray(coords, dtype=float),
        radii=np.full(n, 1.7),
        elements=np.asarray(elements, dtype=object),
        atom_names=np.asarray([f"{e}{i}" for i, e in enumerate(elements)], dtype=object),
        resnames=np.asarray([resname] * n, dtype=object),
        chain_ids=np.asarray(["A"] * n, dtype=object),
        resseqs=np.asarray([1] * n),
        icodes=np.asarray([""] * n, dtype=object),
        model_ids=np.asarray([0] * n),
        occupancies=np.ones(n),
        bfactors=np.zeros(n),
        is_polymer=np.full(n, is_polymer),
        is_solvent=np.zeros(n, dtype=bool),
        residue_index=np.zeros(n, dtype=int),
        residue_keys=[(0, "A", 1, "")],
    )


def _inositol(n_phosphates: int = 0, resname: str = "LIG"):
    """Build an inositol ring with ``n_phosphates`` phosphate groups.

    Each ring carbon carries an oxygen placed radially outward; the first
    ``n_phosphates`` of those oxygens additionally carry a phosphorus.
    """
    ring = _ring_coords()
    coords = list(ring)
    elements = ["C"] * 6

    unit = ring / np.linalg.norm(ring, axis=1, keepdims=True)
    oxygen_coords = ring + unit * 1.43
    for index in range(6):
        coords.append(oxygen_coords[index])
        elements.append("O")

    for index in range(n_phosphates):
        coords.append(ring[index] + unit[index] * (1.43 + 1.60))
        elements.append("P")

    return np.asarray(coords), elements, resname


class TestRingChemistry:
    def test_bare_cyclohexane_is_not_an_inositol(self):
        """A carbon ring with no oxygens must be rejected."""
        ring = _ring_coords()
        arrays = _build_arrays(ring, ["C"] * 6)
        assert detect_inositol_residues(arrays, require_phosphate=False) == []

    def test_hydroxylated_ring_with_no_phosphate_is_insp0(self):
        coords, elements, resname = _inositol(n_phosphates=0)
        arrays = _build_arrays(coords, elements, resname)

        # Not an inositol *phosphate*, so the default call returns nothing.
        assert detect_inositol_residues(arrays) == []

        relaxed = detect_inositol_residues(arrays, require_phosphate=False)
        assert len(relaxed) == 1
        assert relaxed[0].series == "InsP0"
        assert not relaxed[0].is_phosphorylated

    @pytest.mark.parametrize("n_phosphates", [1, 2, 3, 4, 5, 6])
    def test_series_follows_the_phosphate_count(self, n_phosphates):
        coords, elements, resname = _inositol(n_phosphates=n_phosphates)
        arrays = _build_arrays(coords, elements, resname)

        found = detect_inositol_residues(arrays)
        assert len(found) == 1
        assert found[0].n_phosphorus == n_phosphates
        assert found[0].series == f"InsP{n_phosphates}"
        assert found[0].is_phosphorylated

    def test_identity_does_not_depend_on_the_residue_name(self):
        """The whole point: an unknown identifier is still identified."""
        coords, elements, _ = _inositol(n_phosphates=6)
        for resname in ("IHP", "ZZZ", "UNL", "1AB"):
            arrays = _build_arrays(coords, elements, resname)
            found = detect_inositol_residues(arrays)
            assert len(found) == 1, resname
            assert found[0].series == "InsP6", resname
            assert found[0].comp_id == resname


class TestSpeciesSelection:
    def test_phosphorylated_copy_is_ranked_before_free_inositol(self):
        """A structure holding both must not be measured on the free inositol.

        This is the defect the identifier whitelist allowed: it treated INS
        (myo-inositol, no phosphate) as an IP ligand, and burial is measured on
        the most buried matching copy, so a free inositol could displace the
        real ligand.
        """
        free_coords, free_elements, _ = _inositol(n_phosphates=0)
        bound_coords, bound_elements, _ = _inositol(n_phosphates=6)
        bound_coords = bound_coords + np.array([30.0, 0.0, 0.0])

        n_free = len(free_elements)
        n_bound = len(bound_elements)
        coords = np.vstack([free_coords, bound_coords])
        elements = free_elements + bound_elements
        n = n_free + n_bound

        arrays = StructureArrays(
            coords=coords,
            radii=np.full(n, 1.7),
            elements=np.asarray(elements, dtype=object),
            atom_names=np.asarray([f"A{i}" for i in range(n)], dtype=object),
            resnames=np.asarray(["INS"] * n_free + ["IHP"] * n_bound, dtype=object),
            chain_ids=np.asarray(["A"] * n, dtype=object),
            resseqs=np.asarray([1] * n_free + [2] * n_bound),
            icodes=np.asarray([""] * n, dtype=object),
            model_ids=np.asarray([0] * n),
            occupancies=np.ones(n),
            bfactors=np.zeros(n),
            is_polymer=np.zeros(n, dtype=bool),
            is_solvent=np.zeros(n, dtype=bool),
            residue_index=np.asarray([0] * n_free + [1] * n_bound),
            residue_keys=[(0, "A", 1, ""), (0, "A", 2, "")],
        )

        # Requiring phosphate excludes the free inositol outright.
        strict = detect_inositol_residues(arrays)
        assert [r.comp_id for r in strict] == ["IHP"]

        # Even when both are returned, the phosphorylated one ranks first.
        both = detect_inositol_residues(arrays, require_phosphate=False)
        assert [r.series for r in both] == ["InsP6", "InsP0"]
        assert summarise_detection(both) == {"InsP6": 1, "InsP0": 1}
        assert detected_comp_ids(both) == ("IHP", "INS")

    def test_unbonded_phosphate_is_not_credited_to_the_ligand(self):
        """A phosphorus not bonded through a ring oxygen must not be counted."""
        coords, elements, resname = _inositol(n_phosphates=0)
        coords = np.vstack([coords, np.array([8.0, 0.0, 0.0])])
        elements = elements + ["P"]
        arrays = _build_arrays(coords, elements, resname)

        found = detect_inositol_residues(arrays, require_phosphate=False)
        assert len(found) == 1
        assert found[0].n_phosphorus == 0


class TestStructuralContext:
    def test_polymer_residues_are_ignored_by_default(self):
        coords, elements, resname = _inositol(n_phosphates=6)
        arrays = _build_arrays(coords, elements, resname, is_polymer=True)
        assert detect_inositol_residues(arrays) == []
        assert len(detect_inositol_residues(arrays, include_polymer=True)) == 1

    def test_a_five_membered_ring_is_rejected(self):
        """Ribose-like five-membered rings must not pass as inositol."""
        ring = _ring_coords(n=5)
        unit = ring / np.linalg.norm(ring, axis=1, keepdims=True)
        coords = np.vstack([ring, ring + unit * 1.43])
        elements = ["C"] * 5 + ["O"] * 5
        arrays = _build_arrays(coords, elements)
        assert detect_inositol_residues(arrays, require_phosphate=False) == []
