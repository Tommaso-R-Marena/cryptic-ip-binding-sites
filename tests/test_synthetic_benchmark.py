"""End-to-end verification against synthetic structures with known ground truth.

These tests are the offline substitute for a live PDB: they assert that the
measurement code recovers answers that follow from how the structures were
built. They validate the machinery, not the biology.
"""

from pathlib import Path

import numpy as np
import pytest

from cryptic_ip.analysis.features import FEATURE_NAMES, PocketFeatureExtractor
from cryptic_ip.analysis.structure_arrays import load_structure_arrays
from cryptic_ip.testing.synthetic import (
    SyntheticStructureSpec,
    build_synthetic_benchmark,
    inositol_hexakisphosphate_coords,
    write_synthetic_structure,
)
from cryptic_ip.validation.burial_metrics import compute_burial_metrics

# Fewer SASA points keeps the suite fast; the contrasts asserted here are far
# larger than the sampling error at this resolution.
SASA_POINTS = 128


@pytest.fixture(scope="module")
def buried_structure(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("synthetic_buried")
    return write_synthetic_structure(
        SyntheticStructureSpec(name="BURIED", site_offset=0.0, seed=1), out
    )


@pytest.fixture(scope="module")
def surface_structure(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("synthetic_surface")
    return write_synthetic_structure(
        SyntheticStructureSpec(
            name="SURFACE",
            site_offset=24.0,
            protein_radius=22.0,
            seed=2,
            expected_burial_class="surface",
        ),
        out,
    )


def test_ligand_geometry_matches_inositol_hexakisphosphate_composition():
    atoms = inositol_hexakisphosphate_coords([0.0, 0.0, 0.0])
    elements = [element for _, element, _ in atoms]
    # C6 O24 P6: the composition of the real chemical component.
    assert len(atoms) == 36
    assert elements.count("C") == 6
    assert elements.count("P") == 6
    assert elements.count("O") == 24


def test_ligand_bond_lengths_are_physical():
    atoms = inositol_hexakisphosphate_coords([0.0, 0.0, 0.0])
    by_name = {name: coord for name, _, coord in atoms}
    # C-O ester bond and O-P bond.
    assert np.linalg.norm(by_name["O1"] - by_name["C1"]) == pytest.approx(1.43, abs=0.05)
    assert np.linalg.norm(by_name["P1"] - by_name["O1"]) == pytest.approx(1.60, abs=0.05)


def test_generator_is_deterministic(tmp_path: Path):
    spec = SyntheticStructureSpec(name="DET", seed=99)
    first = write_synthetic_structure(spec, tmp_path / "a").read_text()
    second = write_synthetic_structure(spec, tmp_path / "b").read_text()
    assert first == second


def test_buried_ligand_is_measured_as_cryptic(buried_structure: Path):
    metrics = compute_burial_metrics(buried_structure, n_points=SASA_POINTS)
    assert metrics.burial_class == "cryptic"
    assert metrics.relative_sasa < 0.05
    assert metrics.enclosure > 0.95


def test_surface_ligand_is_measured_as_surface(surface_structure: Path):
    metrics = compute_burial_metrics(surface_structure, n_points=SASA_POINTS)
    assert metrics.burial_class == "surface"
    assert metrics.relative_sasa > 0.25
    # A site on the exterior sees solvent across a large share of directions.
    assert metrics.enclosure < 0.75


def test_buried_and_surface_sites_are_clearly_separated(
    buried_structure: Path, surface_structure: Path
):
    buried = compute_burial_metrics(buried_structure, n_points=SASA_POINTS)
    surface = compute_burial_metrics(surface_structure, n_points=SASA_POINTS)
    assert surface.relative_sasa - buried.relative_sasa > 0.3
    assert buried.enclosure > surface.enclosure


def test_burial_depth_of_a_core_site_approaches_the_protein_radius(buried_structure: Path):
    """A ligand at the centre of a 22 A sphere must be about 22 A below the surface."""
    metrics = compute_burial_metrics(buried_structure, n_points=SASA_POINTS)
    assert metrics.burial_depth > 15.0


def test_relative_sasa_is_independent_of_ligand_copy_count(tmp_path: Path):
    """The regression that broke the original dataset: copy count changed the class.

    Absolute summed SASA scales with the number of ligand copies. The relative
    measure must not, because it is normalised per copy.
    """
    spec = SyntheticStructureSpec(name="ONE", site_offset=0.0, seed=4)
    path = write_synthetic_structure(spec, tmp_path)

    # Append a second, independent copy far from the first, on its own chain.
    lines = path.read_text().splitlines()
    body = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
    extra = []
    serial = 90000
    for name, element, coord in inositol_hexakisphosphate_coords([200.0, 0.0, 0.0], seed=7):
        atom_name = name if len(name) >= 4 else f" {name:<3s}"
        extra.append(
            f"HETATM{serial:5d} {atom_name:<4s} IHP Z 999    "
            f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00 50.00          {element:>2s}"
        )
        serial += 1
    two_copy = tmp_path / "TWO.pdb"
    two_copy.write_text("\n".join(body + extra + ["END"]) + "\n", encoding="utf-8")

    single = compute_burial_metrics(path, n_points=SASA_POINTS)
    double = compute_burial_metrics(two_copy, n_points=SASA_POINTS)

    assert double.n_instances == 2
    # The buried copy is reported, and its relative burial is unchanged.
    assert double.burial_class == single.burial_class
    assert double.relative_sasa == pytest.approx(single.relative_sasa, abs=0.02)


def test_feature_extraction_produces_the_full_descriptor_suite(buried_structure: Path):
    arrays = load_structure_arrays(buried_structure)
    extractor = PocketFeatureExtractor(arrays, n_points=SASA_POINTS)
    ligand = arrays.coords[arrays.mask_resnames(["IHP"])]
    features = extractor.extract(1, ligand.mean(axis=0), alpha_sphere_coords=ligand)

    assert set(features.features) >= set(FEATURE_NAMES)
    # A buried, basic-lined site must show the expected qualitative signature.
    assert features.features["enclosure"] > 0.9
    assert features.features["n_basic_residues"] > 0
    assert features.features["coulomb_potential_kt"] > 0
    assert np.isfinite(features.features["hull_volume"])


def test_coulomb_surrogate_is_positive_for_basic_and_negative_for_acidic(tmp_path: Path):
    """The electrostatic surrogate must track the sign of the lining charge."""
    basic = write_synthetic_structure(
        SyntheticStructureSpec(name="BASIC", n_basic_lining=40, seed=5), tmp_path / "basic"
    )
    acidic = write_synthetic_structure(
        SyntheticStructureSpec(
            name="ACIDIC", n_basic_lining=40, acidic_lining=True, seed=5
        ),
        tmp_path / "acidic",
    )

    def potential_at_site(path: Path) -> float:
        arrays = load_structure_arrays(path)
        extractor = PocketFeatureExtractor(arrays, n_points=64)
        centre = arrays.coords[arrays.mask_resnames(["IHP"])].mean(axis=0)
        return extractor.coulomb_potential(centre)

    assert potential_at_site(basic) > 0
    assert potential_at_site(acidic) < 0


def test_benchmark_ground_truth_matches_measurement(tmp_path: Path):
    """Every generated structure must measure as the class it was built to be."""
    bench = build_synthetic_benchmark(tmp_path, n_buried=2, n_surface=2, n_decoy=1)
    assert len(bench.paths) == 5

    for path, spec in zip(bench.paths, bench.specs):
        metrics = compute_burial_metrics(path, n_points=SASA_POINTS)
        if not spec.include_ligand:
            assert metrics.burial_class == "unknown"
            assert metrics.n_instances == 0
            continue
        assert metrics.burial_class == spec.expected_burial_class, (
            f"{spec.name}: expected {spec.expected_burial_class}, "
            f"measured {metrics.burial_class} (relative SASA {metrics.relative_sasa:.3f})"
        )
