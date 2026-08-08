"""Burial metric tests against the real ADAR2 structure.

These run only where the crystal structure is available. The synthetic
counterparts in ``test_synthetic_benchmark.py`` cover the same code with ground
truth known by construction and run everywhere.
"""

import numpy as np
import pytest

from cryptic_ip.validation.burial_metrics import (
    CRYPTIC_RELATIVE_SASA_MAX,
    SEMI_CRYPTIC_RELATIVE_SASA_MAX,
    classify_burial,
    classify_burial_relative,
    compute_burial_metrics,
)
from cryptic_ip.validation.control_scoring import positive_passed, validation_score

# Sampling density for the real structure; high enough that the quantities
# asserted below are stable well inside their margins.
SASA_POINTS = 256


@pytest.fixture(scope="module")
def adar2_pdb(adar2_structure_path):
    return adar2_structure_path


@pytest.fixture(scope="module")
def adar2_metrics(adar2_pdb):
    return compute_burial_metrics(adar2_pdb, n_points=SASA_POINTS)


def test_classify_burial_absolute_thresholds():
    """The legacy absolute-SASA interface keeps its documented behaviour."""
    assert classify_burial(0.0) == "cryptic"
    assert classify_burial(30.0) == "semi_cryptic"
    assert classify_burial(80.0) == "surface"
    assert classify_burial(None) == "unknown"


def test_classify_burial_relative_uses_the_more_exposed_measure():
    """A ligand buried to the ring but with exposed phosphates is not cryptic."""
    assert classify_burial_relative(0.01, relative_phosphate_sasa=0.01) == "cryptic"
    assert classify_burial_relative(0.01, relative_phosphate_sasa=0.40) == "surface"
    assert classify_burial_relative(None) == "unknown"
    assert classify_burial_relative(0.01, is_probable_artifact=True) == "crystal_artifact"


def test_adar2_ligand_is_substantially_buried(adar2_metrics):
    """ADAR2's InsP6 must measure as far more buried than solvent-exposed.

    Macbeth et al. (Science 309:1534, 2005) describe the ligand as encapsulated
    with only a narrow window to the exterior, so most of its surface must be
    occluded. This asserts that claim directly rather than asserting a specific
    class label, because the class boundary is a calibrated parameter and this
    fact is not.
    """
    metrics = adar2_metrics
    assert metrics.relative_sasa is not None
    print(
        f"\nADAR2 burial: relative_sasa={metrics.relative_sasa:.4f} "
        f"relative_phosphate_sasa={metrics.relative_phosphate_sasa:.4f} "
        f"depth={metrics.burial_depth:.2f} enclosure={metrics.enclosure:.3f} "
        f"basic_residues={metrics.n_basic_residues} "
        f"class={metrics.burial_class} copies={metrics.n_instances}"
    )
    # Most of the ligand surface is occluded by protein.
    assert metrics.relative_sasa < 0.25
    assert metrics.enclosure > 0.80
    # A buried structural cofactor is coordinated by a basic cluster.
    assert metrics.n_basic_residues >= 4


def test_adar2_is_classified_as_buried_not_surface(adar2_metrics):
    """The class must land on the buried side of the scale.

    ``cryptic`` versus ``semi_cryptic`` depends on where the relative-SASA
    boundary sits, and that boundary is calibrated from the control panel with
    ``scripts/calibrate_controls.py``. What must hold regardless of where it is
    set is that the paradigm buried site is not called ``surface``.
    """
    assert adar2_metrics.burial_class in {"cryptic", "semi_cryptic"}, (
        f"ADAR2 measured relative SASA {adar2_metrics.relative_sasa:.4f} "
        f"(cryptic <= {CRYPTIC_RELATIVE_SASA_MAX}, "
        f"semi-cryptic <= {SEMI_CRYPTIC_RELATIVE_SASA_MAX})"
    )


def test_adar2_reports_per_copy_measurements(adar2_metrics):
    """Every ligand copy is measured separately, not summed into one number."""
    metrics = adar2_metrics
    assert metrics.n_instances >= 1
    assert len(metrics.instances) == metrics.n_instances
    for instance in metrics.instances:
        assert 0.0 <= instance.relative_sasa <= 1.0 or np.isnan(instance.relative_sasa)
        assert instance.sasa_isolated > 0
    # Copies are ordered most-buried first.
    ratios = [
        inst.relative_sasa for inst in metrics.instances if np.isfinite(inst.relative_sasa)
    ]
    assert ratios == sorted(ratios)


def test_adar2_delta_sasa_is_the_surface_buried_on_binding(adar2_metrics):
    """delta_sasa must be surface *lost*, not a restatement of the complex value."""
    metrics = adar2_metrics
    best = metrics.instances[0]
    assert metrics.delta_sasa == pytest.approx(best.sasa_isolated - best.sasa_complex)
    assert metrics.delta_sasa > 0


def test_semi_cryptic_positive_criteria():
    val = validation_score("positive", 0.75, 36.0, 8)
    assert positive_passed(val, 36.0, 8, 0.75, burial_class="semi_cryptic")
