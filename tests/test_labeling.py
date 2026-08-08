"""Tests for pocket labelling from observed ligand positions.

The central regression: a pocket must be labelled by how much of the ligand it
actually contains, and pockets that partially overlap must be excluded rather
than silently called negative.
"""

import numpy as np
import pytest

from cryptic_ip.analysis.labeling import (
    LigandSite,
    PocketAssignment,
    PocketLabel,
    assign_pocket_labels,
    ligand_overlap_fraction,
    summarise_labels,
)


def _ligand(offset=(0.0, 0.0, 0.0), n: int = 20, spread: float = 3.0) -> LigandSite:
    rng = np.random.default_rng(0)
    coords = rng.uniform(-spread, spread, size=(n, 3)) + np.asarray(offset, dtype=float)
    return LigandSite(instance_id="IHP_A_501", comp_id="IHP", coords=coords)


def test_overlap_is_one_when_the_pocket_covers_the_ligand():
    site = _ligand()
    overlap, min_distance = ligand_overlap_fraction(site.coords, site.coords)
    assert overlap == pytest.approx(1.0)
    assert min_distance == pytest.approx(0.0)


def test_overlap_is_zero_for_a_distant_pocket():
    site = _ligand()
    far = site.coords + np.array([100.0, 0.0, 0.0])
    overlap, min_distance = ligand_overlap_fraction(site.coords, far)
    assert overlap == pytest.approx(0.0)
    assert min_distance > 50.0


def test_overlap_handles_empty_inputs():
    overlap, distance = ligand_overlap_fraction(np.empty((0, 3)), np.zeros((3, 3)))
    assert overlap == 0.0
    assert np.isnan(distance)


def test_covering_pocket_is_positive_and_distant_pocket_is_negative():
    site = _ligand()
    pockets = [
        (1, site.coords.mean(axis=0), site.coords),
        (2, np.array([100.0, 0.0, 0.0]), site.coords + np.array([100.0, 0.0, 0.0])),
    ]
    assignments = assign_pocket_labels(pockets, [site])
    by_id = {a.pocket_id: a for a in assignments}
    assert by_id[1].label is PocketLabel.POSITIVE
    assert by_id[2].label is PocketLabel.NEGATIVE
    assert by_id[1].rank_within_structure == 0


def test_partial_overlap_is_ambiguous_and_excluded_from_training():
    """A pocket holding some of the ligand is neither a site nor a clean decoy."""
    coords = np.column_stack(
        [np.linspace(0.0, 40.0, 20), np.zeros(20), np.zeros(20)]
    )
    site = LigandSite(instance_id="IHP_A_1", comp_id="IHP", coords=coords)
    # Pocket points cover only the first few ligand atoms.
    pocket_points = coords[:3]
    assignments = assign_pocket_labels([(1, pocket_points.mean(axis=0), pocket_points)], [site])
    assert assignments[0].label is PocketLabel.AMBIGUOUS
    assert not assignments[0].label.is_trainable
    assert "partial overlap" in assignments[0].reason


def test_structure_without_ligand_yields_only_negatives():
    pockets = [(1, np.zeros(3), np.zeros((4, 3))), (2, np.ones(3), np.ones((4, 3)))]
    assignments = assign_pocket_labels(pockets, [])
    assert all(a.label is PocketLabel.NEGATIVE for a in assignments)
    assert all("no inositol phosphate" in a.reason for a in assignments)


def test_require_cryptic_excludes_surface_sites_rather_than_negating_them():
    """A surface IP site is not evidence that the pocket cannot bind IP."""
    site = _ligand()
    site.burial_class = "surface"
    assignments = assign_pocket_labels(
        [(1, site.coords.mean(axis=0), site.coords)], [site], require_cryptic=True
    )
    assert assignments[0].label is PocketLabel.AMBIGUOUS

    site.burial_class = "cryptic"
    assignments = assign_pocket_labels(
        [(1, site.coords.mean(axis=0), site.coords)], [site], require_cryptic=True
    )
    assert assignments[0].label is PocketLabel.POSITIVE


def test_inconsistent_thresholds_are_rejected():
    with pytest.raises(ValueError):
        assign_pocket_labels([], [], positive_overlap=0.1, negative_overlap=0.2)


def test_best_matching_ligand_copy_is_reported():
    near = _ligand()
    far = LigandSite(
        instance_id="IHP_B_502", comp_id="IHP", coords=near.coords + np.array([80.0, 0, 0])
    )
    assignments = assign_pocket_labels(
        [(1, near.coords.mean(axis=0), near.coords)], [far, near]
    )
    assert assignments[0].matched_instance_id == "IHP_A_501"


def test_summary_counts_and_detector_recall():
    assignments = {
        "AAAA": [
            PocketAssignment(1, PocketLabel.POSITIVE, 0.9),
            PocketAssignment(2, PocketLabel.NEGATIVE, 0.0),
        ],
        "BBBB": [PocketAssignment(1, PocketLabel.NEGATIVE, 0.0)],
        "CCCC": [PocketAssignment(1, PocketLabel.AMBIGUOUS, 0.2)],
    }
    summary = summarise_labels(assignments, structures_with_ligand=["AAAA", "BBBB", "CCCC"])
    assert summary.n_positive == 1
    assert summary.n_negative == 2
    assert summary.n_ambiguous == 1
    assert summary.n_structures == 3
    # Only AAAA produced a matching pocket, so recall over ligand-bearing
    # structures is one in three.
    assert summary.site_recall == pytest.approx(1 / 3)
    assert summary.positive_rate == pytest.approx(1 / 3)
