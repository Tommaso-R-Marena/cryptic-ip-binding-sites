"""Tests for the architecture-homology metric.

Three things are checked that nothing else can check: that the three distances
really are metrics (identity, symmetry, triangle inequality), that they are
invariant to the things architecture should not depend on, and that the
fold-family result survives the controls which would otherwise explain it.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import architecture as ar
from cryptic_ip.detect import geometry_fixtures as gf
from cryptic_ip.detect import homology as hm


def random_linings(count=10, seed=5):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(count):
        n = int(rng.integers(8, 16))
        span = int(rng.integers(3 * n, 8 * n))
        out.append(sorted(rng.choice(np.arange(span), size=n, replace=False).tolist()))
    return out


# ------------------------------------------------------------ metric axioms

@pytest.mark.parametrize("metric", sorted(hm.METRICS))
def test_a_lining_is_zero_distance_from_itself(metric):
    for lining in random_linings(6):
        assert hm.METRICS[metric](lining, lining) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("metric", sorted(hm.METRICS))
def test_the_distance_is_symmetric(metric):
    linings = random_linings(6)
    for a in linings[:3]:
        for b in linings[3:]:
            assert hm.METRICS[metric](a, b) == pytest.approx(hm.METRICS[metric](b, a))


@pytest.mark.parametrize("metric", sorted(hm.METRICS))
def test_the_triangle_inequality_holds(metric):
    """A combination of metrics is a metric, so this must hold for all three."""
    linings = random_linings(9, seed=7)
    function = hm.METRICS[metric]
    for a in linings[:3]:
        for b in linings[3:6]:
            for c in linings[6:]:
                assert function(a, c) <= function(a, b) + function(b, c) + 1e-9


@pytest.mark.parametrize("metric", sorted(hm.METRICS))
def test_the_distance_is_non_negative(metric):
    linings = random_linings(6)
    for a in linings[:3]:
        for b in linings[3:]:
            assert hm.METRICS[metric](a, b) >= 0.0


# -------------------------------------------------------------- invariance

@pytest.mark.parametrize("metric", sorted(hm.METRICS))
def test_sliding_a_lining_along_the_chain_changes_nothing(metric):
    """Architecture is about spacing, so a constant offset must not move it.

    For the register half this is exact rather than incidental: shifting every
    position by c multiplies the resultant by a unit phase, which the modulus
    discards.
    """
    lining = random_linings(1)[0]
    for shift in (1, 17, 500):
        moved = [p + shift for p in lining]
        assert hm.METRICS[metric](lining, moved) == pytest.approx(0.0, abs=1e-9)


def test_the_spectrum_is_translation_invariant_pointwise():
    lining = random_linings(1)[0]
    first = hm.coherence_spectrum(lining)
    second = hm.coherence_spectrum([p + 97 for p in lining])
    assert np.allclose(first, second, atol=1e-12)


def test_duplicate_positions_do_not_change_a_lining():
    lining = [3, 7, 10, 14, 17, 21]
    assert hm.topology_distance(lining, lining + [7, 14]) == pytest.approx(0.0)


# ----------------------------------------------------------- the topology

def test_the_diagram_is_one_point_per_gap_born_at_zero():
    lining = [1, 4, 8, 13]
    diagram = hm.persistence_diagram(lining)
    assert diagram.shape == (3, 2)
    assert np.all(diagram[:, 0] == 0.0)
    assert sorted(diagram[:, 1].tolist()) == [3.0, 4.0, 5.0]


def test_the_wasserstein_distance_uses_the_diagonal_when_it_is_cheaper():
    """Sorted matching is only optimal at equal cardinality; this is why.

    Measured on random pairs: routing points to the diagonal beat sorted
    matching in nearly half of equal-cardinality cases, by up to ten residues.
    """
    rng = np.random.default_rng(3)
    cheaper = 0
    for _ in range(30):
        n = int(rng.integers(5, 13))
        a = sorted(rng.choice(np.arange(90), size=n, replace=False).tolist())
        b = sorted(rng.choice(np.arange(90), size=n, replace=False).tolist())
        sorted_cost = float(np.abs(np.sort(ar.gaps(a)) - np.sort(ar.gaps(b))).sum())
        exact = hm.topology_distance(a, b)
        assert exact <= sorted_cost + 1e-9
        cheaper += exact < sorted_cost - 1e-9
    assert cheaper > 0


def test_diagrams_of_different_size_are_comparable():
    short = [0, 4, 8]
    long = [0, 4, 8, 12, 16, 20, 24]
    assert hm.topology_distance(short, long) > 0.0
    assert np.isfinite(hm.topology_distance(short, long))


# ----------------------------------------------------------- the register

def test_a_flat_spectrum_standardises_to_zero_without_dividing_by_zero():
    assert np.allclose(hm.standardise(np.ones(50)), 0.0)


def test_the_register_distance_is_bounded_by_two():
    linings = random_linings(8, seed=11)
    for a in linings[:4]:
        for b in linings[4:]:
            assert 0.0 <= hm.register_distance(a, b) <= 2.0 + 1e-9


def test_two_faces_of_the_same_period_are_nearer_than_two_different_periods():
    helical_a = list(np.cumsum([3, 4] * 7))
    helical_b = list(np.cumsum([4, 3] * 7))
    extended = list(range(0, 30, 2))
    assert hm.register_distance(helical_a, helical_b) < hm.register_distance(helical_a, extended)


# ------------------------------------------------------------- the matrix

def test_the_distance_matrix_is_symmetric_with_a_zero_diagonal():
    matrix = hm.distance_matrix(random_linings(6), "register")
    assert np.allclose(matrix, matrix.T)
    assert np.allclose(np.diag(matrix), 0.0)


def test_an_unknown_metric_is_refused():
    with pytest.raises(ar.ArchitectureError):
        hm.distance_matrix(random_linings(4), "vibes")


def test_a_bad_topology_scale_is_refused():
    with pytest.raises(ar.ArchitectureError):
        hm.architecture_distance([0, 3, 7], [0, 4, 8], topology_scale=0.0)


# ------------------------------------------------------------ retrieval

def fold_library():
    """Linings from two folds, matched on size, element length and linker."""
    def atoms(fixture):
        return {int(r): {str(n): c for n, c in
                         zip(fixture.atom_names[fixture.residue_index == r],
                             fixture.coords[fixture.residue_index == r])}
                for r in np.unique(fixture.residue_index)}

    per, linker, size = 18, 10, 16
    renumber = lambda pos: sorted(p + linker * (p // per) for p in pos)  # noqa: E731
    linings, labels = [], []
    for count in (4, 5, 6, 7):
        barrel = atoms(gf.helix_barrel(count, per, 9.0))
        inward = sorted(r for r, d in barrel.items()
                        if "CA" in d and "CB" in d
                        and np.hypot(*d["CB"][:2]) < np.hypot(*d["CA"][:2]))
        sheet = atoms(gf.strand_sheet(count, per, 4.8))
        face = sorted(r for r, d in sheet.items()
                      if "CA" in d and "CB" in d and d["CB"][1] > d["CA"][1])
        for positions, label in ((inward, "barrel"), (face, "sheet")):
            if len(positions) >= size:
                linings.append(renumber(positions[:size]))
                labels.append(label)
    return linings, labels


def test_the_register_metric_retrieves_the_fold():
    linings, labels = fold_library()
    assert len(set(labels)) == 2
    result = hm.retrieval(linings, labels, "register")
    assert result.precision_at_one == 1.0
    assert result.separation > 1.0


def test_the_size_control_cannot_retrieve_a_size_matched_library():
    """If it could, the real metrics would have proved nothing."""
    linings, labels = fold_library()
    assert len({len(set(p)) for p in linings}) == 1, "the library must be size-matched"
    hm.METRICS.setdefault("size_only", lambda a, b: abs(len(set(a)) - len(set(b))))
    result = hm.retrieval(linings, labels, "size_only")
    assert result.separation == pytest.approx(0.0, abs=1e-9)


def test_retrieval_needs_one_label_per_lining():
    linings, labels = fold_library()
    with pytest.raises(ar.ArchitectureError):
        hm.retrieval(linings, labels[:-1], "register")


def test_retrieval_needs_three_linings():
    with pytest.raises(ar.ArchitectureError):
        hm.retrieval([[0, 3, 7], [0, 4, 8]], ["a", "b"], "register")


def test_the_permutation_p_is_small_when_retrieval_is_perfect():
    linings, labels = fold_library()
    assert hm.label_permutation_p(linings, labels, "register", n_draws=500) < 0.05


def test_the_permutation_p_is_large_for_shuffled_labels():
    """The null must be able to return a null."""
    linings, _labels = fold_library()
    rng = np.random.default_rng(2)
    scrambled = list(rng.permutation(["a", "b"] * (len(linings) // 2)))
    scrambled += ["a"] * (len(linings) - len(scrambled))
    assert hm.label_permutation_p(linings, scrambled, "register", n_draws=500) > 0.05


# ------------------------------------------------------------ clustering

def test_single_linkage_splits_the_two_folds():
    linings, labels = fold_library()
    assignment = hm.single_linkage(hm.distance_matrix(linings, "register"), 2)
    assert hm.cluster_purity(assignment, labels) == 1.0


def test_asking_for_one_cluster_returns_one():
    linings, _ = fold_library()
    assert len(set(hm.single_linkage(hm.distance_matrix(linings, "register"), 1))) == 1


def test_asking_for_every_item_as_its_own_cluster_works():
    linings, _ = fold_library()
    assignment = hm.single_linkage(hm.distance_matrix(linings, "register"), len(linings))
    assert len(set(assignment)) == len(linings)


def test_an_impossible_cluster_count_is_refused():
    linings, _ = fold_library()
    with pytest.raises(ar.ArchitectureError):
        hm.single_linkage(hm.distance_matrix(linings, "register"), len(linings) + 1)


def test_purity_needs_one_label_per_item():
    with pytest.raises(ar.ArchitectureError):
        hm.cluster_purity([0, 0, 1], ["a", "b"])


# ------------------------------------------------------- graded and lookup

def test_the_parameter_correlation_returns_a_spearman_and_a_pair_count():
    linings, _ = fold_library()
    out = hm.parameter_correlation(linings, list(range(len(linings))), "register")
    assert -1.0 <= out["spearman"] <= 1.0
    assert out["n_pairs"] == len(linings) * (len(linings) - 1) / 2


def test_the_parameter_correlation_needs_one_value_per_lining():
    linings, _ = fold_library()
    with pytest.raises(ar.ArchitectureError):
        hm.parameter_correlation(linings, [1.0], "register")


def test_nearest_neighbours_come_back_closest_first():
    linings, labels = fold_library()
    library = {f"{label}-{index}": positions
               for index, (positions, label) in enumerate(zip(linings, labels))}
    query = linings[0]
    hits = hm.nearest_neighbours(query, library, "register", k=4)
    assert [d for _, d in hits] == sorted(d for _, d in hits)
    assert hits[0][0].startswith(labels[0])


# ------------------------------------------------------ the sheet fixture

def test_the_sheet_fixture_has_ideal_covalent_geometry():
    report = gf.bond_length_report(gf.strand_sheet(4, 12, 4.8))
    assert report["CA_CB_mean"] == pytest.approx(1.530, abs=1e-3)
    assert report["CA_CB_max_deviation"] < 1e-6


def test_a_sheet_face_is_spaced_two_apart():
    """Which is what makes it a different architecture from a helix barrel."""
    fixture = gf.strand_sheet(4, 14, 4.8)
    face = []
    for residue in np.unique(fixture.residue_index):
        sel = fixture.residue_index == residue
        atoms = {str(n): c for n, c in zip(fixture.atom_names[sel], fixture.coords[sel])}
        if "CA" in atoms and "CB" in atoms and atoms["CB"][1] > atoms["CA"][1]:
            face.append(int(residue))
    gaps = ar.gaps(sorted(face))
    assert float(np.mean(gaps == 2)) > 0.5


def test_a_sheet_needs_two_strands():
    with pytest.raises(gf.FixtureError):
        gf.strand_sheet(1, 12, 4.8)
