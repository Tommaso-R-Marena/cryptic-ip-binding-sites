"""Tests for the lining-architecture invariants.

Organised around the module's own claims. Two identities are checked against
independent implementations rather than asserted; two degeneracies that would
have produced confident nonsense are checked to still be there, so that a
future change cannot quietly reintroduce the statistic they rule out; and the
classifier is checked against architectures that are known from a 3-D
construction rather than from a sequence rule.
"""

from __future__ import annotations

import re
from itertools import permutations
from pathlib import Path

import numpy as np
import pytest

from cryptic_ip.detect import architecture as ar
from cryptic_ip.detect import geometry_fixtures as gf

PLAN = Path(__file__).resolve().parents[1] / "docs" / "ARCHITECTURE_PLAN.md"

PER, NHEL, LINK = 18, 6, 10


# --------------------------------------------------------------- fixtures

def barrel_face():
    """Residues of a 6-helix barrel whose CB points at the channel axis.

    Ground truth from the construction: a residue lines the channel exactly
    when its side chain points inward. No sequence rule is used to find them.
    """
    fx = gf.helix_barrel(NHEL, PER, 9.0)
    inward = []
    for residue in np.unique(fx.residue_index):
        sel = fx.residue_index == residue
        atoms = {n: c for n, c in zip(fx.atom_names[sel], fx.coords[sel])}
        if "CA" in atoms and "CB" in atoms:
            if np.hypot(*atoms["CB"][:2]) < np.hypot(*atoms["CA"][:2]):
                inward.append(int(residue))
    return sorted(inward)


def linkered(positions, linker=LINK):
    """Renumber as a real barrel: elements joined by loops, geometry unchanged."""
    return sorted(int(p) + linker * (int(p) // PER) for p in positions)


def single_linkage_deaths(positions):
    """H0 death scales from explicit single-linkage, as an independent check."""
    clusters = [[p] for p in sorted(positions)]
    deaths = []
    while len(clusters) > 1:
        best, index = float("inf"), 0
        for i in range(len(clusters) - 1):
            distance = min(abs(a - b) for a in clusters[i] for b in clusters[i + 1])
            if distance < best:
                best, index = distance, i
        deaths.append(best)
        clusters[index] = clusters[index] + clusters[index + 1]
        del clusters[index + 1]
    return sorted(deaths)


# ------------------------------------------------------------- the identity

def test_the_barcode_is_the_sorted_gap_multiset():
    """The theorem the whole module rests on, against single-linkage."""
    rng = np.random.default_rng(3)
    for _ in range(60):
        n = int(rng.integers(3, 15))
        pos = sorted(rng.choice(np.arange(80), size=n, replace=False).tolist())
        assert single_linkage_deaths(pos) == sorted(ar.h0_barcode(pos).tolist())


def test_components_at_a_scale_come_straight_off_the_barcode():
    pos = [1, 2, 3, 10, 11, 30]
    assert sorted(ar.gaps(pos).tolist()) == [1, 1, 1, 7, 19]
    # At scale 4 the gaps of 7 and 19 survive, so three runs.
    assert ar.segmentation(pos, 4) == 3
    # At scale 7 only the 19 survives.
    assert ar.segmentation(pos, 7) == 2
    assert ar.segmentation(pos, 19) == 1
    assert ar.segmentation(pos, 100) == 1


def test_the_segmentation_curve_is_monotone_and_ends_at_one():
    curve = ar.segmentation_curve([1, 2, 3, 10, 11, 30])
    values = [curve[s] for s in sorted(curve)]
    assert values == sorted(values, reverse=True)
    assert values[-1] == 1


def test_gaps_need_two_distinct_positions():
    with pytest.raises(ar.ArchitectureError):
        ar.gaps([7])
    with pytest.raises(ar.ArchitectureError):
        ar.gaps([7, 7, 7])


def test_fragmentation_is_zero_for_one_run_and_one_for_all_isolated():
    assert ar.fragmentation(list(range(12))) == 0.0
    assert ar.fragmentation([0, 20, 40, 60, 80]) == 1.0


# ------------------------------------------------------------- occupancy

def test_occupancy_of_a_contiguous_run_is_one():
    assert ar.occupancy(list(range(20))) == pytest.approx(1.0)


def test_occupancy_falls_as_a_lining_spreads():
    dense = ar.occupancy(sorted(range(0, 30, 2)))
    sparse = ar.occupancy(sorted(range(0, 300, 20)))
    assert dense > sparse


def test_occupancy_matches_the_closed_form():
    pos = [3, 9, 14, 22, 40]
    assert ar.occupancy(pos) == pytest.approx(np.log(5) / np.log(40 - 3 + 1))


# ------------------------------------------------------------ gap entropy

def test_equal_gaps_give_maximal_regularity():
    assert ar.gap_entropy(list(range(14))) == pytest.approx(1.0)
    assert ar.gap_entropy(list(range(0, 28, 2))) == pytest.approx(1.0)


def test_one_dominating_gap_lowers_regularity():
    clustered = ar.gap_entropy([0, 1, 2, 3, 4, 60, 61, 62, 63, 64])
    assert clustered < 0.5


# ------------------------------------------------------- multiset orderings

def test_distinct_orderings_is_the_multinomial():
    assert ar.distinct_orderings([0, 3, 7, 10, 14, 17, 21, 24, 28]) == 70
    assert ar.distinct_orderings(list(range(10))) == 1


def test_multiset_permutations_enumerates_each_ordering_once():
    for gapset in ([3, 4, 3, 4], [1, 1, 2], [2, 2, 2], [1, 2, 3, 4]):
        produced = list(ar.multiset_permutations(gapset))
        assert len(produced) == len(set(produced))
        assert set(produced) == set(permutations(sorted(gapset)))


def test_multiset_permutations_count_matches_distinct_orderings():
    pos = [0, 3, 7, 10, 14, 17, 21, 24, 28]
    assert len(list(ar.multiset_permutations(ar.gaps(pos)))) == ar.distinct_orderings(pos)


# ------------------------------------------------------- phase statistics

def test_the_resultant_is_bounded_and_maximal_when_phases_agree():
    assert 0.0 <= ar.phase_resultant([0, 4, 8, 12], 4.0) <= 1.0
    assert ar.phase_resultant([0, 4, 8, 12], 4.0) == pytest.approx(1.0)


def test_the_resultant_rejects_an_impossible_period():
    with pytest.raises(ar.ArchitectureError):
        ar.phase_resultant([0, 3, 7], 1.0)


def test_the_vectorised_resultant_matches_the_scalar_one():
    pos = [1, 5, 8, 12, 15, 19]
    grid = np.linspace(2.0, 12.0, 37)
    vector = ar.resultants(np.asarray(pos, dtype=float), grid)
    for index, period in enumerate(grid):
        assert vector[index] == pytest.approx(ar.phase_resultant(pos, float(period)))


def test_the_vectorised_resultant_handles_a_stack_of_linings():
    stack = np.array([[0.0, 3, 7, 10], [0.0, 4, 8, 12]])
    grid = np.array([3.5, 4.0])
    out = ar.resultants(stack, grid)
    assert out.shape == (2, 2)
    assert out[1, 1] == pytest.approx(1.0)


# ----------------------------------------------------- the two degeneracies

def test_the_axial_resultant_is_identically_one_at_period_two():
    """The reason the axial statistic is not used, pinned so it cannot return.

    Standard directional statistics would use the doubled-angle resultant for
    data with two modes 180 degrees apart. On integer positions the doubled
    angle at period 2 is 2*pi*p, zero for every integer, so it is exactly 1 for
    any lining whatsoever.
    """
    rng = np.random.default_rng(4)
    for _ in range(50):
        n = int(rng.integers(5, 25))
        pos = np.sort(rng.choice(np.arange(300), size=n, replace=False)).astype(float)
        angles = 2.0 * np.pi * pos / 2.0
        assert abs(abs(np.exp(2j * angles).mean()) - 1.0) < 1e-12


def test_the_first_resultant_is_degenerate_only_for_uniform_parity():
    """Which is the strand-face property it should detect, so it is kept."""
    assert ar.phase_resultant(list(range(0, 28, 2)), 2.0) == pytest.approx(1.0)
    assert ar.phase_resultant(sorted(set(list(range(0, 28, 2)) + [9])), 2.0) < 1.0


# ------------------------------------------------------ the peak-width law

def test_probe_density_grows_with_the_span():
    narrow = ar.band_probes(ar.HELIX_BAND, 20)
    wide = ar.band_probes(ar.HELIX_BAND, 400)
    assert narrow < wide
    assert narrow >= ar.MIN_BAND_PROBES
    assert wide <= ar.MAX_BAND_PROBES


def test_probe_density_is_bounded_at_both_ends():
    assert ar.band_probes(ar.HELIX_BAND, 1) == ar.MIN_BAND_PROBES
    assert ar.band_probes(ar.HELIX_BAND, 10 ** 6) == ar.MAX_BAND_PROBES


def test_a_fixed_thirteen_probe_scan_loses_a_real_face():
    """The measurement that made probe density a derived quantity."""
    lining_span, period = 200, 3.47
    pos = np.unique(np.round(np.arange(0, lining_span + 1, period)).astype(int)).astype(float)
    fine = float(ar.resultants(pos, np.linspace(3.0, 4.2, 2401)).max())
    coarse = float(ar.resultants(pos, np.linspace(3.0, 4.2, 13)).max())
    assert fine - coarse > 0.2
    # And the derived density recovers nearly all of it.
    derived = float(ar.resultants(pos, ar.band_grid(ar.HELIX_BAND, lining_span)).max())
    assert fine - derived < 0.05


def test_a_band_needs_a_sane_definition():
    with pytest.raises(ar.ArchitectureError):
        ar.band_probes((4.2, 3.0), 50)
    with pytest.raises(ar.ArchitectureError):
        ar.band_grid(ar.HELIX_BAND, 50, probes=1)


# ------------------------------------------------------- dominant period

def test_the_dominant_period_prefers_the_fundamental_over_a_subharmonic():
    """A 3/4-spaced face peaks equally at 3.5, 7.0 and 2.8; the band must not decide."""
    pos = list(np.cumsum([3, 4] * 6))
    dominant = ar.dominant_period(pos)
    assert dominant["in_helix_band"] == 1.0
    assert dominant["period"] == pytest.approx(3.5, abs=0.2)


def test_the_dominant_period_of_an_alternating_lining_is_in_the_strand_band():
    dominant = ar.dominant_period(list(range(0, 28, 2)))
    assert dominant["in_strand_band"] == 1.0
    assert dominant["period"] == pytest.approx(2.0, abs=0.1)


# --------------------------------------------------------- the two nulls

def test_register_is_exact_for_a_small_gap_multiset():
    pos = list(np.cumsum([3, 4] * 5))
    result = ar.register(pos, ar.HELIX_BAND)
    assert result.exact
    assert result.n_draws == result.n_distinct_orderings
    # An exact p is a multiple of 1/k.
    assert result.p_one_sided * result.n_distinct_orderings == pytest.approx(
        round(result.p_one_sided * result.n_distinct_orderings))


def test_register_p_cannot_fall_below_the_orderings_floor():
    pos = list(np.cumsum([3, 4] * 5))
    result = ar.register(pos, ar.HELIX_BAND)
    assert result.p_one_sided >= result.minimum_attainable_p - 1e-12


def test_register_refuses_a_lining_too_short_for_a_null():
    with pytest.raises(ar.ArchitectureError):
        ar.register([1, 4, 8, 11], ar.HELIX_BAND)


def test_the_verdict_comes_from_the_p_value_and_not_the_z():
    """A z off a 15-ordering null returned +2.60 where p could not beat 0.067."""
    pos = [3, 6, 7, 10, 13, 14, 17]
    assert ar.distinct_orderings(pos) < ar.MIN_ORDERINGS
    arch = ar.describe_architecture(pos, n_draws=200)
    assert arch.register_verdict == "unresolved"
    assert arch.register_reason and "attainable p" in arch.register_reason


def test_the_topology_null_finds_a_face_that_register_cannot():
    """The measurement that made topology the primary test.

    A helical face with gaps drawn at random from {3, 4} has no register at
    all: permuting those gaps gives another helical face.
    """
    rng = np.random.default_rng(11)
    pos = list(np.cumsum(rng.choice([3, 4], size=14)))
    topology = ar.topology(pos, ar.HELIX_FACE_GAPS, n_draws=800)
    register = ar.register(pos, ar.HELIX_BAND, n_draws=800)
    assert topology.observed == pytest.approx(1.0)
    assert topology.p_one_sided < 0.01
    assert register.p_one_sided > 0.05


def test_the_topology_null_refuses_a_contiguous_lining():
    with pytest.raises(ar.ArchitectureError):
        ar.topology(list(range(14)), ar.HELIX_FACE_GAPS, n_draws=50)


def test_the_topology_null_needs_three_positions():
    with pytest.raises(ar.ArchitectureError):
        ar.topology([1, 9], ar.HELIX_FACE_GAPS, n_draws=50)


def test_the_topology_null_preserves_n_and_span():
    """Otherwise it would be testing size rather than spacing."""
    pos = [5, 9, 12, 20, 31, 35, 40, 48]
    result = ar.topology(pos, ar.HELIX_FACE_GAPS, n_draws=200)
    assert 0.0 <= result.null_mean <= 1.0
    assert result.observed == ar.face_gap_fraction(pos, ar.HELIX_FACE_GAPS)


def test_face_gap_fraction_counts_the_right_widths():
    assert ar.face_gap_fraction([0, 3, 7, 10], ar.HELIX_FACE_GAPS) == pytest.approx(1.0)
    assert ar.face_gap_fraction([0, 2, 4, 6], ar.STRAND_FACE_GAPS) == pytest.approx(1.0)
    assert ar.face_gap_fraction([0, 2, 4, 6], ar.HELIX_FACE_GAPS) == 0.0


# ----------------------------------------------------------- multiplicity

def test_holm_is_applied_over_the_tests_that_ran():
    arch = ar.describe_architecture(linkered(barrel_face()), n_draws=400)
    raw, corrected = arch.family, arch.holm
    assert set(raw) == set(corrected)
    for name in raw:
        assert corrected[name] >= raw[name] - 1e-12
    assert all(p <= 1.0 for p in corrected.values())


def test_the_verdicts_read_the_corrected_p():
    arch = ar.describe_architecture(linkered(barrel_face()), n_draws=400)
    corrected = arch.holm
    expected = min(corrected.values())
    assert arch.signature_verdict == ("clear" if expected <= ar.CLEAR_P
                                      else "weak" if expected <= ar.WEAK_P else "none")


def test_a_lining_with_no_runnable_test_is_unresolved_not_none():
    """They are different claims and must not be pooled."""
    arch = ar.describe_architecture(list(range(14)), n_draws=100)
    assert arch.architecture_class == "unresolved"
    assert arch.signature_verdict == "unresolved"
    assert arch.topology_reason and arch.register_reason


# ------------------------------------------- recovery of known architectures

def test_the_barrel_face_is_recovered_as_a_helical_mosaic():
    """The headline check: architecture read from positions alone, truth from 3-D."""
    arch = ar.describe_architecture(linkered(barrel_face()), n_draws=800)
    assert arch.signature_verdict == "clear"
    assert arch.element == "helix_face"
    assert arch.architecture_class == "mosaic"
    assert arch.segments == NHEL
    assert arch.dominant["period"] == pytest.approx(3.6, abs=0.3)


def test_the_same_barrel_without_linkers_is_one_run():
    """Its helices are covalently consecutive, so it really is a single run."""
    arch = ar.describe_architecture(barrel_face(), n_draws=800)
    assert arch.segments == 1
    assert arch.architecture_class == "phase_locked"
    assert arch.element == "helix_face"


def test_a_jittered_helical_face_is_still_called_a_helical_face():
    rng = np.random.default_rng(11)
    arch = ar.describe_architecture(list(np.cumsum(rng.choice([3, 4], size=14))), n_draws=800)
    assert arch.element == "helix_face"
    assert arch.architecture_class == "phase_locked"
    assert arch.register_verdict == "none"


def test_a_strand_face_is_called_a_strand_face():
    arch = ar.describe_architecture(sorted(set(list(range(0, 30, 2)) + [9])), n_draws=800)
    assert arch.element == "strand_face"


def test_a_random_lining_is_distributed_with_no_signature():
    rng = np.random.default_rng(5)
    pos = sorted(rng.choice(np.arange(200), 14, replace=False).tolist())
    arch = ar.describe_architecture(pos, n_draws=800)
    assert arch.signature_verdict == "none"
    assert arch.element == "none"
    assert arch.architecture_class == "distributed"


def test_a_hairpin_is_segmental():
    arch = ar.describe_architecture(list(range(10, 17)) + list(range(40, 47)), n_draws=400)
    assert arch.architecture_class == "segmental"
    assert not arch.face_like


@pytest.mark.parametrize("seed", range(6))
def test_random_linings_are_rarely_called_structured(seed):
    """A light calibration check; the full measurement is in the probe script."""
    rng = np.random.default_rng(900 + seed)
    hits = 0
    for index in range(12):
        n = int(rng.integers(9, 20))
        lining_span = int(rng.integers(int(2.5 * n), 6 * n))
        pos = sorted(rng.choice(np.arange(lining_span), size=n, replace=False).tolist())
        arch = ar.describe_architecture(pos, n_draws=300, seed=index)
        hits += arch.signature_verdict in ("weak", "clear")
    assert hits <= 4, f"{hits} of 12 random linings called structured"


# ---------------------------------------------------------- serialisation

def test_the_architecture_serialises_both_halves():
    arch = ar.describe_architecture(linkered(barrel_face()), n_draws=300)
    out = arch.as_dict()
    assert set(out) >= {
        "n_positions", "span", "segments", "fragmentation", "occupancy",
        "gap_entropy", "dominant_period", "helix_topology", "strand_topology",
        "helix_register", "strand_register", "topology_verdict",
        "register_verdict", "signature_verdict", "element",
        "architecture_class", "p_raw", "p_holm",
    }
    assert out["class_meaning"] == ar.CLASSES[out["architecture_class"]]


def test_every_class_has_a_meaning():
    for name in ("phase_locked", "mosaic", "segmental", "distributed", "unresolved"):
        assert name in ar.CLASSES


# ------------------------------------------- parity with the pre-registration

def test_the_plan_lists_exactly_the_code_s_classes():
    """A pre-registration the code has drifted from is worse than none."""
    text = PLAN.read_text()
    named = set(re.findall(r"\*\*([a-z-]+)\*\*\s*\|", text))
    expected = {"phase-locked", "mosaic", "segmental", "distributed"}
    assert expected <= named, named
    for name in ar.CLASSES:
        assert name.replace("_", "-") in text


def test_the_plan_is_dated_before_any_lining_is_read():
    text = PLAN.read_text()
    assert "before any annotated lining has been read" in text
    assert "2026-10-12" in text


def test_the_plan_records_the_constants_the_code_uses():
    text = PLAN.read_text()
    assert str(ar.MIN_ORDERINGS) in text
    assert "MIN_ORDERINGS" in text
    assert f"{ar.PEAK_WIDTH_FACTOR}" in text
