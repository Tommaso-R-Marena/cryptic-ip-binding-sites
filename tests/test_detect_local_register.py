"""Tests for the locally phase-locked register statistic.

Two of these carry the whole validity argument for the pruned variant and are
the reason the others exist at all:

* the statistic is invariant to sliding a lining along the chain, because the
  resultant's modulus is, and
* the window selection depends only on the **gap multiset** and not on its
  order — which is what makes selecting window sizes from the null free rather
  than a tuning that would invalidate the permutation p-value.

If either stops holding, the pruned statistic is no longer a valid test and the
numbers in `results/architecture/ARCHITECTURE.md` section 10 mean nothing.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import architecture as arch

HELIX = [4, 8, 11, 15, 18, 22, 25, 29, 32, 36, 39, 43]
CONTIGUOUS = list(range(20, 34))


# ------------------------------------------------------------- windows

def test_windows_enumerates_every_offset_of_every_size():
    assert arch.windows(8, (6,)) == [(6, 0), (6, 1), (6, 2)]


def test_a_window_larger_than_the_lining_is_dropped_not_clipped():
    """Clipping would make the statistic depend on n in a way the null cannot copy."""
    assert arch.windows(7, (6, 8, 10)) == [(6, 0), (6, 1)]


def test_a_lining_shorter_than_every_window_enumerates_nothing():
    assert arch.windows(4, (6, 8, 10)) == []


def test_a_window_of_one_position_is_refused():
    with pytest.raises(arch.ArchitectureError):
        arch.windows(10, (1,))


def test_best_window_reports_where_the_window_is():
    best = arch.best_window(HELIX, arch.HELIX_BAND, sizes=(6,))
    assert best["size"] == 6
    assert 0 <= best["start"] <= len(HELIX) - 6
    assert best["first_position"] == HELIX[best["start"]]
    assert best["last_position"] == HELIX[best["start"] + 5]
    assert 0.0 <= best["resultant"] <= 1.0


def test_the_best_window_is_at_least_as_good_as_any_single_window():
    best = arch.best_window(HELIX, arch.HELIX_BAND, sizes=(6, 8))
    for size, start in arch.windows(len(HELIX), (6, 8)):
        window = HELIX[start:start + size]
        value, _ = arch.coherence(window, arch.HELIX_BAND,
                                  probes=arch.band_probes(arch.HELIX_BAND,
                                                          HELIX[-1] - HELIX[0]))
        assert value <= best["resultant"] + 1e-9


# ------------------------------------- the two invariances that matter

def test_sliding_the_lining_along_the_chain_changes_nothing():
    """The resultant's modulus is translation invariant, so the statistic is."""
    here = arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200)
    there = arch.local_register([p + 137 for p in HELIX], arch.HELIX_BAND, n_draws=200)
    assert here.observed == pytest.approx(there.observed, abs=1e-12)
    assert here.p_one_sided == pytest.approx(there.p_one_sided)
    assert here.window_size == there.window_size


def test_the_window_selection_sees_only_the_gap_multiset():
    """The validity argument: the pilot cannot depend on the observed order.

    Reordering the same gaps gives a different lining with the same multiset. If
    the selection changed, it would be a function of the observation and the
    permutation p-value would no longer be exact.
    """
    gaps = arch.gaps(HELIX)
    shuffled = [HELIX[0]]
    for gap in np.random.default_rng(3).permutation(gaps):
        shuffled.append(shuffled[-1] + int(gap))
    assert sorted(arch.gaps(shuffled)) == sorted(gaps)
    assert shuffled != HELIX
    first, means = arch.pilot_sizes(HELIX, arch.HELIX_BAND)
    second, other = arch.pilot_sizes(shuffled, arch.HELIX_BAND)
    assert first == second
    assert means == other


def test_the_pilot_is_not_drawn_from_the_null_sample():
    """Its seed is its own, so the selection is fixed before the null is drawn."""
    result = arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200,
                                 selection="unsaturated", pilot_seed=99)
    assert result.pilot_seed == 99
    unpruned = arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200)
    assert unpruned.pilot_seed is None
    assert unpruned.pilot_mean_by_size == {}


# ------------------------------------------------------------- pruning

def test_pruning_drops_exactly_the_saturated_sizes():
    kept, means = arch.pilot_sizes(HELIX, arch.HELIX_BAND)
    expected = tuple(s for s, mean in sorted(means.items()) if mean < arch.SATURATED_NULL)
    assert kept == expected
    result = arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200, selection="unsaturated")
    assert result.sizes == kept
    assert set(result.pruned_out) == set(means) - set(kept)


def test_pruning_reports_the_pilot_means_it_decided_on():
    result = arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200, selection="unsaturated")
    assert set(result.pilot_mean_by_size) == set(arch.WINDOW_SIZES)
    assert all(0.0 <= v <= 1.0 for v in result.pilot_mean_by_size.values())


def test_the_saturation_threshold_is_the_one_fixed_in_advance():
    """Amendment 3 fixed it before any pruned statistic existed; 4 inherited it."""
    assert arch.SATURATED_NULL == 0.9
    assert arch.WINDOW_SIZES == (6, 8, 10)


def test_a_lining_whose_every_window_saturates_is_unresolved_not_significant():
    """Monkeypatched threshold, because the point is the branch and not the data."""
    lining = arch.local_architecture(HELIX, n_draws=200, selection="unsaturated")
    assert lining.local_p is not None
    saved = arch.SATURATED_NULL
    try:
        arch.SATURATED_NULL = 0.0          # nothing can be below this
        with pytest.raises(arch.ArchitectureError):
            arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200, selection="unsaturated")
        blocked = arch.local_architecture(HELIX, n_draws=200, selection="unsaturated")
    finally:
        arch.SATURATED_NULL = saved
    assert blocked.local_p is None
    assert blocked.verdict == "unresolved"
    assert "saturates" in (blocked.reason or "")


# --------------------------------------------------- unresolved, not p = 1

def test_a_contiguous_lining_is_unresolved_rather_than_null():
    lining = arch.local_architecture(CONTIGUOUS, n_draws=200)
    assert lining.local_p is None
    assert lining.verdict == "unresolved"
    assert "one ordering" in (lining.reason or "")


def test_too_few_positions_is_a_stated_reason():
    lining = arch.local_architecture([3, 7, 11], n_draws=200)
    assert lining.local_p is None
    assert "positions" in (lining.reason or "")


def test_local_register_refuses_a_lining_below_the_floor():
    with pytest.raises(arch.ArchitectureError):
        arch.local_register([3, 7, 11], arch.HELIX_BAND)


# ----------------------------------------------- two detectors, two families

def test_the_local_family_is_disjoint_from_the_pooled_one():
    """The amendment requires it: folding these in would move published verdicts."""
    pooled = arch.describe_architecture(HELIX, n_draws=200)
    local = arch.local_architecture(HELIX, n_draws=200)
    assert set(pooled.family) & set(local.family) == set()
    assert set(local.family) == {"helix_local", "strand_local"}


def test_holm_over_the_local_family_is_never_smaller_than_the_raw_p():
    local = arch.local_architecture(HELIX, n_draws=200)
    corrected = local.holm
    for name, raw in local.family.items():
        assert corrected[name] >= raw - 1e-12


def test_the_local_p_is_the_best_corrected_p():
    local = arch.local_architecture(HELIX, n_draws=200)
    assert local.local_p == pytest.approx(min(local.holm.values()))


# ----------------------------------------------------------- calibration

def test_the_local_statistic_does_not_cry_wolf_on_random_linings():
    """A valid permutation test, so the rate at 0.05 has to stay near it."""
    hits = 0
    measured = 0
    for seed in range(24):
        rng = np.random.default_rng(500 + seed)
        positions = sorted(int(v) for v in rng.choice(220, size=14, replace=False))
        lining = arch.local_architecture(positions, n_draws=200, seed=seed)
        if lining.local_p is None:
            continue
        measured += 1
        hits += int(lining.local_p <= 0.05)
    assert measured >= 20
    assert hits / measured <= 0.15, f"{hits}/{measured} random linings called at 0.05"


def test_a_result_serialises_to_json_safe_types():
    import json
    result = arch.local_architecture(HELIX, n_draws=200, selection="unsaturated")
    text = json.dumps(result.as_dict())
    assert "pilot_mean_by_size" in text
    assert json.loads(text)["verdict"] in {"clear", "weak", "none", "unresolved"}


# ------------------------------------------------- the rank, amendment 5

def test_the_rank_rule_always_selects_exactly_one_size():
    kept, means = arch.pilot_sizes(HELIX, arch.HELIX_BAND, rule="best")
    assert len(kept) == 1
    assert means[kept[0]] == min(means.values())


def test_the_rank_rule_sees_only_the_gap_multiset_too():
    """The validity argument has to cover the rank, not only the threshold."""
    gaps = arch.gaps(HELIX)
    shuffled = [HELIX[0]]
    for gap in np.random.default_rng(11).permutation(gaps):
        shuffled.append(shuffled[-1] + int(gap))
    assert shuffled != HELIX
    first, _ = arch.pilot_sizes(HELIX, arch.HELIX_BAND, rule="best")
    second, _ = arch.pilot_sizes(shuffled, arch.HELIX_BAND, rule="best")
    assert first == second


def test_the_rank_fires_where_the_threshold_need_not():
    """A rank cannot fail to select, which is the whole point of amendment 5."""
    saved = arch.SATURATED_NULL
    try:
        arch.SATURATED_NULL = 0.0
        with pytest.raises(arch.ArchitectureError):
            arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200,
                                selection="unsaturated")
        ranked = arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200, selection="best")
    finally:
        arch.SATURATED_NULL = saved
    assert len(ranked.sizes) == 1


def test_an_unknown_selection_rule_is_refused():
    with pytest.raises(arch.ArchitectureError):
        arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200, selection="whatever")
    with pytest.raises(arch.ArchitectureError):
        arch.pilot_sizes(HELIX, arch.HELIX_BAND, rule="whatever")


def test_the_three_selections_are_nested_by_construction():
    """all contains unsaturated contains best, so they cannot cross by accident."""
    every, _ = tuple(arch.WINDOW_SIZES), None
    unsaturated, _means = arch.pilot_sizes(HELIX, arch.HELIX_BAND, rule="unsaturated")
    best, _again = arch.pilot_sizes(HELIX, arch.HELIX_BAND, rule="best")
    assert set(unsaturated) <= set(every)
    assert set(best) <= set(every)


def test_the_selection_is_recorded_on_the_result():
    for rule in arch.SELECTIONS:
        result = arch.local_register(HELIX, arch.HELIX_BAND, n_draws=200, selection=rule)
        assert result.selection == rule
        assert result.as_dict()["selection"] == rule
