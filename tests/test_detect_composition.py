"""Tests for the composition half of the decomposition (study Q10).

Two fences carry the study, and both are attacked here: every statistic in the
registry must genuinely ignore the order of a window's letters, and
`composition_z` must refuse every one of them rather than return the silent
zero that an order-permuting null produces for a constant.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import composition as cp
from cryptic_ip.detect import sequence as sq

CHAIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"
SHELL = "RKRKRKHHWY"


# ------------------------------------------------- the registry means what it says

@pytest.mark.parametrize("statistic", sorted(cp.COMPOSITION_STATISTICS))
def test_every_registered_statistic_ignores_the_order_of_its_letters(statistic):
    """The claim that makes the gap-closed lining window sound for this study."""
    background = cp.background_of(CHAIN)
    for window in (CHAIN[:40], SHELL, "ACDEFGHIKLMNPQRSTVWY"):
        assert cp.is_order_invariant(statistic, window, background=background), \
            f"{statistic} reads order on {window!r} and must not be in this registry"


def test_the_two_registries_are_the_same_set():
    """`sequence.py` spells the names out to avoid importing this module."""
    assert sq.ORDER_INVARIANT_STATISTICS == cp.COMPOSITION_STATISTICS


def test_no_name_belongs_to_both_halves_of_the_decomposition():
    assert not (cp.COMPOSITION_STATISTICS & frozenset(sq.STATISTICS))


@pytest.mark.parametrize("statistic", sorted(cp.COMPOSITION_STATISTICS))
def test_composition_z_refuses_an_order_invariant_statistic(statistic):
    """Not a zero - a refusal. A zero would read as 'no ordering signal'."""
    with pytest.raises(sq.SequenceError, match="letters alone"):
        sq.composition_z(CHAIN[:40], statistic, n_shuffles=8)


def test_the_order_statistics_still_pass_the_guard():
    """The guard must not catch what it is not for."""
    result = sq.composition_z("LKKLLKKLLKKLLKKLLKKL", "helix_charge_periodicity",
                              n_shuffles=200, seed=1)
    assert result.z > 2.0


# ----------------------------------------------------------- the statistics

def test_the_basic_fraction_counts_arginine_and_lysine_only():
    assert cp.basic_fraction("RKRK") == pytest.approx(1.0)
    assert cp.basic_fraction("AAAA") == pytest.approx(0.0)
    assert cp.basic_fraction("RKAA") == pytest.approx(0.5)


def test_histidine_is_not_counted_as_basic():
    """It is fractionally protonated at pH 7.4, and `CHARGE_OF` gives it 0.1."""
    assert cp.basic_fraction("HHHH") == pytest.approx(0.0)
    assert sq.CHARGE_OF["H"] == pytest.approx(0.1)


def test_the_divergence_is_zero_against_a_window_s_own_background():
    """A window compared with itself has no compositional difference."""
    window = CHAIN[:40]
    assert cp.composition_divergence(window, cp.background_of(window)) == \
        pytest.approx(0.0, abs=1e-9)


def test_the_divergence_is_bounded_and_positive_for_a_skewed_window():
    value = cp.composition_divergence("RRRRRRRRRR", cp.background_of(CHAIN))
    assert 0.0 < value <= 1.0


def test_the_divergence_is_symmetric():
    one, two = cp.distribution(SHELL), cp.distribution(CHAIN)
    assert cp.jensen_shannon(one, two) == pytest.approx(cp.jensen_shannon(two, one))


def test_the_nonbasic_divergence_cannot_see_arginine_or_lysine():
    """Its whole purpose: what is unusual besides being basic."""
    background = cp.background_of(CHAIN)
    plain = "ACDEFGHILM"
    assert cp.nonbasic_divergence(plain, background) == pytest.approx(
        cp.nonbasic_divergence(plain + "RKRKRK", background), abs=0.03)
    # And the omnibus statistic, which can see them, does move.
    assert cp.composition_divergence(plain, background) != pytest.approx(
        cp.composition_divergence(plain + "RKRKRK", background), abs=0.03)


def test_a_basic_window_separates_from_an_acidic_one_on_net_charge():
    assert cp.net_charge_per_residue("RRRKKK") == pytest.approx(1.0)
    assert cp.net_charge_per_residue("DDDEEE") == pytest.approx(-1.0)


def test_mean_hydropathy_is_on_the_kyte_doolittle_scale():
    assert cp.mean_hydropathy("IIII") == pytest.approx(sq.HYDROPATHY["I"])
    assert cp.mean_hydropathy("RRRR") == pytest.approx(sq.HYDROPATHY["R"])


def test_a_short_window_s_divergence_is_finite_without_any_smoothing():
    """Jensen-Shannon needs no pseudocount, so no free parameter is inside it.

    A two-residue window populates two of twenty letters. A one-sided
    Kullback-Leibler divergence would be infinite here and would need a
    constant chosen to tame it.
    """
    value = cp.composition_divergence("RK", cp.background_of(CHAIN))
    assert np.isfinite(value) and 0.0 < value <= 1.0
    assert not hasattr(cp, "PSEUDOCOUNT")


def test_a_letter_absent_from_the_window_does_not_blow_the_divergence_up():
    window = "AAAAAAAAAA"                      # nineteen letters at probability zero
    assert np.isfinite(cp.composition_divergence(window, cp.background_of(CHAIN)))
    assert np.isfinite(cp.nonbasic_divergence(window, cp.background_of(CHAIN)))


# ------------------------------------------------------------ loud about gaps

def test_an_empty_window_is_refused():
    with pytest.raises(sq.SequenceError):
        cp.letter_counts("")


def test_a_window_of_only_non_standard_letters_is_refused():
    with pytest.raises(sq.SequenceError):
        cp.letter_counts("XXXX")


def test_an_unknown_statistic_names_the_ones_that_exist():
    with pytest.raises(sq.SequenceError, match="known:"):
        cp.evaluate("not_a_statistic", CHAIN[:20], cp.background_of(CHAIN))


def test_an_empty_background_is_refused_rather_than_divided_by():
    with pytest.raises(sq.SequenceError):
        cp.composition_divergence(CHAIN[:20], {})


def test_distributions_of_different_lengths_are_refused():
    with pytest.raises(sq.SequenceError):
        cp.jensen_shannon(np.ones(20) / 20, np.ones(18) / 18)


def test_is_order_invariant_catches_a_statistic_that_reads_order(monkeypatch):
    """The checker has to be able to fail, or the registry fence means nothing."""
    monkeypatch.setitem(cp.SIMPLE_STATISTICS, "_probe", sq.helix_charge_periodicity)
    assert not cp.is_order_invariant("_probe", CHAIN[:40])
    # And it still passes the one beside it, so it is not failing everything.
    assert cp.is_order_invariant("basic_fraction", CHAIN[:40])
