"""Tests for the pocket-to-window layer and its two refusals.

Two things are checked here that no amount of testing `sequence.py` would
catch: that a positional statistic cannot be computed on a gap-closed window,
and that the pre-registered test family in the code still matches the one in
`docs/SEQUENCE_PLAN.md`. The second is a parity test against the document
itself, because a pre-registration that the code has quietly drifted away from
is worse than none.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from cryptic_ip.detect import signature as sg

PLAN = Path(__file__).resolve().parents[1] / "docs" / "SEQUENCE_PLAN.md"

#: A chain with basic residues on a roughly helical spacing, so the windows
#: built from it are measurable and have something in them.
CHAIN = "MKVARLSTGGRKLAAVRRTGGKLSAARVRKLGGSARVLKTRGASRVLKTGGASRVL"

#: Ten lining positions three to four apart: a plausible helix-face lining.
LINING = [5, 9, 12, 16, 19, 23, 26, 30, 33, 37]


def windows():
    return sg.site_windows(CHAIN, LINING, offset=1)


# ----------------------------------------------------------- window building

def test_the_lining_window_is_the_lining_residues_in_order():
    built = windows()
    assert built.lining == "".join(CHAIN[p - 1] for p in LINING)
    assert len(built.lining) == len(LINING)


def test_the_local_window_is_the_contiguous_stretch_they_span():
    built = windows()
    assert built.local == CHAIN[LINING[0] - 1:LINING[-1]]
    assert built.local_start == LINING[0]
    assert built.local_end == LINING[-1]
    assert not built.truncated


def test_the_lining_window_is_shorter_than_the_local_one_it_came_from():
    built = windows()
    assert len(built.lining) < len(built.local)


def test_positions_are_deduplicated_and_sorted():
    built = sg.site_windows(CHAIN, [12, 5, 9, 5, 16, 19, 23, 26, 30, 33], offset=1)
    assert built.positions == (5, 9, 12, 16, 19, 23, 26, 30, 33)


def test_an_offset_shifts_the_numbering_and_not_the_window():
    """A chain that does not start at residue 1 must give the same window."""
    plain = windows()
    shifted = sg.site_windows(CHAIN, [p + 100 for p in LINING], offset=101)
    assert shifted.lining == plain.lining
    assert shifted.local == plain.local


def test_a_position_outside_the_chain_raises_rather_than_being_clamped():
    """Getting the offset wrong would otherwise be invisible in the output."""
    with pytest.raises(sg.WindowError) as caught:
        sg.site_windows(CHAIN, [500, 501, 502], offset=1)
    assert "offset" in str(caught.value)


def test_a_pocket_needs_two_positions():
    with pytest.raises(sg.WindowError):
        sg.site_windows(CHAIN, [5], offset=1)


def test_an_empty_chain_raises():
    with pytest.raises(sg.WindowError):
        sg.site_windows("", LINING, offset=1)


def test_a_wide_lining_truncates_its_local_window_and_says_so():
    built = sg.site_windows(CHAIN, [2, 10, 20, 30, 40, 50], offset=1)
    assert built.truncated
    assert len(built.local) <= sg.LOCAL_CAP
    assert built.as_dict()["truncated"] is True


def test_a_truncated_window_stays_inside_the_chain():
    built = sg.site_windows(CHAIN, [2, 10, 20, 30, 40, 50], offset=1)
    assert built.local_start >= 1
    assert built.local_end <= len(CHAIN)
    assert built.local in CHAIN


def test_nonstandard_residues_are_counted_not_bucketed():
    chain = CHAIN[:10] + "XX" + CHAIN[12:]
    built = sg.site_windows(chain, [11, 12, 16, 19, 23, 26, 30, 33, 37, 40], offset=1)
    assert built.n_nonstandard == 2
    assert "X" not in built.lining
    assert built.nonstandard_fraction == pytest.approx(0.2)


# --------------------------------------------------------------- the floors

def test_a_measurable_pocket_says_so():
    assert sg.is_measurable(windows()) == (True, "measurable")


def test_a_short_lining_is_excluded_with_its_reason():
    built = sg.site_windows(CHAIN, [5, 9, 12, 16], offset=1)
    ok, reason = sg.is_measurable(built)
    assert not ok
    assert "lining window" in reason


def test_a_lining_with_too_many_nonstandard_residues_is_excluded():
    chain = "XXX" + CHAIN
    positions = [1, 2, 3] + [p + 3 for p in LINING] + [45, 48, 52, 55, 58]
    built = sg.site_windows(chain, positions, offset=1)
    assert built.nonstandard_fraction > sg.MAX_NONSTANDARD_FRACTION
    assert len(built.lining) >= sg.MIN_WINDOW
    ok, reason = sg.is_measurable(built)
    assert not ok
    assert "non-standard" in reason


def test_the_nonstandard_fraction_is_reported_ahead_of_the_length_it_caused():
    """A lining shortened by dropped residues belongs in the right column."""
    chain = "XXXXXX" + CHAIN
    built = sg.site_windows(chain, [1, 2, 3, 4, 5, 6, 10, 14, 18], offset=1)
    assert len(built.lining) < sg.MIN_WINDOW
    ok, reason = sg.is_measurable(built)
    assert not ok
    assert "non-standard" in reason


def test_the_exclusion_reason_is_carried_into_the_bundle():
    built = sg.site_windows(CHAIN, [5, 9, 12, 16], offset=1)
    out = sg.describe_site(built, n_shuffles=20)
    assert out["measurable"] is False
    assert out["exclusion_reason"]
    assert "local" not in out


# ------------------------------------------------------------- the fence

def test_a_positional_statistic_is_refused_on_the_gap_closed_lining_window():
    """The rule the plan writes down, enforced where it can be enforced.

    Every function in `sequence.py` accepts either window because both are
    strings, so the refusal has to live here.
    """
    built = windows()
    for statistic in sg.POSITIONAL_STATISTICS:
        with pytest.raises(sg.WindowError) as caught:
            sg.statistic_for_window(built, "lining", statistic, n_shuffles=10)
        assert "positional" in str(caught.value)


def test_the_same_statistic_is_allowed_on_the_local_window():
    built = windows()
    result = sg.statistic_for_window(
        built, "local", "helix_charge_periodicity", n_shuffles=50)
    assert result.statistic == "helix_charge_periodicity"


def test_a_nonpositional_statistic_is_allowed_on_either_window():
    built = windows()
    for window in ("lining", "local"):
        assert sg.statistic_for_window(
            built, window, "lz_normalised", n_shuffles=50).n_shuffles == 50


def test_an_unknown_window_is_refused():
    with pytest.raises(sg.WindowError):
        sg.statistic_for_window(windows(), "global", "entropy", n_shuffles=10)


def test_the_lining_bundle_carries_no_periodicity_at_all():
    """Not a zero, which would read as 'no helical face': nothing."""
    out = sg.describe_site(windows(), n_shuffles=50)
    lining = out["lining"]
    for key in lining["raw"]:
        assert not key.startswith(("helix", "strand", "dominant")), key
    for key in lining["nulls"]:
        assert key not in sg.POSITIONAL_STATISTICS, key


def test_the_local_bundle_does_carry_periodicity():
    out = sg.describe_site(windows(), n_shuffles=50)
    assert "helix_face" in out["local"]["raw"]
    assert "dominant_charge_period" in out["local"]["raw"]


def test_the_spacing_bundle_sees_the_helical_gaps():
    out = sg.describe_site(windows(), n_shuffles=20)
    assert out["spacing"]["fraction_helical"] > 0.5
    assert out["spacing"]["n_positions"] == float(len(LINING))


# ------------------------------------------- parity with the pre-registration

def plan_rows():
    """The plan's primary-test table, as (number, window, statistic)."""
    rows = []
    for line in PLAN.read_text().splitlines():
        match = re.match(r"^\|\s*(\d+)\s*\|(.+)$", line)
        if not match:
            continue
        cells = [c.strip() for c in match.group(2).split("|")]
        ticked = [re.findall(r"`([^`]+)`", cell) for cell in cells]
        rows.append((int(match.group(1)), ticked[0][0], ticked[1][0]))
    return rows


def test_the_plan_still_lists_six_primary_tests():
    assert len(plan_rows()) == 6
    assert len(sg.PRIMARY_TESTS) == 6


def test_the_code_and_the_plan_agree_on_every_primary_test():
    """A pre-registration the code has drifted from is worse than none."""
    from_plan = plan_rows()
    from_code = [(t.number, t.window, t.statistic) for t in sg.PRIMARY_TESTS]
    assert from_code == from_plan


def test_the_primary_tests_are_numbered_one_to_six_in_order():
    assert [t.number for t in sg.PRIMARY_TESTS] == [1, 2, 3, 4, 5, 6]


def test_no_primary_test_puts_a_positional_statistic_on_the_lining_window():
    for test in sg.PRIMARY_TESTS:
        if test.window == "lining":
            assert test.statistic not in sg.POSITIONAL_STATISTICS, test


def test_every_primary_window_is_one_the_plan_defines():
    for test in sg.PRIMARY_TESTS:
        assert test.window in sg.WINDOWS


def test_primary_statistics_lists_the_family_in_plan_order():
    assert sg.primary_statistics() == [(t.window, t.statistic) for t in sg.PRIMARY_TESTS]


def test_the_plan_exists_and_is_dated_before_the_study_runs():
    text = PLAN.read_text()
    assert "before any annotated sequence has been read" in text
    assert "2026-10-12" in text


# --------------------------------------------------------------- the refusal

def test_nothing_is_established_yet():
    """The gate's current state, pinned so adding to it is a deliberate act."""
    assert sg.ESTABLISHED_SIGNATURES == ()


def test_the_detector_gets_no_sequence_term_while_nothing_is_established():
    out = sg.describe_site(windows(), n_shuffles=50)
    value, reason = sg.sequence_contribution(out)
    assert value is None
    assert "study Q6" in reason


def test_an_abstention_is_not_a_zero():
    """`monotone.py` treats None as an abstention; a 0.0 would push a score."""
    value, _reason = sg.sequence_contribution(sg.describe_site(windows(), n_shuffles=20))
    assert value is not 0.0  # noqa: F632 - identity is the point: None, not zero
    assert value is None


def test_an_established_signature_contributes_its_scaled_z():
    out = sg.describe_site(windows(), n_shuffles=100)
    z = out["local"]["nulls"]["entropy_rate_k2"]["z"]
    value, reason = sg.sequence_contribution(
        out, established=("local:entropy_rate_k2",))
    assert value == pytest.approx(
        max(-1.0, min(1.0, z / sg.Z_FOR_FULL_CREDIT)))
    assert "1 established" in reason


def test_a_contribution_is_clipped_to_one():
    site = {"local": {"nulls": {"entropy": {"z": 40.0}}}}
    value, _reason = sg.sequence_contribution(site, established=("local:entropy",))
    assert value == 1.0
    site = {"local": {"nulls": {"entropy": {"z": -40.0}}}}
    value, _reason = sg.sequence_contribution(site, established=("local:entropy",))
    assert value == -1.0


def test_full_credit_needs_twice_the_composition_threshold():
    """The scale is chosen from the threshold it is judged by, not fitted."""
    assert sg.Z_FOR_FULL_CREDIT == 4.0
    site = {"local": {"nulls": {"entropy": {"z": 2.0}}}}
    value, _reason = sg.sequence_contribution(site, established=("local:entropy",))
    assert value == pytest.approx(0.5)


def test_an_established_name_must_say_which_window():
    """The same statistic on the other window is a test nobody registered."""
    out = sg.describe_site(windows(), n_shuffles=20)
    with pytest.raises(sg.WindowError):
        sg.sequence_contribution(out, established=("entropy_rate_k2",))


def test_a_pocket_missing_the_established_feature_abstains():
    value, reason = sg.sequence_contribution(
        {"local": {"nulls": {}}}, established=("local:entropy_rate_k2",))
    assert value is None
    assert "measurable" in reason


def test_every_statistic_is_classified_as_positional_or_not():
    """A statistic added to `sequence.STATISTICS` must be classified.

    Otherwise a new positional measure would default to "safe on the gap-closed
    lining window", which is the one thing this fence exists to prevent.
    """
    from cryptic_ip.detect import sequence as sq

    classified = sg.POSITIONAL_STATISTICS | sg.NON_POSITIONAL_STATISTICS
    assert set(sq.STATISTICS) == classified
    assert not (sg.POSITIONAL_STATISTICS & sg.NON_POSITIONAL_STATISTICS)


def test_the_two_classes_agree_with_what_the_statistics_actually_use():
    """A positional statistic must change when the window is reordered.

    The classification is checked against behaviour rather than trusted: a
    statistic that survives a reordering unchanged is not positional.
    """
    import numpy as np

    from cryptic_ip.detect import sequence as sq

    window = "RAASRVGRGSARVARVTGRAGARSTRTGTRVARGVG"
    letters = list(window)
    np.random.default_rng(5).shuffle(letters)
    reordered = "".join(letters)
    for name in sg.POSITIONAL_STATISTICS:
        assert sq.STATISTICS[name](window) != pytest.approx(
            sq.STATISTICS[name](reordered)), name
    # Entropy is the one that provably cannot see an ordering at all.
    assert sq.STATISTICS["entropy"](window) == pytest.approx(
        sq.STATISTICS["entropy"](reordered))


def test_periodicity_power_is_invariant_under_reversal():
    """Reversal is not a test of positionality: it is a symmetry of the fit.

    A least-squares sinusoid fitted to a reversed series has the same amplitude
    and the opposite phase, so the explained variance is unchanged. The first
    version of the test above used a reversal and passed every statistic as
    non-positional for that reason.
    """
    from cryptic_ip.detect import sequence as sq

    window = "RAASRVGRGSARVARVTGRAGARSTRTGTRVARGVG"
    for name in sg.POSITIONAL_STATISTICS:
        assert sq.STATISTICS[name](window) == pytest.approx(
            sq.STATISTICS[name](window[::-1])), name


def test_the_hydropathy_statistics_are_classified_as_positional():
    """They are periodicity measures, so the gap-closed window must refuse them.

    This is the fence working as designed: three statistics were added to
    `sequence.STATISTICS` for study Q9 and CI failed until they were
    classified, which is exactly the failure the fence exists to force.
    """
    for name in ("helix_hydropathy_periodicity", "strand_hydropathy_periodicity",
                 "dominant_hydropathy_period_power"):
        assert name in sg.POSITIONAL_STATISTICS
        assert name not in sg.NON_POSITIONAL_STATISTICS


def test_a_hydropathy_statistic_is_refused_on_the_gap_closed_window():
    windows = sg.SiteWindows(
        lining="LKKLLKKLLKKLLKKL", local="LKKLLKKLLKKLLKKLLKKL",
        positions=(10, 11, 12, 13, 14, 15, 16, 17),
        local_start=10, local_end=29, n_nonstandard=0, truncated=False)
    with pytest.raises(sg.WindowError):
        sg.statistic_for_window(windows, "lining", "helix_hydropathy_periodicity")
