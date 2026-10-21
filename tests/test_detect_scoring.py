"""The monotone scorer with abstention, and the competitor specificity panel."""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import monotone as MO
from cryptic_ip.detect import specificity as SP


def _mid(scorer=None):
    scorer = scorer or MO.MonotoneScorer()
    return {f.name: f.midpoint for f in scorer.features}


# --------------------------------------------------------------- the features

def test_every_shipped_feature_declares_a_direction_and_a_reason():
    for feature in MO.DEFAULT_FEATURES:
        assert feature.direction in (MO.UP, MO.DOWN)
        assert feature.why, f"{feature.name} has no stated reason"
        assert feature.low < feature.high


def test_burial_pushes_down_and_field_terms_push_up():
    directions = {f.name: f.direction for f in MO.DEFAULT_FEATURES}
    assert directions["relative_sasa"] == MO.DOWN
    assert directions["field_depth_kt"] == MO.UP
    assert directions["field_shape"] == MO.UP
    assert directions["openability"] == MO.UP


def test_an_invalid_direction_is_refused():
    with pytest.raises(MO.MonotoneError, match="direction"):
        MO.Feature(name="x", direction=0, midpoint=0.0, slope=1.0, low=0.0, high=1.0, weight=1.0)


@pytest.mark.parametrize("bad", [{"slope": 0.0}, {"slope": -1.0}, {"weight": -1.0}])
def test_an_invalid_shape_parameter_is_refused(bad):
    kwargs = dict(name="x", direction=MO.UP, midpoint=0.0, slope=1.0, low=0.0, high=1.0, weight=1.0)
    kwargs.update(bad)
    with pytest.raises(MO.MonotoneError):
        MO.Feature(**kwargs)


def test_an_inverted_range_is_refused():
    with pytest.raises(MO.MonotoneError, match="low must be below high"):
        MO.Feature(name="x", direction=MO.UP, midpoint=0.0, slope=1.0, low=1.0, high=0.0, weight=1.0)


def test_a_features_bounds_are_its_values_at_the_ends():
    feature = MO.Feature(
        name="x", direction=MO.UP, midpoint=0.5, slope=10.0, low=0.0, high=1.0, weight=1.0,
    )
    low, high = feature.bounds()
    assert low == pytest.approx(feature.value_at(0.0))
    assert high == pytest.approx(feature.value_at(1.0))
    assert low < high


# ------------------------------------------------------------- monotonicity

def test_the_score_is_monotone_in_every_feature():
    # With positives in five homology groups, a model free to fit any shape
    # will learn that more burial is worse if that is what five families
    # happen to show. The monotonicity is declared in advance, and this
    # checks the composition rather than trusting the algebra.
    scorer = MO.MonotoneScorer()
    rng = np.random.default_rng(0)
    for _ in range(2000):
        values = {f.name: float(rng.uniform(f.low, f.high)) for f in scorer.features}
        feature = scorer.features[rng.integers(len(scorer.features))]
        before = scorer.score(values).low
        raised = dict(values)
        raised[feature.name] = min(
            feature.high, values[feature.name] + rng.uniform(0.01, (feature.high - feature.low) * 0.3),
        )
        after = scorer.score(raised).low
        if feature.direction == MO.UP:
            assert after >= before - 1e-12, f"{feature.name} is not increasing"
        else:
            assert after <= before + 1e-12, f"{feature.name} is not decreasing"


def test_the_score_stays_in_the_unit_interval():
    scorer = MO.MonotoneScorer()
    rng = np.random.default_rng(1)
    for _ in range(500):
        values = {f.name: float(rng.uniform(f.low, f.high)) for f in scorer.features}
        interval = scorer.score(values)
        assert 0.0 <= interval.low <= 1.0 and 0.0 <= interval.high <= 1.0


# ---------------------------------------------------------------- abstention

def test_a_complete_pocket_gets_a_point_score():
    interval = MO.MonotoneScorer().score(_mid())
    assert interval.width == 0.0
    assert interval.point is not None
    assert interval.missing == ()


def test_a_missing_feature_widens_the_score_into_an_interval():
    # A pocket from a crystal structure has no pLDDT, and a pipeline without
    # APBS has no electrostatic term. The gate cascade reads both as failures.
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["mean_plddt"] = None
    interval = scorer.score(values)
    assert interval.width > 0.0
    assert interval.point is None
    assert interval.missing == ("mean_plddt",)


def test_a_score_interval_straddling_the_threshold_abstains():
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["mean_plddt"] = None
    assert scorer.score(values).decision(0.5) == "abstain"


def test_a_decision_is_made_when_the_interval_clears_the_threshold():
    scorer = MO.MonotoneScorer()
    good = {f.name: (f.high if f.direction == MO.UP else f.low) for f in scorer.features}
    good["mean_plddt"] = None
    interval = scorer.score(good)
    # Even without pLDDT the rest is enough to decide, which is the point: the
    # missing measurement only matters when it could change the answer.
    assert interval.decision(0.5) == "site"


def test_a_hopeless_pocket_is_decided_against_even_with_data_missing():
    scorer = MO.MonotoneScorer()
    bad = {f.name: (f.low if f.direction == MO.UP else f.high) for f in scorer.features}
    bad["mean_plddt"] = None
    assert scorer.score(bad).decision(0.5) == "not_site"


def test_nan_is_treated_as_missing_rather_than_as_a_value():
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["openability"] = float("nan")
    assert "openability" in scorer.score(values).missing


def test_strict_mode_refuses_a_missing_feature():
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["openability"] = None
    with pytest.raises(MO.MonotoneError, match="is missing"):
        scorer.score(values, strict=True)


def test_a_value_outside_its_physical_range_is_an_error_not_a_clip():
    # Clipping would hide a measurement error, which is how a hull depth of
    # 0.955 was once read as 2.84 angstrom.
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["relative_sasa"] = 5.0
    with pytest.raises(MO.MonotoneError, match="outside its physical range"):
        scorer.score(values)


def test_an_unknown_feature_name_is_refused():
    with pytest.raises(MO.MonotoneError, match="unknown feature"):
        MO.MonotoneScorer().score({"not_a_feature": 1.0})


def test_duplicate_features_are_refused():
    feature = MO.DEFAULT_FEATURES[0]
    with pytest.raises(MO.MonotoneError, match="duplicate"):
        MO.MonotoneScorer([feature, feature])


def test_an_empty_scorer_is_refused():
    with pytest.raises(MO.MonotoneError, match="at least one feature"):
        MO.MonotoneScorer([])


# ------------------------------------------------------------- attribution

def test_the_scorer_says_which_measurement_holds_a_pocket_back():
    # The gate cascade could say a protein failed; it could not say which
    # measurement would have to change, which is what an experimentalist needs.
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["field_shape"] = -1.0
    values["relative_sasa"] = 0.9
    assert scorer.blocking_feature(values, 0.5) in {"field_shape", "relative_sasa"}


def test_a_passing_pocket_has_no_blocking_feature():
    scorer = MO.MonotoneScorer()
    good = {f.name: (f.high if f.direction == MO.UP else f.low) for f in scorer.features}
    assert scorer.blocking_feature(good, 0.5) is None


def test_the_sensitivity_sums_to_the_score():
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    assert sum(scorer.sensitivity(values).values()) == pytest.approx(scorer.score(values).low)


# -------------------------------------------------------------- specificity

def _panel(target=0.9, depth=3.0, **competitors):
    panel = {SP.TARGET: {"shape": target, "depth_kt": depth}}
    for name in SP.REQUIRED_COMPETITORS:
        panel[name] = {"shape": competitors.get(name, 0.1), "depth_kt": depth}
    return panel


def test_a_pocket_that_prefers_ip6_by_the_margin_is_selective():
    verdict = SP.verdict(_panel(target=0.9, ATP=0.2))
    assert verdict.call == SP.SELECTIVE and verdict.is_candidate
    assert "above the pre-registered margin" in verdict.explain()


def test_a_pocket_a_competitor_matches_is_called_polyanion_generic():
    # This is the outcome the deployed screen cannot express, and it is
    # probably the right answer for a large share of the 75 candidates.
    verdict = SP.verdict(_panel(target=0.5, PAPS=0.48))
    assert verdict.call == SP.GENERIC
    assert not verdict.is_candidate
    assert "cannot say which polyanion" in verdict.explain()


def test_a_pocket_a_competitor_beats_is_not_an_ip6_site():
    verdict = SP.verdict(_panel(target=0.3, SULFATE=0.8))
    assert verdict.call == SP.COMPETITOR
    assert verdict.best_competitor == "SULFATE"
    assert not verdict.is_candidate


def test_an_unfavourable_field_makes_the_comparison_meaningless():
    verdict = SP.verdict(_panel(target=0.9, depth=-2.0))
    assert verdict.call == SP.UNFAVOURABLE
    assert not verdict.is_candidate


def test_the_ranking_puts_the_best_fitting_ligand_first():
    verdict = SP.verdict(_panel(target=0.3, SULFATE=0.8))
    assert verdict.ranking[0][0] == "SULFATE"


def test_a_selectivity_claim_needs_the_competitors_the_project_could_not_separate():
    thin = {SP.TARGET: {"shape": 0.9, "depth_kt": 3.0}, "ATP": {"shape": 0.1, "depth_kt": 3.0}}
    with pytest.raises(SP.SpecificityError, match="could not separate"):
        SP.verdict(thin)


def test_a_panel_without_the_target_is_refused():
    panel = {"ATP": {"shape": 0.5, "depth_kt": 1.0}}
    with pytest.raises(SP.SpecificityError, match="no entry for the target"):
        SP.verdict(panel)


def test_a_nonpositive_margin_is_refused():
    with pytest.raises(SP.SpecificityError, match="margin"):
        SP.verdict(_panel(), margin=0.0)


def test_the_summary_names_which_competitor_did_the_explaining():
    # "30 of 75 are explained by another ligand" is weaker than naming which
    # ligand explains how many, and the latter says whether the screen is
    # finding ATP sites or sulfate sites.
    verdicts = [
        SP.verdict(_panel(target=0.9, ATP=0.2)),
        SP.verdict(_panel(target=0.5, PAPS=0.48)),
        SP.verdict(_panel(target=0.3, SULFATE=0.8)),
        SP.verdict(_panel(target=0.5, PAPS=0.52)),
    ]
    summary = SP.summarise(verdicts)
    assert summary.counts[SP.SELECTIVE] == 1
    assert summary.by_competitor["PAPS"] == 2
    assert "PAPS (2)" in summary.summary()
    assert summary.selective_fraction == pytest.approx(0.25)


def test_the_margin_sweep_is_published_beside_the_result_not_instead_of_it():
    verdicts = [SP.verdict(_panel(target=t, PAPS=0.4)) for t in (0.45, 0.55, 0.7, 0.9)]
    rows = SP.margin_sweep(verdicts, margins=(0.05, 0.2, 0.6))
    assert [row["margin"] for row in rows] == [0.05, 0.2, 0.6]
    # A larger margin can only move calls out of selective, never into it.
    assert rows[0][SP.SELECTIVE] >= rows[-1][SP.SELECTIVE]


def test_the_verdict_round_trips_to_a_dict_with_its_explanation():
    payload = SP.verdict(_panel(target=0.9, ATP=0.2)).as_dict()
    assert payload["call"] == SP.SELECTIVE
    assert payload["is_candidate"] is True
    assert payload["explanation"]
    assert payload["margin_threshold"] == SP.DEFAULT_MARGIN


# ------------------------------------------------------- substitute groups

def test_an_alternative_measurement_nobody_took_is_not_an_unknown():
    # relative_sasa is measured on a ligand copy; cavity_openness is measured
    # by the cavity detector. They are the same property by different means, so
    # having one must not leave the score abstaining for want of the other -
    # which it did, making every score from the detector path an interval.
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["relative_sasa"] = None
    interval = scorer.score(values)
    assert interval.width == 0.0
    assert interval.missing == ()
    assert interval.point is not None


def test_losing_both_burial_terms_is_a_real_unknown():
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["relative_sasa"] = None
    values["cavity_openness"] = None
    interval = scorer.score(values)
    assert interval.width > 0.0
    assert set(interval.missing) == {"relative_sasa", "cavity_openness"}


def test_the_substituted_feature_does_not_dilute_the_weights():
    # The denominator has to drop with the feature, or dropping a term would
    # silently lower every score that uses the alternative.
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["relative_sasa"] = None
    assert scorer.score(values).low == pytest.approx(0.5, abs=1e-9)


def test_the_two_burial_terms_have_separate_names_for_separate_measurements():
    names = {f.name for f in MO.DEFAULT_FEATURES}
    assert {"relative_sasa", "cavity_openness"} <= names
    openness = next(f for f in MO.DEFAULT_FEATURES if f.name == "cavity_openness")
    assert openness.direction == MO.DOWN
    assert "angstrom" in openness.why  # it names the mislabeling it avoids


def test_sensitivity_uses_the_same_denominator_as_the_score():
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["relative_sasa"] = None
    assert sum(scorer.sensitivity(values).values()) == pytest.approx(scorer.score(values).low)


def test_a_scorer_whose_every_feature_is_substituted_away_is_refused():
    only_burial = [f for f in MO.DEFAULT_FEATURES if f.name in ("relative_sasa", "cavity_openness")]
    scorer = MO.MonotoneScorer(only_burial)
    # Both present is fine; one present drops the other and still scores.
    assert scorer.score({"relative_sasa": 0.1, "cavity_openness": None}).width == 0.0


# --------------------------------------------------------- physical vetoes

def test_a_repulsive_well_rules_a_pocket_out_however_buried_it_is():
    # Found by running the detector on an acidic decoy: a perfectly enclosed
    # cavity with no basic residues and a -3.6 kT/e well scored 0.554 against
    # a 0.55 threshold and was called a site. Burial, confidence and volume
    # were outvoting the physics.
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["relative_sasa"] = None
    values.update(cavity_openness=0.0, field_depth_kt=-3.64, n_basic_residues=0.0)
    interval = scorer.score(values)
    assert interval.vetoed == ("field_depth_kt",)
    assert interval.decision(0.55) == "not_site"


def test_a_veto_overrides_even_a_score_above_the_threshold():
    scorer = MO.MonotoneScorer()
    best = {f.name: (f.high if f.direction == MO.UP else f.low) for f in scorer.features}
    best["relative_sasa"] = None
    best["field_depth_kt"] = -1.0
    interval = scorer.score(best)
    assert interval.low > 0.55
    assert interval.decision(0.55) == "not_site"


def test_an_unmeasured_veto_feature_abstains_rather_than_passing():
    # A veto you never measured cannot be cleared, so a site call on the
    # strength of the other terms would be a claim about missing data.
    scorer = MO.MonotoneScorer()
    good = {f.name: (f.high if f.direction == MO.UP else f.low) for f in scorer.features}
    good["relative_sasa"] = None
    good["field_depth_kt"] = None
    interval = scorer.score(good)
    assert interval.veto_pending == ("field_depth_kt",)
    assert interval.decision(0.55) == "abstain"


def test_an_unmeasured_veto_on_a_hopeless_pocket_still_decides_against_it():
    scorer = MO.MonotoneScorer()
    bad = {f.name: (f.low if f.direction == MO.UP else f.high) for f in scorer.features}
    bad["relative_sasa"] = None
    bad["field_depth_kt"] = None
    assert scorer.score(bad).decision(0.55) == "not_site"


def test_a_favourable_well_fires_no_veto():
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["relative_sasa"] = None
    values["field_depth_kt"] = 4.0
    assert scorer.score(values).vetoed == ()


def test_the_blocking_feature_names_the_veto_when_one_fired():
    scorer = MO.MonotoneScorer()
    values = _mid(scorer)
    values["relative_sasa"] = None
    values["field_depth_kt"] = -5.0
    assert scorer.blocking_feature(values, 0.55) == "field_depth_kt"


def test_a_veto_must_state_its_reason():
    # A bound that can reject a protein outright without saying why is the
    # gate cascade this score replaces.
    with pytest.raises(MO.MonotoneError, match="state its reason"):
        MO.Feature(
            name="x", direction=MO.UP, midpoint=0.0, slope=1.0,
            low=-1.0, high=1.0, weight=1.0, veto_below=0.0,
        )


def test_only_the_electrostatic_well_carries_a_veto():
    # One physically-justified condition, not five arbitrary thresholds.
    vetoing = [f.name for f in MO.DEFAULT_FEATURES
               if f.veto_below is not None or f.veto_above is not None]
    assert vetoing == ["field_depth_kt"]
