"""Positive-unlabelled learning, conformal sets, and the power calculation.

The estimator is tested by handing it data whose answer is known: a synthetic
population with a chosen label frequency, and - more decisively - the exact
non-traditional classifier, which isolates the estimator from any model's
quality.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import pu as PU


def _population(n_groups=60, per_group=30, seed=1):
    """Groups of correlated members, like homology groups of one protein."""
    rng = np.random.default_rng(seed)
    groups, features, truth = [], [], []
    for g in range(n_groups):
        for _ in range(per_group):
            x = rng.normal(size=2)
            groups.append(f"g{g}")
            features.append(x)
            truth.append(rng.random() < 1.0 / (1.0 + np.exp(-2.0 * x[0])))
    return np.asarray(features), np.asarray(truth), groups, rng


# ------------------------------------------------------------ leakage audit

def test_a_clean_split_shares_no_group():
    groups = ["a", "a", "b", "b", "c", "c"]
    audit = PU.audit_split(groups, [0, 1, 2, 3], [4, 5])
    assert audit.clean and audit.train_groups == 2 and audit.test_groups == 1


def test_a_split_that_shares_a_group_is_reported_as_leaking():
    groups = ["a", "a", "b", "b"]
    audit = PU.audit_split(groups, [0, 2], [1, 3])
    assert not audit.clean
    assert set(audit.shared_groups) == {"a", "b"}


def test_group_folds_never_split_a_group():
    groups = [f"g{i // 7}" for i in range(70)]
    for train, test in PU.group_folds(groups, 5, seed=0):
        assert PU.audit_split(groups, train, test).clean


def test_group_folds_cover_every_row_exactly_once():
    groups = [f"g{i // 7}" for i in range(70)]
    folds = PU.group_folds(groups, 5, seed=0)
    seen = sorted(i for _, test in folds for i in test)
    assert seen == list(range(70))


def test_folds_spread_a_dominant_group_rather_than_failing():
    # One group holding most of the rows is this project's actual situation:
    # 87% of buried-site positives are in a single homology group.
    groups = ["big"] * 80 + [f"g{i}" for i in range(10)]
    folds = PU.group_folds(groups, 5, seed=0)
    sizes = sorted(len(test) for _, test in folds)
    assert sizes[-1] == 80  # the big group cannot be split, and is not
    assert all(PU.audit_split(groups, tr, te).clean for tr, te in folds)


def test_fewer_groups_than_folds_is_refused():
    with pytest.raises(PU.PULearningError, match="cannot be split"):
        PU.group_folds(["a", "a", "b"], 5, seed=0)


def test_one_fold_is_refused():
    with pytest.raises(PU.PULearningError, match="at least 2 folds"):
        PU.group_folds(["a", "b"], 1, seed=0)


# -------------------------------------------------------- label frequency

def test_the_quantile_estimator_beats_the_mean_on_overlapping_classes():
    # Elkan and Noto's mean estimator is biased low by exactly the class
    # overlap, and the posterior divides by c - so a c that is too small
    # inflates P(positive). The bias is in the dangerous direction, which is
    # why the default is a quantile.
    features, truth, groups, rng = _population()
    classifier = PU.logistic_classifier(l2=0.01, iterations=3000, learning_rate=0.5)
    labelled = truth & (rng.random(len(truth)) < 0.5)
    quantile = PU.fit_pu(features, labelled, groups, classifier, n_folds=5).frequency.c
    mean = PU.fit_pu(features, labelled, groups, classifier, n_folds=5, method="mean").frequency.c
    assert abs(quantile - 0.5) < abs(mean - 0.5)


def test_with_the_exact_classifier_the_quantile_estimator_is_nearly_unbiased():
    # This isolates the estimator from any model: it is handed the true
    # g(x) = c * P(y=1|x), so whatever bias remains is the estimator's own.
    rng = np.random.default_rng(3)
    n = 40000
    x = rng.normal(size=n)
    p = 1.0 / (1.0 + np.exp(-2.0 * x))
    y = rng.random(n) < p
    groups = [f"g{i % 50}" for i in range(n)]
    for true_c in (0.2, 0.5, 0.9):
        labelled = y & (rng.random(n) < true_c)
        estimate = PU.estimate_label_frequency(true_c * p, labelled, groups).c
        assert 0.92 < estimate / true_c < 1.08, f"c={true_c} recovered as {estimate}"


def test_the_estimate_is_averaged_over_groups_not_rows():
    # One family holding most of the positives must not set c on its own.
    scores = np.concatenate([np.full(100, 0.9), np.full(2, 0.1)])
    labelled = np.ones(102, dtype=bool)
    groups = ["big"] * 100 + ["small", "small"]
    estimate = PU.estimate_label_frequency(scores, labelled, groups, method="mean")
    assert estimate.c == pytest.approx(0.5)      # (0.9 + 0.1) / 2, not 0.884
    assert estimate.n_groups == 2


def test_the_group_spread_is_reported_as_the_error_bar():
    scores = np.array([0.9, 0.9, 0.1, 0.1])
    estimate = PU.estimate_label_frequency(scores, np.ones(4, bool), ["a", "a", "b", "b"])
    assert estimate.group_spread > 0.5
    assert estimate.thinnest_group == 2


def test_five_groups_is_the_floor_for_the_estimate_being_evidence():
    scores = np.full(4, 0.5)
    estimate = PU.estimate_label_frequency(scores, np.ones(4, bool), ["a", "b", "c", "d"])
    assert estimate.n_groups == 4 and not estimate.is_evidence
    bigger = PU.estimate_label_frequency(
        np.full(5, 0.5), np.ones(5, bool), ["a", "b", "c", "d", "e"],
    )
    assert bigger.is_evidence


def test_no_labelled_positives_is_refused():
    with pytest.raises(PU.PULearningError, match="no labelled positives"):
        PU.estimate_label_frequency([0.5, 0.5], [False, False], ["a", "b"])


def test_scores_outside_zero_to_one_are_refused():
    with pytest.raises(PU.PULearningError, match="probabilities"):
        PU.estimate_label_frequency([1.5], [True], ["a"])


def test_an_unknown_method_is_refused():
    with pytest.raises(PU.PULearningError, match="method must be"):
        PU.estimate_label_frequency([0.5], [True], ["a"], method="eyeball")


# ------------------------------------------------------------- the posterior

def test_the_posterior_divides_by_the_label_frequency():
    frequency = PU.LabelFrequency(c=0.5, per_group={}, per_group_n={}, n_groups=5, method="mean")
    assert PU.pu_posterior([0.1, 0.25], frequency).tolist() == [0.2, 0.5]


def test_the_posterior_saturates_rather_than_exceeding_one():
    # c is estimated, so g(x)/c can pass 1. Clipping is visible; a probability
    # above one that a later step divides by is not.
    frequency = PU.LabelFrequency(c=0.5, per_group={}, per_group_n={}, n_groups=5, method="mean")
    assert PU.pu_posterior([0.9], frequency)[0] == 1.0


# -------------------------------------------------------------- conformal

def test_the_conformal_threshold_uses_the_finite_sample_correction():
    # With few calibration points the corrected quantile and the plain one
    # differ, and the difference is the coverage being quietly lower.
    scores = np.linspace(0.0, 1.0, 20)
    corrected = PU.conformal_threshold(scores, alpha=0.1)
    plain = float(np.quantile(scores, 0.9))
    assert corrected >= plain


def test_too_few_calibration_points_for_the_confidence_is_refused():
    with pytest.raises(PU.PULearningError, match="cannot support"):
        PU.conformal_threshold([0.5] * 5, alpha=0.01)


def test_conformal_sets_can_hold_both_labels_which_is_the_method_abstaining():
    rng = np.random.default_rng(0)
    calibration = rng.uniform(0.4, 0.6, 60)
    groups = [f"g{i % 10}" for i in range(60)]
    sets = PU.conformal_sets([0.5], calibration, groups, alpha=0.1)
    assert sets.contains_positive[0] and sets.contains_negative[0]
    assert sets.sizes[0] == 2
    assert sets.fraction_decided == 0.0


def test_a_confident_row_gets_a_single_label():
    rng = np.random.default_rng(0)
    calibration = rng.uniform(0.8, 0.95, 60)
    groups = [f"g{i % 10}" for i in range(60)]
    sets = PU.conformal_sets([0.99], calibration, groups, alpha=0.1)
    assert sets.sizes[0] == 1 and sets.contains_positive[0]


def test_conformal_coverage_is_close_to_its_target():
    # The guarantee is the point of the method, so it is checked rather than
    # cited: held-out positives must fall in their set about 1-alpha of the time.
    rng = np.random.default_rng(5)
    scores = rng.beta(5, 2, 4000)
    groups = [f"g{i % 40}" for i in range(4000)]
    calibration, test = scores[:2000], scores[2000:]
    sets = PU.conformal_sets(test, calibration, groups[:2000], alpha=0.1)
    coverage = float(np.mean(sets.contains_positive))
    assert 0.85 < coverage < 0.95


def test_calibrating_on_too_few_groups_is_refused():
    with pytest.raises(PU.PULearningError, match="not evidence"):
        PU.conformal_sets([0.5], [0.5] * 20, ["a"] * 10 + ["b"] * 10, alpha=0.1)


def test_an_alpha_outside_the_unit_interval_is_refused():
    with pytest.raises(PU.PULearningError, match="alpha"):
        PU.conformal_threshold([0.1, 0.2], alpha=1.5)


# ----------------------------------------------------- minimum detectable

def test_the_minimum_detectable_auc_falls_as_groups_are_added():
    values = [PU.minimum_detectable_auc(g, 100) for g in (5, 10, 20, 40, 100)]
    assert values == sorted(values, reverse=True)


def test_five_groups_can_only_detect_a_very_large_effect():
    # This is the number that reframes the project's 0.499 AUC. With
    # buried-site positives in five homology groups, nothing below about 0.87
    # is detectable, so a null result is a statement about power rather than
    # about biology.
    assert PU.minimum_detectable_auc(5, 100) > 0.8


def test_a_single_group_cannot_detect_anything():
    assert PU.minimum_detectable_auc(1, 100) == 1.0


def test_an_empty_class_is_refused():
    with pytest.raises(PU.PULearningError, match="at least one group"):
        PU.minimum_detectable_auc(0, 10)


def test_an_unimplemented_power_is_refused_rather_than_approximated():
    with pytest.raises(PU.PULearningError, match="only power=0.8"):
        PU.minimum_detectable_auc(10, 10, power=0.9)


# -------------------------------------------------------------- group AUC

def test_the_group_auc_separates_groups_not_rows():
    # A row-level AUC over 134 entries in 5 families reports the separability
    # of five things with the confidence of 134.
    scores = [0.9] * 50 + [0.1] * 2
    positive = [True] * 50 + [False] * 2
    groups = ["pos"] * 50 + ["neg1", "neg2"]
    assert PU.group_auc(scores, positive, groups) == pytest.approx(1.0)


def test_the_group_auc_is_half_when_the_groups_are_indistinguishable():
    scores = [0.5, 0.5, 0.5, 0.5]
    assert PU.group_auc(scores, [True, True, False, False], ["a", "a", "b", "b"]) == pytest.approx(0.5)


def test_a_group_holding_both_classes_is_refused():
    with pytest.raises(PU.PULearningError, match="both classes"):
        PU.group_auc([0.1, 0.9], [True, False], ["a", "a"])


def test_a_group_auc_needs_both_classes():
    with pytest.raises(PU.PULearningError, match="both classes"):
        PU.group_auc([0.1, 0.9], [True, True], ["a", "b"])


# ------------------------------------------------------------------ the fit

def test_the_fit_scores_every_row_out_of_fold():
    features, truth, groups, rng = _population(n_groups=30, per_group=10)
    labelled = truth & (rng.random(len(truth)) < 0.5)
    fit = PU.fit_pu(features, labelled, groups, PU.logistic_classifier(), n_folds=5)
    assert not np.any(np.isnan(fit.out_of_fold_scores))
    assert fit.leakage_clean
    assert fit.posterior.shape == (len(features),)


def test_the_fit_records_the_power_it_had_before_any_result_is_read():
    features, truth, groups, rng = _population(n_groups=30, per_group=10)
    labelled = truth & (rng.random(len(truth)) < 0.5)
    fit = PU.fit_pu(features, labelled, groups, PU.logistic_classifier(), n_folds=5)
    assert fit.minimum_detectable_auc is not None
    assert 0.5 < fit.minimum_detectable_auc <= 1.0


def test_the_fit_is_reproducible_from_its_seed():
    features, truth, groups, rng = _population(n_groups=30, per_group=10)
    labelled = truth & (rng.random(len(truth)) < 0.5)
    a = PU.fit_pu(features, labelled, groups, PU.logistic_classifier(), seed=7)
    b = PU.fit_pu(features, labelled, groups, PU.logistic_classifier(), seed=7)
    assert np.allclose(a.out_of_fold_scores, b.out_of_fold_scores)


def test_a_classifier_returning_the_wrong_shape_is_refused():
    features, truth, groups, rng = _population(n_groups=20, per_group=5)
    labelled = truth & (rng.random(len(truth)) < 0.5)

    def broken(train_x, train_y, test_x):
        return np.zeros(3)

    with pytest.raises(PU.PULearningError, match="one score per test row"):
        PU.fit_pu(features, labelled, groups, broken)


def test_the_fit_round_trips_to_a_dict():
    features, truth, groups, rng = _population(n_groups=30, per_group=10)
    labelled = truth & (rng.random(len(truth)) < 0.5)
    payload = PU.fit_pu(features, labelled, groups, PU.logistic_classifier()).as_dict()
    assert payload["leakage_clean"] is True
    assert "minimum_detectable_auc" in payload
    assert payload["label_frequency"]["method"] == "quantile"
