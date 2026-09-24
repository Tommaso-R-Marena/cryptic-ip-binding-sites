"""Tests for the statistical utilities added for methodological rigour.

Each assertion checks a published property of the estimator, not a value copied
from the implementation's own output.
"""

import numpy as np
import pytest

from cryptic_ip.analysis.statistical_validation import StatisticalValidation as SV


def test_wilson_interval_stays_inside_the_unit_interval_at_zero_successes():
    """The Wald interval goes negative here; Wilson must not."""
    low, high = SV.wilson_interval(0, 1000)
    assert low == pytest.approx(0.0)
    assert 0.0 < high < 0.01


def test_wilson_interval_brackets_the_observed_proportion():
    low, high = SV.wilson_interval(25, 1000)
    assert low < 0.025 < high


def test_wilson_interval_narrows_as_the_sample_grows():
    narrow = SV.wilson_interval(50, 5000)
    wide = SV.wilson_interval(5, 500)
    assert (narrow[1] - narrow[0]) < (wide[1] - wide[0])


def test_wilson_interval_rejects_impossible_counts():
    with pytest.raises(ValueError):
        SV.wilson_interval(11, 10)
    assert np.isnan(SV.wilson_interval(0, 0)[0])


def test_clopper_pearson_is_conservative_relative_to_wilson():
    """The exact interval must be at least as wide as the score interval."""
    exact = SV.clopper_pearson_interval(5, 200)
    wilson = SV.wilson_interval(5, 200)
    assert (exact[1] - exact[0]) >= (wilson[1] - wilson[0]) - 1e-12


def test_clopper_pearson_boundaries_are_exact_at_the_extremes():
    assert SV.clopper_pearson_interval(0, 20)[0] == 0.0
    assert SV.clopper_pearson_interval(20, 20)[1] == 1.0


def test_holm_is_more_conservative_than_benjamini_hochberg():
    """Family-wise control must reject no more than FDR control."""
    p_values = [0.001, 0.008, 0.02, 0.04, 0.3, 0.7]
    holm = SV.holm_bonferroni(p_values)
    bh = SV.benjamini_hochberg(p_values)
    assert holm["reject"].sum() <= bh["reject"].sum()
    assert (holm["holm_adjusted"] >= np.asarray(p_values) - 1e-12).all()


def test_holm_adjusted_values_are_monotone_in_rank():
    holm = SV.holm_bonferroni([0.01, 0.02, 0.03, 0.5])
    ordered = holm.sort_values("p_value")["holm_adjusted"].to_numpy()
    assert np.all(np.diff(ordered) >= -1e-12)


def test_holm_matches_statsmodels_reference():
    p_values = [0.011, 0.02, 0.03, 0.5, 0.9]
    holm = SV.holm_bonferroni(p_values)["holm_adjusted"].to_numpy()
    statsmodels = pytest.importorskip("statsmodels.stats.multitest")
    expected = statsmodels.multipletests(p_values, method="holm")[1]
    assert np.allclose(holm, expected)


def test_holm_handles_empty_input():
    assert SV.holm_bonferroni([]).empty


def test_cliffs_delta_is_one_for_completely_separated_groups():
    assert SV.cliffs_delta([10, 11, 12], [1, 2, 3]) == pytest.approx(1.0)
    assert SV.cliffs_delta([1, 2, 3], [10, 11, 12]) == pytest.approx(-1.0)


def test_cliffs_delta_is_zero_for_identical_distributions():
    values = [1.0, 2.0, 3.0, 4.0]
    assert SV.cliffs_delta(values, values) == pytest.approx(0.0)


def test_cliffs_delta_of_empty_group_is_nan():
    assert np.isnan(SV.cliffs_delta([], [1.0, 2.0]))


def test_hedges_g_shrinks_cohens_d_for_small_samples():
    a = [10.0, 11.0, 12.0, 11.5]
    b = [1.0, 2.0, 1.5, 2.5]
    d = SV.cohens_d(a, b)
    g = SV.hedges_g(a, b)
    assert 0 < g < d


def test_methods_report_selects_a_nonparametric_test_when_normality_fails():
    """The reported test must follow the assumption check, not ignore it."""
    rng = np.random.default_rng(0)
    skewed_a = rng.exponential(1.0, 60) + 10.0
    skewed_b = rng.exponential(1.0, 60)
    report = SV.methods_report(skewed_a, skewed_b)
    assert "Mann-Whitney U test" in report
    assert "Cliff's delta" in report
    assert "Welch's two-sample t-test" not in report


def test_methods_report_uses_welch_when_normality_holds():
    rng = np.random.default_rng(1)
    report = SV.methods_report(rng.normal(5.0, 1.0, 60), rng.normal(0.0, 1.0, 60))
    assert "Welch's two-sample t-test" in report
    assert "Hedges' g" in report


def test_methods_report_requires_enough_observations():
    with pytest.raises(ValueError):
        SV.methods_report([1.0, 2.0], [3.0, 4.0])


def test_benjamini_hochberg_matches_statsmodels_reference():
    p_values = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205]
    ours = SV.benjamini_hochberg(p_values)["fdr_q_value"].to_numpy()
    statsmodels = pytest.importorskip("statsmodels.stats.multitest")
    expected = statsmodels.multipletests(p_values, method="fdr_bh")[1]
    assert np.allclose(ours, expected)


def test_required_sample_size_grows_as_the_effect_shrinks():
    small_effect = SV.required_sample_size_two_proportion(0.010, 0.012)
    large_effect = SV.required_sample_size_two_proportion(0.010, 0.050)
    assert small_effect > large_effect


def test_wilson_interval_coverage_is_near_nominal():
    """Empirical coverage check: the interval should cover the truth ~95 % of the time."""
    rng = np.random.default_rng(7)
    truth, trials = 0.01, 2000
    covered = 0
    for _ in range(400):
        successes = int(rng.binomial(trials, truth))
        low, high = SV.wilson_interval(successes, trials)
        covered += int(low <= truth <= high)
    assert covered / 400 > 0.90


def test_permutation_enrichment_returns_calibrated_pvalues_under_the_null():
    """With no real enrichment, p-values must not be systematically small."""
    rng = np.random.default_rng(3)
    categories = np.array(["a", "b", "c"] * 40)
    hits = rng.random(120) < 0.2
    result = SV.permutation_enrichment_test(
        categories, hits, n_permutations=400, random_state=1
    )
    assert (result["p_value"] > 0.01).all()
    assert set(result["category"]) == {"a", "b", "c"}
