"""Tests for statistical validation workflows."""

import numpy as np

from cryptic_ip.analysis.statistical_validation import StatisticalValidation


def test_roc_and_pr_bootstrap_outputs():
    """Bootstrap ROC/PR analyses should return sensible bounds and AUCs."""
    y_true = np.array([0, 0, 0, 1, 1, 1, 0, 1, 0, 1])
    y_score = np.array([0.1, 0.2, 0.35, 0.7, 0.8, 0.9, 0.4, 0.75, 0.15, 0.6])

    roc_result = StatisticalValidation.roc_with_bootstrap_ci(
        y_true,
        y_score,
        n_bootstrap=200,
        random_state=7,
    )
    pr_result = StatisticalValidation.precision_recall_with_bootstrap_ci(
        y_true,
        y_score,
        n_bootstrap=200,
        random_state=7,
    )

    assert 0 <= roc_result.auc_value <= 1
    assert 0 <= pr_result.auc_value <= 1
    assert roc_result.auc_ci[0] <= roc_result.auc_ci[1]
    assert pr_result.auc_ci[0] <= pr_result.auc_ci[1]


def test_fdr_and_power_and_effect_size_workflows():
    """FDR correction, power analysis, and effect size estimates should be stable."""
    fdr = StatisticalValidation.benjamini_hochberg([0.001, 0.01, 0.04, 0.2, 0.8], alpha=0.05)
    assert "fdr_q_value" in fdr.columns
    assert fdr["reject"].sum() >= 1

    power_df = StatisticalValidation.pairwise_hit_rate_power(
        {"yeast": 0.03, "human": 0.015, "ddisc": 0.025},
        alpha=0.05,
        power=0.8,
    )
    assert len(power_df) == 3
    assert (power_df["required_n_per_group"] > 0).all()

    effect = StatisticalValidation.cohens_d([0.3, 0.35, 0.33, 0.31], [0.22, 0.2, 0.24, 0.18])
    assert effect > 0


def test_permutation_and_methods_report():
    """Permutation enrichment and methods reporting should produce usable outputs."""
    categories = np.array(["kinase", "kinase", "tf", "tf", "enzyme", "enzyme", "enzyme"])
    hit_mask = np.array([1, 1, 0, 1, 0, 0, 1], dtype=bool)

    perm = StatisticalValidation.permutation_enrichment_test(
        categories,
        hit_mask,
        n_permutations=300,
        random_state=11,
    )
    assert set(["category", "p_value", "fdr_q_value"]).issubset(perm.columns)

    report = StatisticalValidation.methods_report(
        [0.5, 0.6, 0.58, 0.61, 0.63],
        [0.45, 0.4, 0.42, 0.44, 0.39],
    )
    assert "Shapiro-Wilk" in report
    assert "Cohen's d" in report
