"""Statistical validation utilities for proteome-wide cryptic site screening."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def _roc_curve(y_true: np.ndarray, y_score: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    order = np.argsort(-y_score)
    y = y_true[order]
    positives = np.sum(y == 1)
    negatives = np.sum(y == 0)

    tps = np.cumsum(y == 1)
    fps = np.cumsum(y == 0)

    tpr = np.concatenate(([0.0], tps / positives, [1.0]))
    fpr = np.concatenate(([0.0], fps / negatives, [1.0]))
    return fpr, tpr


def _roc_auc_score(y_true: np.ndarray, y_score: np.ndarray) -> float:
    fpr, tpr = _roc_curve(y_true, y_score)
    return float(np.trapz(tpr, fpr))


def _precision_recall_curve(y_true: np.ndarray, y_score: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    order = np.argsort(-y_score)
    y = y_true[order]
    tp = np.cumsum(y == 1)
    fp = np.cumsum(y == 0)
    total_pos = np.sum(y_true == 1)

    precision = tp / np.maximum(tp + fp, 1)
    recall = tp / total_pos

    precision = np.concatenate(([1.0], precision, [y_true.mean()]))
    recall = np.concatenate(([0.0], recall, [1.0]))
    return precision, recall


def _average_precision_score(y_true: np.ndarray, y_score: np.ndarray) -> float:
    precision, recall = _precision_recall_curve(y_true, y_score)
    sort_idx = np.argsort(recall)
    return float(np.trapz(precision[sort_idx], recall[sort_idx]))


@dataclass
class BootstrapCurveResult:
    """Container for bootstrapped curve metrics and confidence intervals."""

    x: np.ndarray
    y: np.ndarray
    y_lower: np.ndarray
    y_upper: np.ndarray
    auc_value: float
    auc_ci: Tuple[float, float]


class StatisticalValidation:
    """High-level statistical workflows for validating cryptic IP site predictions."""

    @staticmethod
    def _rng(seed: Optional[int]) -> np.random.Generator:
        return np.random.default_rng(seed)

    @staticmethod
    def _check_binary(y_true: Sequence[int]) -> np.ndarray:
        y_arr = np.asarray(y_true)
        uniques = np.unique(y_arr)
        if not np.array_equal(np.sort(uniques), np.array([0, 1])):
            raise ValueError("y_true must contain binary labels encoded as 0/1")
        return y_arr

    @staticmethod
    def roc_with_bootstrap_ci(
        y_true: Sequence[int],
        y_score: Sequence[float],
        n_bootstrap: int = 2000,
        ci: float = 0.95,
        random_state: Optional[int] = None,
        grid_size: int = 200,
    ) -> BootstrapCurveResult:
        """Compute ROC curve and bootstrapped confidence intervals."""
        y_arr = StatisticalValidation._check_binary(y_true)
        score_arr = np.asarray(y_score, dtype=float)

        fpr, tpr = _roc_curve(y_arr, score_arr)
        roc_auc = _roc_auc_score(y_arr, score_arr)

        x_grid = np.linspace(0, 1, grid_size)
        interp_tprs: List[np.ndarray] = []
        boot_aucs: List[float] = []
        rng = StatisticalValidation._rng(random_state)

        for _ in range(n_bootstrap):
            sample_idx = rng.integers(0, len(y_arr), size=len(y_arr))
            boot_y = y_arr[sample_idx]
            if len(np.unique(boot_y)) < 2:
                continue
            boot_score = score_arr[sample_idx]
            b_fpr, b_tpr = _roc_curve(boot_y, boot_score)
            interp_tprs.append(np.interp(x_grid, b_fpr, b_tpr))
            boot_aucs.append(_roc_auc_score(boot_y, boot_score))

        if not interp_tprs:
            raise RuntimeError("Bootstrapping failed: no valid resamples with both classes.")

        interp_tprs_arr = np.asarray(interp_tprs)
        alpha = (1 - ci) / 2
        y_lower = np.quantile(interp_tprs_arr, alpha, axis=0)
        y_upper = np.quantile(interp_tprs_arr, 1 - alpha, axis=0)
        auc_ci = (
            float(np.quantile(boot_aucs, alpha)),
            float(np.quantile(boot_aucs, 1 - alpha)),
        )

        return BootstrapCurveResult(
            x=x_grid,
            y=np.interp(x_grid, fpr, tpr),
            y_lower=y_lower,
            y_upper=y_upper,
            auc_value=float(roc_auc),
            auc_ci=auc_ci,
        )

    @staticmethod
    def precision_recall_with_bootstrap_ci(
        y_true: Sequence[int],
        y_score: Sequence[float],
        n_bootstrap: int = 2000,
        ci: float = 0.95,
        random_state: Optional[int] = None,
        grid_size: int = 200,
    ) -> BootstrapCurveResult:
        """Compute precision-recall curve and bootstrapped confidence intervals."""
        y_arr = StatisticalValidation._check_binary(y_true)
        score_arr = np.asarray(y_score, dtype=float)

        precision, recall = _precision_recall_curve(y_arr, score_arr)
        pr_auc = _average_precision_score(y_arr, score_arr)

        x_grid = np.linspace(0, 1, grid_size)
        interp_precisions: List[np.ndarray] = []
        boot_aucs: List[float] = []
        rng = StatisticalValidation._rng(random_state)

        for _ in range(n_bootstrap):
            sample_idx = rng.integers(0, len(y_arr), size=len(y_arr))
            boot_y = y_arr[sample_idx]
            if len(np.unique(boot_y)) < 2:
                continue
            boot_score = score_arr[sample_idx]
            b_precision, b_recall = _precision_recall_curve(boot_y, boot_score)
            sort_idx = np.argsort(b_recall)
            interp_precisions.append(np.interp(x_grid, b_recall[sort_idx], b_precision[sort_idx]))
            boot_aucs.append(_average_precision_score(boot_y, boot_score))

        if not interp_precisions:
            raise RuntimeError("Bootstrapping failed: no valid resamples with both classes.")

        interp_precisions_arr = np.asarray(interp_precisions)
        alpha = (1 - ci) / 2
        y_lower = np.quantile(interp_precisions_arr, alpha, axis=0)
        y_upper = np.quantile(interp_precisions_arr, 1 - alpha, axis=0)
        auc_ci = (
            float(np.quantile(boot_aucs, alpha)),
            float(np.quantile(boot_aucs, 1 - alpha)),
        )

        sort_idx = np.argsort(recall)
        return BootstrapCurveResult(
            x=x_grid,
            y=np.interp(x_grid, recall[sort_idx], precision[sort_idx]),
            y_lower=y_lower,
            y_upper=y_upper,
            auc_value=float(pr_auc),
            auc_ci=auc_ci,
        )

    @staticmethod
    def benjamini_hochberg(p_values: Sequence[float], alpha: float = 0.05) -> pd.DataFrame:
        """Apply Benjamini-Hochberg false discovery rate correction."""
        p_arr = np.asarray(p_values, dtype=float)
        if np.any((p_arr < 0) | (p_arr > 1)):
            raise ValueError("All p-values must be in [0, 1].")

        m = len(p_arr)
        order = np.argsort(p_arr)
        ranked_p = p_arr[order]
        ranks = np.arange(1, m + 1)

        adjusted = ranked_p * m / ranks
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        adjusted = np.clip(adjusted, 0, 1)

        reject_ranked = ranked_p <= (ranks / m) * alpha
        if np.any(reject_ranked):
            max_i = np.where(reject_ranked)[0].max()
            reject_ranked = np.arange(m) <= max_i

        reject = np.zeros(m, dtype=bool)
        reject[order] = reject_ranked

        adj_original = np.zeros(m)
        adj_original[order] = adjusted

        return pd.DataFrame({"p_value": p_arr, "fdr_q_value": adj_original, "reject": reject})

    @staticmethod
    def required_sample_size_two_proportion(
        p1: float,
        p2: float,
        alpha: float = 0.05,
        power: float = 0.8,
        two_sided: bool = True,
    ) -> int:
        """Estimate per-group sample size for detecting difference in hit rates."""
        if not (0 < p1 < 1 and 0 < p2 < 1):
            raise ValueError("p1 and p2 must be between 0 and 1.")
        if p1 == p2:
            raise ValueError("p1 and p2 must differ for power analysis.")

        delta = abs(p1 - p2)
        pbar = (p1 + p2) / 2
        z_alpha = stats.norm.ppf(1 - alpha / 2) if two_sided else stats.norm.ppf(1 - alpha)
        z_beta = stats.norm.ppf(power)

        numerator = (
            z_alpha * np.sqrt(2 * pbar * (1 - pbar))
            + z_beta * np.sqrt(p1 * (1 - p1) + p2 * (1 - p2))
        ) ** 2
        n = numerator / (delta**2)
        return int(np.ceil(n))

    @staticmethod
    def pairwise_hit_rate_power(
        hit_rates: Mapping[str, float],
        alpha: float = 0.05,
        power: float = 0.8,
    ) -> pd.DataFrame:
        """Calculate pairwise required sample sizes for organism-level hit-rate comparisons."""
        rows = []
        for org1, org2 in combinations(hit_rates.keys(), 2):
            p1 = hit_rates[org1]
            p2 = hit_rates[org2]
            required = StatisticalValidation.required_sample_size_two_proportion(
                p1=p1,
                p2=p2,
                alpha=alpha,
                power=power,
            )
            rows.append(
                {
                    "organism_a": org1,
                    "organism_b": org2,
                    "hit_rate_a": p1,
                    "hit_rate_b": p2,
                    "required_n_per_group": required,
                }
            )

        return pd.DataFrame(rows)

    @staticmethod
    def permutation_enrichment_test(
        categories: Sequence[str],
        hit_mask: Sequence[bool],
        n_permutations: int = 5000,
        random_state: Optional[int] = None,
    ) -> pd.DataFrame:
        """Permutation test for enrichment of hits within functional categories."""
        category_arr = np.asarray(categories)
        hit_arr = np.asarray(hit_mask, dtype=bool)

        if len(category_arr) != len(hit_arr):
            raise ValueError("categories and hit_mask must have same length")

        observed = []
        rng = StatisticalValidation._rng(random_state)
        overall_rate = hit_arr.mean()

        for cat in np.unique(category_arr):
            cat_mask = category_arr == cat
            observed_rate = hit_arr[cat_mask].mean()
            observed_delta = observed_rate - overall_rate

            perm_deltas = np.zeros(n_permutations)
            for idx in range(n_permutations):
                perm_hits = rng.permutation(hit_arr)
                perm_deltas[idx] = perm_hits[cat_mask].mean() - perm_hits.mean()

            p_val = (np.sum(np.abs(perm_deltas) >= abs(observed_delta)) + 1) / (n_permutations + 1)
            observed.append(
                {
                    "category": cat,
                    "observed_hit_rate": observed_rate,
                    "overall_hit_rate": overall_rate,
                    "delta_hit_rate": observed_delta,
                    "p_value": p_val,
                }
            )

        results = pd.DataFrame(observed)
        fdr = StatisticalValidation.benjamini_hochberg(results["p_value"].values)
        results["fdr_q_value"] = fdr["fdr_q_value"].values
        results["reject_fdr_0_05"] = fdr["reject"].values
        return results.sort_values("p_value").reset_index(drop=True)

    @staticmethod
    def cohens_d(group_a: Sequence[float], group_b: Sequence[float]) -> float:
        """Compute Cohen's d effect size for two independent groups."""
        a = np.asarray(group_a, dtype=float)
        b = np.asarray(group_b, dtype=float)
        if len(a) < 2 or len(b) < 2:
            raise ValueError("Each group must contain at least two observations.")

        pooled_sd = np.sqrt(
            ((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1))
            / (len(a) + len(b) - 2)
        )
        if pooled_sd == 0:
            raise ValueError("Pooled standard deviation is zero; Cohen's d undefined.")
        return float((np.mean(a) - np.mean(b)) / pooled_sd)

    @staticmethod
    def publication_style() -> None:
        """Apply publication-ready style inspired by Nature/Science figure conventions."""
        plt.style.use("default")
        plt.rcParams.update(
            {
                "font.family": "sans-serif",
                "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                "font.size": 8,
                "axes.labelsize": 8,
                "axes.titlesize": 9,
                "axes.linewidth": 0.8,
                "xtick.labelsize": 7,
                "ytick.labelsize": 7,
                "xtick.major.width": 0.8,
                "ytick.major.width": 0.8,
                "legend.fontsize": 7,
                "figure.dpi": 300,
                "savefig.dpi": 600,
                "savefig.bbox": "tight",
                "axes.spines.top": False,
                "axes.spines.right": False,
            }
        )

    @staticmethod
    def plot_roc_pr_panels(
        roc_result: BootstrapCurveResult,
        pr_result: BootstrapCurveResult,
        output_path: Optional[str] = None,
    ) -> plt.Figure:
        """Create publication-quality ROC and PR figure panels."""
        StatisticalValidation.publication_style()
        fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.0))

        ax = axes[0]
        ax.plot(roc_result.x, roc_result.y, color="#1b9e77", lw=1.4, label=f"AUROC={roc_result.auc_value:.3f}")
        ax.fill_between(roc_result.x, roc_result.y_lower, roc_result.y_upper, color="#1b9e77", alpha=0.2)
        ax.plot([0, 1], [0, 1], "--", color="0.5", lw=1)
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_title("ROC curve")
        ax.legend(frameon=False, loc="lower right")

        ax = axes[1]
        ax.plot(pr_result.x, pr_result.y, color="#d95f02", lw=1.4, label=f"AP={pr_result.auc_value:.3f}")
        ax.fill_between(pr_result.x, pr_result.y_lower, pr_result.y_upper, color="#d95f02", alpha=0.2)
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_title("Precision-recall curve")
        ax.legend(frameon=False, loc="lower left")

        fig.tight_layout()
        if output_path:
            fig.savefig(output_path)
        return fig

    @staticmethod
    def methods_report(
        group_a: Sequence[float],
        group_b: Sequence[float],
        alpha: float = 0.05,
    ) -> str:
        """Generate methods-style statistical report including assumption checks."""
        a = np.asarray(group_a, dtype=float)
        b = np.asarray(group_b, dtype=float)

        shapiro_a = stats.shapiro(a)
        shapiro_b = stats.shapiro(b)
        levene_res = stats.levene(a, b, center="median")
        t_res = stats.ttest_ind(a, b, equal_var=levene_res.pvalue >= alpha)
        effect = StatisticalValidation.cohens_d(a, b)

        normality_ok = shapiro_a.pvalue >= alpha and shapiro_b.pvalue >= alpha
        variance_ok = levene_res.pvalue >= alpha

        assumptions = [
            f"Normality assessed by Shapiro-Wilk: group A p={shapiro_a.pvalue:.4g}, "
            f"group B p={shapiro_b.pvalue:.4g} ({'met' if normality_ok else 'violated'}).",
            f"Homoscedasticity assessed by Levene's test: p={levene_res.pvalue:.4g} "
            f"({'met' if variance_ok else 'violated'}).",
            "Independence of observations assumed by design (distinct proteins screened per organism).",
        ]

        lines = [
            "Statistical analyses were performed using SciPy and scikit-learn in Python.",
            *assumptions,
            (
                "Two-sample t-test was used for comparative proteomic summary metrics "
                f"(t={t_res.statistic:.3f}, p={t_res.pvalue:.4g}, "
                f"two-sided, alpha={alpha})."
            ),
            f"Effect size was reported as Cohen's d={effect:.3f}.",
            "Multiple testing correction used Benjamini-Hochberg FDR control.",
            "Classification performance was quantified via AUROC and average precision with "
            "95% bootstrap confidence intervals.",
            "Functional enrichment significance was evaluated using two-sided permutation tests.",
        ]

        return "\n".join(lines)
