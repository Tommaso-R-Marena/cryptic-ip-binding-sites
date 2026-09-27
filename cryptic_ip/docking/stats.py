"""Estimands and group-bootstrap intervals for docking outcomes.

Two estimands for any per-copy outcome ``y``:

* **per copy** - every copy weighted equally;
* **per group** - every homology group weighted equally (the mean of within-group
  means), so a family with many deposited entries (ADAR2) counts once.

Intervals resample whole groups (``protocol.group_bootstrap_weights``), never copies:
a resample weights each copy by how often its group was drawn. Under that weighting
the per-group estimand is ``sum_g c_g m_g / sum_g c_g`` with ``c_g`` the draw count.
"""

from __future__ import annotations

from typing import Callable, Dict, Optional, Sequence

import numpy as np
import pandas as pd

from ..benchmark import protocol

N_BOOTSTRAP = 2000
MIN_GROUPS = 5


def _codes(groups: Sequence) -> np.ndarray:
    return pd.factorize(pd.Series(np.asarray(groups)).astype(str))[0]


def per_copy_mean(y: np.ndarray, w: np.ndarray) -> float:
    total = float(np.sum(w))
    return float(np.sum(w * y) / total) if total > 0 else float("nan")


def per_group_mean(y: np.ndarray, codes: np.ndarray, w: np.ndarray) -> float:
    """Mean of within-group means, groups weighted by ``w`` (constant within a group)."""
    sizes = np.bincount(codes)
    share = w / sizes[codes]  # each row carries c_g / n_g, so a group sums to c_g
    total = float(np.sum(share))
    return float(np.sum(share * y) / total) if total > 0 else float("nan")


def estimate(stat: Callable[[np.ndarray], float], groups: Sequence, *, null: float = 0.0,
             n_bootstrap: int = N_BOOTSTRAP, seed: int = 20260926) -> Dict[str, float]:
    e = protocol.bootstrap_statistic(stat, np.asarray(groups), null=null, n_bootstrap=n_bootstrap, seed=seed)
    return e.as_dict()


def mean_estimates(y: Sequence[float], groups: Sequence, *, null: float = 0.5, n_bootstrap: int = N_BOOTSTRAP,
                   seed: int = 20260926) -> Dict[str, object]:
    """Per-copy and per-group means of ``y`` with group-bootstrap intervals and counts."""
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(y)
    y, groups = y[keep], np.asarray(groups)[keep]
    n_groups = int(pd.Series(groups).nunique())
    out: Dict[str, object] = {"copies": int(len(y)), "groups": n_groups,
                              "evidence": n_groups >= MIN_GROUPS}
    if len(y) == 0:
        return out
    codes = _codes(groups)
    out["per_copy"] = estimate(lambda w: per_copy_mean(y, w), groups, null=null, n_bootstrap=n_bootstrap, seed=seed)
    out["per_group"] = estimate(lambda w: per_group_mean(y, codes, w), groups, null=null,
                                n_bootstrap=n_bootstrap, seed=seed)
    return out


def difference_estimates(y: Sequence[float], in_a: Sequence[bool], in_b: Sequence[bool], groups: Sequence, *,
                         n_bootstrap: int = N_BOOTSTRAP, seed: int = 20260926) -> Dict[str, object]:
    """mean(a) - mean(b), both estimands, both strata resampled together."""
    y = np.asarray(y, dtype=float)
    a, b = np.asarray(in_a, bool) & np.isfinite(y), np.asarray(in_b, bool) & np.isfinite(y)
    groups = np.asarray(groups)
    keep = a | b
    y, a, b, groups = y[keep], a[keep], b[keep], groups[keep]
    ga, gb = int(pd.Series(groups[a]).nunique()), int(pd.Series(groups[b]).nunique())
    out: Dict[str, object] = {"copies": [int(a.sum()), int(b.sum())], "groups": [ga, gb],
                              "evidence": min(ga, gb) >= MIN_GROUPS}
    if a.sum() == 0 or b.sum() == 0:
        return out
    ca, cb = _codes(groups[a]), _codes(groups[b])

    def copy_diff(w):
        return per_copy_mean(y[a], w[a]) - per_copy_mean(y[b], w[b])

    def group_diff(w):
        return per_group_mean(y[a], ca, w[a]) - per_group_mean(y[b], cb, w[b])

    out["per_copy"] = estimate(copy_diff, groups, n_bootstrap=n_bootstrap, seed=seed)
    out["per_group"] = estimate(group_diff, groups, n_bootstrap=n_bootstrap, seed=seed)
    return out


def paired_difference(y1: Sequence[float], y0: Sequence[float], groups: Sequence, *,
                      n_bootstrap: int = N_BOOTSTRAP, seed: int = 20260926) -> Dict[str, object]:
    """Mean of ``y1 - y0`` over copies with both, both estimands."""
    d = np.asarray(y1, dtype=float) - np.asarray(y0, dtype=float)
    return mean_estimates(d, groups, null=0.0, n_bootstrap=n_bootstrap, seed=seed)


def auc_estimate(labels: Sequence[int], scores: Sequence[float], groups: Sequence, *, null: float = 0.5,
                 n_bootstrap: int = N_BOOTSTRAP, seed: int = 20260926) -> Dict[str, object]:
    """ROC-AUC (higher score ranks higher) with a group bootstrap."""
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    keep = np.isfinite(scores)
    labels, scores, groups = labels[keep], scores[keep], np.asarray(groups)[keep]
    out: Dict[str, object] = {"positives": int(labels.sum()), "negatives": int((1 - labels).sum()),
                              "groups": int(pd.Series(groups).nunique())}
    out["evidence"] = out["groups"] >= MIN_GROUPS
    if labels.sum() == 0 or labels.sum() == len(labels):
        return out
    ranked = protocol.Ranked(labels, scores)
    out["roc_auc"] = estimate(ranked.roc_auc, groups, null=null, n_bootstrap=n_bootstrap, seed=seed)
    return out


def label(estimate_: Optional[Dict[str, float]], *, good: str, bad: str, threshold: float = 0.5,
          bad_upper: Optional[float] = None) -> str:
    """Three-way label: ``good`` if low >= threshold, ``bad`` if high < (bad_upper or threshold)."""
    if not estimate_ or not np.isfinite(estimate_.get("low", np.nan)):
        return "not evaluable"
    if estimate_["low"] >= threshold:
        return good
    if estimate_["high"] < (threshold if bad_upper is None else bad_upper):
        return bad
    return "inconclusive"
