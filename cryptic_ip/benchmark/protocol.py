"""The evaluation protocol of docs/ANALYSIS_PLAN.md, with its leak guards.

Each guard is enforced in code rather than by convention:

* every split is checked for a group on both sides (:func:`assert_disjoint`);
* descriptors that can carry the ligand back in (crystallographic B-factors,
  exposed as ``plddt_*``) are never offered to a model (:data:`EXCLUDED_FEATURES`);
* the model family is chosen **inside** the inner loop, jointly with its
  hyperparameters, so no outer-fold score is used to pick what is reported;
* the decision threshold and the probability calibrator are fitted on the
  chosen configuration's **inner** out-of-fold predictions and applied
  unchanged to the outer test fold;
* the temporal holdout is group-disjoint from development by construction and
  checked again;
* uncertainty comes from resampling whole homology groups, never pockets.

The two arms of a paired comparison draw the same folds and the same candidate
configurations from the same seeds, so they differ only in their descriptors.
"""

from __future__ import annotations

import logging
import math
import warnings
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from cryptic_ip.analysis.features import FEATURE_NAMES

LOGGER = logging.getLogger(__name__)

#: Descriptors never offered to a benchmark model. On a PDB entry the
#: ``plddt_*`` columns are crystallographic B-factors: residues ordered by the
#: bound ligand have low B-factors, so the removed ligand leaks back in - and on
#: an AlphaFold model the same columns mean pLDDT. APBS is not run on the
#: benchmark, so ``electrostatic_potential`` is empty.
EXCLUDED_FEATURES: Tuple[str, ...] = (
    "plddt_mean",
    "plddt_min",
    "plddt_fraction_above_cutoff",
    "electrostatic_potential",
)
BENCHMARK_FEATURES: Tuple[str, ...] = tuple(f for f in FEATURE_NAMES if f not in EXCLUDED_FEATURES)

#: The descriptor sets compared by the primary hypotheses.
ARMS: Dict[str, Tuple[str, ...]] = {
    "full": BENCHMARK_FEATURES,
    "no_hull_depth": tuple(f for f in BENCHMARK_FEATURES if f != "hull_depth"),
}

#: A linear floor, a bagged and a boosted ensemble.
FAMILIES: Tuple[str, ...] = ("logistic_regression", "extra_trees", "hist_gradient_boosting")


class LeakError(AssertionError):
    """A split, holdout or feature set would let information cross into the test data."""


def assert_disjoint(train_groups: Sequence[Any], test_groups: Sequence[Any], what: str = "split") -> None:
    """Raise :class:`LeakError` if any group is on both sides."""
    shared = set(map(str, train_groups)) & set(map(str, test_groups))
    if shared:
        raise LeakError(f"{what}: {len(shared)} group(s) on both sides, e.g. {sorted(shared)[:5]}")


def assert_allowed_features(features: Sequence[str]) -> None:
    """Raise :class:`LeakError` if an excluded descriptor would reach a model."""
    leaked = [f for f in features if f in EXCLUDED_FEATURES]
    if leaked:
        raise LeakError(f"excluded descriptors offered to a model: {leaked}")
    unknown = [f for f in features if f not in FEATURE_NAMES]
    if unknown:
        raise LeakError(f"non-descriptor columns offered to a model: {unknown}")


# --------------------------------------------------------------- candidates
@dataclass(frozen=True)
class Candidate:
    """One model family with one hyperparameter setting."""

    family: str
    params: Tuple[Tuple[str, Any], ...]

    def as_dict(self) -> Dict[str, Any]:
        return {"family": self.family, **{k.replace("classifier__", ""): v for k, v in self.params}}


def benchmark_specs() -> Dict[str, Any]:
    """Model specifications for the benchmark families, with the plan's search spaces.

    Extra trees are held to 300 trees: the benchmark fits each candidate many
    times, and more trees change the ranking of pockets very little.
    """
    from sklearn.ensemble import ExtraTreesClassifier

    from cryptic_ip.analysis.ml_classifier import ModelSpec, default_model_specs

    specs = {spec.name: spec for spec in default_model_specs(include_xgboost=False)}
    specs["extra_trees"] = ModelSpec(
        name="extra_trees",
        build=lambda seed: ExtraTreesClassifier(
            n_estimators=300, random_state=seed, n_jobs=-1, class_weight="balanced_subsample"
        ),
        param_distributions={
            "classifier__max_depth": [None, 12, 20],
            "classifier__min_samples_leaf": [1, 2, 4, 8],
            "classifier__max_features": ["sqrt", 0.5],
        },
    )
    return {name: specs[name] for name in FAMILIES}


def draw_candidates(n_draws: int, seed: int, specs: Optional[Mapping[str, Any]] = None) -> List[Candidate]:
    """``n_draws`` hyperparameter settings per family, reproducibly from ``seed``."""
    from sklearn.model_selection import ParameterSampler

    specs = specs if specs is not None else benchmark_specs()
    out: List[Candidate] = []
    for family, spec in specs.items():
        space = dict(spec.param_distributions)
        n_space = int(np.prod([len(v) for v in space.values()])) if space else 1
        draws = ParameterSampler(space, n_iter=min(n_draws, n_space), random_state=seed) if space else [{}]
        for params in draws:
            out.append(Candidate(family, tuple(sorted(params.items()))))
    return out


def fit_candidate(candidate: Candidate, X: pd.DataFrame, y: np.ndarray, seed: int, specs: Mapping[str, Any]):
    """Fit one candidate on ``(X, y)``."""
    from cryptic_ip.analysis.ml_classifier import build_pipeline

    pipeline = build_pipeline(specs[candidate.family], seed)
    pipeline.set_params(**dict(candidate.params))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pipeline.fit(X, y)
    return pipeline


def _scores(model, X: pd.DataFrame) -> np.ndarray:
    return model.predict_proba(X)[:, 1]


class SplitError(ValueError):
    """No grouped split gives every training fold both classes."""


#: Seeded fold assignments tried at each fold count before trying fewer folds.
SPLIT_ATTEMPTS = 50


def _group_splits(y: np.ndarray, groups: np.ndarray, n_splits: int, seed: int) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Grouped, stratified folds in which every training fold holds both classes.

    With few positive groups - or one large group holding most positives - a
    single shuffled assignment can leave a training fold without positives.
    Seeded reassignments are tried, then fewer folds. The result depends only
    on the labels, the groups and the seed, so the paired arms of a comparison
    receive identical folds.
    """
    from sklearn.model_selection import StratifiedGroupKFold

    positive_groups = len(np.unique(groups[y == 1]))
    negative_groups = len(np.unique(groups[y == 0]))
    k_max = min(n_splits, positive_groups, negative_groups)
    if k_max < 2:
        raise SplitError(f"cannot split: {positive_groups} positive and {negative_groups} negative groups")
    for k in range(k_max, 1, -1):
        for attempt in range(SPLIT_ATTEMPTS):
            splitter = StratifiedGroupKFold(n_splits=k, shuffle=True, random_state=seed + 7919 * attempt)
            splits = list(splitter.split(np.zeros(len(y)), y, groups=groups))
            if all(len(np.unique(y[train])) == 2 for train, _ in splits):
                for train, test in splits:
                    assert_disjoint(groups[train], groups[test])
                if k < n_splits:
                    LOGGER.warning("using %d folds instead of %d: no valid %d-fold split", k, n_splits, k + 1)
                return splits
    raise SplitError(
        f"no grouped split of {positive_groups} positive groups gives every training fold both classes"
    )


def _average_precision(y: np.ndarray, s: np.ndarray, w: Optional[np.ndarray] = None) -> float:
    from sklearn.metrics import average_precision_score

    if len(np.unique(y)) < 2:
        return float("nan")
    return float(average_precision_score(y, s, sample_weight=w))


@dataclass
class Selection:
    """The inner loop's choice for one training set."""

    candidate: Candidate
    inner_ap: float
    inner_oof: np.ndarray
    threshold: float
    calibrator: Optional[Any]
    all_inner_ap: Dict[str, float] = field(default_factory=dict)


def select_candidate(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    candidates: Sequence[Candidate],
    *,
    n_inner: int,
    seed: int,
    specs: Mapping[str, Any],
) -> Selection:
    """Choose family and hyperparameters jointly by pooled inner out-of-fold AP.

    The chosen candidate's inner out-of-fold scores then fix the MCC-optimal
    threshold and a Platt calibrator for this training set.
    """
    from sklearn.linear_model import LogisticRegression

    from cryptic_ip.analysis.ml_classifier import select_threshold

    splits = _group_splits(y, groups, n_inner, seed)
    best: Optional[Tuple[float, int, np.ndarray]] = None
    scores: Dict[str, float] = {}
    for index, candidate in enumerate(candidates):
        oof = np.full(len(y), np.nan)
        for train, valid in splits:
            model = fit_candidate(candidate, X.iloc[train], y[train], seed, specs)
            oof[valid] = _scores(model, X.iloc[valid])
        ap = _average_precision(y, oof)
        scores[str(candidate.as_dict())] = ap
        # Ties go to the earlier candidate, so the choice is deterministic.
        if np.isfinite(ap) and (best is None or ap > best[0]):
            best = (ap, index, oof)
    if best is None:
        raise ValueError("no candidate produced a finite inner average precision")
    ap, index, oof = best
    threshold = float(select_threshold(y, oof, objective="mcc"))
    calibrator = None
    if len(np.unique(y)) == 2:
        calibrator = LogisticRegression(C=1e6, max_iter=1000).fit(oof.reshape(-1, 1), y)
    return Selection(candidates[index], float(ap), oof, threshold, calibrator, scores)


def _calibrate(selection: Selection, scores: np.ndarray) -> np.ndarray:
    if selection.calibrator is None:
        return np.full(len(scores), np.nan)
    return selection.calibrator.predict_proba(np.asarray(scores).reshape(-1, 1))[:, 1]


def _rule_threshold(y: np.ndarray, rule: np.ndarray) -> float:
    from cryptic_ip.analysis.ml_classifier import select_threshold

    return float(select_threshold(y, rule, objective="mcc"))


def run_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    rule_scores: np.ndarray,
    *,
    n_outer: int = 5,
    n_inner: int = 3,
    n_draws: int = 10,
    fold_seed: int = 0,
    model_seed: int = 0,
    candidate_seed: int = 7,
    specs: Optional[Mapping[str, Any]] = None,
    progress: Optional[Callable[[str], None]] = None,
) -> pd.DataFrame:
    """Nested grouped cross-validation; one out-of-fold prediction per row.

    Returns:
        A frame aligned with ``X``: ``fold``, ``score`` (uncalibrated),
        ``calibrated``, ``threshold``, ``family``, ``candidate``,
        ``rule_threshold``.
    """
    assert_allowed_features(list(X.columns))
    specs = specs if specs is not None else benchmark_specs()
    candidates = draw_candidates(n_draws, candidate_seed, specs)
    out = pd.DataFrame(
        {
            "fold": -1,
            "score": np.nan,
            "calibrated": np.nan,
            "threshold": np.nan,
            "family": "",
            "candidate": "",
            "rule_threshold": np.nan,
        },
        index=X.index,
    )
    for fold, (train, test) in enumerate(_group_splits(y, groups, n_outer, fold_seed)):
        selection = select_candidate(
            X.iloc[train], y[train], groups[train], candidates,
            n_inner=n_inner, seed=fold_seed * 1000 + fold + 1, specs=specs,
        )
        model = fit_candidate(selection.candidate, X.iloc[train], y[train], model_seed, specs)
        scores = _scores(model, X.iloc[test])
        rows = out.index[test]
        out.loc[rows, "fold"] = fold
        out.loc[rows, "score"] = scores
        out.loc[rows, "calibrated"] = _calibrate(selection, scores)
        out.loc[rows, "threshold"] = selection.threshold
        out.loc[rows, "family"] = selection.candidate.family
        out.loc[rows, "candidate"] = str(selection.candidate.as_dict())
        out.loc[rows, "rule_threshold"] = _rule_threshold(y[train], rule_scores[train])
        if progress:
            progress(
                f"fold {fold}: {selection.candidate.as_dict()} inner AP {selection.inner_ap:.3f}"
            )
    if (out["fold"] < 0).any():
        raise RuntimeError("some rows received no out-of-fold prediction")
    return out


def run_locked(
    X_dev: pd.DataFrame,
    y_dev: np.ndarray,
    groups_dev: np.ndarray,
    X_hold: pd.DataFrame,
    groups_hold: np.ndarray,
    rule_dev: np.ndarray,
    *,
    n_inner: int = 3,
    n_draws: int = 10,
    fold_seed: int = 0,
    model_seed: int = 0,
    candidate_seed: int = 7,
    specs: Optional[Mapping[str, Any]] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Fit the locked model on the development set and score the holdout once."""
    assert_allowed_features(list(X_dev.columns))
    assert_disjoint(groups_dev, groups_hold, "holdout")
    specs = specs if specs is not None else benchmark_specs()
    candidates = draw_candidates(n_draws, candidate_seed, specs)
    selection = select_candidate(
        X_dev, y_dev, groups_dev, candidates, n_inner=n_inner, seed=fold_seed * 1000 + 999, specs=specs
    )
    model = fit_candidate(selection.candidate, X_dev, y_dev, model_seed, specs)
    scores = _scores(model, X_hold)
    frame = pd.DataFrame(
        {
            "fold": -1,
            "score": scores,
            "calibrated": _calibrate(selection, scores),
            "threshold": selection.threshold,
            "family": selection.candidate.family,
            "candidate": str(selection.candidate.as_dict()),
            "rule_threshold": _rule_threshold(y_dev, rule_dev),
        },
        index=X_hold.index,
    )
    return frame, {"candidate": selection.candidate.as_dict(), "inner_ap": selection.inner_ap}


def temporal_holdout(
    entries: pd.DataFrame,
    *,
    group_column: str,
    date_column: str = "release_date",
    fraction: float = 0.20,
) -> set:
    """Entries in the latest ``fraction`` of groups, by each group's earliest release.

    Chosen from dates and groups alone - no label, no model. Ties in date are
    broken by group name, so the choice is deterministic.
    """
    table = entries[["pdb_id", group_column, date_column]].copy()
    table[date_column] = pd.to_datetime(table[date_column], errors="coerce")
    if table[date_column].isna().any():
        missing = table.loc[table[date_column].isna(), "pdb_id"].tolist()
        raise ValueError(f"release date missing for {len(missing)} entries, e.g. {missing[:5]}")
    first = table.groupby(group_column)[date_column].min().reset_index()
    first = first.sort_values([date_column, group_column]).reset_index(drop=True)
    n_hold = int(math.ceil(fraction * len(first)))
    held = set(first[group_column].iloc[len(first) - n_hold:])
    chosen = set(table.loc[table[group_column].isin(held), "pdb_id"].str.upper())
    assert_disjoint(
        table.loc[~table[group_column].isin(held), group_column],
        table.loc[table[group_column].isin(held), group_column],
        "holdout",
    )
    return chosen


# ------------------------------------------------------------ inference
class Ranked:
    """Scores ranked once, so weighted ROC-AUC and average precision cost O(n).

    A group bootstrap recomputes each metric thousands of times on tens of
    thousands of pockets; re-sorting every time dominates. Tied scores form one
    block, handled as scikit-learn does (a tie counts one half in ROC-AUC; average
    precision steps once per distinct threshold), and results match
    ``roc_auc_score`` / ``average_precision_score`` with ``sample_weight``.
    """

    def __init__(self, y: Sequence[int], scores: Sequence[float]) -> None:
        y = np.asarray(y, dtype=int)
        scores = np.asarray(scores, dtype=float)
        if not np.isfinite(scores).all():
            raise ValueError("scores must be finite")
        _, self._block = np.unique(-scores, return_inverse=True)  # block 0 = highest score
        self._n_blocks = int(self._block.max()) + 1 if len(scores) else 0
        self._pos = (y == 1).astype(float)
        self._neg = 1.0 - self._pos

    def _sums(self, w: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        tp = np.bincount(self._block, weights=w * self._pos, minlength=self._n_blocks)
        fp = np.bincount(self._block, weights=w * self._neg, minlength=self._n_blocks)
        return tp, fp

    def roc_auc(self, w: np.ndarray) -> float:
        tp, fp = self._sums(np.asarray(w, dtype=float))
        p, n = tp.sum(), fp.sum()
        if p <= 0 or n <= 0:
            return float("nan")
        negatives_below = n - np.cumsum(fp)
        return float(np.sum(tp * (negatives_below + 0.5 * fp)) / (p * n))

    def average_precision(self, w: np.ndarray) -> float:
        tp, fp = self._sums(np.asarray(w, dtype=float))
        p = tp.sum()
        if p <= 0 or fp.sum() <= 0:
            return float("nan")
        ctp, cfp = np.cumsum(tp), np.cumsum(fp)
        precision = np.divide(ctp, ctp + cfp, out=np.zeros_like(ctp), where=(ctp + cfp) > 0)
        return float(np.sum(tp / p * precision))


METRICS: Dict[str, Callable[["Ranked", np.ndarray], float]] = {
    "roc_auc": lambda ranked, w: ranked.roc_auc(w),
    "pr_auc": lambda ranked, w: ranked.average_precision(w),
}


def group_bootstrap_weights(groups: np.ndarray, n: int, seed: int):
    """Yield ``n`` row-weight vectors, each row weighted by how often its group was drawn.

    Groups are drawn with replacement, as many as there are groups. Generated
    one resample at a time: a dense ``(n, rows)`` matrix would not fit in memory
    for a benchmark of tens of thousands of pockets.
    """
    codes, uniques = pd.factorize(pd.Series(np.asarray(groups)).astype(str))
    rng = np.random.default_rng(seed)
    for _ in range(n):
        counts = np.bincount(rng.integers(0, len(uniques), size=len(uniques)), minlength=len(uniques))
        yield counts[codes].astype(float)


@dataclass
class Estimate:
    """A point estimate with its group-bootstrap interval and two-sided p-value against ``null``."""

    point: float
    low: float
    high: float
    p_value: float
    n_bootstrap: int

    def as_dict(self) -> Dict[str, float]:
        return {"point": self.point, "low": self.low, "high": self.high, "p_value": self.p_value}


def bootstrap_statistic(
    statistic: Callable[[np.ndarray], float],
    groups: np.ndarray,
    *,
    null: float = 0.0,
    n_bootstrap: int = 2000,
    seed: int = 20260923,
) -> Estimate:
    """Percentile interval and bootstrap p-value for ``statistic(weights)``.

    ``statistic`` receives per-row weights (all ones for the point estimate),
    so paired quantities - two arms scored on the same rows - are resampled
    together.
    """
    point = statistic(np.ones(len(groups)))
    weights = group_bootstrap_weights(groups, n_bootstrap, seed)
    values = np.array([statistic(w) for w in weights])
    values = values[np.isfinite(values)]
    if values.size == 0:
        return Estimate(point, float("nan"), float("nan"), float("nan"), 0)
    low, high = np.percentile(values, [2.5, 97.5])
    below = float(np.mean(values <= null))
    above = float(np.mean(values >= null))
    p = min(1.0, 2.0 * min(below, above))
    return Estimate(float(point), float(low), float(high), p, int(values.size))


def holm(p_values: Mapping[str, float]) -> Dict[str, float]:
    """Holm-adjusted p-values."""
    ordered = sorted(p_values.items(), key=lambda kv: kv[1])
    m = len(ordered)
    adjusted: Dict[str, float] = {}
    running = 0.0
    for rank, (name, p) in enumerate(ordered):
        running = max(running, min(1.0, (m - rank) * p))
        adjusted[name] = running
    return adjusted


def decide(
    development: Mapping[str, Estimate],
    adjusted_p: float,
    holdout: Optional[Estimate],
    holdout_positive_groups: int,
    *,
    equivalence_margin: float = 0.01,
    min_holdout_groups: int = 10,
) -> str:
    """The plan's decision rule (section 7) for one hypothesis.

    Args:
        development: Paired ROC-AUC difference under each grouping
            (``"sequence"`` and ``"strict"``).
        adjusted_p: Holm-adjusted p-value under the sequence grouping.
        holdout: Holdout difference, or ``None``.
        holdout_positive_groups: Positive groups in the holdout.
    """
    groupings = [development[g] for g in ("sequence", "strict") if g in development]
    if len(groupings) < 2 or any(not np.isfinite(e.low) for e in groupings):
        return "inconclusive"
    holdout_powered = holdout is not None and holdout_positive_groups >= min_holdout_groups
    if (
        all(e.low > 0 for e in groupings)
        and adjusted_p < 0.05
        and (not holdout_powered or (holdout.point > 0 and holdout.low > 0))
    ):
        return "supported"
    if all(e.high < 0 for e in groupings) or all(
        -equivalence_margin < e.low and e.high < equivalence_margin for e in groupings
    ):
        return "refuted"
    return "inconclusive"
