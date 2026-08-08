"""Machine-learning models for cryptic IP-binding pocket classification.

Evaluation design
-----------------
The previous implementation reported performance that could not be trusted, for
three separate reasons. All three are fixed here.

**1. Selection bias.** Hyperparameters were chosen by ``GridSearchCV`` over a
cross-validation split, and then the *same* split was reused to produce the
"validation" probabilities. Scoring a model on the folds that selected it is
optimistic by construction. Here, hyperparameter search runs in an **inner** loop
nested inside an **outer** loop that the search never sees, so every reported
number comes from data untouched by model selection.

**2. Leakage across related structures.** Splits were stratified but not grouped.
Pockets from the same PDB entry - and entries of the same protein - landed on both
sides of a split, so a model could memorise a protein rather than learn a site.
Every split here is grouped, and the group key is caller-supplied so it can be a
sequence cluster rather than merely an entry identifier.

**3. Miscalibrated probabilities.** Raw tree-ensemble scores are not probabilities;
using them at a 0.5 threshold, or comparing them to a rule-based score, is
meaningless. Each outer-fold model is calibrated on a held-out, group-disjoint
slice of its own training data, and calibration quality is reported (Brier score
and expected calibration error) alongside discrimination.

Two further points of rigour:

* **Class imbalance is handled by weighting, not resampling.** Oversampling
  minority pockets across a grouped split duplicates rows from the same protein
  on both sides of the calibration boundary; class weights avoid that entirely.
* **Screening-relevant metrics are reported.** A proteome screen inspects the top
  fraction of a ranked list, so enrichment factor and precision at *k* describe
  utility far better than accuracy at an arbitrary threshold.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .features import FEATURE_NAMES, LEGACY_FEATURE_NAMES

LOGGER = logging.getLogger(__name__)

try:
    import joblib
    from sklearn.base import BaseEstimator, clone
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.ensemble import (
        ExtraTreesClassifier,
        HistGradientBoostingClassifier,
        RandomForestClassifier,
    )
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        average_precision_score,
        balanced_accuracy_score,
        brier_score_loss,
        f1_score,
        matthews_corrcoef,
        precision_recall_curve,
        precision_score,
        recall_score,
        roc_auc_score,
        roc_curve,
    )
    from sklearn.model_selection import (
        GroupShuffleSplit,
        RandomizedSearchCV,
        StratifiedGroupKFold,
    )
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    SKLEARN_AVAILABLE = True
except ImportError:  # pragma: no cover - import guard for optional dependency
    joblib = None
    BaseEstimator = Any
    SKLEARN_AVAILABLE = False


def _ensure_sklearn_installed() -> None:
    """Raise a helpful error when scikit-learn is unavailable."""
    if not SKLEARN_AVAILABLE:
        raise ImportError(
            "scikit-learn and joblib are required for ML classification. "
            "Install dependencies with `pip install scikit-learn joblib`."
        )


#: Backwards-compatible six-feature schema. Retained so models serialised by
#: earlier versions load and score, and so existing imports keep working.
FEATURE_COLUMNS: Tuple[str, ...] = LEGACY_FEATURE_NAMES

#: Full descriptor schema used by newly trained models.
DEFAULT_FEATURE_COLUMNS: Tuple[str, ...] = FEATURE_NAMES

#: Fractions of a ranked list at which screening metrics are reported.
TOP_FRACTIONS: Tuple[float, ...] = (0.01, 0.05, 0.10)


# --------------------------------------------------------------------- metrics


def expected_calibration_error(
    y_true: np.ndarray, y_prob: np.ndarray, *, n_bins: int = 10
) -> float:
    """Expected calibration error over equal-width probability bins.

    ECE is the mean absolute gap between predicted probability and observed
    frequency, weighted by bin population. A model can discriminate perfectly
    (AUROC 1.0) and still be badly calibrated, which matters here because
    downstream filters apply probability thresholds.

    Args:
        y_true: Binary labels.
        y_prob: Predicted probabilities.
        n_bins: Number of equal-width bins.

    Returns:
        ECE in ``[0, 1]``; ``nan`` for empty input.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_prob = np.asarray(y_prob, dtype=float)
    if y_true.size == 0:
        return float("nan")
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    total = 0.0
    for lower, upper in zip(edges[:-1], edges[1:]):
        mask = (y_prob > lower) & (y_prob <= upper) if lower > 0 else (y_prob <= upper)
        if not np.any(mask):
            continue
        total += np.mean(mask) * abs(float(np.mean(y_true[mask]) - np.mean(y_prob[mask])))
    return float(total)


def enrichment_factor(
    y_true: np.ndarray, y_score: np.ndarray, *, fraction: float = 0.01
) -> float:
    """Enrichment of positives in the top ``fraction`` of a ranked list.

    An enrichment factor of 10 at 1 % means the top 1 % of ranked pockets contains
    ten times the positive density of the full set. This is the metric that
    predicts how much wet-lab effort a screen saves, which accuracy does not.

    Args:
        y_true: Binary labels.
        y_score: Ranking scores, higher is better.
        fraction: Top fraction to inspect.

    Returns:
        The enrichment factor, or ``nan`` when undefined.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    n = y_true.size
    if n == 0 or not 0 < fraction <= 1:
        return float("nan")
    baseline = float(np.mean(y_true))
    if baseline <= 0:
        return float("nan")
    k = max(1, int(np.ceil(fraction * n)))
    top = np.argsort(-y_score)[:k]
    return float(np.mean(y_true[top]) / baseline)


def precision_at_k(y_true: np.ndarray, y_score: np.ndarray, *, fraction: float = 0.01) -> float:
    """Precision within the top ``fraction`` of a ranked list.

    Args:
        y_true: Binary labels.
        y_score: Ranking scores, higher is better.
        fraction: Top fraction to inspect.

    Returns:
        Precision in ``[0, 1]``, or ``nan`` for empty input.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    n = y_true.size
    if n == 0 or not 0 < fraction <= 1:
        return float("nan")
    k = max(1, int(np.ceil(fraction * n)))
    top = np.argsort(-y_score)[:k]
    return float(np.mean(y_true[top]))


def classification_metrics(
    y_true: Sequence[int],
    y_prob: Sequence[float],
    *,
    threshold: float = 0.5,
) -> Dict[str, float]:
    """Compute the full metric panel for one set of predictions.

    Args:
        y_true: Binary labels.
        y_prob: Predicted probabilities.
        threshold: Decision threshold for the hard-label metrics.

    Returns:
        Mapping of metric name to value. Metrics that are undefined for the input
        (e.g. AUROC with a single class present) are ``nan`` rather than an error,
        so a degenerate fold does not abort a whole evaluation.
    """
    _ensure_sklearn_installed()
    y_true_arr = np.asarray(list(y_true), dtype=int)
    y_prob_arr = np.asarray(list(y_prob), dtype=float)
    predicted = (y_prob_arr >= float(threshold)).astype(int)

    both_classes = len(np.unique(y_true_arr)) == 2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        metrics: Dict[str, float] = {
            "n": float(y_true_arr.size),
            "n_positive": float(int(np.sum(y_true_arr == 1))),
            "prevalence": float(np.mean(y_true_arr)) if y_true_arr.size else float("nan"),
            "roc_auc": float(roc_auc_score(y_true_arr, y_prob_arr)) if both_classes else float("nan"),
            "pr_auc": (
                float(average_precision_score(y_true_arr, y_prob_arr))
                if both_classes
                else float("nan")
            ),
            "mcc": float(matthews_corrcoef(y_true_arr, predicted)) if both_classes else float("nan"),
            "f1": float(f1_score(y_true_arr, predicted, zero_division=0)),
            "precision": float(precision_score(y_true_arr, predicted, zero_division=0)),
            "recall": float(recall_score(y_true_arr, predicted, zero_division=0)),
            "balanced_accuracy": (
                float(balanced_accuracy_score(y_true_arr, predicted))
                if both_classes
                else float("nan")
            ),
            "brier": float(brier_score_loss(y_true_arr, y_prob_arr)) if y_true_arr.size else np.nan,
            "ece": expected_calibration_error(y_true_arr, y_prob_arr),
            "threshold": float(threshold),
        }
    for fraction in TOP_FRACTIONS:
        tag = f"{int(fraction * 100)}pct"
        metrics[f"enrichment_at_{tag}"] = enrichment_factor(
            y_true_arr, y_prob_arr, fraction=fraction
        )
        metrics[f"precision_at_{tag}"] = precision_at_k(y_true_arr, y_prob_arr, fraction=fraction)
    return metrics


def grouped_bootstrap_ci(
    y_true: Sequence[int],
    y_prob: Sequence[float],
    groups: Sequence[Any],
    *,
    metric: Callable[[np.ndarray, np.ndarray], float],
    n_bootstrap: int = 2000,
    alpha: float = 0.95,
    random_state: int = 42,
) -> Tuple[float, float, float]:
    """Bootstrap a metric by resampling **groups**, not rows.

    Resampling individual pockets treats several pockets from one protein as
    independent observations. They are not: they share a fold, a sequence and a
    structure, so row-level bootstrap understates uncertainty - often severely
    when a few large entries dominate. Resampling whole groups respects the
    dependence structure and yields honest intervals.

    Args:
        y_true: Binary labels.
        y_prob: Predicted probabilities.
        groups: Group identifier per row.
        metric: Callable ``(y_true, y_prob) -> float``.
        n_bootstrap: Number of resamples.
        alpha: Confidence level.
        random_state: Seed.

    Returns:
        ``(point_estimate, lower, upper)``; bounds are ``nan`` if no resample
        contained both classes.
    """
    y_true_arr = np.asarray(list(y_true), dtype=int)
    y_prob_arr = np.asarray(list(y_prob), dtype=float)
    group_arr = np.asarray(list(groups), dtype=object)

    point = float(metric(y_true_arr, y_prob_arr))
    unique_groups = np.unique(group_arr)
    if unique_groups.size < 2:
        return point, float("nan"), float("nan")

    index_by_group = {group: np.flatnonzero(group_arr == group) for group in unique_groups}
    rng = np.random.default_rng(random_state)
    values: List[float] = []
    for _ in range(int(n_bootstrap)):
        sampled = rng.choice(unique_groups, size=unique_groups.size, replace=True)
        indices = np.concatenate([index_by_group[group] for group in sampled])
        if len(np.unique(y_true_arr[indices])) < 2:
            continue
        try:
            values.append(float(metric(y_true_arr[indices], y_prob_arr[indices])))
        except ValueError:
            continue

    if not values:
        return point, float("nan"), float("nan")
    lower_q = (1.0 - alpha) / 2.0
    return point, float(np.quantile(values, lower_q)), float(np.quantile(values, 1.0 - lower_q))


def delong_roc_test(
    y_true: Sequence[int], scores_a: Sequence[float], scores_b: Sequence[float]
) -> Tuple[float, float]:
    """DeLong test for two correlated ROC curves on the same samples.

    Comparing two models by their bootstrap confidence intervals is the wrong
    test: the models score the *same* samples, so their errors are correlated and
    overlapping intervals do not imply the difference is insignificant. DeLong's
    method uses that correlation explicitly (DeLong et al., *Biometrics*
    44:837-845, 1988).

    Args:
        y_true: Binary labels.
        scores_a: Scores from the first model.
        scores_b: Scores from the second model.

    Returns:
        ``(auc_difference, two_sided_p_value)``.
    """
    from scipy import stats

    y = np.asarray(list(y_true), dtype=int)
    a = np.asarray(list(scores_a), dtype=float)
    b = np.asarray(list(scores_b), dtype=float)
    positives = y == 1
    negatives = ~positives
    m, n = int(np.sum(positives)), int(np.sum(negatives))
    if m == 0 or n == 0:
        return float("nan"), float("nan")

    def structural_components(scores: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
        pos, neg = scores[positives], scores[negatives]
        # Midrank kernel: 1 for a win, 0.5 for a tie.
        comparison = (pos[:, None] > neg[None, :]).astype(float)
        comparison += 0.5 * (pos[:, None] == neg[None, :])
        return comparison.mean(axis=1), comparison.mean(axis=0), float(comparison.mean())

    v10_a, v01_a, auc_a = structural_components(a)
    v10_b, v01_b, auc_b = structural_components(b)

    s10 = np.cov(np.vstack((v10_a, v10_b)), ddof=1) if m > 1 else np.zeros((2, 2))
    s01 = np.cov(np.vstack((v01_a, v01_b)), ddof=1) if n > 1 else np.zeros((2, 2))
    covariance = s10 / m + s01 / n

    difference = auc_a - auc_b
    variance = float(covariance[0, 0] + covariance[1, 1] - 2.0 * covariance[0, 1])
    if variance <= 0:
        return float(difference), 1.0 if difference == 0 else 0.0
    z = difference / np.sqrt(variance)
    return float(difference), float(2.0 * (1.0 - stats.norm.cdf(abs(z))))


def select_threshold(
    y_true: Sequence[int],
    y_prob: Sequence[float],
    *,
    objective: str = "mcc",
    beta: float = 1.0,
    min_precision: Optional[float] = None,
) -> float:
    """Choose a decision threshold from out-of-fold predictions.

    The threshold is a fitted parameter, so choosing it on the data used to report
    performance would bias that report. Callers must pass **out-of-fold**
    predictions; the returned threshold is then applied unchanged to held-out data.

    Args:
        y_true: Binary labels.
        y_prob: Predicted probabilities.
        objective: ``"mcc"``, ``"f1"``, ``"fbeta"`` or ``"youden"``.
        beta: Beta for the ``"fbeta"`` objective; ``beta < 1`` favours precision.
        min_precision: When set, only thresholds achieving at least this precision
            are considered - useful when follow-up experiments are expensive.

    Returns:
        The selected threshold in ``(0, 1)``.
    """
    _ensure_sklearn_installed()
    y = np.asarray(list(y_true), dtype=int)
    p = np.asarray(list(y_prob), dtype=float)
    if y.size == 0 or len(np.unique(y)) < 2:
        return 0.5

    candidates = np.unique(np.clip(p, 1e-6, 1 - 1e-6))
    if candidates.size > 512:
        candidates = np.quantile(candidates, np.linspace(0.0, 1.0, 512))

    best_threshold, best_value = 0.5, -np.inf
    for threshold in candidates:
        predicted = (p >= threshold).astype(int)
        if min_precision is not None:
            precision = precision_score(y, predicted, zero_division=0)
            if precision < min_precision:
                continue
        if objective == "mcc":
            value = matthews_corrcoef(y, predicted)
        elif objective == "f1":
            value = f1_score(y, predicted, zero_division=0)
        elif objective == "fbeta":
            precision = precision_score(y, predicted, zero_division=0)
            recall = recall_score(y, predicted, zero_division=0)
            denominator = beta**2 * precision + recall
            value = (
                (1 + beta**2) * precision * recall / denominator if denominator > 0 else 0.0
            )
        elif objective == "youden":
            recall = recall_score(y, predicted, zero_division=0)
            specificity = recall_score(1 - y, 1 - predicted, zero_division=0)
            value = recall + specificity - 1.0
        else:
            raise ValueError(f"Unknown objective: {objective!r}")
        if value > best_value:
            best_threshold, best_value = float(threshold), float(value)
    return best_threshold


# ---------------------------------------------------------------- model zoo


@dataclass
class ModelSpec:
    """A candidate model: an estimator factory plus its search space.

    Attributes:
        name: Identifier used in reports and artefacts.
        build: Callable returning a fresh, unfitted estimator.
        param_distributions: Search space keyed by pipeline parameter name.
        needs_scaling: Whether the estimator requires standardised features.
        description: One-line rationale for including this model.
    """

    name: str
    build: Callable[[int], Any]
    param_distributions: Dict[str, Sequence[Any]] = field(default_factory=dict)
    needs_scaling: bool = False
    description: str = ""


def default_model_specs(*, include_xgboost: bool = True) -> List[ModelSpec]:
    """Return the candidate models compared during training.

    The set spans model families deliberately. A regularised linear model is
    included as an honest floor: if it matches the ensembles, the extra capacity
    is not buying anything and should not be claimed. Tree ensembles capture the
    conjunctive structure of the hypothesis (buried **and** basic **and** the
    right size), which a linear model cannot express.

    Args:
        include_xgboost: Include XGBoost when the package is importable.

    Returns:
        Model specifications.
    """
    _ensure_sklearn_installed()
    specs: List[ModelSpec] = [
        ModelSpec(
            name="logistic_regression",
            build=lambda seed: LogisticRegression(
                max_iter=5000, class_weight="balanced", random_state=seed
            ),
            param_distributions={
                "classifier__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0],
                "classifier__penalty": ["l2"],
            },
            needs_scaling=True,
            description="Regularised linear baseline; the floor any ensemble must beat.",
        ),
        ModelSpec(
            name="random_forest",
            build=lambda seed: RandomForestClassifier(
                random_state=seed, n_jobs=-1, class_weight="balanced_subsample"
            ),
            param_distributions={
                "classifier__n_estimators": [300, 600, 1000],
                "classifier__max_depth": [None, 6, 10, 16],
                "classifier__min_samples_leaf": [1, 2, 4, 8],
                "classifier__max_features": ["sqrt", "log2", 0.5],
            },
            description="Bagged trees; robust to irrelevant features and mixed scales.",
        ),
        ModelSpec(
            name="extra_trees",
            build=lambda seed: ExtraTreesClassifier(
                random_state=seed, n_jobs=-1, class_weight="balanced_subsample"
            ),
            param_distributions={
                "classifier__n_estimators": [300, 600, 1000],
                "classifier__max_depth": [None, 8, 14],
                "classifier__min_samples_leaf": [1, 2, 4],
                "classifier__max_features": ["sqrt", 0.5],
            },
            description="Extremely randomised trees; lower variance than a random forest.",
        ),
        ModelSpec(
            name="hist_gradient_boosting",
            build=lambda seed: HistGradientBoostingClassifier(
                random_state=seed, class_weight="balanced"
            ),
            param_distributions={
                "classifier__max_iter": [200, 400, 800],
                "classifier__learning_rate": [0.02, 0.05, 0.1],
                "classifier__max_leaf_nodes": [15, 31, 63],
                "classifier__min_samples_leaf": [5, 10, 20],
                "classifier__l2_regularization": [0.0, 0.1, 1.0],
            },
            description="Gradient boosting with native missing-value handling.",
        ),
    ]

    if include_xgboost:
        try:
            from xgboost import XGBClassifier

            specs.append(
                ModelSpec(
                    name="xgboost",
                    build=lambda seed: XGBClassifier(
                        random_state=seed,
                        eval_metric="aucpr",
                        tree_method="hist",
                        n_jobs=-1,
                    ),
                    param_distributions={
                        "classifier__n_estimators": [300, 600],
                        "classifier__max_depth": [3, 5, 7],
                        "classifier__learning_rate": [0.02, 0.05, 0.1],
                        "classifier__subsample": [0.7, 0.9, 1.0],
                        "classifier__colsample_bytree": [0.6, 0.8, 1.0],
                        "classifier__min_child_weight": [1, 5, 10],
                        "classifier__reg_lambda": [1.0, 5.0],
                    },
                    description="Gradient boosting; usually the strongest single model here.",
                )
            )
        except ImportError:  # pragma: no cover - optional dependency
            LOGGER.info("xgboost unavailable; continuing without it")
    return specs


def build_pipeline(spec: ModelSpec, random_state: int) -> "Pipeline":
    """Assemble the preprocessing and estimator pipeline for a model spec.

    Missing values are median-imputed **with an indicator column**. The indicator
    matters: a missing APBS potential is informative (the calculation failed or
    was skipped), and silently imputing it away would discard that signal while
    pretending a value was measured.

    Args:
        spec: Model specification.
        random_state: Seed for the estimator.

    Returns:
        An unfitted pipeline.
    """
    _ensure_sklearn_installed()
    # keep_empty_features keeps the column geometry stable when a descriptor is
    # missing for every row (e.g. APBS never run), so a model trained with the
    # column present can still score data where it is absent.
    steps: List[Tuple[str, Any]] = [
        (
            "imputer",
            SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
        )
    ]
    if spec.needs_scaling:
        steps.append(("scaler", StandardScaler()))
    steps.append(("classifier", spec.build(random_state)))
    return Pipeline(steps=steps)


# --------------------------------------------------------- training results


@dataclass
class FoldResult:
    """Metrics and predictions from one outer cross-validation fold."""

    fold: int
    n_train: int
    n_test: int
    n_train_positive: int
    n_test_positive: int
    best_params: Dict[str, Any]
    metrics: Dict[str, float]
    test_indices: List[int] = field(default_factory=list)
    test_probabilities: List[float] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serialisable representation without raw predictions."""
        payload = asdict(self)
        payload.pop("test_indices", None)
        payload.pop("test_probabilities", None)
        return payload


@dataclass
class TrainingResults:
    """Nested cross-validation results for one model.

    Attributes:
        model_name: Model identifier.
        feature_names: Features the model consumes, in order.
        folds: Per-fold results from the outer loop.
        oof_probabilities: Out-of-fold probability for every training row.
        oof_labels: Labels aligned with :attr:`oof_probabilities`.
        oof_groups: Group identifier aligned with :attr:`oof_probabilities`.
        metrics: Pooled out-of-fold metrics at the selected threshold.
        metric_cis: Group-bootstrap confidence intervals for key metrics.
        selected_threshold: Threshold chosen on out-of-fold predictions.
        best_params: Hyperparameters of the final model refit on all data.
        n_positive: Positive rows in the training set.
        n_negative: Negative rows in the training set.
        n_groups: Distinct groups in the training set.
    """

    model_name: str
    feature_names: Tuple[str, ...]
    folds: List[FoldResult] = field(default_factory=list)
    oof_probabilities: np.ndarray = field(default_factory=lambda: np.empty(0))
    oof_labels: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    oof_groups: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=object))
    metrics: Dict[str, float] = field(default_factory=dict)
    metric_cis: Dict[str, Tuple[float, float, float]] = field(default_factory=dict)
    selected_threshold: float = 0.5
    best_params: Dict[str, Any] = field(default_factory=dict)
    n_positive: int = 0
    n_negative: int = 0
    n_groups: int = 0

    # Legacy attribute names kept so older reporting code keeps working.
    @property
    def roc_auc(self) -> float:
        """Pooled out-of-fold AUROC."""
        return float(self.metrics.get("roc_auc", float("nan")))

    @property
    def pr_auc(self) -> float:
        """Pooled out-of-fold average precision."""
        return float(self.metrics.get("pr_auc", float("nan")))

    @property
    def roc_auc_ci(self) -> Tuple[float, float]:
        """Confidence interval for the pooled AUROC."""
        _, lower, upper = self.metric_cis.get("roc_auc", (np.nan, np.nan, np.nan))
        return (float(lower), float(upper))

    @property
    def pr_auc_ci(self) -> Tuple[float, float]:
        """Confidence interval for the pooled average precision."""
        _, lower, upper = self.metric_cis.get("pr_auc", (np.nan, np.nan, np.nan))
        return (float(lower), float(upper))

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serialisable summary."""
        return {
            "model_name": self.model_name,
            "feature_names": list(self.feature_names),
            "n_positive": self.n_positive,
            "n_negative": self.n_negative,
            "n_groups": self.n_groups,
            "selected_threshold": self.selected_threshold,
            "best_params": {k: _jsonable(v) for k, v in self.best_params.items()},
            "metrics": self.metrics,
            "metric_cis": {k: list(v) for k, v in self.metric_cis.items()},
            "folds": [fold.to_dict() for fold in self.folds],
        }


def _jsonable(value: Any) -> Any:
    """Coerce a value into something ``json.dumps`` accepts."""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.ndarray, list, tuple)):
        return [_jsonable(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


# ------------------------------------------------------------- the classifier


class CrypticSiteMLClassifier:
    """Nested-cross-validated, calibrated classifier for cryptic IP-binding sites.

    Args:
        model_type: Model name from :func:`default_model_specs`.
        random_state: Seed controlling every stochastic step.
        n_splits: Outer cross-validation folds.
        inner_splits: Inner folds used for hyperparameter search.
        n_search_iter: Randomised-search budget per inner loop.
        feature_names: Feature schema the model consumes.
        calibrate: Fit a probability calibrator on a group-disjoint slice.
        calibration_fraction: Share of each training fold reserved for calibration.
        scoring: Inner-loop selection metric. Average precision is the default
            because the positive class is rare, and AUROC is insensitive to
            precision changes in that regime.
    """

    def __init__(
        self,
        model_type: str = "random_forest",
        random_state: int = 42,
        n_splits: int = 5,
        *,
        inner_splits: int = 3,
        n_search_iter: int = 40,
        feature_names: Sequence[str] = DEFAULT_FEATURE_COLUMNS,
        calibrate: bool = True,
        calibration_fraction: float = 0.25,
        scoring: str = "average_precision",
    ) -> None:
        _ensure_sklearn_installed()
        self.model_type = model_type
        self.random_state = int(random_state)
        self.n_splits = int(n_splits)
        self.inner_splits = int(inner_splits)
        self.n_search_iter = int(n_search_iter)
        self.feature_names: Tuple[str, ...] = tuple(feature_names)
        self.calibrate = bool(calibrate)
        self.calibration_fraction = float(calibration_fraction)
        self.scoring = scoring

        self.spec = self._resolve_spec(model_type)
        self.pipeline = build_pipeline(self.spec, self.random_state)
        self.best_estimator_: Optional[Any] = None
        self.results_: Optional[TrainingResults] = None
        self.decision_threshold_: float = 0.5

    @staticmethod
    def _resolve_spec(model_type: str) -> ModelSpec:
        """Look up a model specification by name."""
        for spec in default_model_specs():
            if spec.name == model_type:
                return spec
        available = ", ".join(spec.name for spec in default_model_specs())
        raise ValueError(f"Unknown model_type {model_type!r}. Available: {available}")

    # ------------------------------------------------------------ validation

    def _validate_inputs(
        self,
        features: pd.DataFrame,
        labels: Iterable[int],
        groups: Optional[Sequence[Any]] = None,
    ) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
        """Check the feature matrix, labels and groups, returning aligned arrays.

        Args:
            features: Feature table.
            labels: Binary labels.
            groups: Group identifiers; defaults to one group per row, which
                disables grouping and is only appropriate for genuinely
                independent samples.

        Returns:
            ``(X, y, groups)``.

        Raises:
            ValueError: On missing columns, non-binary labels, or too few
                examples per class.
        """
        missing = [name for name in self.feature_names if name not in features.columns]
        if missing:
            raise ValueError(
                f"Missing feature columns: {missing}. "
                f"Expected schema: {list(self.feature_names)}"
            )

        y = np.asarray(list(labels), dtype=int)
        if y.size != len(features):
            raise ValueError(
                f"Label count ({y.size}) does not match feature rows ({len(features)})"
            )
        unique = np.unique(y)
        if not np.array_equal(np.sort(unique), np.array([0, 1])):
            raise ValueError(f"Binary 0/1 labels are required; found {unique.tolist()}")

        counts = np.bincount(y, minlength=2)
        if int(np.min(counts)) < 2:
            raise ValueError(
                f"At least 2 examples per class are required; class counts are {counts.tolist()}"
            )

        if groups is None:
            LOGGER.warning(
                "No groups supplied; every row is treated as an independent group. "
                "Pass groups (PDB entry or, better, sequence cluster) to prevent "
                "leakage between related structures."
            )
            group_arr = np.arange(len(features), dtype=object)
        else:
            group_arr = np.asarray(list(groups), dtype=object)
            if group_arr.size != len(features):
                raise ValueError("groups must align with features")

        return features.loc[:, list(self.feature_names)].copy(), y, group_arr

    def _effective_splits(self, y: np.ndarray, groups: np.ndarray, requested: int) -> int:
        """Reduce the fold count when the data cannot support the request.

        A grouped split needs at least one group of each class per fold. When
        positives are concentrated in few groups, asking for more folds than
        positive groups produces folds with no positives and meaningless metrics.
        """
        positive_groups = len(np.unique(groups[y == 1]))
        negative_groups = len(np.unique(groups[y == 0]))
        return max(2, min(int(requested), positive_groups, negative_groups))

    # -------------------------------------------------------------- training

    def fit(
        self,
        features: pd.DataFrame,
        labels: Iterable[int],
        *,
        groups: Optional[Sequence[Any]] = None,
        class_weight: Optional[str] = None,
        n_bootstrap: int = 1000,
        threshold_objective: str = "mcc",
    ) -> TrainingResults:
        """Run nested cross-validation, then refit the final model on all data.

        The outer loop yields unbiased performance estimates; the final model is
        refit on everything (with hyperparameters re-searched on the full set) for
        deployment, which is standard practice and does not affect the reported
        numbers.

        Args:
            features: Feature table.
            labels: Binary labels.
            groups: Group identifiers used by every split.
            class_weight: Deprecated; class weighting is part of the model specs.
                Accepted for backwards compatibility and ignored with a warning.
            n_bootstrap: Group-bootstrap resamples for confidence intervals.
            threshold_objective: Objective passed to :func:`select_threshold`.

        Returns:
            The nested cross-validation results.
        """
        if class_weight is not None:
            LOGGER.warning(
                "class_weight=%r is ignored: class weighting is configured per model "
                "specification so that every candidate handles imbalance consistently.",
                class_weight,
            )

        X, y, group_arr = self._validate_inputs(features, labels, groups)
        outer_splits = self._effective_splits(y, group_arr, self.n_splits)
        outer = StratifiedGroupKFold(
            n_splits=outer_splits, shuffle=True, random_state=self.random_state
        )

        oof_probabilities = np.full(len(X), np.nan, dtype=float)
        folds: List[FoldResult] = []

        for fold_index, (train_idx, test_idx) in enumerate(outer.split(X, y, groups=group_arr)):
            X_train, y_train = X.iloc[train_idx], y[train_idx]
            X_test, y_test = X.iloc[test_idx], y[test_idx]
            groups_train = group_arr[train_idx]

            estimator, best_params = self._search_and_fit(X_train, y_train, groups_train)
            probabilities = estimator.predict_proba(X_test)[:, 1]
            oof_probabilities[test_idx] = probabilities

            folds.append(
                FoldResult(
                    fold=fold_index,
                    n_train=int(len(train_idx)),
                    n_test=int(len(test_idx)),
                    n_train_positive=int(np.sum(y_train == 1)),
                    n_test_positive=int(np.sum(y_test == 1)),
                    best_params={k: _jsonable(v) for k, v in best_params.items()},
                    metrics=classification_metrics(y_test, probabilities),
                    test_indices=[int(i) for i in test_idx],
                    test_probabilities=[float(p) for p in probabilities],
                )
            )

        evaluated = ~np.isnan(oof_probabilities)
        threshold = select_threshold(
            y[evaluated], oof_probabilities[evaluated], objective=threshold_objective
        )
        self.decision_threshold_ = threshold

        metrics = classification_metrics(
            y[evaluated], oof_probabilities[evaluated], threshold=threshold
        )
        metric_cis = {
            "roc_auc": grouped_bootstrap_ci(
                y[evaluated],
                oof_probabilities[evaluated],
                group_arr[evaluated],
                metric=lambda a, b: float(roc_auc_score(a, b)),
                n_bootstrap=n_bootstrap,
                random_state=self.random_state,
            ),
            "pr_auc": grouped_bootstrap_ci(
                y[evaluated],
                oof_probabilities[evaluated],
                group_arr[evaluated],
                metric=lambda a, b: float(average_precision_score(a, b)),
                n_bootstrap=n_bootstrap,
                random_state=self.random_state,
            ),
            "mcc": grouped_bootstrap_ci(
                y[evaluated],
                oof_probabilities[evaluated],
                group_arr[evaluated],
                metric=lambda a, b: float(matthews_corrcoef(a, (b >= threshold).astype(int))),
                n_bootstrap=n_bootstrap,
                random_state=self.random_state,
            ),
        }

        # Final deployment model: search once more on the full data set.
        self.best_estimator_, final_params = self._search_and_fit(X, y, group_arr)

        results = TrainingResults(
            model_name=self.model_type,
            feature_names=self.feature_names,
            folds=folds,
            oof_probabilities=oof_probabilities[evaluated],
            oof_labels=y[evaluated],
            oof_groups=group_arr[evaluated],
            metrics=metrics,
            metric_cis=metric_cis,
            selected_threshold=threshold,
            best_params={k: _jsonable(v) for k, v in final_params.items()},
            n_positive=int(np.sum(y == 1)),
            n_negative=int(np.sum(y == 0)),
            n_groups=int(len(np.unique(group_arr))),
        )
        self.results_ = results
        return results

    def _search_and_fit(
        self, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray
    ) -> Tuple[Any, Dict[str, Any]]:
        """Search hyperparameters with grouped inner CV, then fit and calibrate.

        Args:
            X: Feature table for this training fold.
            y: Labels.
            groups: Group identifiers.

        Returns:
            ``(fitted_estimator, best_params)``.
        """
        inner_splits = self._effective_splits(y, groups, self.inner_splits)
        inner = StratifiedGroupKFold(
            n_splits=inner_splits, shuffle=True, random_state=self.random_state
        )
        pipeline = build_pipeline(self.spec, self.random_state)

        if self.spec.param_distributions:
            search = RandomizedSearchCV(
                estimator=pipeline,
                param_distributions=dict(self.spec.param_distributions),
                n_iter=self.n_search_iter,
                scoring=self.scoring,
                cv=inner.split(X, y, groups=groups),
                random_state=self.random_state,
                refit=False,
                n_jobs=1,
                error_score=np.nan,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                search.fit(X, y, groups=groups)
            best_params = dict(search.best_params_)
            pipeline = build_pipeline(self.spec, self.random_state)
            pipeline.set_params(**best_params)
        else:
            best_params = {}

        estimator = self._fit_with_calibration(pipeline, X, y, groups)
        return estimator, best_params

    def _fit_with_calibration(
        self, pipeline: Any, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray
    ) -> Any:
        """Fit the pipeline, calibrating on a group-disjoint held-out slice.

        Calibrating on data the model was fitted on would map an overfitted,
        over-confident score onto an equally over-confident probability. Holding
        out whole groups keeps the calibration set independent in the same sense
        the evaluation is.

        Args:
            pipeline: Configured, unfitted pipeline.
            X: Feature table.
            y: Labels.
            groups: Group identifiers.

        Returns:
            A fitted estimator exposing ``predict_proba``.
        """
        if not self.calibrate:
            pipeline.fit(X, y)
            return pipeline

        n_groups = len(np.unique(groups))
        if n_groups < 4 or int(np.sum(y == 1)) < 8:
            # Too little data to spare a calibration split without destabilising
            # the fit; an uncalibrated model is preferable to a calibrator fitted
            # on a handful of points.
            LOGGER.info(
                "Skipping calibration: %d groups and %d positives are too few.",
                n_groups,
                int(np.sum(y == 1)),
            )
            pipeline.fit(X, y)
            return pipeline

        splitter = GroupShuffleSplit(
            n_splits=1, test_size=self.calibration_fraction, random_state=self.random_state
        )
        fit_idx, calib_idx = next(splitter.split(X, y, groups=groups))
        if len(np.unique(y[calib_idx])) < 2 or len(np.unique(y[fit_idx])) < 2:
            LOGGER.info("Skipping calibration: a split lacked both classes.")
            pipeline.fit(X, y)
            return pipeline

        base = clone(pipeline)
        base.fit(X.iloc[fit_idx], y[fit_idx])
        # Isotonic regression needs a few hundred points to beat Platt scaling;
        # below that it overfits the calibration set.
        method = "isotonic" if int(np.sum(y[calib_idx] == 1)) >= 50 else "sigmoid"
        calibrated = CalibratedClassifierCV(base, method=method, cv="prefit")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            calibrated.fit(X.iloc[calib_idx], y[calib_idx])
        return calibrated

    # ------------------------------------------------------------- inference

    def predict_proba(self, features: pd.DataFrame) -> np.ndarray:
        """Predict cryptic-site probabilities.

        Args:
            features: Feature table containing at least the model's schema.

        Returns:
            Probability of the positive class per row.

        Raises:
            RuntimeError: If the model has not been fitted or loaded.
        """
        if self.best_estimator_ is None:
            raise RuntimeError("Model is not trained. Run fit() first.")
        missing = [name for name in self.feature_names if name not in features.columns]
        if missing:
            raise ValueError(f"Missing feature columns for inference: {missing}")
        X = features.loc[:, list(self.feature_names)].copy()
        return self.best_estimator_.predict_proba(X)[:, 1]

    def predict(self, features: pd.DataFrame) -> np.ndarray:
        """Predict hard labels using the fitted decision threshold."""
        return (self.predict_proba(features) >= self.decision_threshold_).astype(int)

    def permutation_importance_(
        self,
        features: pd.DataFrame,
        labels: Sequence[int],
        *,
        n_repeats: int = 20,
        scoring: str = "average_precision",
    ) -> pd.DataFrame:
        """Permutation importance of each feature.

        Permutation importance is reported alongside SHAP because the two answer
        different questions: SHAP attributes individual predictions, while
        permutation importance measures how much a feature contributes to the
        metric that matters. Correlated features share credit in both, so neither
        should be read as a causal claim.

        Args:
            features: Feature table.
            labels: Binary labels.
            n_repeats: Permutations per feature.
            scoring: Metric whose degradation is measured.

        Returns:
            A frame of ``feature``, ``importance_mean`` and ``importance_std``,
            sorted by descending importance.
        """
        from sklearn.inspection import permutation_importance

        if self.best_estimator_ is None:
            raise RuntimeError("Model is not trained. Run fit() first.")
        X = features.loc[:, list(self.feature_names)].copy()
        y = np.asarray(list(labels), dtype=int)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = permutation_importance(
                self.best_estimator_,
                X,
                y,
                n_repeats=int(n_repeats),
                random_state=self.random_state,
                scoring=scoring,
                n_jobs=1,
            )
        return pd.DataFrame(
            {
                "feature": list(self.feature_names),
                "importance_mean": result.importances_mean,
                "importance_std": result.importances_std,
            }
        ).sort_values("importance_mean", ascending=False, ignore_index=True)

    def shap_values(self, features: pd.DataFrame) -> pd.DataFrame:
        """Mean absolute SHAP value per feature.

        Args:
            features: Feature table.

        Returns:
            A frame of ``feature`` and ``mean_abs_shap``, sorted descending.

        Raises:
            ImportError: If SHAP is not installed.
            RuntimeError: If the model has not been fitted.
        """
        if self.best_estimator_ is None:
            raise RuntimeError("Model is not trained. Run fit() first.")
        try:
            import shap
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError("SHAP is required for feature attribution.") from exc

        X = features.loc[:, list(self.feature_names)].copy()
        estimator = self.best_estimator_
        # Unwrap the calibration layer: SHAP explains the underlying model.
        if hasattr(estimator, "calibrated_classifiers_"):
            inner = estimator.calibrated_classifiers_[0]
            estimator = getattr(inner, "estimator", getattr(inner, "base_estimator", estimator))

        transformed = X
        model = estimator
        if hasattr(estimator, "named_steps"):
            transformed = estimator.named_steps["imputer"].transform(X)
            if "scaler" in estimator.named_steps:
                transformed = estimator.named_steps["scaler"].transform(transformed)
            model = estimator.named_steps["classifier"]

        try:
            explainer = shap.TreeExplainer(model)
        except Exception:  # noqa: BLE001 - non-tree models need the generic path
            explainer = shap.Explainer(model, transformed)
        values = explainer.shap_values(transformed) if hasattr(explainer, "shap_values") else (
            explainer(transformed).values
        )
        if isinstance(values, list):
            values = values[-1]
        values = np.asarray(values)
        if values.ndim == 3:
            # (n_samples, n_features, n_classes) -> positive class
            values = values[:, :, -1]

        importance = np.abs(values).mean(axis=0)
        # The imputer's indicator columns extend the matrix beyond the schema.
        names = list(self.feature_names)
        if importance.shape[0] > len(names):
            names = names + [
                f"missing_indicator_{i}" for i in range(importance.shape[0] - len(names))
            ]
        return pd.DataFrame(
            {"feature": names[: importance.shape[0]], "mean_abs_shap": importance}
        ).sort_values("mean_abs_shap", ascending=False, ignore_index=True)

    def compute_curves(
        self, features: pd.DataFrame, labels: Iterable[int]
    ) -> Dict[str, np.ndarray]:
        """Compute ROC and precision-recall curve coordinates.

        Prefers the stored out-of-fold predictions, because curves drawn from
        predictions on training data are optimistic. Falls back to scoring the
        supplied rows when no cross-validation results are stored.

        Args:
            features: Feature table.
            labels: Binary labels.

        Returns:
            Curve coordinates and areas.
        """
        _ensure_sklearn_installed()
        if self.results_ is not None and self.results_.oof_probabilities.size:
            y = self.results_.oof_labels
            probabilities = self.results_.oof_probabilities
            source = "out_of_fold"
        else:
            X, y, _ = self._validate_inputs(features, labels, None)
            probabilities = self.predict_proba(X)
            source = "in_sample"

        fpr, tpr, _ = roc_curve(y, probabilities)
        precision, recall, _ = precision_recall_curve(y, probabilities)
        return {
            "fpr": fpr,
            "tpr": tpr,
            "roc_auc": np.array([roc_auc_score(y, probabilities)]),
            "precision": precision,
            "recall": recall,
            "pr_auc": np.array([average_precision_score(y, probabilities)]),
            "source": np.array([source]),
        }

    # ------------------------------------------------------------ persistence

    def save(self, path: str) -> None:
        """Serialise the fitted model with its schema and threshold.

        Args:
            path: Destination file path.

        Raises:
            RuntimeError: If the model has not been fitted.
        """
        if self.best_estimator_ is None:
            raise RuntimeError("Model is not trained. Run fit() first.")
        payload = {
            "format_version": 2,
            "model_type": self.model_type,
            "random_state": self.random_state,
            "n_splits": self.n_splits,
            "feature_names": list(self.feature_names),
            "decision_threshold": float(self.decision_threshold_),
            "calibrated": bool(self.calibrate),
            "model": self.best_estimator_,
            "metrics": self.results_.metrics if self.results_ else {},
        }
        joblib.dump(payload, path)

    @classmethod
    def load(cls, path: str) -> "CrypticSiteMLClassifier":
        """Load a serialised model, including artefacts from format version 1.

        Args:
            path: Model file path.

        Returns:
            The loaded classifier.
        """
        _ensure_sklearn_installed()
        payload = joblib.load(path)
        feature_names = payload.get("feature_names") or list(FEATURE_COLUMNS)
        model_type = payload.get("model_type", "random_forest")
        try:
            classifier = cls(
                model_type=model_type,
                random_state=int(payload.get("random_state", 42)),
                n_splits=int(payload.get("n_splits", 5)),
                feature_names=feature_names,
            )
        except ValueError:
            # A model trained under an older, differently named zoo entry.
            classifier = cls(
                model_type="random_forest",
                random_state=int(payload.get("random_state", 42)),
                feature_names=feature_names,
            )
            classifier.model_type = model_type
        classifier.best_estimator_ = payload["model"]
        classifier.decision_threshold_ = float(payload.get("decision_threshold", 0.5))
        return classifier


# ---------------------------------------------------------------- comparison


def compare_models(
    features: pd.DataFrame,
    labels: Sequence[int],
    groups: Sequence[Any],
    *,
    model_names: Optional[Sequence[str]] = None,
    feature_names: Sequence[str] = DEFAULT_FEATURE_COLUMNS,
    random_state: int = 42,
    n_splits: int = 5,
    inner_splits: int = 3,
    n_search_iter: int = 30,
    n_bootstrap: int = 1000,
) -> Dict[str, TrainingResults]:
    """Train and nested-cross-validate several models on identical splits.

    Identical seeds and identical grouped splits are used for every candidate, so
    the comparison is paired and a DeLong test on the pooled out-of-fold
    predictions is valid.

    Args:
        features: Feature table.
        labels: Binary labels.
        groups: Group identifiers.
        model_names: Models to compare; defaults to the whole zoo.
        feature_names: Feature schema.
        random_state: Seed shared by every model.
        n_splits: Outer folds.
        inner_splits: Inner folds.
        n_search_iter: Randomised-search budget.
        n_bootstrap: Group-bootstrap resamples.

    Returns:
        Results keyed by model name.
    """
    _ensure_sklearn_installed()
    wanted = list(model_names) if model_names else [spec.name for spec in default_model_specs()]
    out: Dict[str, TrainingResults] = {}
    for name in wanted:
        LOGGER.info("Training %s", name)
        classifier = CrypticSiteMLClassifier(
            model_type=name,
            random_state=random_state,
            n_splits=n_splits,
            inner_splits=inner_splits,
            n_search_iter=n_search_iter,
            feature_names=feature_names,
        )
        try:
            out[name] = classifier.fit(
                features, labels, groups=groups, n_bootstrap=n_bootstrap
            )
            out[name].estimator = classifier  # type: ignore[attr-defined]
        except Exception as exc:  # noqa: BLE001 - one failed model must not abort
            LOGGER.warning("Model %s failed to train: %s", name, exc)
    return out


def model_comparison_table(results: Mapping[str, TrainingResults]) -> pd.DataFrame:
    """Summarise several models' out-of-fold metrics in one table.

    Args:
        results: Results keyed by model name.

    Returns:
        One row per model, sorted by descending average precision.
    """
    rows: List[Dict[str, Any]] = []
    for name, result in results.items():
        _, pr_low, pr_high = result.metric_cis.get("pr_auc", (np.nan, np.nan, np.nan))
        _, roc_low, roc_high = result.metric_cis.get("roc_auc", (np.nan, np.nan, np.nan))
        rows.append(
            {
                "model": name,
                "roc_auc": result.metrics.get("roc_auc", np.nan),
                "roc_auc_ci_low": roc_low,
                "roc_auc_ci_high": roc_high,
                "pr_auc": result.metrics.get("pr_auc", np.nan),
                "pr_auc_ci_low": pr_low,
                "pr_auc_ci_high": pr_high,
                "mcc": result.metrics.get("mcc", np.nan),
                "f1": result.metrics.get("f1", np.nan),
                "brier": result.metrics.get("brier", np.nan),
                "ece": result.metrics.get("ece", np.nan),
                "enrichment_at_1pct": result.metrics.get("enrichment_at_1pct", np.nan),
                "precision_at_1pct": result.metrics.get("precision_at_1pct", np.nan),
                "threshold": result.selected_threshold,
                "n_positive": result.n_positive,
                "n_negative": result.n_negative,
                "n_groups": result.n_groups,
            }
        )
    return pd.DataFrame(rows).sort_values("pr_auc", ascending=False, ignore_index=True)


def pairwise_delong(results: Mapping[str, TrainingResults]) -> pd.DataFrame:
    """DeLong comparisons between every pair of models on shared out-of-fold rows.

    Args:
        results: Results keyed by model name.

    Returns:
        A frame of pairwise AUROC differences and p-values, empty when fewer than
        two models share an aligned out-of-fold set.
    """
    from itertools import combinations

    names = [
        name
        for name, result in results.items()
        if result.oof_probabilities.size and result.oof_labels.size
    ]
    rows: List[Dict[str, Any]] = []
    for a, b in combinations(names, 2):
        ra, rb = results[a], results[b]
        if ra.oof_labels.size != rb.oof_labels.size or not np.array_equal(
            ra.oof_labels, rb.oof_labels
        ):
            LOGGER.warning("Skipping DeLong for %s vs %s: out-of-fold rows differ", a, b)
            continue
        difference, p_value = delong_roc_test(ra.oof_labels, ra.oof_probabilities, rb.oof_probabilities)
        rows.append(
            {
                "model_a": a,
                "model_b": b,
                "auc_a": ra.metrics.get("roc_auc", np.nan),
                "auc_b": rb.metrics.get("roc_auc", np.nan),
                "auc_difference": difference,
                "p_value": p_value,
            }
        )
    return pd.DataFrame(rows)


class MLPocketScorer:
    """Adapter presenting a trained classifier through the scorer interface.

    Args:
        classifier: A fitted :class:`CrypticSiteMLClassifier`.

    Raises:
        ValueError: If the classifier is not fitted.
    """

    def __init__(self, classifier: CrypticSiteMLClassifier) -> None:
        if classifier.best_estimator_ is None:
            raise ValueError("Classifier must be trained before using MLPocketScorer.")
        self.classifier = classifier

    @property
    def feature_names(self) -> Tuple[str, ...]:
        """Feature schema the wrapped model consumes."""
        return self.classifier.feature_names

    def calculate_composite_score(
        self,
        volume: float,
        depth: float,
        sasa: float,
        basic_count: int,
        potential: Optional[float] = None,
        plddt_confidence: float = float("nan"),
        **extra: float,
    ) -> float:
        """Score a single pocket described by the legacy argument set.

        Any descriptor outside the legacy six may be supplied through ``extra``;
        anything the model expects but does not receive is passed as ``nan`` and
        handled by the pipeline's imputer.

        Args:
            volume: Pocket volume (Å³).
            depth: Geometric burial depth (Å).
            sasa: Mean lining-residue SASA (Å²).
            basic_count: Basic residues in the lining shell.
            potential: Electrostatic potential (kT/e).
            plddt_confidence: Mean pLDDT of lining residues.
            **extra: Additional named descriptors.

        Returns:
            Probability of a cryptic IP-binding site.
        """
        row: Dict[str, float] = {name: float("nan") for name in self.feature_names}
        supplied = {
            "pocket_volume": volume,
            "pocket_depth": depth,
            "burial_depth": depth,
            "sasa": sasa,
            "sasa_mean": sasa,
            "n_basic_residues": float(basic_count),
            "electrostatic_potential": potential if potential is not None else float("nan"),
            "plddt_confidence": plddt_confidence,
            "plddt_mean": plddt_confidence,
        }
        supplied.update({key: float(value) for key, value in extra.items()})
        for name, value in supplied.items():
            if name in row:
                row[name] = float("nan") if value is None else float(value)
        return float(self.classifier.predict_proba(pd.DataFrame([row]))[0])

    def calculate_composite_scores(self, samples: pd.DataFrame) -> np.ndarray:
        """Vectorised scoring over many pockets.

        Args:
            samples: Table containing the model's feature columns. Missing
                columns are added as ``nan`` so callers built for the legacy
                schema keep working against a richer model.

        Returns:
            Probabilities, one per row.
        """
        frame = samples.copy()
        for name in self.feature_names:
            if name not in frame.columns:
                frame[name] = np.nan
        return self.classifier.predict_proba(frame)

    def classify_site(self, score: float) -> str:
        """Map a probability onto a confidence label.

        The boundaries follow the model's own fitted decision threshold rather
        than fixed constants, so the labels track the operating point that was
        actually selected.

        Args:
            score: Predicted probability.

        Returns:
            A human-readable confidence label.
        """
        threshold = self.classifier.decision_threshold_
        if score >= max(threshold, 0.5) + 0.5 * (1.0 - max(threshold, 0.5)):
            return "High confidence cryptic IP site"
        if score >= threshold:
            return "Moderate confidence candidate"
        if score >= 0.5 * threshold:
            return "Low confidence - manual inspection recommended"
        return "Unlikely cryptic IP site"
