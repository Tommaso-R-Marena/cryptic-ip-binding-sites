"""Positive-unlabelled learning and conformal sets for a screen with no negatives.

Why this exists
---------------
The project's machine-learning attempt reported a held-out random-forest
ROC-AUC of 0.499 against a threshold-score baseline of 0.556, and concluded
that threshold scoring stays the deployment mode. The stated cause was severe
pocket-level class imbalance: a handful of known buried sites against millions
of ordinary pockets.

That diagnosis is right about the symptom and wrong about the disease. The
training set was built by treating every pocket that is not a known site as a
negative. Almost none of those are known negatives - they are *unlabelled*.
Nobody has shown that the 78,920th pocket does not bind an inositol phosphate;
it has simply never been co-crystallised with one. A classifier trained with
unlabelled data as negatives learns the sampling of the PDB, and the PDB's
sampling of IP6 complexes is 134 InsP6 entries in 5 strict homology groups with
one group holding 87% of them.

So the correct frame is positive-unlabelled learning: a labelled positive set,
an unlabelled pool, and an unknown class prior. Under the standard
selected-completely-at-random assumption, a classifier trained to separate
labelled positives from the unlabelled pool estimates ``P(labelled | x)``,
which is ``c * P(positive | x)`` for a constant ``c`` - the label frequency.
Elkan and Noto's result is that ``c`` can be estimated from the classifier's
own output on held-out positives, which turns an uncalibrated score into a
probability.

What this module adds on top
----------------------------
1. ``c`` is estimated **per homology group**, held out as a whole, because the
   one thing this dataset guarantees is that copies of the same protein are not
   independent. A label-frequency estimate that resamples copies would be
   estimated on 134 entries' worth of apparent information and 5 groups' worth
   of real information.
2. Every prediction is returned as a **conformal set**, not a point estimate.
   Split-conformal prediction with group-wise exchangeability gives a set of
   labels that contains the truth with a stated frequency; on this dataset the
   honest output for most pockets is "either", and a method that cannot say so
   will say "positive" instead.
3. A **leakage audit** that fails loudly. The group-aware split is the single
   assumption the whole evaluation rests on, so it is checked rather than
   trusted, and the check runs on every fit.
4. A **minimum detectable effect** computed from the group count before any
   model is fitted, so the study knows in advance whether it could see
   anything.

What this cannot fix
--------------------
No estimator recovers information the data does not hold. With buried-site
positives in 5 homology groups, one holding 87%, the effective sample size for
anything that generalises across families is about five. This module makes that
visible instead of letting a 0.5 AUC look like a modelling failure, and the
pre-registered decision rule in docs/DETECTOR_PLAN.md stops the study when the
interval is too wide to decide anything.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

#: Below this many groups in either class, no estimate is evidence. The same
#: floor the rest of the project uses, restated here so the modelling code
#: cannot quietly adopt a looser one.
MIN_GROUPS = 5

#: Default miscoverage for conformal sets: 0.1 gives sets that contain the
#: truth 90% of the time, under exchangeability.
DEFAULT_ALPHA = 0.1


class PULearningError(ValueError):
    """A dataset or split from which no honest estimate can be made."""


# --------------------------------------------------------------- leakage audit

@dataclass(frozen=True)
class LeakageAudit:
    """Whether a train/test split keeps homology groups whole."""

    n_train: int
    n_test: int
    shared_groups: Tuple[object, ...]
    train_groups: int
    test_groups: int

    @property
    def clean(self) -> bool:
        return not self.shared_groups

    def as_dict(self) -> Dict[str, object]:
        return {
            "n_train": self.n_train,
            "n_test": self.n_test,
            "shared_groups": [str(g) for g in self.shared_groups],
            "train_groups": self.train_groups,
            "test_groups": self.test_groups,
            "clean": self.clean,
        }


def audit_split(groups: Sequence[object], train: Sequence[int], test: Sequence[int]) -> LeakageAudit:
    """Check that no homology group appears on both sides of a split."""
    groups = list(groups)
    train_groups = {groups[i] for i in train}
    test_groups = {groups[i] for i in test}
    shared = tuple(sorted(train_groups & test_groups, key=str))
    return LeakageAudit(
        n_train=len(train),
        n_test=len(test),
        shared_groups=shared,
        train_groups=len(train_groups),
        test_groups=len(test_groups),
    )


def group_folds(groups: Sequence[object], n_folds: int, *, seed: int) -> List[Tuple[List[int], List[int]]]:
    """Folds that split whole groups, balanced by group size.

    Groups are assigned to folds largest-first, each going to the fold holding
    fewest rows so far. With one group holding 87% of the positives, a random
    assignment leaves most folds with no positive at all; this is the greedy
    bin-packing that at least spreads what there is.
    """
    groups = list(groups)
    if n_folds < 2:
        raise PULearningError(f"need at least 2 folds, got {n_folds}")
    unique: Dict[object, List[int]] = {}
    for index, group in enumerate(groups):
        unique.setdefault(group, []).append(index)
    if len(unique) < n_folds:
        raise PULearningError(
            f"{len(unique)} homology group(s) cannot be split into {n_folds} folds without "
            "putting the same group on both sides, which is the leak this split exists to prevent"
        )
    rng = np.random.default_rng(seed)
    order = sorted(unique.items(), key=lambda item: (-len(item[1]), str(item[0])))
    buckets: List[List[int]] = [[] for _ in range(n_folds)]
    for _, members in order:
        target = min(range(n_folds), key=lambda k: (len(buckets[k]), rng.random()))
        buckets[target].extend(members)
    folds: List[Tuple[List[int], List[int]]] = []
    for k in range(n_folds):
        test = sorted(buckets[k])
        train = sorted(i for j, bucket in enumerate(buckets) if j != k for i in bucket)
        audit = audit_split(groups, train, test)
        if not audit.clean:
            raise PULearningError(f"fold {k} leaks groups {audit.shared_groups}")
        folds.append((train, test))
    return folds


# ----------------------------------------------------------- label frequency

#: How the label frequency is estimated. Elkan and Noto's own first estimator -
#: the mean of g(x) over labelled positives - is biased low by exactly the
#: class overlap, and that bias matters in the dangerous direction: the
#: posterior divides by c, so a c that is too small inflates P(positive).
#:
#: Measured with the estimator isolated from any classifier, by handing it the
#: exact g(x) = c * P(y=1|x) on 200,000 synthetic rows. Each entry is the
#: estimate divided by the truth, so 1.00 is unbiased:
#:
#:   true c    mean    q90    q95    q99    max
#:      0.2    0.70   0.96   0.98   0.99   1.00
#:      0.5    0.69   0.96   0.98   0.99   1.00
#:      0.9    0.70   0.96   0.98   0.99   1.00
#:
#: The bias in the mean is 30% and does not depend on c, because it is
#: E[P(y=1|x) | y=1], a property of how much the classes overlap. A high
#: quantile removes almost all of it. q95 is the default rather than max or
#: q99 because those rest on one or two rows, and with 134 entries in 5
#: homology groups there are not enough labelled positives per group to spend
#: on an extreme quantile.
LABEL_FREQUENCY_METHODS = ("quantile", "mean")

#: The quantile used by the default method, from the table above.
DEFAULT_LABEL_QUANTILE = 0.95


@dataclass(frozen=True)
class LabelFrequency:
    """The Elkan-Noto label frequency ``c`` and what it was estimated on."""

    c: float
    per_group: Dict[object, float]
    per_group_n: Dict[object, int]
    n_groups: int
    method: str
    quantile: Optional[float] = None

    @property
    def is_evidence(self) -> bool:
        """Whether enough groups contributed for this to mean anything."""
        return self.n_groups >= MIN_GROUPS

    @property
    def thinnest_group(self) -> int:
        """Labelled positives in the group that had fewest of them."""
        return min(self.per_group_n.values()) if self.per_group_n else 0

    @property
    def group_spread(self) -> float:
        """Standard deviation of the per-group estimates: the honest error bar.

        With five groups this is a standard deviation of five numbers, and
        reporting ``c`` without it would imply a precision that 134 entries
        suggest and 5 groups do not support.
        """
        values = list(self.per_group.values())
        return float(np.std(values, ddof=1)) if len(values) > 1 else float("nan")

    def as_dict(self) -> Dict[str, object]:
        return {
            "c": float(self.c),
            "n_groups": self.n_groups,
            "method": self.method,
            "quantile": self.quantile,
            "is_evidence": self.is_evidence,
            "group_spread": self.group_spread,
            "thinnest_group": self.thinnest_group,
            "per_group": {str(k): float(v) for k, v in self.per_group.items()},
            "per_group_n": {str(k): int(v) for k, v in self.per_group_n.items()},
        }


def estimate_label_frequency(
    scores: Sequence[float],
    labelled: Sequence[bool],
    groups: Sequence[object],
    *,
    method: str = "quantile",
    quantile: float = DEFAULT_LABEL_QUANTILE,
) -> LabelFrequency:
    """Estimate ``c = P(labelled | positive)``, holding homology groups whole.

    ``scores`` are a non-traditional classifier's outputs - the probability
    that a row is *labelled*, not that it is positive. The estimate is only
    unbiased if those scores were produced out-of-fold; the caller is
    responsible for that and :func:`fit_pu` does it.

    The statistic is computed within each homology group and then averaged
    across groups, so one family holding 87% of the positives cannot set ``c``
    on its own. See ``LABEL_FREQUENCY_METHODS`` for why the default is a
    quantile rather than a mean.
    """
    scores = np.asarray(scores, dtype=float)
    labelled = np.asarray(labelled, dtype=bool)
    groups = list(groups)
    if not (len(scores) == len(labelled) == len(groups)):
        raise PULearningError("scores, labels and groups must be the same length")
    if not np.any(labelled):
        raise PULearningError("no labelled positives: c is not estimable")
    if np.any((scores < 0) | (scores > 1)):
        raise PULearningError("scores must be probabilities in [0, 1]")
    if method not in LABEL_FREQUENCY_METHODS:
        raise PULearningError(f"method must be one of {LABEL_FREQUENCY_METHODS}, got {method!r}")
    if not 0.0 < quantile <= 1.0:
        raise PULearningError(f"quantile must be in (0, 1], got {quantile}")

    per_group: Dict[object, float] = {}
    per_group_n: Dict[object, int] = {}
    for group in sorted({groups[i] for i in np.nonzero(labelled)[0]}, key=str):
        mask = np.array([g == group for g in groups]) & labelled
        values = scores[mask]
        per_group_n[group] = int(values.size)
        if method == "mean":
            per_group[group] = float(np.mean(values))
        else:
            per_group[group] = float(np.quantile(values, quantile))
    c = float(np.mean(list(per_group.values())))
    if c <= 0.0:
        raise PULearningError(
            "the classifier gives labelled positives zero probability of being labelled, "
            "so c is zero and no probability can be recovered"
        )
    return LabelFrequency(
        c=c,
        per_group=per_group,
        per_group_n=per_group_n,
        n_groups=len(per_group),
        method=method,
        quantile=quantile if method == "quantile" else None,
    )


def pu_posterior(scores: Sequence[float], frequency: LabelFrequency) -> np.ndarray:
    """``P(positive | x) = g(x) / c``, clipped to [0, 1].

    The clip is not cosmetic. ``c`` is estimated, so ``g(x)/c`` can exceed one
    for a row the classifier is more confident about than any held-out
    positive; clipping is the only honest thing to do with that, and it is
    visible in the output as a saturated probability rather than as a number
    above one that some later step divides by.
    """
    scores = np.asarray(scores, dtype=float)
    if frequency.c <= 0:
        raise PULearningError("label frequency must be positive")
    return np.clip(scores / frequency.c, 0.0, 1.0)


# ------------------------------------------------------------------ conformal

@dataclass(frozen=True)
class ConformalSets:
    """Set-valued predictions with a stated coverage target."""

    alpha: float
    threshold: float
    contains_positive: np.ndarray
    contains_negative: np.ndarray
    n_calibration: int
    n_calibration_groups: int

    @property
    def sizes(self) -> np.ndarray:
        return self.contains_positive.astype(int) + self.contains_negative.astype(int)

    @property
    def fraction_decided(self) -> float:
        """How often the set holds exactly one label - an actual decision."""
        return float(np.mean(self.sizes == 1))

    @property
    def fraction_empty(self) -> float:
        """An empty set means the row looks like neither class."""
        return float(np.mean(self.sizes == 0))

    def as_dict(self) -> Dict[str, object]:
        return {
            "alpha": self.alpha,
            "threshold": self.threshold,
            "n_calibration": self.n_calibration,
            "n_calibration_groups": self.n_calibration_groups,
            "fraction_decided": self.fraction_decided,
            "fraction_empty": self.fraction_empty,
            "mean_set_size": float(np.mean(self.sizes)),
        }


def conformal_threshold(
    calibration_scores: Sequence[float],
    *,
    alpha: float = DEFAULT_ALPHA,
) -> float:
    """The split-conformal threshold at miscoverage ``alpha``.

    Uses the finite-sample corrected quantile - ``ceil((n+1)(1-alpha))/n`` -
    rather than the plain empirical one, because with a few dozen calibration
    points the difference between them is the difference between the stated
    coverage and a coverage that is quietly lower.
    """
    scores = np.asarray(calibration_scores, dtype=float)
    if scores.size == 0:
        raise PULearningError("no calibration scores")
    if not 0.0 < alpha < 1.0:
        raise PULearningError(f"alpha must be in (0, 1), got {alpha}")
    n = scores.size
    rank = int(np.ceil((n + 1) * (1.0 - alpha)))
    if rank > n:
        raise PULearningError(
            f"{n} calibration points cannot support alpha={alpha}: the corrected quantile "
            f"needs rank {rank}. Lower the confidence or calibrate on more groups."
        )
    return float(np.sort(scores)[rank - 1])


def conformal_sets(
    test_scores: Sequence[float],
    calibration_scores: Sequence[float],
    calibration_groups: Sequence[object],
    *,
    alpha: float = DEFAULT_ALPHA,
) -> ConformalSets:
    """Label sets for ``test_scores``, calibrated on held-out positives.

    The nonconformity score is ``1 - P(positive | x)`` on known positives, so a
    test row is admitted to the positive class when its probability is at least
    as high as the alpha-quantile of the calibration positives. A row can end
    up in both classes, which is the method saying it does not know.
    """
    calibration = np.asarray(calibration_scores, dtype=float)
    groups = list(calibration_groups)
    if len(groups) != calibration.size:
        raise PULearningError("one group per calibration score")
    n_groups = len({*groups})
    if n_groups < MIN_GROUPS:
        raise PULearningError(
            f"conformal calibration on {n_groups} homology group(s) is not evidence; "
            f"the pre-registered floor is {MIN_GROUPS}. Exchangeability across copies of one "
            "protein is not exchangeability across proteins."
        )
    threshold = conformal_threshold(1.0 - calibration, alpha=alpha)
    test = np.asarray(test_scores, dtype=float)
    # Positive is in the set when its nonconformity does not exceed the
    # threshold; negative is in the set when the complementary score does not.
    contains_positive = (1.0 - test) <= threshold
    contains_negative = test <= threshold
    return ConformalSets(
        alpha=float(alpha),
        threshold=float(threshold),
        contains_positive=contains_positive,
        contains_negative=contains_negative,
        n_calibration=int(calibration.size),
        n_calibration_groups=n_groups,
    )


# ------------------------------------------------------- minimum detectable

def minimum_detectable_auc(n_positive_groups: int, n_negative_groups: int, *, power: float = 0.8) -> float:
    """The smallest AUC a group-level comparison of this size could detect.

    Computed from the group counts, not the row counts, and computed *before*
    fitting anything. With buried-site positives in 5 homology groups this
    number is large, and knowing that in advance is the difference between
    reporting "no signal" and reporting "no power".

    Uses the Hanley-McNeil variance of an AUC under the normal approximation,
    solved for the effect size that a two-sided 0.05 test would reject with
    the stated power. It is an approximation and it is quoted to two decimals
    for that reason.
    """
    if n_positive_groups < 1 or n_negative_groups < 1:
        raise PULearningError("both classes need at least one group")
    z_alpha, z_beta = 1.959963985, 0.8416212336  # two-sided 0.05, power 0.8
    if not 0.5 < power < 1.0:
        raise PULearningError(f"power must be in (0.5, 1), got {power}")
    if power != 0.8:
        # Only the common case is tabulated; anything else is refused rather
        # than silently computed with the wrong constant.
        raise PULearningError("only power=0.8 is implemented")
    m, n = float(n_positive_groups), float(n_negative_groups)
    # Variance of the AUC at the null (0.5) is (m + n + 1) / (12 * m * n).
    null_se = np.sqrt((m + n + 1.0) / (12.0 * m * n))
    delta = (z_alpha + z_beta) * null_se
    return float(min(0.5 + delta, 1.0))


# ------------------------------------------------------------------- the fit

@dataclass
class PUFit:
    """A fitted positive-unlabelled model and everything needed to judge it."""

    frequency: LabelFrequency
    out_of_fold_scores: np.ndarray
    posterior: np.ndarray
    audits: List[LeakageAudit] = field(default_factory=list)
    minimum_detectable_auc: Optional[float] = None
    n_folds: int = 0
    seed: int = 0

    @property
    def leakage_clean(self) -> bool:
        return all(audit.clean for audit in self.audits)

    def as_dict(self) -> Dict[str, object]:
        return {
            "label_frequency": self.frequency.as_dict(),
            "n_folds": self.n_folds,
            "seed": self.seed,
            "leakage_clean": self.leakage_clean,
            "minimum_detectable_auc": self.minimum_detectable_auc,
            "audits": [audit.as_dict() for audit in self.audits],
            "posterior_summary": {
                "mean": float(np.mean(self.posterior)),
                "max": float(np.max(self.posterior)),
                "saturated_fraction": float(np.mean(self.posterior >= 1.0)),
            },
        }


#: A non-traditional classifier: fit on (features, is_labelled) and return
#: P(labelled | x) for the rows it is given. Anything with this shape works,
#: which keeps the statistics here independent of the model choice.
Classifier = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]


def fit_pu(
    features: np.ndarray,
    labelled: Sequence[bool],
    groups: Sequence[object],
    classifier: Classifier,
    *,
    n_folds: int = 5,
    seed: int = 20261010,
    method: str = "quantile",
    quantile: float = DEFAULT_LABEL_QUANTILE,
) -> PUFit:
    """Fit a positive-unlabelled model with group-held-out scoring.

    ``classifier(train_features, train_labelled, test_features)`` returns
    ``P(labelled | x)`` for the test rows. Every score used to estimate the
    label frequency is produced by a model that never saw that row's homology
    group, which is the condition the Elkan-Noto estimator needs and the thing
    the original 0.499 experiment could not have had with copies split across
    folds.
    """
    features = np.asarray(features, dtype=float)
    labelled_arr = np.asarray(labelled, dtype=bool)
    groups = list(groups)
    if features.ndim != 2:
        raise PULearningError("features must be a 2-d array")
    if not (len(features) == labelled_arr.size == len(groups)):
        raise PULearningError("features, labels and groups must agree in length")
    if not np.any(labelled_arr):
        raise PULearningError("no labelled positives to learn from")

    positive_groups = len({groups[i] for i in np.nonzero(labelled_arr)[0]})
    unlabelled_groups = len({groups[i] for i in np.nonzero(~labelled_arr)[0]})

    folds = group_folds(groups, n_folds, seed=seed)
    out_of_fold = np.full(len(features), np.nan)
    audits: List[LeakageAudit] = []
    for train, test in folds:
        audit = audit_split(groups, train, test)
        audits.append(audit)
        if not audit.clean:
            raise PULearningError(f"split leaks groups {audit.shared_groups}")
        predicted = np.asarray(
            classifier(features[train], labelled_arr[train], features[test]), dtype=float,
        )
        if predicted.shape != (len(test),):
            raise PULearningError("the classifier must return one score per test row")
        out_of_fold[test] = predicted
    if np.any(np.isnan(out_of_fold)):
        raise PULearningError("the folds did not cover every row")

    frequency = estimate_label_frequency(
        out_of_fold, labelled_arr, groups, method=method, quantile=quantile,
    )
    posterior = pu_posterior(out_of_fold, frequency)
    return PUFit(
        frequency=frequency,
        out_of_fold_scores=out_of_fold,
        posterior=posterior,
        audits=audits,
        minimum_detectable_auc=minimum_detectable_auc(
            max(positive_groups, 1), max(unlabelled_groups, 1),
        ),
        n_folds=n_folds,
        seed=seed,
    )


def group_auc(scores: Sequence[float], positive: Sequence[bool], groups: Sequence[object]) -> float:
    """ROC-AUC computed on per-group mean scores, not per-row.

    A row-level AUC over 134 entries in 5 families reports the separability of
    5 things with the confidence of 134. Averaging within a group first is the
    cheapest way to stop that.
    """
    scores = np.asarray(scores, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    groups = list(groups)
    if not (scores.size == positive.size == len(groups)):
        raise PULearningError("scores, labels and groups must be the same length")
    per_group: Dict[object, List[float]] = {}
    per_label: Dict[object, bool] = {}
    for score, label, group in zip(scores, positive, groups):
        per_group.setdefault(group, []).append(float(score))
        if group in per_label and per_label[group] != bool(label):
            raise PULearningError(
                f"group {group!r} holds both classes, so a group-level AUC is not defined for it"
            )
        per_label[group] = bool(label)
    pos = [np.mean(v) for g, v in per_group.items() if per_label[g]]
    neg = [np.mean(v) for g, v in per_group.items() if not per_label[g]]
    if not pos or not neg:
        raise PULearningError("a group-level AUC needs groups of both classes")
    wins = 0.0
    for a in pos:
        for b in neg:
            wins += 1.0 if a > b else (0.5 if a == b else 0.0)
    return float(wins / (len(pos) * len(neg)))


def logistic_classifier(
    l2: float = 1.0,
    iterations: int = 200,
    learning_rate: float = 0.1,
) -> Classifier:
    """A small, dependency-free non-traditional classifier for the PU fit.

    Plain regularised logistic regression by gradient descent on standardised
    features. It is here so the statistics above can be exercised and tested
    without pulling in a model choice: anything matching :data:`Classifier`
    can replace it, and a real run should use whatever the study pre-registers.
    """

    def fit(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray) -> np.ndarray:
        mean = train_x.mean(axis=0)
        scale = train_x.std(axis=0)
        scale[scale == 0] = 1.0
        x = (train_x - mean) / scale
        y = train_y.astype(float)
        weights = np.zeros(x.shape[1])
        bias = 0.0
        for _ in range(iterations):
            logit = x @ weights + bias
            prediction = 1.0 / (1.0 + np.exp(-np.clip(logit, -60, 60)))
            error = prediction - y
            weights -= learning_rate * (x.T @ error / len(x) + l2 * weights / len(x))
            bias -= learning_rate * float(np.mean(error))
        z = (test_x - mean) / scale
        logit = z @ weights + bias
        return 1.0 / (1.0 + np.exp(-np.clip(logit, -60, 60)))

    return fit
