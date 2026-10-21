"""A monotone additive score that abstains instead of failing on missing input.

Why this exists
---------------
The deployed screen calls a protein a hit when its best confident pocket
clears every one of five hard gates: composite score >= 0.75, lining SASA
<= 10 A^2, at least 4 basic residues, volume 300-1600 A^3, mean pLDDT >= 70.
The project's own known-binder check shows what that costs. ADAR2 and ADAR1
pass; ADAT1, PDS5B and HDAC1 are blocked at the score gate and HDAC3 at the
volume gate. All four are annotated binders. A cascade of hard gates is a
conjunction, so its false-negative rate is the union of five false-negative
rates, and nothing about the cascade says which gate did the blocking in a way
that a reader can weigh.

Two things are wrong and both are fixable without fitting anything new.

**A missing input is not a failing input.** If a pocket has no pLDDT because
it came from a crystal structure, or no electrostatic term because APBS was
not installed - which this pipeline degrades past with a warning - the gate
reads a missing value as a failure. The honest answer is "unknown". This module
computes the score over the *interval* a missing feature could take and
abstains when that interval spans the decision boundary, so a protein is never
rejected for a measurement nobody made.

**Monotonicity should be guaranteed, not hoped for.** With buried-site
positives in five homology groups, a model free to fit any shape will learn
that more burial is worse if that is what five families happen to show. Every
feature here is given a direction in advance, the shape functions are monotone
by construction, the weights are non-negative, and
``test_the_score_is_monotone_in_every_feature`` checks the composition over
random inputs rather than trusting the algebra.

What this is not
----------------
It is not a fitted model and it does not learn. It is the interpretable
baseline the project already had, made monotone by construction and made
honest about missing data. The learned part lives in
``cryptic_ip.detect.pu``, and the pre-registered comparison between them is in
docs/DETECTOR_PLAN.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

#: A feature's direction: +1 when a larger value makes a site more likely,
#: -1 when it makes it less likely. Declared in advance, which is the point.
UP, DOWN = 1, -1

#: Groups of features that measure the same property by different means. When
#: any member of a group is present, the absent members are dropped from the
#: score rather than widening its interval: an alternative measurement nobody
#: needed to take is not an unknown. Without this, a pocket whose burial came
#: from the cavity detector would still be penalised for lacking the
#: ligand-copy relative SASA, and every score would abstain.
SUBSTITUTE_GROUPS: Tuple[Tuple[str, ...], ...] = (
    ("relative_sasa", "cavity_openness"),
)


def _is_given(value) -> bool:
    """Whether a feature value counts as measured."""
    if value is None:
        return False
    try:
        return not np.isnan(float(value))
    except (TypeError, ValueError):
        return False


class MonotoneError(ValueError):
    """A specification or input from which no score can be computed."""


@dataclass(frozen=True)
class Feature:
    """One term of the score: a monotone ramp over a bounded input.

    ``low`` and ``high`` are the range the feature can physically take, and
    they are what the interval arithmetic uses when the value is missing. They
    are not clipping bounds for present values - a value outside them is a
    measurement error and is reported as such.
    """

    name: str
    direction: int
    midpoint: float
    slope: float
    low: float
    high: float
    weight: float
    units: str = ""
    why: str = ""
    veto_below: Optional[float] = None
    veto_above: Optional[float] = None
    veto_why: str = ""

    def __post_init__(self) -> None:
        if self.direction not in (UP, DOWN):
            raise MonotoneError(f"{self.name}: direction must be +1 or -1, got {self.direction}")
        if self.slope <= 0:
            raise MonotoneError(f"{self.name}: slope must be positive, got {self.slope}")
        if self.weight < 0:
            raise MonotoneError(f"{self.name}: weight cannot be negative, got {self.weight}")
        if not self.low < self.high:
            raise MonotoneError(f"{self.name}: low must be below high, got {self.low} and {self.high}")
        if (self.veto_below is not None or self.veto_above is not None) and not self.veto_why:
            raise MonotoneError(
                f"{self.name}: a veto has to state its reason. A bound that can reject a protein "
                "outright without saying why is the gate cascade this score replaces."
            )

    def value_at(self, x: float) -> float:
        """The feature's contribution in [0, 1], monotone in ``x`` by construction."""
        signed = self.direction * self.slope * (float(x) - self.midpoint)
        return float(1.0 / (1.0 + np.exp(-np.clip(signed, -60.0, 60.0))))

    def vetoes(self, x: float) -> bool:
        """Whether ``x`` puts this feature outside physical admissibility.

        A veto is not a tuned gate. It is a statement that no amount of the
        other terms can compensate - a pocket whose electrostatic well repels
        the ligand is not a site for it however buried it is. Each veto has to
        carry its reason, and a missing veto feature abstains rather than
        passing.
        """
        value = float(x)
        if self.veto_below is not None and value < self.veto_below:
            return True
        if self.veto_above is not None and value > self.veto_above:
            return True
        return False

    def bounds(self) -> Tuple[float, float]:
        """The contribution's range over the feature's whole physical range.

        Because the ramp is monotone, the extremes are at the ends - which is
        what makes the interval arithmetic exact rather than a bound.
        """
        a, b = self.value_at(self.low), self.value_at(self.high)
        return (min(a, b), max(a, b))

    def as_dict(self) -> Dict[str, object]:
        return {
            "name": self.name, "direction": self.direction, "midpoint": self.midpoint,
            "slope": self.slope, "low": self.low, "high": self.high,
            "weight": self.weight, "units": self.units, "why": self.why,
        }


#: The specification. Every midpoint is the value at which the term is 0.5 and
#: carries the observation that motivates it; every direction is a claim about
#: physics that the fitting cannot overturn.
DEFAULT_FEATURES: Tuple[Feature, ...] = (
    Feature(
        name="relative_sasa", direction=DOWN, midpoint=0.12, slope=25.0,
        low=0.0, high=1.0, weight=1.0, units="fraction",
        why="the project's own burial cutoff; the survey found burial is a continuum, so this "
            "is a ramp rather than the cut it replaced",
    ),
    Feature(
        name="cavity_openness", direction=DOWN, midpoint=0.30, slope=12.0,
        low=0.0, high=1.0, weight=1.0, units="fraction",
        why="burial measured by the cavity detector itself: one minus the mean fraction of the "
            "seven scan directions in which the cavity is enclosed. A different quantity from "
            "relative_sasa, which is measured on a ligand copy, and given its own name for that "
            "reason - a burial term computed one way and named after another is how this "
            "pipeline came to read a hydrophobic density as a depth in angstrom",
    ),
    Feature(
        name="field_depth_kt", direction=UP, midpoint=1.0, slope=1.0,
        low=-20.0, high=20.0, weight=1.5, units="kT/e",
        why="a -9 polyanion needs a positive well, and this is the term the basic-residue count "
            "was standing in for",
        veto_below=0.0,
        veto_why="an electrostatic well that repels the ligand cannot hold it, however buried "
                 "the pocket is. Without this, a perfectly enclosed acidic cavity with no basic "
                 "residues scored 0.554 against a 0.55 threshold and was called a site - burial "
                 "and confidence outvoting the physics",
    ),
    Feature(
        name="field_shape", direction=UP, midpoint=0.5, slope=6.0,
        low=-1.0, high=1.0, weight=1.5, units="correlation",
        why="whether the positive charge forms a shell at the phosphate radius rather than a "
            "surface patch; this is what separates an IP6 site from an ATP site",
    ),
    Feature(
        name="openability", direction=UP, midpoint=0.3, slope=8.0,
        low=0.0, high=1.0, weight=1.0, units="overlap",
        why="a cryptic site is closed in the apo model, so being open is the wrong thing to "
            "require; this asks whether a soft mode opens it",
    ),
    Feature(
        name="n_basic_residues", direction=UP, midpoint=3.5, slope=1.2,
        low=0.0, high=30.0, weight=0.5, units="count",
        why="kept at a low weight because it is the term that lets a surface patch score well; "
            "it is a necessary condition, not a discriminating one",
    ),
    Feature(
        name="pocket_volume", direction=UP, midpoint=300.0, slope=0.01,
        low=0.0, high=5000.0, weight=0.5, units="A^3",
        why="an IP6 needs room; the upper gate that blocked HDAC3 is deliberately absent, "
            "because a large pocket is not evidence against a site within it",
    ),
    Feature(
        name="mean_plddt", direction=UP, midpoint=70.0, slope=0.15,
        low=0.0, high=100.0, weight=0.75, units="pLDDT",
        why="a pocket whose lining is a low-confidence prediction is evidence about the "
            "predictor, not the protein; missing for crystal structures, which is what the "
            "abstention is for",
    ),
)


@dataclass(frozen=True)
class ScoreInterval:
    """A score that may be an interval, because some input was missing."""

    low: float
    high: float
    present: Tuple[str, ...]
    missing: Tuple[str, ...]
    contributions: Dict[str, float]
    vetoed: Tuple[str, ...] = ()
    veto_pending: Tuple[str, ...] = ()

    @property
    def point(self) -> Optional[float]:
        """The score, when nothing was missing."""
        return float(self.low) if not self.missing else None

    @property
    def width(self) -> float:
        return float(self.high - self.low)

    def decision(self, threshold: float) -> str:
        """``site``, ``not_site`` or ``abstain``.

        A fired veto is decisive: no weighted sum overrides a physical
        impossibility. A veto feature that was never measured cannot be
        cleared, so the best available answer is abstention rather than a
        site call on the strength of the other terms.

        Otherwise abstention is the whole point: a pocket whose score interval
        straddles the threshold has not been measured well enough to call, and
        calling it either way is a claim about data nobody collected.
        """
        if self.vetoed:
            return "not_site"
        if self.veto_pending:
            return "abstain" if self.high >= threshold else "not_site"
        if self.low >= threshold:
            return "site"
        if self.high < threshold:
            return "not_site"
        return "abstain"

    @property
    def veto_reasons(self) -> Tuple[str, ...]:
        """Why each fired veto fired, for the report."""
        return self.vetoed

    def as_dict(self) -> Dict[str, object]:
        return {
            "low": float(self.low), "high": float(self.high), "width": self.width,
            "point": self.point,
            "present": list(self.present), "missing": list(self.missing),
            "vetoed": list(self.vetoed), "veto_pending": list(self.veto_pending),
            "contributions": {k: float(v) for k, v in self.contributions.items()},
        }


class MonotoneScorer:
    """A weighted sum of monotone ramps, with exact intervals for missing input."""

    def __init__(
        self,
        features: Sequence[Feature] = DEFAULT_FEATURES,
        substitute_groups: Sequence[Sequence[str]] = SUBSTITUTE_GROUPS,
    ) -> None:
        if not features:
            raise MonotoneError("a scorer needs at least one feature")
        names = [f.name for f in features]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise MonotoneError(f"duplicate feature(s) {duplicates}")
        total = float(sum(f.weight for f in features))
        if total <= 0:
            raise MonotoneError("feature weights must sum to a positive value")
        self.features: Tuple[Feature, ...] = tuple(features)
        self.total_weight = total
        # Either burial term satisfies the score; they are separate names
        # because they are separate measurements.
        self.burial_terms = tuple(
            name for name in ("relative_sasa", "cavity_openness") if name in names
        )
        self.substitute_groups = tuple(tuple(group) for group in substitute_groups)

    @property
    def names(self) -> Tuple[str, ...]:
        return tuple(f.name for f in self.features)

    def score(self, values: Mapping[str, Optional[float]], *, strict: bool = False) -> ScoreInterval:
        """Score one pocket. Missing or ``None`` features widen the interval.

        With ``strict`` a missing feature is an error instead, which is what a
        benchmark wants when every input is supposed to be present.
        """
        unknown = set(values) - set(self.names)
        if unknown:
            raise MonotoneError(f"unknown feature(s) {sorted(unknown)}; known: {list(self.names)}")
        # Within a substitute group, an absent member is irrelevant rather than
        # unknown as long as a sibling was measured.
        covered: set = set()
        for group in self.substitute_groups:
            members = [name for name in group if name in self.names]
            if any(_is_given(values.get(name)) for name in members):
                covered.update(name for name in members if not _is_given(values.get(name)))

        low = high = 0.0
        weight_used = 0.0
        present: List[str] = []
        missing: List[str] = []
        vetoed: List[str] = []
        veto_pending: List[str] = []
        contributions: Dict[str, float] = {}
        for feature in self.features:
            if feature.name in covered:
                continue
            weight_used += feature.weight
            raw = values.get(feature.name)
            if raw is None or (isinstance(raw, float) and np.isnan(raw)):
                if strict:
                    raise MonotoneError(f"feature {feature.name!r} is missing")
                bound_low, bound_high = feature.bounds()
                low += feature.weight * bound_low
                high += feature.weight * bound_high
                missing.append(feature.name)
                if feature.veto_below is not None or feature.veto_above is not None:
                    veto_pending.append(feature.name)
                continue
            value = float(raw)
            if not feature.low - 1e-9 <= value <= feature.high + 1e-9:
                raise MonotoneError(
                    f"feature {feature.name!r} is {value}{feature.units and ' ' + feature.units}, "
                    f"outside its physical range [{feature.low}, {feature.high}]. That is a "
                    "measurement error rather than an extreme pocket, and clipping it would hide it."
                )
            if feature.vetoes(value):
                vetoed.append(feature.name)
            contribution = feature.value_at(value)
            contributions[feature.name] = contribution
            low += feature.weight * contribution
            high += feature.weight * contribution
            present.append(feature.name)
        if weight_used <= 0:
            raise MonotoneError("every feature was substituted away; nothing left to score")
        return ScoreInterval(
            low=low / weight_used,
            high=high / weight_used,
            present=tuple(present),
            missing=tuple(missing),
            contributions=contributions,
            vetoed=tuple(vetoed),
            veto_pending=tuple(veto_pending),
        )

    def sensitivity(self, values: Mapping[str, Optional[float]]) -> Dict[str, float]:
        """How much each present feature's term contributes to the score.

        A ranking that cannot say why a pocket ranks where it does is not usable
        by anyone deciding what to test at a bench.
        """
        interval = self.score(values)
        used = sum(
            f.weight for f in self.features
            if f.name in interval.contributions or f.name in interval.missing
        )
        used = used or self.total_weight
        return {
            name: float(interval.contributions[name] * feature.weight / used)
            for feature, name in ((f, f.name) for f in self.features)
            if name in interval.contributions
        }

    def blocking_feature(self, values: Mapping[str, Optional[float]], threshold: float) -> Optional[str]:
        """Which present feature most holds a below-threshold pocket back.

        The gate cascade could say a protein failed; it could not say which
        measurement would have to change, which is what an experimentalist
        needs. This answers that by finding the feature whose term is furthest
        below its own best case.
        """
        interval = self.score(values)
        if interval.vetoed:
            return interval.vetoed[0]
        if interval.decision(threshold) == "site":
            return None
        shortfalls: Dict[str, float] = {}
        for feature in self.features:
            if feature.name not in interval.contributions:
                continue
            best = max(feature.bounds())
            shortfalls[feature.name] = (best - interval.contributions[feature.name]) * feature.weight
        if not shortfalls:
            return None
        return max(shortfalls.items(), key=lambda item: item[1])[0]

    def as_dict(self) -> Dict[str, object]:
        return {
            "features": [f.as_dict() for f in self.features],
            "total_weight": self.total_weight,
        }
