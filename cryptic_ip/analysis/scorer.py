"""Transparent rule-based scoring for cryptic IP binding sites.

Role of this scorer
-------------------
The rule-based score is the interpretable baseline. It encodes the published
criteria for a buried inositol phosphate site directly, so a result can be
explained without reference to a fitted model, and so the machine-learned model
has something meaningful to be compared against. It is *not* a probability, and
:class:`CalibratedRuleScorer` exists to convert it into one when labelled data is
available.

Three changes make it defensible:

**Smooth component functions.** The original components were step functions -
a basic-residue count of 4 scored 0.8 and a count of 3 scored 0.4. Discontinuities
that large make the score unstable under any measurement noise, and they are not
justified by anything in the underlying biology. Every component here is a smooth
monotone function of its input, so a small change in a measurement produces a
small change in the score.

**Correct inputs.** ``depth`` means the geometric burial depth - distance from the
pocket centre to the nearest solvent-exposed atom - not fpocket's mean local
hydrophobic density, which is what the previous pipeline actually passed.

**Explicit, documented parameters.** Every midpoint and slope is a named constant
with the observation that motivates it, and weights are overridable, so the
scoring function can be recalibrated on controls rather than being folded into
the code.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Mapping, Optional, Sequence

import numpy as np

LOGGER = logging.getLogger(__name__)


def _logistic(x: float, midpoint: float, slope: float) -> float:
    """Logistic ramp: 0.5 at ``midpoint``, rising with positive ``slope``.

    Args:
        x: Input value.
        midpoint: Value at which the output is 0.5.
        slope: Steepness; larger values approach a step function.

    Returns:
        A value in ``(0, 1)``.
    """
    # Clip the exponent to avoid overflow warnings at extreme inputs.
    exponent = float(np.clip(-slope * (x - midpoint), -60.0, 60.0))
    return float(1.0 / (1.0 + np.exp(exponent)))


@dataclass
class ScoringParameters:
    """Tunable parameters of the rule-based score.

    Attributes:
        weights: Component weights; normalised to sum to one at use time.
        volume_optimum_low: Lower edge of the ideal cavity volume (Å³). InsP3
            occupies roughly 300 Å³ and InsP6 roughly 600 Å³ including its
            hydration shell.
        volume_optimum_high: Upper edge of the ideal cavity volume (Å³).
        volume_tolerance: Width over which the volume score decays outside the
            optimum (Å³).
        depth_midpoint: Burial depth scoring 0.5 (Å). The published criterion is
            "deeper than 15 Å"; a midpoint of 12 Å with a moderate slope reaches
            ~0.8 at 15 Å rather than switching abruptly there.
        depth_slope: Steepness of the depth ramp (per Å).
        sasa_midpoint: Mean lining-residue SASA scoring 0.5 (Å²).
        sasa_slope: Steepness of the SASA ramp (per Å²); negative direction is
            applied internally since lower SASA is better.
        basic_midpoint: Basic-residue count scoring 0.5. ADAR2 coordinates its
            InsP6 with six basic residues, so the ramp centres just below that.
        basic_slope: Steepness of the basic-residue ramp.
        potential_midpoint: Electrostatic potential scoring 0.5 (kT/e).
        potential_slope: Steepness of the electrostatic ramp (per kT/e).
        enclosure_midpoint: Enclosure fraction scoring 0.5.
        enclosure_slope: Steepness of the enclosure ramp.
        missing_component_score: Score assigned to a component whose input is
            unavailable. Neutral by design, so a missing measurement neither
            rewards nor penalises a pocket.
    """

    weights: Dict[str, float] = field(
        default_factory=lambda: {
            "volume": 0.10,
            "depth": 0.22,
            "sasa": 0.25,
            "basic_residues": 0.20,
            "electrostatics": 0.10,
            "enclosure": 0.13,
        }
    )
    volume_optimum_low: float = 300.0
    volume_optimum_high: float = 800.0
    volume_tolerance: float = 400.0
    depth_midpoint: float = 12.0
    depth_slope: float = 0.45
    sasa_midpoint: float = 20.0
    sasa_slope: float = 0.12
    basic_midpoint: float = 3.5
    basic_slope: float = 1.1
    potential_midpoint: float = 3.0
    potential_slope: float = 0.5
    enclosure_midpoint: float = 0.75
    enclosure_slope: float = 12.0
    missing_component_score: float = 0.5

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serialisable representation."""
        return asdict(self)


class PocketScorer:
    """Score pockets against the cryptic IP-binding site criteria.

    Criteria, and the component that expresses each:

    1. Deeply buried below the surface - ``depth`` (geometric burial depth).
    2. Little solvent exposure - ``sasa`` (mean lining-residue SASA).
    3. Enclosed rather than merely concave - ``enclosure``.
    4. A cluster of basic residues - ``basic_residues``.
    5. Positive electrostatic potential - ``electrostatics``.
    6. A cavity sized for InsP3-InsP6 - ``volume``.

    Args:
        weights: Optional component weights, merged over the defaults and
            normalised to sum to one.
        parameters: Optional full parameter set.
    """

    def __init__(
        self,
        weights: Optional[Mapping[str, float]] = None,
        parameters: Optional[ScoringParameters] = None,
    ) -> None:
        self.parameters = parameters or ScoringParameters()
        if weights:
            merged = dict(self.parameters.weights)
            merged.update({str(k): float(v) for k, v in weights.items()})
            self.parameters.weights = merged
        self.weights = self._normalise_weights(self.parameters.weights)

    @staticmethod
    def _normalise_weights(weights: Mapping[str, float]) -> Dict[str, float]:
        """Normalise weights to sum to one.

        Args:
            weights: Raw weights.

        Returns:
            Normalised weights.

        Raises:
            ValueError: If the weights are empty, negative, or sum to zero.
        """
        if not weights:
            raise ValueError("At least one scoring weight is required")
        if any(value < 0 for value in weights.values()):
            raise ValueError("Scoring weights must be non-negative")
        total = float(sum(weights.values()))
        if total <= 0:
            raise ValueError("Scoring weights must sum to a positive value")
        return {key: float(value) / total for key, value in weights.items()}

    # ------------------------------------------------------------ components

    def score_volume(self, volume: Optional[float]) -> float:
        """Score cavity volume against the InsP3-InsP6 range.

        Full credit inside the optimum, with a smooth decay outside it. Decay is
        gentler above the range than below: a cavity larger than the ligand can
        still host it (often with ordered water), whereas one that is too small
        cannot.

        Args:
            volume: Cavity volume in Å³.

        Returns:
            A score in ``[0, 1]``.
        """
        if volume is None or not np.isfinite(volume):
            return self.parameters.missing_component_score
        volume = float(volume)
        low, high = self.parameters.volume_optimum_low, self.parameters.volume_optimum_high
        if low <= volume <= high:
            return 1.0
        if volume < low:
            return float(max(0.0, volume / low) ** 1.5)
        return float(np.exp(-(volume - high) / self.parameters.volume_tolerance))

    def score_depth(self, depth: Optional[float]) -> float:
        """Score geometric burial depth.

        Args:
            depth: Distance from the pocket centre to the nearest solvent-exposed
                atom, in Å.

        Returns:
            A score in ``[0, 1]``.
        """
        if depth is None or not np.isfinite(depth):
            return self.parameters.missing_component_score
        return _logistic(float(depth), self.parameters.depth_midpoint, self.parameters.depth_slope)

    def score_sasa(self, sasa: Optional[float]) -> float:
        """Score solvent accessibility; lower is better.

        Args:
            sasa: Mean SASA of pocket-lining residues in Å².

        Returns:
            A score in ``[0, 1]``.
        """
        if sasa is None or not np.isfinite(sasa):
            return self.parameters.missing_component_score
        return 1.0 - _logistic(
            float(sasa), self.parameters.sasa_midpoint, self.parameters.sasa_slope
        )

    def score_basic_residues(self, count: Optional[float]) -> float:
        """Score the basic-residue cluster.

        Args:
            count: Number of Arg/Lys/His residues lining the pocket.

        Returns:
            A score in ``[0, 1]``.
        """
        if count is None or not np.isfinite(count):
            return self.parameters.missing_component_score
        return _logistic(
            float(count), self.parameters.basic_midpoint, self.parameters.basic_slope
        )

    def score_electrostatics(self, potential: Optional[float]) -> float:
        """Score the electrostatic potential at the pocket centre.

        Args:
            potential: Potential in kT/e; ``None`` when not computed.

        Returns:
            A score in ``[0, 1]``, neutral when the potential is unavailable.
        """
        if potential is None or not np.isfinite(potential):
            return self.parameters.missing_component_score
        return _logistic(
            float(potential), self.parameters.potential_midpoint, self.parameters.potential_slope
        )

    def score_enclosure(self, enclosure: Optional[float]) -> float:
        """Score how completely the site is surrounded by protein.

        Args:
            enclosure: Fraction of directions blocked by protein atoms.

        Returns:
            A score in ``[0, 1]``.
        """
        if enclosure is None or not np.isfinite(enclosure):
            return self.parameters.missing_component_score
        return _logistic(
            float(enclosure),
            self.parameters.enclosure_midpoint,
            self.parameters.enclosure_slope,
        )

    # -------------------------------------------------------------- combined

    def component_scores(
        self,
        volume: Optional[float] = None,
        depth: Optional[float] = None,
        sasa: Optional[float] = None,
        basic_count: Optional[float] = None,
        potential: Optional[float] = None,
        enclosure: Optional[float] = None,
    ) -> Dict[str, float]:
        """Return every component score, for explanation and diagnostics.

        Args:
            volume: Cavity volume (Å³).
            depth: Geometric burial depth (Å).
            sasa: Mean lining-residue SASA (Å²).
            basic_count: Basic residues lining the pocket.
            potential: Electrostatic potential (kT/e).
            enclosure: Enclosure fraction.

        Returns:
            Component name to score.
        """
        return {
            "volume": self.score_volume(volume),
            "depth": self.score_depth(depth),
            "sasa": self.score_sasa(sasa),
            "basic_residues": self.score_basic_residues(basic_count),
            "electrostatics": self.score_electrostatics(potential),
            "enclosure": self.score_enclosure(enclosure),
        }

    def calculate_composite_score(
        self,
        volume: Optional[float] = None,
        depth: Optional[float] = None,
        sasa: Optional[float] = None,
        basic_count: Optional[float] = None,
        potential: Optional[float] = None,
        enclosure: Optional[float] = None,
        **_ignored: Any,
    ) -> float:
        """Weighted combination of the component scores.

        Components absent from the configured weights are excluded and the
        remaining weights are renormalised, so a caller who supplies no enclosure
        weight gets the same scale as one who does.

        Args:
            volume: Cavity volume (Å³).
            depth: Geometric burial depth (Å).
            sasa: Mean lining-residue SASA (Å²).
            basic_count: Basic residues lining the pocket.
            potential: Electrostatic potential (kT/e).
            enclosure: Enclosure fraction.
            **_ignored: Extra descriptors, accepted and ignored so callers can
                pass a full feature row.

        Returns:
            A composite score in ``[0, 1]``.
        """
        components = self.component_scores(
            volume=volume,
            depth=depth,
            sasa=sasa,
            basic_count=basic_count,
            potential=potential,
            enclosure=enclosure,
        )
        active = {name: weight for name, weight in self.weights.items() if name in components}
        total_weight = sum(active.values())
        if total_weight <= 0:
            return 0.0
        return float(
            sum(components[name] * weight for name, weight in active.items()) / total_weight
        )

    def score_frame(self, frame: "Any") -> "Any":
        """Score every row of a feature table.

        Args:
            frame: A DataFrame with descriptor columns. Both the descriptor
                names from :mod:`cryptic_ip.analysis.features` and the legacy
                names are recognised.

        Returns:
            An array of composite scores.
        """
        import pandas as pd

        def pick(row: Mapping[str, Any], *names: str) -> Optional[float]:
            for name in names:
                if name in row:
                    value = row[name]
                    if value is not None and np.isfinite(np.asarray(value, dtype=float)):
                        return float(value)
            return None

        scores = []
        for _, row in pd.DataFrame(frame).iterrows():
            scores.append(
                self.calculate_composite_score(
                    volume=pick(row, "pocket_volume", "volume", "hull_volume"),
                    depth=pick(row, "burial_depth", "pocket_depth", "depth"),
                    sasa=pick(row, "sasa_mean", "sasa"),
                    basic_count=pick(row, "n_basic_residues", "basic_residues"),
                    potential=pick(row, "electrostatic_potential", "coulomb_potential_kt"),
                    enclosure=pick(row, "enclosure"),
                )
            )
        return np.asarray(scores, dtype=float)

    def classify_site(self, score: float) -> str:
        """Map a composite score onto a confidence label.

        Args:
            score: Composite score.

        Returns:
            A human-readable label.
        """
        if score >= 0.75:
            return "High confidence cryptic IP site"
        if score >= 0.60:
            return "Moderate confidence candidate"
        if score >= 0.40:
            return "Low confidence - manual inspection recommended"
        return "Unlikely cryptic IP site"


class CalibratedRuleScorer:
    """Map the rule-based composite score onto a calibrated probability.

    The composite score is a weighted average of bounded components, so it is
    confined to a narrow band around 0.5 and cannot be read as a probability - a
    score of 0.7 does not mean a 70 % chance of a cryptic site. Fitting a
    one-dimensional logistic regression from composite score to label makes the
    output comparable with the machine-learned model's probabilities while
    keeping the underlying rule fully interpretable.

    Args:
        scorer: The underlying rule-based scorer.
    """

    def __init__(self, scorer: Optional[PocketScorer] = None) -> None:
        self.scorer = scorer or PocketScorer()
        self._model: Optional[Any] = None

    def fit(self, frame: "Any", labels: Sequence[int]) -> "CalibratedRuleScorer":
        """Fit the calibration curve.

        Args:
            frame: Feature table.
            labels: Binary labels.

        Returns:
            ``self``, for chaining.
        """
        from sklearn.linear_model import LogisticRegression

        raw = self.scorer.score_frame(frame).reshape(-1, 1)
        y = np.asarray(list(labels), dtype=int)
        self._model = LogisticRegression(class_weight="balanced", max_iter=1000)
        self._model.fit(raw, y)
        return self

    def predict_proba(self, frame: "Any") -> np.ndarray:
        """Return calibrated probabilities.

        Args:
            frame: Feature table.

        Returns:
            Probability of the positive class per row.

        Raises:
            RuntimeError: If the calibrator has not been fitted.
        """
        if self._model is None:
            raise RuntimeError("CalibratedRuleScorer is not fitted. Call fit() first.")
        raw = self.scorer.score_frame(frame).reshape(-1, 1)
        return self._model.predict_proba(raw)[:, 1]

    def calculate_composite_scores(self, samples: "Any") -> np.ndarray:
        """Alias of :meth:`predict_proba` for the scorer interface."""
        return self.predict_proba(samples)

    def classify_site(self, score: float) -> str:
        """Map a calibrated probability onto a confidence label."""
        if score >= 0.80:
            return "High confidence cryptic IP site"
        if score >= 0.60:
            return "Moderate confidence candidate"
        if score >= 0.40:
            return "Low confidence - manual inspection recommended"
        return "Unlikely cryptic IP site"
