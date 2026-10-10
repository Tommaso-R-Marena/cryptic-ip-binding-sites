"""Does this pocket prefer an inositol phosphate, or any polyanion at all?

Why this exists
---------------
This is the project's oldest unresolved problem and the one a reviewer reaches
first. Study C re-weighted the proteome ranking by P(IP | site) and did not
improve it - delta AUC -0.015 [-0.044, 0.015] - and 30 of the 75 candidates
are already flagged as explained by some other phosphate- or sulfate-rich
ligand. A professor's reaction to the candidate list was that many of these are
known to bind highly negatively charged ligands, which is exactly right and is
the screen working as specified rather than failing: it looks for pockets with
the physical character of a polyanion site, and it finds them.

The fix is not a better descriptor of IP6 sites. It is to stop asking an
absolute question. "Could this pocket bind a polyanion" has the answer yes for
every basic pocket in the proteome. The question that discriminates is
comparative: **given that this pocket binds something negative, would it
prefer an inositol phosphate over the competitors?**

So this module scores a pocket against the whole competitor panel and reports
a *margin*, with three outcomes rather than two:

* ``ip6_selective``        - IP6 fits better than every competitor by the margin
* ``polyanion_generic``    - some competitor is within the margin; the pocket is
                             a polyanion site and this method cannot say which
* ``competitor_preferred`` - a competitor fits better by the margin

The middle outcome is the one the deployed screen cannot express, and it is
probably the right answer for a large share of the 75. Reporting it is not a
weaker claim than a candidate list; it is the claim the data supports.

The margin is pre-registered
----------------------------
``DEFAULT_MARGIN`` is fixed before any candidate is scored and may not be
tuned on results, for the same reason the 2 A redocking criterion and the
5-group floor were fixed in advance. It is deliberately large: with a
descriptor nobody has validated, a small margin would manufacture selectivity
out of noise.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np

#: The ligand this project is about. Everything else in the panel is a
#: competitor, including IP3 - a pocket that cannot tell InsP6 from InsP3 is
#: an inositol-phosphate site but not an InsP6 site, and the ground truth is
#: 134 InsP6 entries against 1 InsP5, so that distinction is untested here.
TARGET = "IP6"

#: How much better IP6 has to fit than the best competitor before the pocket is
#: called selective. In units of the shape correlation, which runs -1 to 1.
#: Fixed in advance; see the module docstring.
DEFAULT_MARGIN = 0.15

#: Competitors whose presence in the panel is required for a selectivity claim
#: to mean anything, because they are the ones the project's own specificity
#: work could not separate from IP6.
REQUIRED_COMPETITORS = ("ATP", "PAPS", "UDP_GLUCOSE", "SULFATE", "PHOSPHATE")

SELECTIVE = "ip6_selective"
GENERIC = "polyanion_generic"
COMPETITOR = "competitor_preferred"
UNFAVOURABLE = "no_polyanion_site"


class SpecificityError(ValueError):
    """A panel from which no comparative statement can be made."""


@dataclass(frozen=True)
class SpecificityVerdict:
    """The comparative result for one pocket."""

    call: str
    margin: float
    target_shape: float
    target_depth_kt: float
    best_competitor: str
    best_competitor_shape: float
    ranking: Tuple[Tuple[str, float], ...]
    margin_threshold: float

    @property
    def is_candidate(self) -> bool:
        """Whether this pocket may be put forward as an IP6 candidate at all."""
        return self.call == SELECTIVE

    def explain(self) -> str:
        """One sentence a reader can check against the ranking."""
        if self.call == SELECTIVE:
            return (
                f"IP6 fits this pocket better than {self.best_competitor} by "
                f"{self.margin:.3f}, above the pre-registered margin of {self.margin_threshold}."
            )
        if self.call == GENERIC:
            return (
                f"{self.best_competitor} fits within {abs(self.margin):.3f} of IP6, inside the "
                f"pre-registered margin of {self.margin_threshold}: this is a polyanion pocket "
                "and this method cannot say which polyanion."
            )
        if self.call == COMPETITOR:
            return (
                f"{self.best_competitor} fits better than IP6 by {abs(self.margin):.3f}; on this "
                "evidence the pocket is not an IP6 site."
            )
        return (
            f"The field is unfavourable for a polyanion at all ({self.target_depth_kt:.2f} kT/e), "
            "so no comparison between polyanions is meaningful."
        )

    def as_dict(self) -> Dict[str, object]:
        return {
            "call": self.call,
            "margin": float(self.margin),
            "margin_threshold": float(self.margin_threshold),
            "target_shape": float(self.target_shape),
            "target_depth_kt": float(self.target_depth_kt),
            "best_competitor": self.best_competitor,
            "best_competitor_shape": float(self.best_competitor_shape),
            "ranking": [[name, float(value)] for name, value in self.ranking],
            "is_candidate": self.is_candidate,
            "explanation": self.explain(),
        }


def verdict(
    per_ligand: Mapping[str, Mapping[str, float]],
    *,
    margin: float = DEFAULT_MARGIN,
    target: str = TARGET,
    require: Sequence[str] = REQUIRED_COMPETITORS,
    min_depth_kt: float = 0.0,
) -> SpecificityVerdict:
    """Compare the target against every competitor in ``per_ligand``.

    ``per_ligand`` is what :func:`cryptic_ip.detect.field.describe_pocket`
    produces: a mapping from ligand name to a dict holding ``shape`` and
    ``depth_kt``.
    """
    if target not in per_ligand:
        raise SpecificityError(f"the panel has no entry for the target {target!r}")
    absent = [name for name in require if name not in per_ligand]
    if absent:
        raise SpecificityError(
            f"a selectivity claim needs the competitors this project could not separate from "
            f"{target}; missing {absent}"
        )
    if margin <= 0:
        raise SpecificityError(f"margin must be positive, got {margin}")

    target_entry = per_ligand[target]
    target_shape = float(target_entry["shape"])
    target_depth = float(target_entry["depth_kt"])

    competitors = {
        name: float(entry["shape"]) for name, entry in per_ligand.items() if name != target
    }
    if not competitors:
        raise SpecificityError("a comparative statement needs at least one competitor")
    best_name, best_shape = max(competitors.items(), key=lambda item: item[1])
    gap = target_shape - best_shape

    ranking = tuple(
        sorted(
            [(target, target_shape)] + list(competitors.items()),
            key=lambda item: -item[1],
        )
    )

    if target_depth <= min_depth_kt:
        call = UNFAVOURABLE
    elif gap >= margin:
        call = SELECTIVE
    elif gap <= -margin:
        call = COMPETITOR
    else:
        call = GENERIC
    return SpecificityVerdict(
        call=call,
        margin=float(gap),
        target_shape=target_shape,
        target_depth_kt=target_depth,
        best_competitor=best_name,
        best_competitor_shape=float(best_shape),
        ranking=ranking,
        margin_threshold=float(margin),
    )


@dataclass
class PanelSummary:
    """What a whole candidate list looks like once the margin is applied."""

    counts: Dict[str, int]
    by_competitor: Dict[str, int]
    n_pockets: int

    @property
    def selective_fraction(self) -> float:
        return self.counts.get(SELECTIVE, 0) / self.n_pockets if self.n_pockets else 0.0

    def as_dict(self) -> Dict[str, object]:
        return {
            "n_pockets": self.n_pockets,
            "counts": dict(self.counts),
            "by_competitor": dict(self.by_competitor),
            "selective_fraction": self.selective_fraction,
        }

    def summary(self) -> str:
        lines = [f"{self.n_pockets} pocket(s) scored against the competitor panel:"]
        for call in (SELECTIVE, GENERIC, COMPETITOR, UNFAVOURABLE):
            if self.counts.get(call):
                lines.append(f"  {call}: {self.counts[call]}")
        if self.by_competitor:
            worst = sorted(self.by_competitor.items(), key=lambda item: -item[1])
            named = ", ".join(f"{name} ({count})" for name, count in worst[:5])
            lines.append(f"  competitors that matched or beat IP6: {named}")
        return "\n".join(lines)


def summarise(verdicts: Sequence[SpecificityVerdict]) -> PanelSummary:
    """Aggregate verdicts, keeping which competitor did the explaining.

    The per-competitor tally is the useful part: "30 of 75 are explained by
    another ligand" is a weaker statement than naming which ligand explains how
    many, and the latter is what tells you whether the screen is finding ATP
    sites or sulfate sites.
    """
    counts: Dict[str, int] = {}
    by_competitor: Dict[str, int] = {}
    for item in verdicts:
        counts[item.call] = counts.get(item.call, 0) + 1
        if item.call in (GENERIC, COMPETITOR):
            by_competitor[item.best_competitor] = by_competitor.get(item.best_competitor, 0) + 1
    return PanelSummary(counts=counts, by_competitor=by_competitor, n_pockets=len(verdicts))


def margin_sweep(
    verdicts: Sequence[SpecificityVerdict],
    margins: Sequence[float] = (0.05, 0.10, 0.15, 0.20, 0.30),
) -> List[Dict[str, object]]:
    """How the call counts move with the margin, recomputed from stored gaps.

    Published next to the result, not instead of it. The pre-registered margin
    is the one that decides; the sweep exists so a reader can see how much the
    answer depends on it, which is the check nobody can run on a tuned
    threshold.
    """
    out: List[Dict[str, object]] = []
    gaps = np.array([item.margin for item in verdicts], dtype=float)
    unfavourable = np.array([item.call == UNFAVOURABLE for item in verdicts], dtype=bool)
    for margin in margins:
        selective = int(np.sum((gaps >= margin) & ~unfavourable))
        competitor = int(np.sum((gaps <= -margin) & ~unfavourable))
        generic = int(np.sum((np.abs(gaps) < margin) & ~unfavourable))
        out.append({
            "margin": float(margin),
            SELECTIVE: selective,
            GENERIC: generic,
            COMPETITOR: competitor,
            UNFAVOURABLE: int(np.sum(unfavourable)),
        })
    return out
