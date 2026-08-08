"""Pocket labelling from experimentally observed ligand positions.

The defect this replaces
------------------------
The previous labelling rule was:

1. Classify the whole *structure* from a SASA sum over all ligand copies.
2. If the structure was ``Cryptic`` or ``Semi-cryptic``, label pockets within 8 Å
   of the ligand centroid positive; otherwise label **every** pocket negative.

Both steps were unsound, and the consequences were measurable in the shipped
artefacts (12 190 pockets, 5 positives, AUROC 0.50, AUPRC 0.0002):

* Because the structure-level class came from a copy-summed SASA, nearly every
  entry landed in ``Surface``, so almost every genuine inositol phosphate site in
  the dataset was labelled **negative**. The training signal was inverted for the
  majority of the positive evidence available.
* Centroid-to-centroid distance is a poor overlap test. InsP6 spans about 11 Å,
  so its centroid can sit more than 8 Å from the centre of the very pocket that
  holds it, while an adjacent unrelated pocket can fall inside 8 Å.
* Pockets that neither contain the ligand nor are safely far from it were labelled
  negative, injecting label noise exactly where the decision boundary lies.

What this module does instead
-----------------------------
Labelling is separated into two questions that the previous scheme conflated:

**Question 1 - is this pocket the ligand site?** Decided by *atom-level overlap*:
the fraction of ligand heavy atoms lying within a contact distance of the
pocket's alpha spheres (falling back to lining atoms). Overlap above
``positive_overlap`` is a site; below ``negative_overlap`` is a decoy; in between
is ``ambiguous`` and is **excluded from training** rather than being called a
negative.

**Question 2 - is that site cryptic?** Decided separately from the per-copy
burial measurements in :mod:`cryptic_ip.validation.burial_metrics`.

Keeping the two apart means the IP-binding classifier trains on hundreds of real
positives instead of a handful, and the burial question is answered by a direct
measurement instead of being folded into a class label.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

LOGGER = logging.getLogger(__name__)

#: Heavy-atom contact distance for overlap counting (Å). A ligand atom within
#: this distance of a pocket alpha sphere is considered to occupy the pocket.
OVERLAP_CONTACT_DISTANCE = 4.0

#: Ligand-atom overlap fraction at or above which a pocket is the ligand site.
POSITIVE_OVERLAP = 0.30

#: Ligand-atom overlap fraction at or below which a pocket is a decoy. The gap
#: between this and :data:`POSITIVE_OVERLAP` is the ambiguity zone.
NEGATIVE_OVERLAP = 0.05


class PocketLabel(Enum):
    """Training label for a pocket."""

    NEGATIVE = 0
    POSITIVE = 1
    AMBIGUOUS = -1

    @property
    def is_trainable(self) -> bool:
        """Whether the label may be used for supervised training."""
        return self is not PocketLabel.AMBIGUOUS


@dataclass
class PocketAssignment:
    """Result of matching one pocket against the ligand copies of a structure.

    Attributes:
        pocket_id: Pocket identifier.
        label: Training label.
        overlap_fraction: Best ligand-atom overlap fraction across copies.
        matched_instance_id: Identifier of the best-matching ligand copy.
        matched_comp_id: Component identifier of that copy.
        centroid_distance: Distance from the pocket centre to that copy's centroid.
        min_atom_distance: Closest pocket-to-ligand heavy-atom distance.
        rank_within_structure: 0 for the best-matching pocket of the structure.
        reason: Human-readable explanation of the label.
    """

    pocket_id: int
    label: PocketLabel
    overlap_fraction: float
    matched_instance_id: Optional[str] = None
    matched_comp_id: Optional[str] = None
    centroid_distance: float = float("nan")
    min_atom_distance: float = float("nan")
    rank_within_structure: int = -1
    reason: str = ""

    def to_row(self) -> Dict[str, object]:
        """Return a flat mapping suitable for a DataFrame row."""
        return {
            "pocket_id": self.pocket_id,
            "label": self.label.value,
            "label_name": self.label.name.lower(),
            "overlap_fraction": self.overlap_fraction,
            "matched_instance_id": self.matched_instance_id,
            "matched_comp_id": self.matched_comp_id,
            "ligand_centroid_distance": self.centroid_distance,
            "ligand_min_atom_distance": self.min_atom_distance,
            "site_rank_within_structure": self.rank_within_structure,
            "label_reason": self.reason,
        }


@dataclass
class LigandSite:
    """Coordinates of one ligand copy for labelling purposes.

    Attributes:
        instance_id: Stable identifier, e.g. ``"IHP_A_501"``.
        comp_id: Chemical component identifier.
        coords: Heavy-atom coordinates, shape ``(n, 3)``.
        burial_class: Burial class of the copy, when known.
    """

    instance_id: str
    comp_id: str
    coords: np.ndarray
    burial_class: str = "unknown"

    @property
    def centroid(self) -> np.ndarray:
        """Centroid of the copy's heavy atoms."""
        return np.asarray(self.coords, dtype=float).mean(axis=0)


def ligand_overlap_fraction(
    ligand_coords: np.ndarray,
    pocket_points: np.ndarray,
    *,
    contact_distance: float = OVERLAP_CONTACT_DISTANCE,
) -> Tuple[float, float]:
    """Fraction of ligand atoms in contact with a pocket, and the closest distance.

    Using the fraction of *ligand atoms* covered - rather than a centre-to-centre
    distance - makes the test insensitive to ligand size and to how the detector
    happened to split a large cavity, because it asks the physical question: does
    this pocket enclose the ligand?

    Args:
        ligand_coords: Ligand heavy-atom coordinates, shape ``(n, 3)``.
        pocket_points: Alpha-sphere centres or lining atoms, shape ``(m, 3)``.
        contact_distance: Contact cutoff in Å.

    Returns:
        ``(overlap_fraction, min_distance)``. Empty input gives ``(0.0, nan)``.
    """
    from scipy.spatial import cKDTree

    ligand_coords = np.asarray(ligand_coords, dtype=float)
    pocket_points = np.asarray(pocket_points, dtype=float)
    if ligand_coords.size == 0 or pocket_points.size == 0:
        return 0.0, float("nan")

    distances, _ = cKDTree(pocket_points).query(ligand_coords, k=1)
    return float(np.mean(distances <= float(contact_distance))), float(np.min(distances))


def assign_pocket_labels(
    pockets: Sequence[Tuple[int, np.ndarray, np.ndarray]],
    ligand_sites: Sequence[LigandSite],
    *,
    contact_distance: float = OVERLAP_CONTACT_DISTANCE,
    positive_overlap: float = POSITIVE_OVERLAP,
    negative_overlap: float = NEGATIVE_OVERLAP,
    require_cryptic: bool = False,
) -> List[PocketAssignment]:
    """Label every pocket of one structure against its ligand copies.

    Args:
        pockets: ``(pocket_id, centre, pocket_points)`` per pocket, where
            ``pocket_points`` are alpha-sphere centres or lining atoms.
        ligand_sites: Ligand copies present in the structure. An empty sequence
            means the structure is a decoy: every pocket becomes a negative,
            which is exactly the intended use of IP-free entries.
        contact_distance: Heavy-atom contact cutoff (Å).
        positive_overlap: Overlap fraction at or above which a pocket is positive.
        negative_overlap: Overlap fraction at or below which a pocket is negative.
        require_cryptic: When ``True``, only copies whose ``burial_class`` is
            ``cryptic`` produce positives, and matches to non-cryptic copies
            become ``AMBIGUOUS`` rather than negative - a surface IP site is not
            evidence of *absence* of IP binding, so calling it negative would
            teach the model the wrong thing.

    Returns:
        One :class:`PocketAssignment` per input pocket, in input order.
    """
    if positive_overlap <= negative_overlap:
        raise ValueError("positive_overlap must exceed negative_overlap")

    assignments: List[PocketAssignment] = []
    for pocket_id, centre, points in pockets:
        centre_arr = np.asarray(centre, dtype=float)
        if not ligand_sites:
            assignments.append(
                PocketAssignment(
                    pocket_id=int(pocket_id),
                    label=PocketLabel.NEGATIVE,
                    overlap_fraction=0.0,
                    reason="structure contains no inositol phosphate ligand",
                )
            )
            continue

        best_overlap = -1.0
        best_site: Optional[LigandSite] = None
        best_min_distance = float("nan")
        for site in ligand_sites:
            overlap, min_distance = ligand_overlap_fraction(
                site.coords, points, contact_distance=contact_distance
            )
            if overlap > best_overlap:
                best_overlap = overlap
                best_site = site
                best_min_distance = min_distance

        assert best_site is not None
        centroid_distance = float(np.linalg.norm(best_site.centroid - centre_arr))

        if best_overlap >= positive_overlap:
            if require_cryptic and best_site.burial_class != "cryptic":
                label = PocketLabel.AMBIGUOUS
                reason = (
                    f"ligand site but burial class {best_site.burial_class!r} is not cryptic; "
                    "excluded rather than treated as a non-binding pocket"
                )
            else:
                label = PocketLabel.POSITIVE
                reason = f"{best_overlap:.0%} of ligand atoms within {contact_distance:.1f} A"
        elif best_overlap <= negative_overlap:
            label = PocketLabel.NEGATIVE
            reason = f"only {best_overlap:.0%} ligand-atom overlap"
        else:
            label = PocketLabel.AMBIGUOUS
            reason = (
                f"partial overlap {best_overlap:.0%} between "
                f"{negative_overlap:.0%} and {positive_overlap:.0%}"
            )

        assignments.append(
            PocketAssignment(
                pocket_id=int(pocket_id),
                label=label,
                overlap_fraction=float(best_overlap),
                matched_instance_id=best_site.instance_id,
                matched_comp_id=best_site.comp_id,
                centroid_distance=centroid_distance,
                min_atom_distance=best_min_distance,
                reason=reason,
            )
        )

    # Rank pockets by overlap so the primary site of each structure is identifiable.
    order = sorted(
        range(len(assignments)), key=lambda i: -assignments[i].overlap_fraction
    )
    for rank, index in enumerate(order):
        assignments[index].rank_within_structure = rank
    return assignments


@dataclass
class LabelSummary:
    """Counts describing a labelled dataset.

    Attributes:
        n_pockets: Total pockets considered.
        n_positive: Pockets labelled positive.
        n_negative: Pockets labelled negative.
        n_ambiguous: Pockets excluded as ambiguous.
        n_structures: Structures contributing pockets.
        n_structures_with_site: Structures with at least one positive pocket.
        n_structures_with_ligand_but_no_site: Structures that contain a ligand yet
            yielded no positive pocket - a detector-recall failure that must be
            reported rather than silently absorbed into the negatives.
    """

    n_pockets: int = 0
    n_positive: int = 0
    n_negative: int = 0
    n_ambiguous: int = 0
    n_structures: int = 0
    n_structures_with_site: int = 0
    n_structures_with_ligand_but_no_site: int = 0
    per_structure: Dict[str, Dict[str, int]] = field(default_factory=dict)

    @property
    def positive_rate(self) -> float:
        """Fraction of trainable pockets that are positive."""
        trainable = self.n_positive + self.n_negative
        return self.n_positive / trainable if trainable else float("nan")

    @property
    def site_recall(self) -> float:
        """Fraction of ligand-bearing structures where a pocket matched the site.

        This is the pocket detector's recall on known sites. It bounds the whole
        pipeline: a site the detector never proposes cannot be scored, however
        good the classifier is.
        """
        with_ligand = self.n_structures_with_site + self.n_structures_with_ligand_but_no_site
        return self.n_structures_with_site / with_ligand if with_ligand else float("nan")

    def to_dict(self) -> Dict[str, object]:
        """Return a JSON-serialisable representation."""
        return {
            "n_pockets": self.n_pockets,
            "n_positive": self.n_positive,
            "n_negative": self.n_negative,
            "n_ambiguous": self.n_ambiguous,
            "n_structures": self.n_structures,
            "n_structures_with_site": self.n_structures_with_site,
            "n_structures_with_ligand_but_no_site": self.n_structures_with_ligand_but_no_site,
            "positive_rate": self.positive_rate,
            "site_recall": self.site_recall,
        }


def summarise_labels(
    per_structure_assignments: Dict[str, Sequence[PocketAssignment]],
    structures_with_ligand: Sequence[str] = (),
) -> LabelSummary:
    """Aggregate label counts across structures.

    Args:
        per_structure_assignments: Assignments keyed by structure identifier.
        structures_with_ligand: Identifiers known to contain an inositol
            phosphate, used to compute detector recall on known sites.

    Returns:
        The aggregated summary.
    """
    summary = LabelSummary()
    with_ligand = {str(identifier) for identifier in structures_with_ligand}

    for identifier, assignments in per_structure_assignments.items():
        counts = {"positive": 0, "negative": 0, "ambiguous": 0}
        for assignment in assignments:
            counts[assignment.label.name.lower()] += 1
        summary.n_pockets += len(assignments)
        summary.n_positive += counts["positive"]
        summary.n_negative += counts["negative"]
        summary.n_ambiguous += counts["ambiguous"]
        summary.n_structures += 1
        summary.per_structure[str(identifier)] = counts
        if counts["positive"]:
            summary.n_structures_with_site += 1
        elif str(identifier) in with_ligand:
            summary.n_structures_with_ligand_but_no_site += 1

    return summary
