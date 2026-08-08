"""Select the pocket that actually holds the ligand in a control structure.

Why this exists
---------------
The control validators previously chose the pocket whose **centre lay closest to
the ligand centroid**. On the ADAR2 crystal structure that picks the wrong
pocket: the calibration diagnostic reports the pocket with 100 % ligand-atom
overlap scoring 0.644, while the validator was assessing a different, smaller
pocket scoring 0.432. The gate was therefore grading the pipeline on a pocket
that does not contain the inositol phosphate.

Centre-to-centroid distance fails here for the same reason it fails as a
labelling rule: InsP6 spans roughly 11 A, so a small satellite pocket can have
its centre nearer the ligand centroid than the large cavity that encloses the
ligand does.

This module reuses the atom-level overlap criterion from
:mod:`cryptic_ip.analysis.labeling`, so control validation and training labels
answer "which pocket holds the ligand" the same way.
"""

from __future__ import annotations

import logging
from typing import Any, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..analysis.labeling import OVERLAP_CONTACT_DISTANCE, ligand_overlap_fraction

LOGGER = logging.getLogger(__name__)


def select_ligand_pocket(
    scored: pd.DataFrame,
    analyzer: Any,
    ligand_coords: Optional[np.ndarray],
    *,
    contact_distance: float = OVERLAP_CONTACT_DISTANCE,
    fallback_site_residues: Optional[Sequence[int]] = None,
) -> Tuple[pd.Series, float]:
    """Return the scored pocket with the greatest ligand-atom overlap.

    Args:
        scored: Scored pockets, one row per pocket.
        analyzer: A :class:`~cryptic_ip.analysis.analyzer.ProteinAnalyzer` used to
            retrieve alpha-sphere coordinates and pocket-lining residues.
        ligand_coords: Ligand heavy-atom coordinates, shape ``(n, 3)``. When
            ``None`` the residue-overlap fallback is used.
        contact_distance: Heavy-atom contact cutoff in A.
        fallback_site_residues: Reference site residue numbers, used only when no
            ligand coordinates are available.

    Returns:
        ``(row, overlap_fraction)``. ``overlap_fraction`` is ``nan`` when the
        selection came from the fallback path.

    Raises:
        ValueError: If ``scored`` is empty.
    """
    if scored.empty:
        raise ValueError("No scored pockets available for site selection")

    if ligand_coords is None or len(ligand_coords) == 0:
        return _select_by_residue_overlap(scored, analyzer, fallback_site_residues)

    ligand_coords = np.asarray(ligand_coords, dtype=float)
    best_row = scored.iloc[0]
    best_overlap = -1.0
    best_distance = float("inf")

    for _, row in scored.iterrows():
        pocket_id = int(row["pocket_id"])
        points = _pocket_points(analyzer, row, pocket_id)
        if points is None or len(points) == 0:
            continue
        overlap, min_distance = ligand_overlap_fraction(
            ligand_coords, points, contact_distance=contact_distance
        )
        # Overlap first; the closest approach breaks ties between pockets that
        # both fully cover a small ligand.
        if overlap > best_overlap or (
            overlap == best_overlap and min_distance < best_distance
        ):
            best_overlap = overlap
            best_distance = min_distance
            best_row = row

    if best_overlap <= 0.0:
        LOGGER.warning(
            "No pocket overlaps the ligand; falling back to residue overlap. "
            "This means the pocket detector did not propose the known site."
        )
        return _select_by_residue_overlap(scored, analyzer, fallback_site_residues)

    return best_row, float(best_overlap)


def _pocket_points(analyzer: Any, row: pd.Series, pocket_id: int) -> Optional[np.ndarray]:
    """Return alpha-sphere coordinates for a pocket, falling back to its centre."""
    getter = getattr(analyzer, "_pocket_alpha_spheres", None)
    if getter is not None:
        try:
            spheres = getter(pocket_id)
        except Exception:  # noqa: BLE001 - fall back to the centre
            spheres = None
        if spheres is not None and len(spheres):
            return np.asarray(spheres, dtype=float)

    centre = row.get("center")
    if centre is not None and not isinstance(centre, float):
        return np.asarray(centre, dtype=float).reshape(1, 3)

    pockets = getattr(analyzer, "pockets", None)
    if pockets is None or pockets.empty:
        return None
    match = pockets[pockets["pocket_id"] == pocket_id]
    if match.empty:
        return None
    pocket = match.iloc[0]
    return np.asarray(
        [[pocket["center_x"], pocket["center_y"], pocket["center_z"]]], dtype=float
    )


def _select_by_residue_overlap(
    scored: pd.DataFrame, analyzer: Any, site_residues: Optional[Sequence[int]]
) -> Tuple[pd.Series, float]:
    """Fallback: pick the pocket sharing the most residues with a reference site.

    Used when the structure has no modelled ligand, so overlap cannot be
    measured. Highest composite score wins when no reference site is known.
    """
    if not site_residues:
        best = scored.sort_values("composite_score", ascending=False).iloc[0]
        return best, float("nan")

    wanted = {int(r) for r in site_residues}
    best_row = scored.iloc[0]
    best_key = (-1, -1.0)
    for _, row in scored.iterrows():
        pocket_id = int(row["pocket_id"])
        try:
            residues = set(analyzer.get_pocket_residues(pocket_id, distance_cutoff=8.0))
        except Exception:  # noqa: BLE001
            continue
        key = (len(residues & wanted), float(row["composite_score"]))
        if key > best_key:
            best_key = key
            best_row = row
    return best_row, float("nan")
