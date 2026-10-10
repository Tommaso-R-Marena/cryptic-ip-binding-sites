"""An in-process cavity detector, so a pocket can be found without fpocket.

Why this exists
---------------
Two reasons, one practical and one scientific.

Practically, `cryptic-ip analyze` and the whole docking path raise a
``RuntimeError`` without fpocket on the PATH. fpocket is a conda package, not a
pip one, so every environment that is built from `requirements.txt` alone
cannot detect a pocket at all. A detector in process removes that cliff.

Scientifically, the reason is the one that matters. Finding a *cryptic* pocket
means detecting pockets in many conformers and asking which ones appear and
disappear. That is hundreds of detection runs per protein, and shelling out to
an external binary for each one is the difference between a proteome-scale
ensemble search and a pilot.

The method is LIGSITE's protein-solvent-protein count, which is old, simple
and well understood: lay a grid over the structure, mark the free points, and
for each free point count the directions in which protein is met on both
sides. A point enclosed in many directions is inside a cavity. Nothing here is
novel and that is deliberate - the novelty in this package is what gets built
on top in ``cryptic_ip.detect.ensemble``, and a new cavity detector underneath
it would just be a second thing to validate.

What it is not
--------------
It is not fpocket and it does not reproduce fpocket's volumes. fpocket fits
alpha spheres to a Voronoi tessellation; this counts enclosed grid points.
The two agree about where large cavities are and disagree about their extent,
so a volume from here must not be compared with a volume from there, and
``docs/DETECTOR_PLAN.md`` pre-registers the agreement study rather than
assuming it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage

#: Grid spacing in angstrom. 0.8 is finer than LIGSITE's original 1.0 because
#: an inositol phosphate's phosphates are 2.5 Å apart and a coarse grid merges
#: the sub-pockets between them.
DEFAULT_SPACING = 0.8

#: Solvent probe radius in angstrom: water.
PROBE_RADIUS = 1.4

#: The seven scan directions: three axes and four body diagonals of the cube.
#: LIGSITE uses exactly these; adding the twelve face diagonals changes the
#: counts but not the ranking, and costs twice as much.
SCAN_DIRECTIONS: Tuple[Tuple[int, int, int], ...] = (
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
)

#: How many of the seven directions must be enclosed for a point to count as
#: pocket. 5 of 7 is LIGSITE's own recommendation for buried sites; lower
#: values pick up surface grooves, which is the failure mode this project
#: already has.
DEFAULT_MIN_ENCLOSURE = 5

#: How far to look for the enclosing protein, in angstrom. Beyond this a point
#: is "open" in that direction. 12 Å is about the width of a pocket that could
#: hold an inositol phosphate plus its shell of water.
DEFAULT_SCAN_REACH = 12.0

#: Smallest cavity reported, in cubic angstrom. An IP6 displaces roughly
#: 350 Å^3, so anything under 50 cannot be a site for one and is noise.
MIN_CAVITY_VOLUME = 50.0

#: Default van der Waals radii in angstrom, by element.
VDW_RADII: Dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "MG": 1.73, "ZN": 1.39, "CA": 2.31, "FE": 2.00, "MN": 2.05, "NA": 2.27, "K": 2.75,
}
DEFAULT_RADIUS = 1.70


class CavityError(ValueError):
    """A structure from which no cavity map can be built."""


@dataclass(frozen=True)
class Cavity:
    """One detected cavity."""

    centre: Tuple[float, float, float]
    volume: float
    n_points: int
    mean_enclosure: float
    max_enclosure: int
    depth: float

    def as_dict(self) -> Dict[str, object]:
        return {
            "centre": list(self.centre),
            "volume": float(self.volume),
            "n_points": int(self.n_points),
            "mean_enclosure": float(self.mean_enclosure),
            "max_enclosure": int(self.max_enclosure),
            "depth": float(self.depth),
        }


@dataclass(frozen=True)
class CavityMap:
    """Every cavity in one conformer, with the grid it was found on."""

    cavities: Tuple[Cavity, ...]
    spacing: float
    origin: Tuple[float, float, float]
    shape: Tuple[int, int, int]
    min_enclosure: int

    @property
    def total_volume(self) -> float:
        return float(sum(c.volume for c in self.cavities))

    def largest(self) -> Optional[Cavity]:
        return max(self.cavities, key=lambda c: c.volume) if self.cavities else None

    def as_dict(self) -> Dict[str, object]:
        return {
            "n_cavities": len(self.cavities),
            "total_volume": self.total_volume,
            "spacing": self.spacing,
            "min_enclosure": self.min_enclosure,
            "cavities": [c.as_dict() for c in self.cavities],
        }


def atom_radii(elements: Sequence[str]) -> np.ndarray:
    """Van der Waals radii for ``elements``, defaulting to carbon."""
    return np.array(
        [VDW_RADII.get(str(element).upper(), DEFAULT_RADIUS) for element in elements],
        dtype=float,
    )


def occupancy_grid(
    coords: np.ndarray,
    radii: np.ndarray,
    *,
    spacing: float = DEFAULT_SPACING,
    probe: float = PROBE_RADIUS,
    padding: float = 4.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """``(occupied, origin)``: which grid points lie inside an atom plus probe.

    The grid is axis aligned and its origin is the structure's lower corner
    minus ``padding``, so the whole molecule plus a shell of solvent fits.
    """
    coords = np.asarray(coords, dtype=float)
    radii = np.asarray(radii, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise CavityError("coordinates must have shape (n, 3)")
    if radii.shape != (len(coords),):
        raise CavityError("one radius per atom")
    if spacing <= 0:
        raise CavityError(f"spacing must be positive, got {spacing}")
    if len(coords) == 0:
        raise CavityError("no atoms")

    origin = coords.min(axis=0) - padding
    extent = coords.max(axis=0) + padding - origin
    shape = tuple(int(np.ceil(e / spacing)) + 1 for e in extent)
    if min(shape) < 3:
        raise CavityError(f"the grid {shape} is too small to scan")
    total = int(np.prod(shape))
    if total > 80_000_000:
        raise CavityError(
            f"a {shape} grid is {total / 1e6:.0f} million points; coarsen the spacing or "
            "crop to a region rather than gridding the whole structure"
        )

    occupied = np.zeros(shape, dtype=bool)
    inflated = radii + probe
    # Stamp each atom into the grid over its own bounding box only.
    for index in range(len(coords)):
        reach = inflated[index]
        lower = np.maximum(((coords[index] - reach - origin) / spacing).astype(int), 0)
        upper = np.minimum(
            np.ceil((coords[index] + reach - origin) / spacing).astype(int) + 1,
            np.array(shape),
        )
        if np.any(lower >= upper):
            continue
        axes = [np.arange(lower[d], upper[d]) * spacing + origin[d] for d in range(3)]
        dx = axes[0][:, None, None] - coords[index, 0]
        dy = axes[1][None, :, None] - coords[index, 1]
        dz = axes[2][None, None, :] - coords[index, 2]
        inside = (dx * dx + dy * dy + dz * dz) <= reach * reach
        block = occupied[lower[0]:upper[0], lower[1]:upper[1], lower[2]:upper[2]]
        np.logical_or(block, inside, out=block)
    return occupied, origin


def enclosure_count(
    occupied: np.ndarray,
    *,
    spacing: float = DEFAULT_SPACING,
    reach: float = DEFAULT_SCAN_REACH,
) -> np.ndarray:
    """For each free point, in how many of the seven directions it is enclosed.

    A direction counts when protein is met within ``reach`` going both ways
    along it. This is LIGSITE's protein-solvent-protein event count.
    """
    if occupied.ndim != 3:
        raise CavityError("the occupancy grid must be 3-d")
    steps = max(1, int(np.ceil(reach / spacing)))
    counts = np.zeros(occupied.shape, dtype=np.int8)
    for direction in SCAN_DIRECTIONS:
        forward = np.zeros(occupied.shape, dtype=bool)
        backward = np.zeros(occupied.shape, dtype=bool)
        for step in range(1, steps + 1):
            shift = tuple(step * d for d in direction)
            forward |= _shifted(occupied, shift)
            backward |= _shifted(occupied, tuple(-s for s in shift))
        counts += (forward & backward).astype(np.int8)
    counts[occupied] = 0
    return counts


def _shifted(array: np.ndarray, shift: Sequence[int]) -> np.ndarray:
    """``array`` translated by ``shift``, with False shifted in from outside.

    Padding with False - rather than wrapping, which `np.roll` would do - is
    what makes a point near the box edge read as open in that direction
    instead of being enclosed by the far side of the protein.
    """
    out = np.zeros_like(array)
    # A shift at least as large as the axis moves everything off the grid, so
    # the answer is all-False. Without this the slice pair below becomes
    # inconsistent - one side empty, the other not - and numpy raises. It
    # happens on any structure thinner than the scan reach, such as a slab.
    if any(abs(offset) >= array.shape[axis] for axis, offset in enumerate(shift)):
        return out
    src = [slice(None)] * 3
    dst = [slice(None)] * 3
    for axis, offset in enumerate(shift):
        if offset > 0:
            src[axis] = slice(0, array.shape[axis] - offset)
            dst[axis] = slice(offset, array.shape[axis])
        elif offset < 0:
            src[axis] = slice(-offset, array.shape[axis])
            dst[axis] = slice(0, array.shape[axis] + offset)
    out[tuple(dst)] = array[tuple(src)]
    return out


def detect_cavities(
    coords: np.ndarray,
    elements: Sequence[str],
    *,
    spacing: float = DEFAULT_SPACING,
    min_enclosure: int = DEFAULT_MIN_ENCLOSURE,
    reach: float = DEFAULT_SCAN_REACH,
    min_volume: float = MIN_CAVITY_VOLUME,
) -> CavityMap:
    """Find the cavities of one conformer."""
    radii = atom_radii(elements)
    occupied, origin = occupancy_grid(coords, radii, spacing=spacing)
    counts = enclosure_count(occupied, spacing=spacing, reach=reach)
    pocket = counts >= min_enclosure
    if not np.any(pocket):
        return CavityMap((), spacing, tuple(float(v) for v in origin), occupied.shape, min_enclosure)

    labels, n_labels = ndimage.label(pocket)
    voxel = spacing ** 3
    cavities: List[Cavity] = []
    for label in range(1, n_labels + 1):
        mask = labels == label
        n_points = int(np.sum(mask))
        volume = n_points * voxel
        if volume < min_volume:
            continue
        indices = np.argwhere(mask)
        centre = indices.mean(axis=0) * spacing + origin
        enclosures = counts[mask]
        # Depth: how far the cavity's deepest point is from the nearest open
        # (non-pocket, non-protein) point. A groove has a shallow depth even
        # when it is large.
        distance = ndimage.distance_transform_edt(mask) * spacing
        cavities.append(Cavity(
            centre=(float(centre[0]), float(centre[1]), float(centre[2])),
            volume=float(volume),
            n_points=n_points,
            mean_enclosure=float(np.mean(enclosures)),
            max_enclosure=int(np.max(enclosures)),
            depth=float(np.max(distance)),
        ))
    cavities.sort(key=lambda c: -c.volume)
    return CavityMap(
        tuple(cavities), spacing, tuple(float(v) for v in origin), occupied.shape, min_enclosure,
    )
