"""Vectorised molecular geometry primitives: SASA, burial depth and enclosure.

Why a local implementation
--------------------------
The pipeline needs solvent accessible surface area for *arbitrary atom subsets
evaluated in arbitrary molecular contexts*: the same ligand atoms are measured
inside the protein complex and again in isolation, so that a size-normalised
burial fraction can be formed. Whole-structure SASA helpers cannot express that,
because the atom set that *occludes* and the atom set that is *measured* must be
chosen independently.

The implementation below is a direct, vectorised Shrake-Rupley numerical
integration (Shrake & Rupley, *J Mol Biol* 79:351-371, 1973):

1. Each atom is assigned a van der Waals radius, expanded by the solvent probe
   radius (1.4 Å, a water molecule).
2. A quasi-uniform point set (Fibonacci sphere) is placed on the expanded
   sphere of each atom.
3. A point is buried when it lies inside any neighbouring atom's expanded
   sphere. SASA is the accessible fraction times the sphere area.

The accessible fraction converges as ``O(1/sqrt(n_points))``; the default of 512
points per atom gives per-atom SASA reproducible to well under 1 Å², which is
far below the differences the classifier resolves. Neighbour queries use a
k-d tree, so cost is linear in atom count for typical protein densities.

All functions are pure and deterministic: identical inputs give bit-identical
outputs, with no random number generator involved.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

#: Solvent probe radius in Å (water).
PROBE_RADIUS = 1.4

#: Default number of surface sample points per atom.
DEFAULT_N_POINTS = 512

#: Van der Waals radii in Å (Bondi, *J Phys Chem* 68:441-451, 1964; values for
#: biologically common elements, matching the radii used by FreeSASA and
#: Biopython so results are comparable across tools).
VDW_RADII: Dict[str, float] = {
    "H": 1.20,
    "D": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "NA": 2.27,
    "MG": 1.73,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "K": 2.75,
    "CA": 2.31,
    "MN": 1.61,
    "FE": 1.56,
    "CO": 1.50,
    "NI": 1.63,
    "CU": 1.40,
    "ZN": 1.39,
    "SE": 1.90,
    "BR": 1.85,
    "MO": 1.75,
    "I": 1.98,
    "CD": 1.58,
    "HG": 1.55,
}

#: Radius used for elements missing from :data:`VDW_RADII` (carbon-like).
DEFAULT_VDW_RADIUS = 1.70


def element_radius(element: Optional[str], *, default: float = DEFAULT_VDW_RADIUS) -> float:
    """Return the van der Waals radius for an element symbol.

    Args:
        element: Element symbol, any capitalisation. ``None`` or empty yields
            ``default``.
        default: Radius used for unknown elements.

    Returns:
        Radius in Å.
    """
    if not element:
        return default
    return VDW_RADII.get(str(element).strip().upper(), default)


def fibonacci_sphere(n_points: int = DEFAULT_N_POINTS) -> np.ndarray:
    """Generate ``n_points`` quasi-uniform points on the unit sphere.

    The Fibonacci (golden-angle) spiral gives a far more uniform covering than
    random sampling at the same point count, which is what makes a modest number
    of points sufficient for stable SASA values.

    Args:
        n_points: Number of points; must be positive.

    Returns:
        Array of shape ``(n_points, 3)`` of unit vectors.
    """
    n = int(n_points)
    if n <= 0:
        raise ValueError("n_points must be positive")
    indices = np.arange(n, dtype=float) + 0.5
    # Uniform in z, golden-angle in azimuth.
    z = 1.0 - 2.0 * indices / n
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    theta = np.pi * (1.0 + 5.0**0.5) * indices
    return np.column_stack((radius * np.cos(theta), radius * np.sin(theta), z))


def shrake_rupley_sasa(
    coords: np.ndarray,
    radii: np.ndarray,
    *,
    subset: Optional[Sequence[int]] = None,
    probe_radius: float = PROBE_RADIUS,
    n_points: int = DEFAULT_N_POINTS,
) -> np.ndarray:
    """Compute per-atom SASA by Shrake-Rupley numerical integration.

    Args:
        coords: Atom coordinates, shape ``(n_atoms, 3)``.
        radii: Van der Waals radii, shape ``(n_atoms,)``.
        subset: Indices of atoms to report SASA for. All atoms in ``coords``
            always act as occluders; only ``subset`` atoms are measured. This
            separation is what allows a ligand to be measured inside a complex
            and again in isolation. ``None`` measures every atom.
        probe_radius: Solvent probe radius in Å.
        n_points: Sample points per atom.

    Returns:
        SASA in Å² for each requested atom, ordered as ``subset`` (or as
        ``coords`` when ``subset`` is ``None``).

    Raises:
        ValueError: If ``coords`` and ``radii`` shapes are inconsistent.
    """
    coords = np.asarray(coords, dtype=float)
    radii = np.asarray(radii, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coords must have shape (n_atoms, 3)")
    if radii.shape != (coords.shape[0],):
        raise ValueError("radii must have shape (n_atoms,)")

    n_atoms = coords.shape[0]
    indices = np.arange(n_atoms) if subset is None else np.asarray(subset, dtype=int)
    if n_atoms == 0 or indices.size == 0:
        return np.zeros(indices.size, dtype=float)

    expanded = radii + probe_radius
    sphere = fibonacci_sphere(n_points)
    tree = cKDTree(coords)
    max_expanded = float(np.max(expanded))

    out = np.zeros(indices.size, dtype=float)
    for slot, atom_index in enumerate(indices):
        r_i = expanded[atom_index]
        centre = coords[atom_index]
        # Any atom whose expanded sphere can reach atom i's sphere is a
        # potential occluder.
        neighbours = tree.query_ball_point(centre, r_i + max_expanded)
        neighbours = [j for j in neighbours if j != atom_index]
        area = 4.0 * np.pi * r_i * r_i
        if not neighbours:
            out[slot] = area
            continue

        test_points = centre + r_i * sphere
        # A point is buried if it falls inside any neighbour's expanded sphere.
        buried = np.zeros(test_points.shape[0], dtype=bool)
        for j in neighbours:
            delta = test_points - coords[j]
            buried |= np.einsum("ij,ij->i", delta, delta) < expanded[j] ** 2
            if buried.all():
                break
        out[slot] = area * float(np.count_nonzero(~buried)) / test_points.shape[0]
    return out


def total_sasa(
    coords: np.ndarray,
    radii: np.ndarray,
    *,
    subset: Optional[Sequence[int]] = None,
    probe_radius: float = PROBE_RADIUS,
    n_points: int = DEFAULT_N_POINTS,
) -> float:
    """Return the summed SASA of the requested atoms.

    Args:
        coords: Atom coordinates.
        radii: Van der Waals radii.
        subset: Atoms to measure; all atoms occlude.
        probe_radius: Solvent probe radius in Å.
        n_points: Sample points per atom.

    Returns:
        Total SASA in Å².
    """
    return float(
        np.sum(
            shrake_rupley_sasa(
                coords, radii, subset=subset, probe_radius=probe_radius, n_points=n_points
            )
        )
    )


def burial_depth(
    point: Sequence[float],
    context_coords: np.ndarray,
    context_sasa: np.ndarray,
    *,
    min_atom_sasa: float = 5.0,
) -> float:
    """Distance from a point to the nearest solvent-exposed atom.

    This is the geometric definition of burial depth used throughout the
    pipeline: how far below the molecular surface a site sits.

    The surface must be defined on the **holo** structure (ligand present, waters
    removed). Removing the ligand first would let the solvent probe enter the
    vacated cavity, so cavity-lining atoms would register as "exposed" and a
    fully enclosed site would report a depth near zero - the opposite of the
    truth. Keeping the ligand in place leaves only the true exterior accessible.

    Args:
        point: Query point (e.g. ligand centroid or pocket centre).
        context_coords: Coordinates of the context atoms, shape ``(n, 3)``.
        context_sasa: Per-atom SASA of the same atoms in the same order.
        min_atom_sasa: SASA above which an atom counts as exposed (Å²). The
            default of 5 Å² requires genuine solvent contact rather than a sliver
            of accessibility. A near-zero threshold makes the measure fragile:
            atoms bordering a narrow crevice pick up a fraction of an Å², so a
            fully enclosed site can report a depth of a few Å instead of the
            tens of Å that its geometry implies. On the synthetic closed-shell
            benchmark, thresholds below ~1 Å² give 5 Å for a ligand at the centre
            of a 22 Å sphere, while 5 Å² recovers the correct 22 Å.

    Returns:
        Distance in Å to the nearest exposed atom. ``nan`` when the context has
        no exposed atom (an ill-posed input, not a deep site).
    """
    context_coords = np.asarray(context_coords, dtype=float)
    context_sasa = np.asarray(context_sasa, dtype=float)
    if context_coords.size == 0:
        return float("nan")
    exposed = context_coords[context_sasa > float(min_atom_sasa)]
    if exposed.size == 0:
        return float("nan")
    delta = exposed - np.asarray(point, dtype=float)
    return float(np.sqrt(np.min(np.einsum("ij,ij->i", delta, delta))))


def enclosure_fraction(
    point: Sequence[float],
    context_coords: np.ndarray,
    context_radii: np.ndarray,
    *,
    n_rays: int = 256,
    max_distance: float = 25.0,
) -> float:
    """Fraction of directions from a point that are blocked by protein atoms.

    Enclosure captures "crypticness" in a way SASA cannot. SASA answers *is
    solvent touching this atom*; enclosure answers *is this site surrounded*. A
    site in a deep groove can have near-zero SASA yet remain open on one side,
    whereas the buried inositol phosphate site of ADAR2 is enclosed in nearly
    every direction. The two measures are complementary, and enclosure is
    invariant to ligand size and copy number.

    Args:
        point: Origin of the rays (ligand centroid or pocket centre).
        context_coords: Atom coordinates, shape ``(n, 3)``.
        context_radii: Van der Waals radii of those atoms.
        n_rays: Number of ray directions (Fibonacci sphere).
        max_distance: Maximum distance along a ray to look for a blocker (Å).

    Returns:
        Blocked fraction in ``[0, 1]``; ``nan`` for an empty context.
    """
    context_coords = np.asarray(context_coords, dtype=float)
    context_radii = np.asarray(context_radii, dtype=float)
    if context_coords.size == 0:
        return float("nan")

    origin = np.asarray(point, dtype=float)
    offsets = context_coords - origin
    distances_sq = np.einsum("ij,ij->i", offsets, offsets)
    within = distances_sq <= float(max_distance) ** 2
    if not np.any(within):
        return 0.0
    offsets = offsets[within]
    radii = context_radii[within]

    directions = fibonacci_sphere(n_rays)
    # Projection of each atom vector on each ray direction.
    projection = offsets @ directions.T  # (n_atoms, n_rays)
    perpendicular_sq = np.einsum("ij,ij->i", offsets, offsets)[:, None] - projection**2
    # A ray is blocked when some atom lies ahead of the origin and the ray
    # passes within that atom's van der Waals radius.
    blocked = (projection > 0.0) & (perpendicular_sq <= (radii**2)[:, None])
    return float(np.mean(np.any(blocked, axis=0)))


@dataclass(frozen=True)
class PointCloudShape:
    """Shape descriptors of a point cloud derived from its gyration tensor."""

    radius_of_gyration: float
    asphericity: float
    max_extent: float

    def to_dict(self) -> Dict[str, float]:
        """Return a JSON-serialisable representation."""
        return {
            "radius_of_gyration": self.radius_of_gyration,
            "asphericity": self.asphericity,
            "max_extent": self.max_extent,
        }


def point_cloud_shape(coords: np.ndarray) -> PointCloudShape:
    """Describe the size and anisotropy of a point cloud.

    Shape matters for inositol phosphate sites: the ligand is a compact,
    near-spherical polyanion, so an elongated crevice of the same volume is a
    poorer match than a globular cavity. Asphericity is computed from the
    eigenvalues of the gyration tensor and is ``0`` for an isotropic cloud.

    Args:
        coords: Point coordinates, shape ``(n, 3)``.

    Returns:
        A :class:`PointCloudShape`. Degenerate inputs give ``nan`` fields.
    """
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[0] < 2:
        return PointCloudShape(float("nan"), float("nan"), float("nan"))

    centred = coords - coords.mean(axis=0)
    gyration = (centred.T @ centred) / coords.shape[0]
    eigenvalues = np.sort(np.linalg.eigvalsh(gyration))[::-1]
    trace = float(np.sum(eigenvalues))
    rg = float(np.sqrt(max(trace, 0.0)))
    if trace <= 0:
        asphericity = float("nan")
    else:
        # Normalised asphericity: 0 for a sphere, 1 for a line.
        asphericity = float((eigenvalues[0] - 0.5 * (eigenvalues[1] + eigenvalues[2])) / trace)
    distances = np.sqrt(np.einsum("ij,ij->i", centred, centred))
    return PointCloudShape(rg, asphericity, float(2.0 * np.max(distances)))


def convex_hull_volume(coords: np.ndarray) -> float:
    """Convex hull volume of a point cloud in Å³.

    Args:
        coords: Point coordinates, shape ``(n, 3)``.

    Returns:
        Hull volume, or ``nan`` when the cloud is degenerate (fewer than four
        points or coplanar).
    """
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[0] < 4:
        return float("nan")
    try:
        from scipy.spatial import ConvexHull

        return float(ConvexHull(coords).volume)
    except Exception:  # noqa: BLE001 - QhullError and degenerate inputs
        return float("nan")


def pairwise_min_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Minimum distance between two coordinate sets.

    Args:
        a: First coordinate set, shape ``(n, 3)``.
        b: Second coordinate set, shape ``(m, 3)``.

    Returns:
        Minimum Euclidean distance, or ``nan`` when either set is empty.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    tree = cKDTree(b)
    distances, _ = tree.query(a, k=1)
    return float(np.min(distances))


def count_within(
    query: np.ndarray, targets: np.ndarray, cutoff: float
) -> Tuple[int, np.ndarray]:
    """Count target points within ``cutoff`` of any query point.

    Args:
        query: Query coordinates, shape ``(n, 3)``.
        targets: Target coordinates, shape ``(m, 3)``.
        cutoff: Distance cutoff in Å.

    Returns:
        ``(count, indices)`` where ``indices`` are the matching target rows.
    """
    query = np.asarray(query, dtype=float)
    targets = np.asarray(targets, dtype=float)
    if query.size == 0 or targets.size == 0:
        return 0, np.empty(0, dtype=int)
    tree = cKDTree(targets)
    hits = tree.query_ball_point(query, float(cutoff))
    matched = sorted({index for group in hits for index in group})
    return len(matched), np.asarray(matched, dtype=int)
