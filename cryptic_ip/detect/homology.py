"""A metric on lining architectures, and what it can and cannot relate.

`architecture.py` reads one lining. This asks whether two linings are *the same
kind of thing* - whether a metric on their sequence-space structure puts
linings from the same fold near each other.

What this is not
----------------
This is **architecture homology**, not sequence or structural homology. It does
not detect that two proteins are evolutionarily related, it does not replace
Foldseek or MMseqs2, and it must never be used in place of the project's
`homology_group_strict` for grouping - the statistical independence that
grouping provides comes from sequence and structure comparison, and a metric on
lining spacing provides none of it. What it can do is say that two pockets are
lined in the same *way*, which is a different and weaker claim.

The metric, and why it is the one it is
---------------------------------------
`architecture.py` proves that a lining factorises into its gap multiset (the
complete H0 invariant) and the order of those gaps. A metric has to respect
that factorisation, and each half has a canonical choice.

**Topology.** The gap multiset *is* the H0 persistence diagram, with every
point born at 0 and dying at its gap. So the canonical distance is the one
persistent homology already has: the **Wasserstein distance between persistence
diagrams**, computed exactly here by optimal matching with the diagonal
available for unmatched points. Nothing is invented; the identity in
`architecture.py` is what makes the standard TDA metric applicable to a lining
at all.

**Register.** The coherence spectrum ``w -> R(w)`` is **exactly invariant to
translation along the chain**: shifting every position by c multiplies the
resultant by a unit phase, which the modulus discards. That makes the spectrum
a position-free signature, and comparing *shapes* rather than levels removes
the dependence on lining size, since the expected resultant of a random lining
falls as ``sqrt(pi / 4n)``. The shape distance is the Euclidean distance
between standardised spectra, which equals ``sqrt(2 (1 - r))`` for Pearson r -
a true metric, bounded on [0, 2].

Both halves are metrics, so any non-negative combination of them is a metric,
and `architecture_distance` is one with a stated scale rather than a fitted
weight.

One honest tension
------------------
`architecture.band_probes` derives probe density from the span, because a
coherence peak has half-width ``0.58 w^2 / span``. Comparing spectra requires a
**common** grid, so `SPECTRUM_GRID` is fixed and a long lining is therefore
compared at coarser relative resolution than it would be measured at. That is a
property of comparison, not a bug, and it is why the topology and register
distances are always reported separately as well as combined.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

from .architecture import (
    ArchitectureError,
    SCAN_BAND,
    gaps,
    resultants,
)

#: The common grid every coherence spectrum is sampled on, so that two linings
#: of different span are comparable at all. 241 probes over 2-12 residues is a
#: step of 0.042, which resolves the peak of a lining spanning up to about 200
#: residues; past that the comparison is coarser than the measurement would be.
SPECTRUM_GRID = np.linspace(SCAN_BAND[0], SCAN_BAND[1], 241)

#: Residues of Wasserstein distance counted as one unit of register distance,
#: when the two halves are combined. Four, one face gap: a lining whose gap
#: multiset is a whole face-width away from another's is as different as two
#: linings a quarter of the maximum spectrum distance apart. It is a stated
#: convention for putting two different units on one axis, not a fitted
#: weight, and `topology_distance` and `register_distance` are always available
#: separately so a reader never has to accept it.
TOPOLOGY_SCALE = 4.0

#: Wasserstein order. 1 rather than infinity: the bottleneck distance is
#: decided by a single worst-matched gap, which on a dozen gaps makes the
#: metric a function of one number.
WASSERSTEIN_P = 1


# ------------------------------------------------------------- the topology

def persistence_diagram(positions: Sequence[int]) -> np.ndarray:
    """The H0 diagram of a lining: one point ``(0, gap)`` per gap.

    By the identity in `architecture.py` this is the complete H0 invariant, so
    it is built from the gaps directly rather than by running a filtration.
    """
    return np.column_stack([np.zeros_like(gaps(positions), dtype=float),
                            gaps(positions).astype(float)])


def topology_distance(a: Sequence[int], b: Sequence[int], *, p: int = WASSERSTEIN_P) -> float:
    """Exact Wasserstein-``p`` distance between two linings' H0 diagrams.

    Every point is born at 0, so the ground distance between two points is
    ``|g - h|`` and a point's distance to the diagonal is ``g / 2``. Unmatched
    points go to the diagonal, which is what lets diagrams of different size be
    compared - a lining with more gaps is not automatically far from one with
    fewer.

    Solved as an assignment problem rather than by matching sorted gaps,
    because sorted matching is only optimal when the two diagrams have the same
    cardinality; with the diagonal in play it is not.
    """
    first = gaps(a).astype(float)
    second = gaps(b).astype(float)
    m, n = first.size, second.size
    size = m + n
    big = float(first.sum() + second.sum() + 1.0) * 10.0
    cost = np.full((size, size), big)
    cost[:m, :n] = np.abs(first[:, None] - second[None, :])
    # Each real point may instead go to the diagonal, at half its death.
    for index in range(m):
        cost[index, n + index] = first[index] / 2.0
    for index in range(n):
        cost[m + index, index] = second[index] / 2.0
    cost[m:, n:] = 0.0
    rows, columns = linear_sum_assignment(cost ** p)
    total = float(np.sum(cost[rows, columns] ** p))
    return float(total ** (1.0 / p))


# ------------------------------------------------------------- the register

def coherence_spectrum(positions: Sequence[int], grid: np.ndarray = SPECTRUM_GRID) -> np.ndarray:
    """The lining's resultant at every period on the common grid.

    Exactly invariant to translation along the chain, because shifting all
    positions by c multiplies the resultant by ``exp(2 pi i c / w)`` and the
    modulus discards it.
    """
    ordered = np.asarray(sorted({int(x) for x in positions}), dtype=float)
    if ordered.size < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    return resultants(ordered, np.asarray(grid, dtype=float))


def standardise(spectrum: np.ndarray) -> np.ndarray:
    """Centre and scale a spectrum to its shape, discarding its level.

    The level is a function of lining size - a random lining's expected
    resultant is about ``sqrt(pi / 4n)`` - so comparing levels would compare
    sizes. What carries architecture is where the peaks are.
    """
    values = np.asarray(spectrum, dtype=float)
    centred = values - values.mean()
    scale = float(np.sqrt(np.sum(centred * centred)))
    if scale <= 1e-12:
        return np.zeros_like(centred)
    return centred / scale


def register_distance(a: Sequence[int], b: Sequence[int],
                      grid: np.ndarray = SPECTRUM_GRID) -> float:
    """Shape distance between two coherence spectra, on [0, 2].

    The Euclidean distance between standardised spectra, which equals
    ``sqrt(2 (1 - r))`` for Pearson correlation r. A true metric, and zero
    exactly when one spectrum is a positive affine image of the other.
    """
    first = standardise(coherence_spectrum(a, grid))
    second = standardise(coherence_spectrum(b, grid))
    return float(np.linalg.norm(first - second))


# -------------------------------------------------------------- combined

def architecture_distance(
    a: Sequence[int],
    b: Sequence[int],
    *,
    grid: np.ndarray = SPECTRUM_GRID,
    topology_scale: float = TOPOLOGY_SCALE,
) -> float:
    """Both halves, on one axis, with a stated scale rather than a fitted weight."""
    if topology_scale <= 0.0:
        raise ArchitectureError(f"topology scale must be positive, got {topology_scale}")
    return (register_distance(a, b, grid)
            + topology_distance(a, b) / float(topology_scale))


#: The three distances, by name, so an experiment can report all of them and a
#: reader can see which half carries the signal.
METRICS = {
    "topology": lambda a, b: topology_distance(a, b),
    "register": lambda a, b: register_distance(a, b),
    "architecture": lambda a, b: architecture_distance(a, b),
}


def distance_matrix(linings: Sequence[Sequence[int]], metric: str = "architecture") -> np.ndarray:
    """A symmetric distance matrix over linings, with a zero diagonal."""
    function = METRICS.get(metric)
    if function is None:
        raise ArchitectureError(f"unknown metric {metric!r}; known: {sorted(METRICS)}")
    count = len(linings)
    out = np.zeros((count, count))
    for i in range(count):
        for j in range(i + 1, count):
            value = float(function(linings[i], linings[j]))
            out[i, j] = out[j, i] = value
    return out


# ------------------------------------------------------------- retrieval

@dataclass(frozen=True)
class Retrieval:
    """How well a metric groups linings by a label they were not told."""

    metric: str
    n_items: int
    labels: Tuple[str, ...]
    precision_at_one: float
    mean_rank_of_first_hit: float
    within_mean: float
    between_mean: float
    separation: float

    def as_dict(self) -> Dict[str, object]:
        return {
            "metric": self.metric,
            "n_items": self.n_items,
            "labels": sorted(set(self.labels)),
            "precision_at_one": float(self.precision_at_one),
            "mean_rank_of_first_hit": float(self.mean_rank_of_first_hit),
            "within_mean": float(self.within_mean),
            "between_mean": float(self.between_mean),
            "separation": float(self.separation),
        }


def retrieval(
    linings: Sequence[Sequence[int]],
    labels: Sequence[str],
    metric: str = "architecture",
) -> Retrieval:
    """Is a lining's nearest neighbour from its own family?

    ``precision_at_one`` is the headline: the share of linings whose closest
    other lining carries the same label. ``separation`` is the gap between the
    mean within-label and between-label distance, divided by their pooled
    spread, so it says whether the families are separated or merely ordered.

    The chance level for precision at one is not 0.5 unless the families are
    balanced, so the caller is responsible for balancing them or for reporting
    the label counts, which `as_dict` includes.
    """
    if len(linings) != len(labels):
        raise ArchitectureError("one label per lining")
    if len(linings) < 3:
        raise ArchitectureError("retrieval needs at least three linings")
    matrix = distance_matrix(linings, metric)
    labels = tuple(str(x) for x in labels)
    count = len(linings)
    hits, ranks = 0, []
    within, between = [], []
    for index in range(count):
        order = [j for j in np.argsort(matrix[index]) if j != index]
        hits += labels[order[0]] == labels[index]
        same = [position for position, j in enumerate(order, start=1)
                if labels[j] == labels[index]]
        ranks.append(float(same[0]) if same else float(count))
        for j in range(count):
            if j == index:
                continue
            (within if labels[j] == labels[index] else between).append(matrix[index, j])
    within_mean = float(np.mean(within)) if within else float("nan")
    between_mean = float(np.mean(between)) if between else float("nan")
    pooled = float(np.sqrt((np.var(within) + np.var(between)) / 2.0)) if within and between else 0.0
    separation = 0.0 if pooled <= 1e-12 else (between_mean - within_mean) / pooled
    return Retrieval(
        metric=metric, n_items=count, labels=labels,
        precision_at_one=hits / count,
        mean_rank_of_first_hit=float(np.mean(ranks)),
        within_mean=within_mean, between_mean=between_mean,
        separation=float(separation),
    )


def label_permutation_p(
    linings: Sequence[Sequence[int]],
    labels: Sequence[str],
    metric: str = "architecture",
    *,
    n_draws: int = 2000,
    seed: int = 20261012,
) -> float:
    """Permutation p for precision at one, with the labels shuffled.

    The chance level depends on how many families there are and how big each
    is, so it is measured rather than assumed. The distance matrix is computed
    once and only the labels move, which is what makes this cheap and also what
    makes it the right null: the geometry is held fixed.
    """
    matrix = distance_matrix(linings, metric)
    labels = np.asarray([str(x) for x in labels])
    count = len(labels)
    nearest = np.array([[j for j in np.argsort(matrix[i]) if j != i][0] for i in range(count)])
    observed = float(np.mean(labels[nearest] == labels))
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(int(n_draws)):
        shuffled = rng.permutation(labels)
        if float(np.mean(shuffled[nearest] == shuffled)) >= observed:
            extreme += 1
    return float((extreme + 1) / (int(n_draws) + 1))


# ------------------------------------------------------------- clustering

def single_linkage(matrix: np.ndarray, n_clusters: int) -> List[int]:
    """Single-linkage clusters, in the same spirit as the H0 filtration.

    Single linkage is the clustering whose merge scales *are* the H0 barcode of
    the distance matrix, so it is the clustering consistent with the rest of
    this module rather than an arbitrary pick.
    """
    matrix = np.asarray(matrix, dtype=float)
    count = matrix.shape[0]
    if n_clusters < 1 or n_clusters > count:
        raise ArchitectureError(f"cannot ask for {n_clusters} clusters from {count} items")
    assignment = list(range(count))
    members: Dict[int, List[int]] = {i: [i] for i in range(count)}
    while len(members) > n_clusters:
        best, pair = float("inf"), None
        keys = sorted(members)
        for i_index, i in enumerate(keys):
            for j in keys[i_index + 1:]:
                linkage = min(matrix[a, b] for a in members[i] for b in members[j])
                if linkage < best:
                    best, pair = linkage, (i, j)
        if pair is None:
            break
        i, j = pair
        members[i] = members[i] + members[j]
        del members[j]
    for label, (_key, group) in enumerate(sorted(members.items())):
        for item in group:
            assignment[item] = label
    return assignment


def cluster_purity(assignment: Sequence[int], labels: Sequence[str]) -> float:
    """Share of items in the majority label of their own cluster."""
    if len(assignment) != len(labels):
        raise ArchitectureError("one label per item")
    total = 0
    for cluster in set(assignment):
        members = [str(labels[i]) for i, c in enumerate(assignment) if c == cluster]
        counts: Dict[str, int] = {}
        for label in members:
            counts[label] = counts.get(label, 0) + 1
        total += max(counts.values())
    return total / len(labels)


# -------------------------------------------------- graded relationships

def parameter_correlation(
    linings: Sequence[Sequence[int]],
    parameters: Sequence[float],
    metric: str = "architecture",
) -> Dict[str, float]:
    """Does architecture distance track a continuous structural parameter?

    A family relation is not only "same or different". If the metric is reading
    structure then two barrels differing by one helix should be nearer than two
    differing by four, and this measures that as the Spearman correlation
    between architecture distance and ``|parameter difference|`` over all pairs.
    """
    matrix = distance_matrix(linings, metric)
    values = np.asarray(parameters, dtype=float)
    if values.size != matrix.shape[0]:
        raise ArchitectureError("one parameter per lining")
    rows, columns = np.triu_indices(matrix.shape[0], k=1)
    distances = matrix[rows, columns]
    differences = np.abs(values[rows] - values[columns])

    def rank(array: np.ndarray) -> np.ndarray:
        order = np.argsort(array)
        out = np.empty_like(array, dtype=float)
        out[order] = np.arange(array.size, dtype=float)
        return out

    first, second = rank(distances), rank(differences)
    if first.std() <= 1e-12 or second.std() <= 1e-12:
        return {"spearman": 0.0, "n_pairs": float(distances.size)}
    corr = float(np.corrcoef(first, second)[0, 1])
    return {"spearman": corr, "n_pairs": float(distances.size)}


def nearest_neighbours(
    query: Sequence[int],
    library: Mapping[str, Sequence[int]],
    metric: str = "architecture",
    *,
    k: int = 3,
) -> List[Tuple[str, float]]:
    """The ``k`` nearest linings in a named library, closest first."""
    function = METRICS.get(metric)
    if function is None:
        raise ArchitectureError(f"unknown metric {metric!r}; known: {sorted(METRICS)}")
    scored = [(name, float(function(query, positions))) for name, positions in library.items()]
    return sorted(scored, key=lambda pair: pair[1])[:int(k)]
