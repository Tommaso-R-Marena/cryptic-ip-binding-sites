"""The architecture of a pocket lining in sequence space.

A pocket is the statement that residues which are *far apart along the chain*
are *close together in space*. So the sequence-space object that matters is not
the window's letters - `sequence.py` handles those - but the **set of lining
positions**, and the question is what structure that set has.

The decomposition
-----------------
Everything here rests on one small theorem and the gap it leaves.

**The complete H0 invariant of a finite set of integer positions, under the
gap (Vietoris-Rips) filtration, is its sorted multiset of gaps.** Each point is
born at scale 0 and dies at the gap that merges it into its neighbour, so the
barcode *is* the gap multiset. Verified against explicit single-linkage on 300
random position sets: zero disagreements (`tests/test_detect_architecture.py`).

Two consequences, and they set up the whole module:

1. Computing persistent homology on a lining buys nothing. The gaps are the
   answer, and `segmentation` reads components straight off them.
2. The gap multiset is therefore **exactly the topological part** of a lining's
   structure, and whatever is left over lives in the *order* of the gaps. That
   leftover is real: over the 70 distinct orderings of the gap multiset
   {3,4,3,4,3,4,3,4}, helical-band phase coherence runs from 0.572 to 0.902.
   Same topology, same n, same span, same occupancy, same gap entropy - and a
   0.330 spread in coherence.

So a lining factorises into **topology (the gap multiset) + register (the order
of the gaps)**, and the canonical null for register is permutation of the gap
multiset - the quotient by the topological invariant. That is not an arbitrary
choice of null; it is the only one that holds the topology fixed and varies
nothing else.

Three coordinates
-----------------
* **fragmentation** - how many separate runs the lining has at a face-sized
  gap. Read off the barcode, exact, no estimation.
* **occupancy** - ``log n / log(span + 1)``, how densely the lining fills the
  stretch it spans. This is the full-range box-counting dimension of the set in
  closed form; the fitted slope agrees to about 0.03 for regular sets and
  deviates for clumped ones, which is information `gap_entropy` reports
  directly, so nothing is fitted here.
* **register** - phase coherence against the gap-permutation null, as a
  z-score. Measured on 600 random linings: mean -0.04, sd 1.04, |z| > 2 in
  5.2% of them, so it is calibrated and it is decorrelated from occupancy
  (r = +0.08) and from gap entropy (r = +0.05).

Why raw coherence is not used
-----------------------------
Two reasons, both measured.

**It is confounded by occupancy**, r = -0.65 over random linings: a contiguous
run samples every phase and so scores low, a sparse set scores high for free.
The null fixes this by construction, because permuting gaps cannot change n,
span or the gap distribution.

**And a scanned coherence has no fixed threshold.** The analytic Rayleigh test
gives a period-independent 95th percentile of 0.500 at n = 12. The true null is
period-*dependent* - using 0.500 as the cutoff gives a real false-positive rate
of 0.019 at period 3.5 and 0.102 at period 2, a five-fold spread, because
integer positions on a rational period cannot produce uniform phases. Worse,
the statistic a band scan actually uses is a maximum over probes, and there
0.500 admits **26%** of random linings over the helical band and **92%** over a
full 2-12 scan. Any analytic threshold on a scanned coherence is therefore
useless, and the only honest null resamples the same statistic.

Local register, and why it did not work
---------------------------------------
The pooled resultant is Neyman-Pearson optimal for a face of one period and
**one** unknown phase, and a real helix with backbone scatter has no single
phase: the rise per residue drifts and a bulge resets it. `local_register`
therefore asks the weaker question - whether the lining is **piecewise** phase
locked - by maximising over windows of consecutive lining positions, with the
identical maximum taken on every draw of the same gap-permutation null.
`pilot_sizes` can restrict which window sizes count, chosen **from the null**:
the resultant's modulus is translation invariant, so a window's null depends
only on the gap multiset, which the register null holds fixed, and a rule
stated in terms of it is the same function on the observation and on every
permutation, which leaves the p-value exact.

It does not beat the pooled statistic on any geometry arm, and the diagnostic
is one number: the null's own mean at a six-position window is 0.88 out of a
possible 1. Kept in the module with its proofs, its calibration and its
negative attached, because the negative bounds a whole family of statistics
rather than one implementation - see section 10 of
`results/architecture/ARCHITECTURE.md`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from .sequence import HELIX_BAND, STRAND_BAND, SequenceError

#: Positions closer than this stay on the same face of one secondary-structure
#: element: a helix face repeats every 3.6 residues, so gaps of 3 and 4 keep
#: you on it and 5 leaves it. Used to count runs, and fixed from that geometry
#: rather than tuned.
FACE_GAP = 4

#: The fewest lining positions from which register can be asked at all. Five
#: gaps give 120 orderings, which is enough for a permutation null to have
#: resolution; below that the null is degenerate and the answer is an
#: abstention rather than a number.
MIN_POSITIONS = 6

#: Coherence peak half-width in period, as a fraction of ``w^2 / span``.
#: Measured by bisecting the peak of a perfect face at period 3.6 over spans
#: from 20 to 320 residues: the ratio is 0.556, 0.58, 0.58, 0.59, 0.59 - flat
#: across a sixteen-fold range of span, so the law ``dw = 0.58 w^2 / span``
#: holds and probe spacing follows from it rather than from a guess.
PEAK_WIDTH_FACTOR = 0.58

#: Probes per peak half-width. Two, so a peak cannot fall between probes.
PROBES_PER_WIDTH = 2

#: Floor and ceiling on the probe count, for short and very long linings.
MIN_BAND_PROBES = 13
MAX_BAND_PROBES = 801

#: Draws for the gap-permutation null.
DEFAULT_DRAWS = 2000

#: |z| at which register would be called present *if the z were calibrated*.
#: It is kept as an effect size only. The verdict comes from the permutation
#: p-value, for a measured reason: on a gap multiset with 15 orderings the
#: sampled z reported +2.60 while the smallest p that null can even reach is
#: 0.0667, so the z claimed significance the data could not supply.
REGISTER_Z = 2.0

#: A permutation null over k distinct orderings cannot return a p below about
#: 1/k, so significance at 0.05 is unreachable below 20. Enumerated exactly:
#: 15 orderings floor at p = 0.0667, 35 at 0.0286, 252 at 0.0079. A lining
#: below this floor is `unresolved`, not "no register".
MIN_ORDERINGS = 20

#: Below this many distinct orderings the null is enumerated exactly rather
#: than sampled, which removes Monte-Carlo error from the p-value entirely.
#: Above it, sampling is used and the p-value carries the usual +1.
MAX_EXACT_ORDERINGS = 20000

#: Runs at or above which a lining counts as multi-element. Three, because one
#: secondary-structure element contributes one run and two runs is still a
#: hairpin or a helix read in two pieces.
MULTI_ELEMENT_SEGMENTS = 3

#: Permutation p-value bands for the graded verdict.
CLEAR_P = 0.01
WEAK_P = 0.05


class ArchitectureError(SequenceError):
    """A lining from which an architecture cannot honestly be read."""


# ------------------------------------------------------------ the topology

def gaps(positions: Sequence[int]) -> np.ndarray:
    """The gaps between consecutive lining positions, in order along the chain."""
    ordered = np.asarray(sorted({int(p) for p in positions}), dtype=int)
    if ordered.size < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    return np.diff(ordered)


def h0_barcode(positions: Sequence[int]) -> np.ndarray:
    """Death scales of the H0 persistence barcode, ascending.

    Equal to the sorted gap multiset, by the theorem in this module's
    docstring. Provided under its topological name because that is where the
    quantity comes from, and tested against explicit single-linkage so the
    identity is checked rather than asserted.
    """
    return np.sort(gaps(positions))


def segmentation(positions: Sequence[int], scale: int = FACE_GAP) -> int:
    """Number of runs separated by a gap wider than ``scale``.

    This is the number of H0 components still alive at that filtration scale,
    which the barcode gives directly as ``1 + #{gaps > scale}``.
    """
    return int(1 + np.sum(gaps(positions) > int(scale)))


def segmentation_curve(positions: Sequence[int]) -> Dict[int, int]:
    """Components alive at every scale from 1 to the span: the whole barcode.

    Carries exactly the information in the gap multiset, in the form that reads
    as "how many pieces is this lining, and at what scale do they join".
    """
    g = gaps(positions)
    return {int(s): int(1 + np.sum(g > s)) for s in range(1, int(g.sum()) + 1)}


def fragmentation(positions: Sequence[int], scale: int = FACE_GAP) -> float:
    """Runs per residue, on [0, 1]: 0 is one unbroken run, 1 is all isolated."""
    ordered = sorted({int(p) for p in positions})
    if len(ordered) < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    return float((segmentation(ordered, scale) - 1) / (len(ordered) - 1))


def span(positions: Sequence[int]) -> int:
    ordered = sorted({int(p) for p in positions})
    return int(ordered[-1] - ordered[0])


def occupancy(positions: Sequence[int]) -> float:
    """``log n / log(span + 1)``: how densely the lining fills its own stretch.

    This is the full-range box-counting dimension in closed form. A contiguous
    run gives 1; n positions scattered over a span much larger than n tend to
    0. Nothing is fitted, which matters: a fitted log-log slope over a handful
    of scales on 10 to 20 points is an estimate with no error bar, and the
    closed form is exact.
    """
    ordered = sorted({int(p) for p in positions})
    if len(ordered) < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    extent = ordered[-1] - ordered[0] + 1
    return float(np.log(len(ordered)) / np.log(extent))


def gap_entropy(positions: Sequence[int]) -> float:
    """Normalised entropy of the gap distribution, on [0, 1].

    The persistence entropy of the H0 barcode, by the same identity. It reads
    as *regularity*: 1 when every gap is the same width - a contiguous run and
    an evenly spaced face both score 1 - and falling toward 0 when one gap
    dominates, which is what two distant segments look like.
    """
    g = gaps(positions).astype(float)
    if g.size < 2:
        return 0.0
    p = g / g.sum()
    entropy = float(-np.sum(p * np.log2(p)))
    return float(entropy / np.log2(g.size))


# ------------------------------------------------------------- the register

def band_probes(band: Sequence[float], lining_span: int) -> int:
    """How finely a band must be scanned before a peak can hide between probes.

    The peak half-width is ``0.58 * w^2 / span`` (measured; see
    `PEAK_WIDTH_FACTOR`), narrowest at the band's low edge, and two probes per
    half-width keeps it from being missed.

    This is not a refinement. A fixed thirteen-probe scan of the helical band
    loses up to **0.287** in resultant on a span-200 lining - a perfect face
    reading 0.871 comes back as 0.584 - and at the barrel span of 104 it loses
    up to 0.072. Those are missed faces, not rounding.
    """
    low, high = float(band[0]), float(band[1])
    if not low < high:
        raise ArchitectureError(f"a band needs low below high, got {band}")
    width = PEAK_WIDTH_FACTOR * low * low / max(1, int(lining_span))
    needed = int(np.ceil(PROBES_PER_WIDTH * (high - low) / width)) + 1
    return int(min(MAX_BAND_PROBES, max(MIN_BAND_PROBES, needed)))


def resultants(positions: np.ndarray, periods: np.ndarray) -> np.ndarray:
    """Mean resultant length at every period, vectorised over both axes.

    ``positions`` may be a single lining ``(n,)`` or a stack of them
    ``(draws, n)``; the result is ``(len(periods),)`` or
    ``(draws, len(periods))``.
    """
    pos = np.asarray(positions, dtype=float)
    per = np.asarray(periods, dtype=float)
    angles = 2.0 * np.pi * pos[..., :, None] / per[None, :]
    return np.abs(np.exp(1j * angles).mean(axis=-2))


#: The axial (doubled-angle) resultant, standard in directional statistics for
#: data with two modes 180 degrees apart - which is what a beta strand's two
#: faces are - is **not** used here, and the reason is a provable degeneracy
#: rather than a preference. On integer positions at period 2 the doubled angle
#: is ``2 * pi * p``, identically zero for every integer, so the axial
#: resultant is exactly 1 for *any* lining whatsoever: measured
#: ``max|R2 - 1| = 0`` over 2000 arbitrary integer position sets. At period 4
#: it is exactly 1 for any set of uniform parity. A strand band containing
#: period 2 would therefore call every lining axially coherent. The first
#: resultant has the same degeneracy only at period 2 and only for
#: uniform-parity sets, which is precisely the strand-face property it should
#: detect, so it is kept.
AXIAL_UNUSABLE = (
    "the doubled-angle resultant is identically 1 at period 2 for any integer "
    "position set, so it cannot be used on sequence positions"
)


#: **Why the resultant and not something else.** Model a secondary-structure
#: face as a von Mises preference on phase: a residue at position p lines the
#: pocket with density proportional to ``exp(kappa cos(2 pi p / w - phi))``,
#: where ``phi`` is the element's face phase and is unknown. Marginalising that
#: nuisance over a uniform prior gives
#:
#:     (1 / 2 pi) * integral exp(kappa * sum_p cos(theta_p - phi)) dphi
#:         = I0(kappa * n * Rbar)
#:
#: which is monotone in the resultant. So **the resultant is the Neyman-Pearson
#: optimal statistic for a face of known period and unknown phase** - it is not
#: a reasonable choice among many, it is the choice. Verified numerically
#: against the integral to nine decimal places.
#:
#: The same derivation has a consequence that cost an experiment. For k
#: elements with *independent* face phases the optimum is a sum of per-element
#: terms, ``sum_e log I0(kappa n_e Rbar_e)``, while a resultant pooled over all
#: positions implicitly requires every element's face to be in phase. Measured
#: on six elements of six residues at random independent phases: pooled Rbar
#: 0.327 against a per-element mean of 0.886.
#:
#: That predicts a per-element statistic should beat the pooled one on
#: multi-element linings, and **it does not**: on ring-shaped barrel linings it
#: detected 0 of 12 where the pooled statistic detected 10. The reason is a
#: second-order effect the derivation ignores - a band *scan* over a run of two
#: or three residues can fit almost any spacing, so the per-element statistic
#: saturates and its permutation null saturates with it. The derivation is
#: sound and the estimator built from it is not, at the run lengths a real
#: lining gives. Recorded here because it is the kind of negative that is
#: otherwise repeated.
VON_MISES_OPTIMALITY = (
    "marginalising an unknown face phase over a uniform prior gives "
    "I0(kappa n Rbar), so the resultant is the Neyman-Pearson statistic"
)


def phase_resultant(positions: Sequence[int], period: float) -> float:
    """Circular concentration of the lining's phases at one period, on [0, 1].

    Each position is mapped to an angle ``2*pi*p/period`` and the mean resultant
    length is taken. This asks a question the charge periodogram in
    `sequence.py` cannot: whether the *lining residues themselves* sit on one
    face, with no reference to what those residues are.
    """
    if period <= 1.0:
        raise ArchitectureError(f"period must exceed 1 residue, got {period}")
    ordered = np.asarray(sorted({int(p) for p in positions}), dtype=float)
    if ordered.size < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    return float(abs(np.exp(2j * np.pi * ordered / float(period)).mean()))


def band_grid(band: Sequence[float], lining_span: int, probes: Optional[int] = None) -> np.ndarray:
    """The periods a band is scanned at, spaced from the measured peak width."""
    low, high = float(band[0]), float(band[1])
    if not low < high:
        raise ArchitectureError(f"a band needs low below high, got {band}")
    count = int(probes) if probes is not None else band_probes(band, lining_span)
    if count < 2:
        raise ArchitectureError(f"a band needs at least 2 probes, got {count}")
    return np.linspace(low, high, count)


def coherence(
    positions: Sequence[int],
    band: Sequence[float] = HELIX_BAND,
    *,
    probes: Optional[int] = None,
) -> Tuple[float, float]:
    """The best phase concentration over a band, and the period that achieved it.

    Reported as a pair because the value alone is not interpretable: it is a
    maximum over correlated probes and its null depends on the band scanned.
    `register` is what turns it into a statement, and `dominant_period` is what
    turns it into a structural label - the two are different questions and are
    never read off each other.
    """
    ordered = np.asarray(sorted({int(p) for p in positions}), dtype=float)
    if ordered.size < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    grid = band_grid(band, int(ordered[-1] - ordered[0]), probes)
    values = resultants(ordered, grid)
    best = int(np.argmax(values))
    return float(values[best]), float(grid[best])


#: The widest period worth scanning for a structural label. Past twelve a
#: "period" on a pocket lining is one cycle.
SCAN_BAND = (2.0, 12.0)


def dominant_period(positions: Sequence[int], band: Sequence[float] = SCAN_BAND) -> Dict[str, float]:
    """The period that concentrates the lining's phases best, over a wide scan.

    This is the **structural** label and it is deliberately not the statistic
    the register test uses. The two come apart: a lining with gaps alternating
    3 and 4 - plainly a helical face - peaks at 0.902 at period 3.50, and also
    at 0.902 at period 7.00 and 0.902 at 2.80, its subharmonics. Picking the
    band with the best *z* can therefore return a strand-band period of 2.32
    for a helix, which is statistically correct and structurally wrong. A band
    result only earns a secondary-structure name when the global maximum lies
    inside that band.
    """
    ordered = np.asarray(sorted({int(p) for p in positions}), dtype=float)
    if ordered.size < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    lining_span = int(ordered[-1] - ordered[0])
    grid = band_grid(band, lining_span, None)
    values = resultants(ordered, grid)
    best = int(np.argmax(values))
    period = float(grid[best])
    return {
        "period": period,
        "resultant": float(values[best]),
        "in_helix_band": float(HELIX_BAND[0] <= period <= HELIX_BAND[1]),
        "in_strand_band": float(STRAND_BAND[0] <= period <= STRAND_BAND[1]),
        "probes": float(grid.size),
    }


@dataclass(frozen=True)
class RegisterResult:
    """Phase coherence beside the gap-permutation null it was judged against."""

    band: Tuple[float, float]
    observed: float
    period: float
    null_mean: float
    null_sd: float
    z: float
    p_one_sided: float
    n_draws: int
    n_distinct_orderings: int
    exact: bool = False

    @property
    def minimum_attainable_p(self) -> float:
        """The smallest p this null could return, which is about 1/k."""
        return 1.0 / max(1, self.n_distinct_orderings)

    @property
    def in_register(self) -> bool:
        """Whether the ordering carries coherence the gap multiset does not force.

        Decided by the permutation p-value and not by the z. The z is an effect
        size on a null that may be discrete and skewed; the permutation p-value
        is valid however few orderings there are - it simply cannot go below
        1/k, which is why `MIN_ORDERINGS` exists.
        """
        return self.p_one_sided <= WEAK_P

    @property
    def verdict(self) -> str:
        """Graded: ``none``, ``weak`` or ``clear``.

        A three-level scale rather than a binary, because the honest answer to
        "is there a signature" on a dozen residues is usually "sort of", and a
        p of 0.04 and a p of 0.0005 should not read the same.
        """
        if self.n_distinct_orderings < MIN_ORDERINGS:
            return "unresolved"
        if self.p_one_sided <= CLEAR_P:
            return "clear"
        if self.p_one_sided <= WEAK_P:
            return "weak"
        return "none"

    def as_dict(self) -> Dict[str, object]:
        return {
            "band": [float(self.band[0]), float(self.band[1])],
            "observed": float(self.observed),
            "period": float(self.period),
            "null_mean": float(self.null_mean),
            "null_sd": float(self.null_sd),
            "z": float(self.z),
            "p_one_sided": float(self.p_one_sided),
            "n_draws": int(self.n_draws),
            "n_distinct_orderings": int(self.n_distinct_orderings),
            "exact": bool(self.exact),
            "minimum_attainable_p": float(self.minimum_attainable_p),
            "in_register": self.in_register,
            "verdict": self.verdict,
        }


def distinct_orderings(positions: Sequence[int]) -> int:
    """How many distinct arrangements the gap multiset admits.

    The multinomial over repeated gap widths. Reported because it is the
    resolution limit of the null: a lining whose gaps are all equal has exactly
    one ordering, the null is a point mass, and no register can be detected -
    which is a fact about the lining, not a failure, and must not be returned
    as z = 0 without saying so.
    """
    from math import factorial

    g = gaps(positions)
    _, counts = np.unique(g, return_counts=True)
    total = factorial(int(g.size))
    for count in counts:
        total //= factorial(int(count))
    return int(total)


def multiset_permutations(items: Sequence[int]):
    """Every distinct ordering of a multiset, each exactly once.

    Written out rather than taken from a library because
    ``itertools.permutations`` on ten gaps yields 3.6 million tuples to reach
    252 distinct orderings, and deduplicating that is the slow way to do an
    exact test.
    """
    pool = sorted(int(x) for x in items)
    n = len(pool)

    def walk(remaining: list, built: list):
        if not remaining:
            yield tuple(built)
            return
        seen = None
        for index, value in enumerate(remaining):
            if value == seen:
                continue
            seen = value
            yield from walk(remaining[:index] + remaining[index + 1:], built + [value])

    if n == 0:
        return
    yield from walk(pool, [])


def register(
    positions: Sequence[int],
    band: Sequence[float] = HELIX_BAND,
    *,
    probes: Optional[int] = None,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = 20261012,
) -> RegisterResult:
    """Is the lining more phase-coherent than its own gaps, reordered?

    The null permutes the gap multiset, which holds the lining's topology - and
    with it n, span, occupancy and gap entropy - exactly fixed, and varies only
    the order. By the decomposition in this module's docstring that is the
    whole of what is left, so this is the complete test of register.

    It also computes the null for the *scanned* statistic, which is the only
    way to get a calibrated answer: a fixed 0.500 cutoff on a band maximum
    admits 26% of random linings over the helical band and 92% over a full
    scan.
    """
    ordered = sorted({int(p) for p in positions})
    if len(ordered) < MIN_POSITIONS:
        raise ArchitectureError(
            f"{len(ordered)} positions is too few for a register null; "
            f"{MIN_POSITIONS} are needed and five gaps give 120 orderings"
        )
    g = gaps(ordered)
    start = ordered[0]
    # The grid is fixed by the observed lining's span, and every permutation of
    # its gaps has the same span, so one grid serves the observation and the
    # whole null. The null must be computed on the same scan as the
    # observation: a maximum over probes has a null that depends on the probes.
    grid = band_grid(band, int(ordered[-1] - ordered[0]), probes)
    observed_scan = resultants(np.asarray(ordered, dtype=float), grid)
    best = int(np.argmax(observed_scan))
    observed, period = float(observed_scan[best]), float(grid[best])
    orderings = distinct_orderings(ordered)

    def scan(orders: np.ndarray) -> np.ndarray:
        """Band maxima for a stack of gap orderings, in chunks to bound memory."""
        out = np.empty(orders.shape[0])
        step = max(1, int(4_000_000 // (orders.shape[1] * grid.size + 1)))
        for begin in range(0, orders.shape[0], step):
            block = orders[begin:begin + step]
            built = start + np.cumsum(block, axis=1)
            stacked = np.concatenate(
                [np.full((block.shape[0], 1), float(start)), built], axis=1)
            out[begin:begin + step] = resultants(stacked, grid).max(axis=1)
        return out

    exact = orderings <= MAX_EXACT_ORDERINGS
    if exact:
        # Every ordering once: the p-value is then exact, with no Monte-Carlo
        # error, and the observed ordering is itself one of the draws so the
        # count is at least one and no +1 correction is needed.
        orders = np.asarray(list(multiset_permutations(g)), dtype=float)
        draws = scan(orders)
        p_value = float(np.mean(draws >= observed - 1e-12))
        used = int(draws.size)
    else:
        rng = np.random.default_rng(seed)
        orders = np.asarray([rng.permutation(g) for _ in range(int(n_draws))], dtype=float)
        draws = scan(orders)
        extreme = int(np.sum(draws >= observed - 1e-12))
        p_value = float((extreme + 1) / (int(n_draws) + 1))
        used = int(n_draws)
    mean = float(draws.mean())
    sd = float(draws.std(ddof=1)) if draws.size > 1 else 0.0
    z = 0.0 if sd <= 1e-12 else float((observed - mean) / sd)
    return RegisterResult(
        band=(float(band[0]), float(band[1])),
        observed=observed, period=period,
        null_mean=mean, null_sd=sd, z=z,
        p_one_sided=p_value,
        n_draws=used, n_distinct_orderings=orderings, exact=exact,
    )


# ------------------------------------------- register, but only locally

#: Window sizes, in **lining positions** rather than residues, for the local
#: register statistic. Fixed in `docs/ARCHITECTURE_PLAN_AMENDMENT_3.md` before
#: the statistic was computed on anything, and deliberately not tuned against
#: the result: the lower end is `MIN_POSITIONS`, below which a permutation null
#: has too few orderings, and the upper end is where a window stops being local
#: on the 15-position median real lining.
WINDOW_SIZES = (6, 8, 10)


def windows(n_positions: int, sizes: Sequence[int] = WINDOW_SIZES) -> List[Tuple[int, int]]:
    """Every (size, start) window of consecutive lining positions, by index.

    Windows are runs in the *lining*, not in the sequence, so a window spans
    whatever stretch of chain its positions happen to occupy. Sizes larger than
    the lining are dropped rather than clipped, because a clipped window is a
    different window and would make the statistic depend on n in a way the null
    could not reproduce.
    """
    out: List[Tuple[int, int]] = []
    for size in sorted({int(s) for s in sizes}):
        if size < 2:
            raise ArchitectureError(f"a window needs at least 2 positions, got {size}")
        for start in range(0, int(n_positions) - size + 1):
            out.append((size, start))
    return out


def windowed_maxima(
    stacked: np.ndarray,
    grid: np.ndarray,
    sizes: Sequence[int] = WINDOW_SIZES,
) -> Dict[int, np.ndarray]:
    """Per-window-size maximum resultant for a stack of linings ``(draws, n)``.

    Returned per size rather than already pooled because the per-size null mean
    is the **saturation diagnostic**: a null mean near 1 means the band scan
    fits any spacing at that window length, which is exactly what killed the
    per-element statistic, and it has to be visible rather than inferred.
    """
    pos = np.asarray(stacked, dtype=float)
    if pos.ndim != 2:
        raise ArchitectureError(f"expected a (draws, n) stack, got shape {pos.shape}")
    draws, n = pos.shape
    out: Dict[int, np.ndarray] = {}
    for size in sorted({int(s) for s in sizes}):
        if size > n:
            continue
        best = np.zeros(draws)
        # Chunked over draws on the same element budget `register` uses: the
        # intermediate is (chunk, size, probes) complex and nothing here is
        # allowed to grow with the number of orderings.
        step = max(1, int(4_000_000 // (size * grid.size + 1)))
        for start in range(0, n - size + 1):
            block = pos[:, start:start + size]
            for begin in range(0, draws, step):
                chunk = block[begin:begin + step]
                best[begin:begin + step] = np.maximum(
                    best[begin:begin + step], resultants(chunk, grid).max(axis=1))
        out[size] = best
    if not out:
        raise ArchitectureError(
            f"no window of {tuple(sizes)} fits a lining of {n} positions")
    return out


def best_window(
    positions: Sequence[int],
    band: Sequence[float] = HELIX_BAND,
    *,
    sizes: Sequence[int] = WINDOW_SIZES,
    probes: Optional[int] = None,
) -> Dict[str, float]:
    """The single best window of a lining, and where it is.

    The grid is the **whole lining's** grid, not the window's, and that is
    deliberate. Probe spacing comes from the measured peak-width law
    ``0.58 * omega^2 / span``, so a shorter span needs *fewer* probes; using
    the full lining's grid is therefore at least as fine as any window requires
    and, far more importantly, is identical for the observation and for every
    permutation of the null. A grid that varied with the window would make the
    maximum incomparable across draws.
    """
    ordered = np.asarray(sorted({int(p) for p in positions}), dtype=float)
    if ordered.size < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    grid = band_grid(band, int(ordered[-1] - ordered[0]), probes)
    best = {"resultant": -1.0, "period": float("nan"), "size": 0, "start": -1}
    for size, start in windows(int(ordered.size), sizes):
        values = resultants(ordered[start:start + size], grid)
        index = int(np.argmax(values))
        if float(values[index]) > best["resultant"]:
            best = {"resultant": float(values[index]), "period": float(grid[index]),
                    "size": int(size), "start": int(start)}
    if best["start"] < 0:
        raise ArchitectureError(
            f"no window of {tuple(sizes)} fits a lining of {ordered.size} positions")
    best["first_position"] = float(ordered[best["start"]])
    best["last_position"] = float(ordered[best["start"] + best["size"] - 1])
    return best


#: A null mean at or above this is **saturation**: the statistic's own null
#: already sits where the observation would have to be, so the test has no
#: power left at that window length. Fixed in amendment 3, before any pruned
#: statistic existed, and inherited unchanged by the pruning in amendment 4 -
#: which is the only reason the pruning is not a tuned threshold.
SATURATED_NULL = 0.9


#: How window sizes are chosen. ``all`` keeps every size that fits, as
#: amendment 3 fixed it. ``unsaturated`` keeps the sizes whose null mean is
#: below `SATURATED_NULL` (amendment 4). ``best`` keeps the single size with
#: the **lowest** null mean (amendment 5) - the same selection argument with a
#: rank in place of a threshold, which is what makes it fire on every lining
#: instead of only on the linings that happen to cross 0.9.
SELECTIONS = ("all", "unsaturated", "best")


#: Permutations used to *select* window sizes, drawn with their own seed and
#: independently of the null sample. See `pilot_sizes` for why the separation
#: matters.
PILOT_DRAWS = 200


def pilot_sizes(
    positions: Sequence[int],
    band: Sequence[float] = HELIX_BAND,
    *,
    sizes: Sequence[int] = WINDOW_SIZES,
    probes: Optional[int] = None,
    pilot_draws: int = PILOT_DRAWS,
    pilot_seed: int = 20261013,
    rule: str = "unsaturated",
) -> Tuple[Tuple[int, ...], Dict[int, float]]:
    """Which window sizes are not saturated, decided **from the null alone**.

    Tuning window sizes against the observation would invalidate the
    permutation p-value. Selecting them from the null does not, and the reason
    is exact rather than approximate: the resultant's modulus is invariant to
    translation, so the null distribution of a window's maximum resultant
    depends only on the **gap multiset** - not on its order, and therefore not
    on the observation. The register null holds that multiset fixed, so a rule
    stated in terms of it selects the same window set for the observation and
    for every permutation, and a permutation test whose statistic is the same
    function on every draw is exactly valid.

    The one piece of care needed is that an estimate taken from the null sample
    itself is a statistic of that sample. So the estimate comes from a
    **pilot** sample with its own seed, drawn independently, and the pilot seed
    is part of the statistic's definition. Below `MAX_EXACT_ORDERINGS` the
    enumeration is complete and the selection is deterministic, with no pilot
    involved at all.

    Returns ``(kept sizes, pilot null mean per size)``; an empty tuple of kept
    sizes means every window length saturates on this lining, which is
    *unresolved* and never a p of 1.
    """
    ordered = sorted({int(p) for p in positions})
    if len(ordered) < MIN_POSITIONS:
        raise ArchitectureError(
            f"{len(ordered)} positions is too few for a register null; "
            f"{MIN_POSITIONS} are needed")
    usable = [int(s) for s in sorted({int(s) for s in sizes}) if int(s) <= len(ordered)]
    if not usable:
        raise ArchitectureError(
            f"no window of {tuple(sizes)} fits a lining of {len(ordered)} positions")
    g = gaps(ordered)
    start = ordered[0]
    grid = band_grid(band, int(ordered[-1] - ordered[0]), probes)
    orderings = distinct_orderings(ordered)
    if orderings <= MAX_EXACT_ORDERINGS:
        orders = np.asarray(list(multiset_permutations(g)), dtype=float)
    else:
        rng = np.random.default_rng(pilot_seed)
        orders = np.asarray([rng.permutation(g) for _ in range(int(pilot_draws))],
                            dtype=float)
    built = start + np.cumsum(orders, axis=1)
    stacked = np.concatenate(
        [np.full((orders.shape[0], 1), float(start)), built], axis=1)
    means = {size: float(values.mean())
             for size, values in windowed_maxima(stacked, grid, usable).items()}
    if rule == "unsaturated":
        kept = tuple(size for size in usable if means[size] < SATURATED_NULL)
    elif rule == "best":
        # The lowest null mean is the size with the most room left between the
        # null and 1, which is where power can live. A rank, not a threshold,
        # so it always selects something - and it is still a function of the
        # gap multiset alone, which is the whole of the validity argument.
        kept = (min(usable, key=lambda size: (means[size], size)),)
    else:
        raise ArchitectureError(f"unknown selection rule {rule!r}")
    return kept, means


@dataclass(frozen=True)
class LocalRegisterResult:
    """Local phase coherence beside the same gap-permutation null.

    This is a *separate* test from `RegisterResult` with its own Holm family,
    as the amendment that introduced it requires, so adding it changes no
    published number of the pooled statistic.
    """

    band: Tuple[float, float]
    observed: float
    period: float
    window_size: int
    window_start: int
    window_span: Tuple[int, int]
    null_mean: float
    null_sd: float
    null_mean_by_size: Dict[int, float]
    z: float
    p_one_sided: float
    n_draws: int
    n_distinct_orderings: int
    sizes: Tuple[int, ...]
    exact: bool = False
    selection: str = "all"
    pruned_out: Tuple[int, ...] = ()
    pilot_mean_by_size: Dict[int, float] = field(default_factory=dict)
    pilot_seed: Optional[int] = None

    @property
    def minimum_attainable_p(self) -> float:
        return 1.0 / max(1, self.n_distinct_orderings)

    @property
    def saturated(self) -> Tuple[int, ...]:
        """Window sizes whose null mean exceeds `SATURATED_NULL`.

        A saturated size cannot contribute power: if a random reordering
        already reaches the same coherence, no observation can look unusual
        against it. Reported rather than silently dropped, because dropping it
        after seeing the result is the selection this module exists to avoid.
        """
        return tuple(size for size, mean in sorted(self.null_mean_by_size.items())
                     if mean >= SATURATED_NULL)

    @property
    def verdict(self) -> str:
        if self.n_distinct_orderings < MIN_ORDERINGS:
            return "unresolved"
        if self.p_one_sided <= CLEAR_P:
            return "clear"
        return "weak" if self.p_one_sided <= WEAK_P else "none"

    def as_dict(self) -> Dict[str, object]:
        return {
            "band": [float(self.band[0]), float(self.band[1])],
            "observed": float(self.observed),
            "period": float(self.period),
            "window_size": int(self.window_size),
            "window_start": int(self.window_start),
            "window_span": [int(self.window_span[0]), int(self.window_span[1])],
            "null_mean": float(self.null_mean),
            "null_sd": float(self.null_sd),
            "null_mean_by_size": {int(k): float(v) for k, v in self.null_mean_by_size.items()},
            "z": float(self.z),
            "p_one_sided": float(self.p_one_sided),
            "n_draws": int(self.n_draws),
            "n_distinct_orderings": int(self.n_distinct_orderings),
            "sizes": [int(s) for s in self.sizes],
            "exact": bool(self.exact),
            "selection": str(self.selection),
            "pruned_out": [int(s) for s in self.pruned_out],
            "pilot_mean_by_size": {int(k): float(v) for k, v in self.pilot_mean_by_size.items()},
            "pilot_seed": None if self.pilot_seed is None else int(self.pilot_seed),
            "minimum_attainable_p": float(self.minimum_attainable_p),
            "saturated": [int(s) for s in self.saturated],
            "verdict": self.verdict,
        }


def local_register(
    positions: Sequence[int],
    band: Sequence[float] = HELIX_BAND,
    *,
    sizes: Sequence[int] = WINDOW_SIZES,
    probes: Optional[int] = None,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = 20261012,
    selection: str = "all",
    pilot_draws: int = PILOT_DRAWS,
    pilot_seed: int = 20261013,
) -> LocalRegisterResult:
    """Is any *stretch* of the lining more phase-coherent than its gaps allow?

    The pooled resultant is the Neyman-Pearson statistic for a face of one
    period and one unknown phase. A real helix with backbone scatter does not
    have one phase: the rise per residue drifts, phase error accumulates along
    the chain, and a bulge resets it. This asks the weaker question the
    structure can actually answer - whether the lining is **piecewise** phase
    locked - by maximising over windows of consecutive lining positions.

    The maximum over windows and sizes is taken *inside* the statistic and the
    identical maximum is taken on every permutation of the null, so there is
    nothing to correct for afterwards. The null is the same gap-multiset
    permutation null as `register`, which holds n, span, occupancy, gap entropy
    and the whole topological half exactly fixed.
    """
    ordered = sorted({int(p) for p in positions})
    if len(ordered) < MIN_POSITIONS:
        raise ArchitectureError(
            f"{len(ordered)} positions is too few for a register null; "
            f"{MIN_POSITIONS} are needed and five gaps give 120 orderings"
        )
    usable = [int(s) for s in sorted({int(s) for s in sizes}) if int(s) <= len(ordered)]
    if not usable:
        raise ArchitectureError(
            f"no window of {tuple(sizes)} fits a lining of {len(ordered)} positions")
    if selection not in SELECTIONS:
        raise ArchitectureError(
            f"selection must be one of {SELECTIONS}, got {selection!r}")
    pilot_means: Dict[int, float] = {}
    dropped: Tuple[int, ...] = ()
    if selection != "all":
        kept, pilot_means = pilot_sizes(
            ordered, band, sizes=usable, probes=probes,
            pilot_draws=pilot_draws, pilot_seed=pilot_seed, rule=selection)
        dropped = tuple(s for s in usable if s not in kept)
        if not kept:
            raise ArchitectureError(
                "every window size saturates its own null on this lining "
                f"(pilot means {pilot_means}), so local register is unresolved")
        usable = list(kept)
    g = gaps(ordered)
    start = ordered[0]
    grid = band_grid(band, int(ordered[-1] - ordered[0]), probes)
    best = best_window(ordered, band, sizes=usable, probes=probes)
    observed = float(best["resultant"])
    orderings = distinct_orderings(ordered)

    def build(block: np.ndarray) -> np.ndarray:
        built = start + np.cumsum(block, axis=1)
        return np.concatenate([np.full((block.shape[0], 1), float(start)), built], axis=1)

    exact = orderings <= MAX_EXACT_ORDERINGS
    if exact:
        orders = np.asarray(list(multiset_permutations(g)), dtype=float)
        used = int(orders.shape[0])
    else:
        rng = np.random.default_rng(seed)
        orders = np.asarray([rng.permutation(g) for _ in range(int(n_draws))], dtype=float)
        used = int(n_draws)
    per_size = windowed_maxima(build(orders), grid, usable)
    draws = np.max(np.stack([per_size[s] for s in usable]), axis=0)
    if exact:
        p_value = float(np.mean(draws >= observed - 1e-12))
    else:
        p_value = float((int(np.sum(draws >= observed - 1e-12)) + 1) / (int(n_draws) + 1))
    mean = float(draws.mean())
    sd = float(draws.std(ddof=1)) if draws.size > 1 else 0.0
    return LocalRegisterResult(
        band=(float(band[0]), float(band[1])),
        observed=observed, period=float(best["period"]),
        window_size=int(best["size"]), window_start=int(best["start"]),
        window_span=(int(best["first_position"]), int(best["last_position"])),
        null_mean=mean, null_sd=sd,
        null_mean_by_size={int(s): float(per_size[s].mean()) for s in usable},
        z=0.0 if sd <= 1e-12 else float((observed - mean) / sd),
        p_one_sided=p_value, n_draws=used, n_distinct_orderings=orderings,
        sizes=tuple(usable), exact=exact,
        selection=selection, pruned_out=dropped,
        pilot_mean_by_size=pilot_means,
        pilot_seed=int(pilot_seed) if selection != "all" else None,
    )


@dataclass(frozen=True)
class LocalArchitecture:
    """The local statistic over both bands, with its **own** Holm family.

    Kept separate from `LiningArchitecture` rather than added to it, because
    folding a fifth and sixth test into that family would change every
    published verdict of the pooled statistic. Two detectors, two families.
    """

    helix: Optional[LocalRegisterResult] = None
    strand: Optional[LocalRegisterResult] = None
    reason: Optional[str] = None

    @property
    def family(self) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for name, result in (("helix_local", self.helix), ("strand_local", self.strand)):
            if result is not None:
                out[name] = float(result.p_one_sided)
        return out

    @property
    def holm(self) -> Dict[str, float]:
        """Holm over this statistic's own family, by the same rule as the pooled one."""
        raw = self.family
        if not raw:
            return {}
        items = sorted(raw.items(), key=lambda kv: kv[1])
        size = len(items)
        out: Dict[str, float] = {}
        running = 0.0
        for index, (name, value) in enumerate(items):
            running = max(running, min(1.0, (size - index) * value))
            out[name] = running
        return out

    @property
    def local_p(self) -> Optional[float]:
        corrected = self.holm
        return min(corrected.values()) if corrected else None

    @property
    def verdict(self) -> str:
        p = self.local_p
        if p is None:
            return "unresolved"
        if p <= CLEAR_P:
            return "clear"
        return "weak" if p <= WEAK_P else "none"

    @property
    def saturated(self) -> Tuple[int, ...]:
        """Window sizes saturated in any band that ran."""
        hit: Set[int] = set()
        for result in (self.helix, self.strand):
            if result is not None:
                hit.update(result.saturated)
        return tuple(sorted(hit))

    def as_dict(self) -> Dict[str, object]:
        return {
            "helix": None if self.helix is None else self.helix.as_dict(),
            "strand": None if self.strand is None else self.strand.as_dict(),
            "reason": self.reason,
            "holm": {k: float(v) for k, v in self.holm.items()},
            "local_p": None if self.local_p is None else float(self.local_p),
            "verdict": self.verdict,
            "saturated": [int(s) for s in self.saturated],
        }


def local_architecture(
    positions: Sequence[int],
    *,
    sizes: Sequence[int] = WINDOW_SIZES,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = 20261012,
    selection: str = "all",
) -> LocalArchitecture:
    """The local register statistic over both bands, or a stated reason it could not run.

    The gate is the same one `describe_architecture` applies to the pooled
    statistic, for the same reason: below `MIN_ORDERINGS` distinct gap
    orderings the smallest attainable p is too large for the question to be
    asked, and that is *unresolved*, never a p of 1.
    """
    ordered = sorted({int(p) for p in positions})
    if len(ordered) < MIN_POSITIONS:
        return LocalArchitecture(reason=(
            f"{len(ordered)} positions; {MIN_POSITIONS} are needed for a register null"))
    orderings = distinct_orderings(ordered)
    if orderings < MIN_ORDERINGS:
        return LocalArchitecture(reason=(
            "every gap is the same width, so the gap multiset has one ordering "
            "and register cannot be separated from topology"
            if orderings < 2 else
            f"the gap multiset has {orderings} orderings, so the smallest "
            f"attainable p is {1.0 / orderings:.3f} and significance at "
            f"{WEAK_P} is out of reach"))
    try:
        helix = local_register(ordered, HELIX_BAND, sizes=sizes, n_draws=n_draws,
                               seed=seed, selection=selection)
        strand = local_register(ordered, STRAND_BAND, sizes=sizes, n_draws=n_draws,
                                seed=seed + 1, selection=selection)
    except ArchitectureError as error:
        return LocalArchitecture(reason=str(error))
    return LocalArchitecture(helix=helix, strand=strand)


# ------------------------------------------------- the topology null

#: Gap widths that keep a residue on the same face of an alpha helix. A helix
#: repeats every 3.6 residues, so the next residue on the same face is 3 or 4
#: along. This is the *topological* signature of a helical face.
HELIX_FACE_GAPS = (3, 4)

#: The same for an extended strand, which alternates faces every residue.
STRAND_FACE_GAPS = (2,)


def face_gap_fraction(positions: Sequence[int], widths: Sequence[int] = HELIX_FACE_GAPS) -> float:
    """Fraction of gaps whose width keeps the lining on one face."""
    g = gaps(positions)
    wanted = {int(w) for w in widths}
    return float(np.mean([int(x) in wanted for x in g]))


@dataclass(frozen=True)
class TopologyResult:
    """A face-spacing fraction beside the uniform-in-span null."""

    widths: Tuple[int, ...]
    observed: float
    null_mean: float
    null_sd: float
    z: float
    p_one_sided: float
    n_draws: int

    @property
    def face_like(self) -> bool:
        return self.p_one_sided <= WEAK_P

    @property
    def verdict(self) -> str:
        if self.p_one_sided <= CLEAR_P:
            return "clear"
        return "weak" if self.p_one_sided <= WEAK_P else "none"

    def as_dict(self) -> Dict[str, object]:
        return {
            "widths": list(self.widths),
            "observed": float(self.observed),
            "null_mean": float(self.null_mean),
            "null_sd": float(self.null_sd),
            "z": float(self.z),
            "p_one_sided": float(self.p_one_sided),
            "n_draws": int(self.n_draws),
            "face_like": self.face_like,
            "verdict": self.verdict,
        }


def topology(
    positions: Sequence[int],
    widths: Sequence[int] = HELIX_FACE_GAPS,
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = 20261012,
) -> TopologyResult:
    """Is the lining's *gap multiset* face-like, against a uniform-in-span null?

    The complement of `register`, and the two together are complete. `register`
    permutes the gap multiset and so holds the topology fixed; this resamples
    positions uniformly within the same span and so holds n and span fixed
    while destroying the gap multiset. Between them they decompose the whole of
    a lining's sequence-space structure.

    This is the test that actually detects a secondary-structure face, and that
    is not obvious - it was measured. A helical face whose gaps are drawn at
    random from {3, 4} has **no register**: permuting those gaps produces
    another helical face, so the ordering carries nothing and the register test
    correctly returns p = 0.22. The periodicity of such a lining lives entirely
    in its gap multiset. Reading face-ness off `register` would therefore miss
    every helical face whose spacing is irregular, which is most of them.
    """
    ordered = sorted({int(p) for p in positions})
    if len(ordered) < 3:
        raise ArchitectureError("a topology null needs at least three positions")
    observed = face_gap_fraction(ordered, widths)
    n, lining_span = len(ordered), ordered[-1] - ordered[0]
    if lining_span + 1 <= n:
        raise ArchitectureError(
            "the lining is contiguous, so a uniform-in-span null has one draw "
            "and the gap multiset is forced"
        )
    rng = np.random.default_rng(seed)
    wanted = {int(w) for w in widths}
    interior = np.arange(ordered[0] + 1, ordered[-1])
    draws = np.empty(int(n_draws))
    for index in range(int(n_draws)):
        # Endpoints are held so the span is exactly preserved; the interior
        # positions are redrawn.
        middle = rng.choice(interior, size=n - 2, replace=False)
        rebuilt = np.sort(np.concatenate([[ordered[0], ordered[-1]], middle]))
        g = np.diff(rebuilt)
        draws[index] = float(np.mean([int(x) in wanted for x in g]))
    mean = float(draws.mean())
    sd = float(draws.std(ddof=1))
    extreme = int(np.sum(draws >= observed - 1e-12))
    return TopologyResult(
        widths=tuple(int(w) for w in widths),
        observed=observed, null_mean=mean, null_sd=sd,
        z=0.0 if sd <= 1e-12 else float((observed - mean) / sd),
        p_one_sided=float((extreme + 1) / (int(n_draws) + 1)),
        n_draws=int(n_draws),
    )


# ---------------------------------------------------------- the categories

#: The architecture classes. Two axes: how many elements the lining is drawn
#: from (the run count, which the H0 barcode gives exactly) and whether its gap
#: multiset is a secondary-structure face (the topology test). Register is
#: reported as a third, finer coordinate - order beyond what the gaps force -
#: and is deliberately not a classification axis, because a helical face with
#: irregular 3/4 spacing has no register at all and would be misclassified by
#: it.
CLASSES: Dict[str, str] = {
    "phase_locked": "one or two runs on a secondary-structure face: a single element lines the pocket",
    "mosaic": "three or more runs, each on a face: several elements each contribute one",
    "segmental": "one or two runs, spacing not face-like: a contiguous stretch or loop",
    "distributed": "three or more runs, spacing not face-like: scattered, interface-like",
    "unresolved": "neither half of the decomposition can be tested on this lining",
}

#: `fragmentation` is reported as a continuous coordinate but is **not** the
#: classification axis; see `LiningArchitecture.architecture_class`.


@dataclass(frozen=True)
class LiningArchitecture:
    """One lining's sequence-space architecture: topology, register, class.

    The two tests are complementary halves of a decomposition, not two
    attempts at the same thing. `helix_topology`/`strand_topology` ask whether
    the gap multiset is a face; `helix`/`strand` ask whether the *order* of
    those gaps adds coherence the multiset does not already force.
    """

    n_positions: int
    span: int
    segments: int
    fragmentation: float
    occupancy: float
    gap_entropy: float
    dominant: Dict[str, float]
    helix_topology: Optional[TopologyResult] = None
    strand_topology: Optional[TopologyResult] = None
    helix: Optional[RegisterResult] = None
    strand: Optional[RegisterResult] = None
    unresolved_reason: Optional[str] = None
    topology_reason: Optional[str] = None
    register_reason: Optional[str] = None

    # ------------------------------------------------- the test family

    @property
    def family(self) -> Dict[str, float]:
        """Raw one-sided p-values of every test that actually ran."""
        out: Dict[str, float] = {}
        for name, result in (("helix_topology", self.helix_topology),
                             ("strand_topology", self.strand_topology),
                             ("helix_register", self.helix),
                             ("strand_register", self.strand)):
            if result is not None:
                out[name] = float(result.p_one_sided)
        return out

    @property
    def holm(self) -> Dict[str, float]:
        """Holm-corrected p-values over the tests that ran.

        Four tests are asked of every lining - two face widths, two bands - and
        taking the best of them inflates the false-positive rate from a nominal
        0.05 to a **measured 0.136** on 500 random linings. Each test alone is
        calibrated or conservative (0.028, 0.024, 0.044, 0.044); it is the
        maximum over the family that is not. Holm over the family that actually
        ran brings it back to 0.030 at nominal 0.05 and 0.006 at 0.01, so every
        verdict below reads the corrected p and never the raw one.
        """
        raw = self.family
        if not raw:
            return {}
        items = sorted(raw.items(), key=lambda kv: kv[1])
        size = len(items)
        out: Dict[str, float] = {}
        running = 0.0
        for index, (name, value) in enumerate(items):
            running = max(running, min(1.0, (size - index) * value))
            out[name] = running
        return out

    @staticmethod
    def _grade(p: Optional[float]) -> str:
        if p is None:
            return "unresolved"
        if p <= CLEAR_P:
            return "clear"
        return "weak" if p <= WEAK_P else "none"

    # ----------------------------------------------------------- topology

    @property
    def best_topology(self) -> Optional[TopologyResult]:
        candidates = [r for r in (self.helix_topology, self.strand_topology) if r is not None]
        return min(candidates, key=lambda r: r.p_one_sided) if candidates else None

    @property
    def topology_p(self) -> Optional[float]:
        """The smaller Holm-corrected p of the two face-width tests."""
        corrected = self.holm
        values = [corrected[k] for k in ("helix_topology", "strand_topology") if k in corrected]
        return min(values) if values else None

    @property
    def topology_verdict(self) -> str:
        return self._grade(self.topology_p)

    @property
    def face_like(self) -> bool:
        p = self.topology_p
        return bool(p is not None and p <= WEAK_P)

    # ----------------------------------------------------------- register

    @property
    def best_register(self) -> Optional[RegisterResult]:
        candidates = [r for r in (self.helix, self.strand) if r is not None]
        return max(candidates, key=lambda r: r.z) if candidates else None

    @property
    def register_z(self) -> float:
        best = self.best_register
        return 0.0 if best is None else best.z

    @property
    def register_p(self) -> Optional[float]:
        """The smaller Holm-corrected p of the two band register tests."""
        corrected = self.holm
        values = [corrected[k] for k in ("helix_register", "strand_register") if k in corrected]
        return min(values) if values else None

    @property
    def register_verdict(self) -> str:
        return self._grade(self.register_p)

    # -------------------------------------------------- the graded answer

    @property
    def signature_verdict(self) -> str:
        """``unresolved``, ``none``, ``weak`` or ``clear``: the yes/sort-of/no scale.

        The stronger of the two halves, because either one is a sequence-space
        signature. Both are also reported separately, since they mean different
        things: a face in the gap multiset is a statement about spacing, and
        register is a statement about arrangement.
        """
        corrected = self.holm
        return self._grade(min(corrected.values()) if corrected else None)

    @property
    def element(self) -> str:
        """The secondary-structure name the lining's spacing earns, or ``none``.

        A significant band result does not by itself earn a name: a band
        maximum is confounded by subharmonics - a 3/4-spaced face peaks equally
        at 3.50, 7.00 and 2.80 - so the global maximum over a wide scan must
        also fall inside the band before the name is used.
        """
        if not self.face_like:
            return "none"
        corrected = self.holm
        helix_p = corrected.get("helix_topology", 1.0)
        strand_p = corrected.get("strand_topology", 1.0)
        if helix_p <= strand_p:
            return "helix_face" if self.dominant.get("in_helix_band") else "periodic_other"
        return "strand_face" if self.dominant.get("in_strand_band") else "periodic_other"

    @property
    def architecture_class(self) -> str:
        if self.topology_verdict == "unresolved" and self.register_verdict == "unresolved":
            return "unresolved"
        many = self.segments >= MULTI_ELEMENT_SEGMENTS
        if self.face_like:
            return "mosaic" if many else "phase_locked"
        return "distributed" if many else "segmental"

    def as_dict(self) -> Dict[str, object]:
        return {
            "n_positions": self.n_positions,
            "span": self.span,
            "segments": self.segments,
            "fragmentation": float(self.fragmentation),
            "occupancy": float(self.occupancy),
            "gap_entropy": float(self.gap_entropy),
            "dominant_period": float(self.dominant.get("period", float("nan"))),
            "dominant_resultant": float(self.dominant.get("resultant", float("nan"))),
            "helix_topology": None if self.helix_topology is None else self.helix_topology.as_dict(),
            "strand_topology": None if self.strand_topology is None else self.strand_topology.as_dict(),
            "helix_register": None if self.helix is None else self.helix.as_dict(),
            "strand_register": None if self.strand is None else self.strand.as_dict(),
            "topology_verdict": self.topology_verdict,
            "register_verdict": self.register_verdict,
            "signature_verdict": self.signature_verdict,
            "p_raw": self.family,
            "p_holm": self.holm,
            "element": self.element,
            "architecture_class": self.architecture_class,
            "class_meaning": CLASSES[self.architecture_class],
            "topology_reason": self.topology_reason,
            "register_reason": self.register_reason,
        }


def describe_architecture(
    positions: Sequence[int],
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = 20261012,
    face_gap: int = FACE_GAP,
) -> LiningArchitecture:
    """Read a lining's architecture: both halves of the decomposition, and a class.

    Each half can fail independently and for a stated reason, and a half that
    could not be tested is reported as unresolved rather than as absent
    structure. The two are different claims: a contiguous lining has no
    testable gap multiset *and* no register, while a perfectly regular lining
    has a testable multiset and no register at all.
    """
    ordered = sorted({int(p) for p in positions})
    if len(ordered) < 2:
        raise ArchitectureError("a lining needs at least two distinct positions")
    shared = dict(
        n_positions=len(ordered), span=span(ordered),
        segments=segmentation(ordered, face_gap),
        fragmentation=fragmentation(ordered, face_gap),
        occupancy=occupancy(ordered), gap_entropy=gap_entropy(ordered),
        dominant=dominant_period(ordered),
    )

    helix_topology = strand_topology = None
    topology_reason = None
    try:
        helix_topology = topology(ordered, HELIX_FACE_GAPS, n_draws=n_draws, seed=seed)
        strand_topology = topology(ordered, STRAND_FACE_GAPS, n_draws=n_draws, seed=seed + 2)
    except ArchitectureError as error:
        topology_reason = str(error)

    helix_register = strand_register = None
    register_reason = None
    if len(ordered) < MIN_POSITIONS:
        register_reason = (
            f"{len(ordered)} positions; {MIN_POSITIONS} are needed for a register null"
        )
    else:
        orderings = distinct_orderings(ordered)
        if orderings < MIN_ORDERINGS:
            register_reason = (
                "every gap is the same width, so the gap multiset has one ordering "
                "and register cannot be separated from topology"
                if orderings < 2 else
                f"the gap multiset has {orderings} orderings, so the smallest "
                f"attainable p is {1.0 / orderings:.3f} and significance at "
                f"{WEAK_P} is out of reach"
            )
        else:
            helix_register = register(ordered, HELIX_BAND, n_draws=n_draws, seed=seed)
            strand_register = register(ordered, STRAND_BAND, n_draws=n_draws, seed=seed + 1)

    unresolved = None
    if topology_reason is not None and register_reason is not None:
        unresolved = f"topology: {topology_reason}; register: {register_reason}"
    return LiningArchitecture(
        **shared,
        helix_topology=helix_topology, strand_topology=strand_topology,
        helix=helix_register, strand=strand_register,
        unresolved_reason=unresolved,
        topology_reason=topology_reason, register_reason=register_reason,
    )
