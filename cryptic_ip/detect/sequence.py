"""Information-theoretic descriptors of a pocket's sequence, done honestly.

The question
-----------
Is there a signature in the sequence at an annotated binding site - in its
entropy, its compressibility, its periodicity - that says a pocket of a
particular kind forms there? It is a good question and there are three ways to
answer it wrongly, all of which this module is built around.

**Kolmogorov complexity is uncomputable.** Nothing here computes it and nothing
here is named after it. What is computable is an upper bound: the length of a
description a particular algorithm produces. `lempel_ziv_complexity` is the
exact LZ76 parse count, which needs no compressor and has no header;
`compressed_bits` calls a real compressor and subtracts its measured empty-input
overhead. Both bound K from above up to a constant, and neither is K.

**Plug-in Shannon entropy is badly biased on a window the size of a pocket
lining.** Mean bias over 4000 uniform random windows per row, seed 20261012,
against a uniform source over the 20 amino acids whose true entropy is 4.322
bits:

    window N   plug-in   Miller-Madow   Grassberger
          10    -1.408         -0.900        -0.821
          20    -0.799         -0.372        -0.316
          50    -0.302         -0.050        -0.029
         100    -0.144         -0.008        -0.002

A pocket lining is 10 to 20 residues. At that size plug-in entropy understates
the truth by 0.8 to 1.4 bits, and - this is the part that ruins a comparison -
*the bias depends on N*. Comparing a 12-residue lining with a 20-residue one
using plug-in entropy reports a difference that is pure sample size. Grassberger
is the default and removes most of it, and window length is still controlled
explicitly rather than trusted to the estimator.

**Composition is the confounder that matters.** These sites are arginine- and
lysine-rich; the project established that long ago. Any entropy or complexity
measure of a window is a function of its composition, so a "signature" can be
nothing but "this window has a lot of arginine" wearing a different name. So
the headline statistic in this module is not a raw measure at all: it is
`composition_z`, the measure's z-score against shuffles of the same window.
That null holds both the length and the letter counts fixed, which makes the
z-score free of the small-sample bias above *by construction* and makes it
answer the only interesting version of the question - whether the ordering
carries information that the composition does not.

What a signature would have to survive
--------------------------------------
`docs/SEQUENCE_PLAN.md` pre-registers the study. The short form: a signature is
a finding only if it survives composition-matched controls, is estimated with
whole homology groups resampled, clears a minimum detectable effect computed
before any result is read, and holds after Holm correction across the features
tried. A difference in a raw measure between sites and non-sites is not a
finding, and this module makes it awkward to report one.
"""

from __future__ import annotations

import bz2
import lzma
import zlib
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.special import digamma

LN2 = float(np.log(2.0))

#: The twenty amino acids, as the identity alphabet.
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

#: Reduced alphabets. A signature may live in a coarse-grained view rather than
#: in the 20-letter one: what matters at a polyanion site is where the charges
#: are, not whether a given basic residue is arginine or lysine. Each is a
#: hypothesis about which distinction carries the information, and they are
#: tried in parallel with a Holm correction rather than picked after the fact.
ALPHABETS: Dict[str, Dict[str, str]] = {
    "identity": {letter: letter for letter in AMINO_ACIDS},
    # Charge at pH 7.4. Histidine is its own class because it is the one
    # residue whose charge is genuinely ambiguous there.
    "charge": {
        **{c: "+" for c in "RK"}, "H": "h",
        **{c: "-" for c in "DE"},
        **{c: "0" for c in "ACFGILMNPQSTVWY"},
    },
    # Kyte-Doolittle sign, plus a class for the two that break helices.
    "hydrophobicity": {
        **{c: "H" for c in "AVILMFCW"},
        **{c: "P" for c in "STNQY"},
        **{c: "C" for c in "DEKR"},
        **{c: "S" for c in "GHP"},
    },
    # Volume classes, because a lining's packing is a size question.
    "size": {
        **{c: "s" for c in "AGSCP"},
        **{c: "m" for c in "DNTVIL"},
        **{c: "l" for c in "EQMHKFRYW"},
    },
}

#: Compressor overhead on empty input, in bytes, measured at import. On a
#: twenty-residue string gzip emits 31 bytes of which 20 are header, so a raw
#: compression ratio on a pocket lining is mostly measuring the header. These
#: are subtracted, and LZ76 - which has no header at all - is the default.
COMPRESSOR_OVERHEAD: Dict[str, int] = {
    "zlib": len(zlib.compress(b"")),
    "bz2": len(bz2.compress(b"")),
    "lzma": len(lzma.compress(b"")),
}

#: Helical and extended periodicities, in residues per turn, kept as the band
#: centres. A helix face pointing into a pocket puts its lining residues 3.6
#: apart; a strand alternates, so 2.0.
HELIX_PERIOD = 3.6
STRAND_PERIOD = 2.0

#: Shuffles for the composition-preserving null. 2000 matches the project's
#: bootstrap convention.
DEFAULT_SHUFFLES = 2000

#: Below this many homology groups on either side, a comparison is not evidence
#: whichever way it points. The project's standing floor, named here rather than
#: written as a literal inside `is_evidence`, so a report can state the floor it
#: is being judged against instead of repeating the number.
MIN_GROUPS = 5


class SequenceError(ValueError):
    """A sequence or window from which no descriptor can be computed."""


# --------------------------------------------------------------- alphabets

def recode(sequence: str, alphabet: str = "identity") -> str:
    """Map ``sequence`` into one of the reduced alphabets.

    Anything outside the twenty standard residues - an X, a gap, a selenocysteine
    - is dropped rather than bucketed, because assigning it a class would invent
    information. The number dropped is available from the length difference.
    """
    table = ALPHABETS.get(alphabet)
    if table is None:
        raise SequenceError(f"unknown alphabet {alphabet!r}; known: {sorted(ALPHABETS)}")
    return "".join(table[c] for c in str(sequence).upper() if c in table)


def symbol_counts(sequence: str) -> np.ndarray:
    """Counts of each distinct symbol, descending. Empty for an empty sequence."""
    if not sequence:
        return np.zeros(0, dtype=int)
    _, counts = np.unique(np.frombuffer(sequence.encode("ascii", "ignore"), dtype=np.uint8),
                          return_counts=True)
    return np.sort(counts)[::-1]


# ----------------------------------------------------------------- entropy

def shannon_entropy(sequence: str, *, estimator: str = "grassberger") -> float:
    """Shannon entropy in bits, with a small-sample correction.

    ``plugin`` is the naive frequency estimate and is here only so the bias can
    be shown. ``miller_madow`` adds the first-order correction
    ``(K - 1) / 2N``. ``grassberger`` uses ``psi(N) - sum n_i psi(n_i) / N``,
    which measurement in this module's docstring shows is the best of the three
    at the window sizes a pocket lining actually has.
    """
    counts = symbol_counts(sequence)
    if counts.size == 0:
        raise SequenceError("no symbols to measure")
    if counts.size == 1:
        return 0.0
    total = float(counts.sum())
    if estimator == "plugin":
        p = counts / total
        return float(-np.sum(p * np.log2(p)))
    if estimator == "miller_madow":
        p = counts / total
        return float(-np.sum(p * np.log2(p)) + (counts.size - 1) / (2.0 * total * LN2))
    if estimator == "grassberger":
        return float((digamma(total) - np.sum(counts * digamma(counts)) / total) / LN2)
    raise SequenceError(f"unknown estimator {estimator!r}")


def block_entropy(sequence: str, k: int, *, estimator: str = "grassberger") -> float:
    """Entropy of overlapping ``k``-mers, in bits.

    The number of distinct k-mers grows as the alphabet size to the k, so for
    k above about 2 on a twenty-residue window every k-mer is unique and the
    estimate saturates at log2(N - k + 1). That saturation is not structure and
    `block_entropy_is_saturated` says when it has happened.
    """
    if k < 1:
        raise SequenceError(f"k must be at least 1, got {k}")
    if len(sequence) < k:
        raise SequenceError(f"a {len(sequence)}-symbol sequence has no {k}-mers")
    blocks = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    counts = np.array(sorted((blocks.count(b) for b in set(blocks)), reverse=True), dtype=int)
    total = float(counts.sum())
    if estimator == "grassberger":
        if counts.size == 1:
            return 0.0
        return float((digamma(total) - np.sum(counts * digamma(counts)) / total) / LN2)
    p = counts / total
    naive = float(-np.sum(p * np.log2(p)))
    if estimator == "miller_madow":
        return naive + (counts.size - 1) / (2.0 * total * LN2)
    return naive


def block_entropy_is_saturated(sequence: str, k: int) -> bool:
    """Whether every ``k``-mer in the window is unique, so the estimate is a count.

    Reported rather than hidden: a saturated block entropy carries no
    information about the sequence beyond its length.
    """
    if len(sequence) < k:
        return True
    blocks = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    return len(set(blocks)) == len(blocks)


def entropy_rate(sequence: str, k: int = 2, *, estimator: str = "grassberger") -> float:
    """``H(k) - H(k-1)``: the entropy of the next symbol given the previous ones.

    This is the measure that separates structure from composition. A window
    whose letters are shuffled has the same single-symbol entropy and a higher
    entropy rate, because shuffling destroys whatever the ordering carried.
    """
    if k < 2:
        raise SequenceError("an entropy rate needs k of at least 2")
    return block_entropy(sequence, k, estimator=estimator) - block_entropy(
        sequence, k - 1, estimator=estimator)


def excess_entropy(sequence: str, k_max: int = 3, *, estimator: str = "grassberger") -> float:
    """How much the block entropy exceeds a memoryless source of the same rate.

    Sometimes called predictive information. Positive when the sequence has
    structure a first-order model would miss.
    """
    if k_max < 2:
        raise SequenceError("an excess entropy needs k_max of at least 2")
    rate = entropy_rate(sequence, k_max, estimator=estimator)
    return float(block_entropy(sequence, k_max, estimator=estimator) - k_max * rate)


# -------------------------------------------------------------- complexity

def lempel_ziv_phrases(sequence: str) -> List[str]:
    """The LZ76 exhaustive-history parse of ``sequence``, phrase by phrase.

    Written straight from the definition: each phrase is the shortest string
    starting where the last one ended that has not already appeared as a
    substring of everything scanned before its own final symbol. The tail
    counts as one phrase even when it is not novel, which is the usual
    convention for an incomplete last phrase.

    This is spelled out rather than taken from the circulated Kaspar-Schuster
    pseudocode because the widely copied form of that snippet drops the
    longest-match bookkeeping, which makes it advance one symbol per phrase and
    report nearly ``n`` phrases for any string without an exact early repeat.
    It was used here first and agreed with every short textbook case -
    ``aaaaaaaa`` to 2, ``abababab`` to 3, ``abcdefgh`` to 8 - while scoring a
    random 500-residue window at 485 phrases against the correct 213. The
    textbook cases cannot see that fault, so the parse is returned and tested
    against its own definition instead.
    """
    text = str(sequence)
    n = len(text)
    if n == 0:
        raise SequenceError("no sequence to parse")
    phrases: List[str] = []
    start = 0
    while start < n:
        length = 1
        while start + length <= n and text[start:start + length] in text[:start + length - 1]:
            length += 1
        phrases.append(text[start:start + length])
        start += length
    return phrases


def lempel_ziv_complexity(sequence: str) -> int:
    """The LZ76 parse count: how many new phrases the sequence needs.

    This is the complexity proxy to prefer on short strings. It has no header,
    no dictionary warm-up and no implementation-defined parameters, so unlike a
    compression ratio it measures the string rather than the compressor.

    It is an upper bound on Kolmogorov complexity up to a constant, which is
    all any computable quantity can be.
    """
    return len(lempel_ziv_phrases(sequence))


def lz_normalised(sequence: str) -> float:
    """LZ76 complexity divided by its asymptotic value for a random string.

    The normaliser is ``n / log_a(n)`` with ``a`` the **observed** alphabet
    size, which matters at these lengths: a random 20-residue window contains
    12.7 of the twenty amino acids on average, not 20, so normalising by the
    nominal alphabet would be wrong by that much.

    Even so the expectation is *asymptotic* and is not reached at the length of
    a pocket lining. Measured on 300 uniform random windows each, the mean
    normalised value is 0.88 at N=20, 0.84 at N=36, 0.83 at N=100 and 0.88 at
    N=500 - below 1 and **not monotone in N**, because the observed-alphabet
    correction shrinks as the window grows while the finite-size deficit
    closes. On a binary alphabet, where the classic ``n / log2(n)`` limit
    applies cleanly, it reaches 1.04 by N=4000.

    So "near 1" is not the random baseline at the sizes this project measures,
    and - as with plug-in entropy - the residual depends on N. That makes a raw
    comparison between windows of different length unsafe. Use
    `composition_z`, whose null has the same length and the same letters.
    """
    text = str(sequence)
    n = len(text)
    if n < 2:
        raise SequenceError("too short to normalise")
    alphabet = len(set(text))
    if alphabet < 2:
        return 0.0
    expected = n / (np.log(n) / np.log(alphabet))
    return float(lempel_ziv_complexity(text) / expected)


def compressed_bits(sequence: str, method: str = "zlib") -> float:
    """Compressed size in bits, with the compressor's empty-input header removed.

    Even with the header subtracted this is unreliable on a pocket-sized window
    - a dictionary needs more than twenty symbols to warm up - and
    `COMPRESSOR_OVERHEAD` records how much of the raw output was header. Use
    `lempel_ziv_complexity` unless a whole-protein sequence is being measured.
    """
    if method not in COMPRESSOR_OVERHEAD:
        raise SequenceError(f"unknown method {method!r}; known: {sorted(COMPRESSOR_OVERHEAD)}")
    data = str(sequence).encode("ascii", "ignore")
    if not data:
        raise SequenceError("no sequence to compress")
    compress = {"zlib": zlib.compress, "bz2": bz2.compress, "lzma": lzma.compress}[method]
    return float(max(0, len(compress(data)) - COMPRESSOR_OVERHEAD[method]) * 8)


# ------------------------------------------------------------- periodicity

#: Charge per residue at pH 7.4, for the periodicity measures.
CHARGE_OF: Dict[str, float] = {
    **{c: 1.0 for c in "RK"}, "H": 0.1,
    **{c: -1.0 for c in "DE"},
    **{c: 0.0 for c in "ACFGILMNPQSTVWY"},
}

#: Kyte-Doolittle hydropathy, for the same purpose.
HYDROPATHY: Dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def as_track(sequence: str, track: str = "charge") -> np.ndarray:
    """A numeric series along the sequence: charge or hydropathy."""
    table = {"charge": CHARGE_OF, "hydropathy": HYDROPATHY}.get(track)
    if table is None:
        raise SequenceError(f"unknown track {track!r}")
    values = [table[c] for c in str(sequence).upper() if c in table]
    if not values:
        raise SequenceError("no recognised residues in the sequence")
    return np.asarray(values, dtype=float)


@lru_cache(maxsize=64)
def _grid(low: float, high: float, probes: int) -> Tuple[float, ...]:
    """A period grid as a hashable tuple, so the basis cache keys on it."""
    return tuple(float(value) for value in np.linspace(float(low), float(high), int(probes)))


#: A basis column whose norm falls below this is a quadrature that is
#: numerically absent - which is what happens to the sine at the Nyquist period
#: of 2 - and fitting it would fit noise.
QUADRATURE_FLOOR = 1e-12


@lru_cache(maxsize=512)
def _probe_basis(n: int, periods: Tuple[float, ...]) -> Tuple[np.ndarray, np.ndarray]:
    """Cos/sin bases for a grid of periods on ``n`` points, and their Gram inverses.

    Both depend only on the window length and the period grid, never on the
    series, so they are built once and reused across every draw of an order
    null. That is the whole trick: a 2000-shuffle null on a hundred-probe scan
    used to rebuild a hundred bases two thousand times over.

    The second return value is the pseudo-inverse of each probe's 2x2 Gram
    matrix. ``lstsq`` returns the minimum-norm least-squares solution, which is
    ``pinv(B) @ y``, and ``pinv(B) == pinv(B.T @ B) @ B.T`` for every ``B`` - so
    a pseudo-inverse of the Gram reproduces the fit exactly, rank-deficient
    cases included, at the cost of a 2x2 decomposition instead of an n-by-2 one.
    """
    grid = np.asarray(periods, dtype=float)
    angles = 2.0 * np.pi * np.arange(n, dtype=float)[None, :] / grid[:, None]
    basis = np.stack([np.cos(angles), np.sin(angles)], axis=-1)
    absent = np.linalg.norm(basis, axis=1) <= QUADRATURE_FLOOR
    basis = np.where(absent[:, None, :], 0.0, basis)
    gram = np.einsum("pnk,pnl->pkl", basis, basis)
    return basis, np.linalg.pinv(gram)


def probe_powers(series: np.ndarray, periods: Sequence[float]) -> np.ndarray:
    """Explained variance at every period in a grid, in one pass.

    The vectorised form of :func:`periodicity_power`. Same fit, same
    normalisation, same quadrature drop; it differs only in the order the
    arithmetic happens, so values agree to floating-point noise and
    :mod:`tests.test_detect_sequence` pins that agreement against the
    least-squares reference.
    """
    series = np.asarray(series, dtype=float)
    n = series.size
    if n < 4:
        raise SequenceError("too short for a periodicity estimate")
    grid = tuple(float(period) for period in periods)
    if not grid:
        raise SequenceError("a period grid needs at least one period")
    for period in grid:
        if period <= 1.0:
            raise SequenceError(f"period must exceed 1 residue, got {period}")
    centred = series - series.mean()
    variance = float(centred @ centred)
    if variance == 0.0:
        return np.zeros(len(grid))
    basis, inverse = _probe_basis(n, grid)
    projected = np.einsum("pnk,n->pk", basis, centred)
    coefficients = np.einsum("pkl,pl->pk", inverse, projected)
    # ||B c||**2 == c.T B.T B c == c.T (B.T y) when c solves the normal
    # equations, so the fitted series never has to be formed.
    explained = np.einsum("pk,pk->p", coefficients, projected)
    return np.minimum(1.0, np.maximum(0.0, explained / variance))


def periodicity_power(series: np.ndarray, period: float) -> float:
    """The fraction of ``series`` variance explained by a sinusoid at ``period``.

    A least-squares fit of ``a cos + b sin`` at exactly the period asked for,
    reported as explained variance. Three reasons it is done this way:

    * A Fourier bin is the wrong tool for a non-integer period on a
      twenty-point window - the nearest bin to 3.6 residues per turn is
      several residues wide.
    * Projecting onto the two quadratures and normalising by ``variance * n/2``
      - which is what this did first - is wrong at the Nyquist period of 2,
      where the sine vanishes and only one quadrature contributes. It returned
      2.0 for a pure alternating signal, and an explained-variance fraction
      cannot exceed 1.
    * Least squares handles both cases without a special case, because it
      normalises by the basis it actually fitted.
    """
    return float(probe_powers(series, (float(period),))[0])


#: The alpha-helical band, in residues per turn. A single exact probe is
#: brittle: a strictly period-4 charge pattern is *exactly* orthogonal to a 3.6
#: probe over 36 residues, so the estimator returned 0.000 for a sequence that
#: plainly has a helical face. Real helices are not exactly 3.6 either - 3.5 to
#: 3.7 for an alpha helix, 3.0 for a 3-10 and 4.1 for a pi - so the band is
#: scanned and the best fit kept. The lower edge is 3.0 and not 3.2 because a
#: 3-10 helix is exactly 3.0 residues per turn and the first version excluded
#: it, scoring a perfect period-3 charge pattern at 0.099.
HELIX_BAND = (3.0, 4.2)

#: The extended band. A strand alternates faces, so 2 residues per turn, with a
#: little room for the twist real sheets have.
STRAND_BAND = (1.9, 2.4)

#: How many probes across a band.
BAND_PROBES = 11


def band_periodicity(
    series: np.ndarray,
    band: Sequence[float] = HELIX_BAND,
    *,
    probes: int = BAND_PROBES,
) -> float:
    """The best explained variance over a band of periods.

    Scanning a band rather than one period is not a free parameter dressed up:
    the band is the physical range of the secondary structure being looked for,
    and the alternative - one infinitely sharp probe - provably misses signals
    it should find.
    """
    low, high = float(band[0]), float(band[1])
    if not low < high:
        raise SequenceError(f"a band needs low below high, got {band}")
    if probes < 2:
        raise SequenceError(f"a band needs at least 2 probes, got {probes}")
    return float(probe_powers(series, _grid(low, high, probes)).max())


#: The widest period worth scanning, in residues. Past about twelve a "period"
#: on a twenty-residue window is one cycle and means nothing.
SCAN_BAND = (2.0, 12.0)
SCAN_PROBES = 101


def dominant_period(
    series: np.ndarray,
    band: Sequence[float] = SCAN_BAND,
    *,
    probes: int = SCAN_PROBES,
) -> Dict[str, float]:
    """The period that explains most of ``series``, and how much.

    This exists because band power on its own is confounded by harmonics. A
    sparse charge pattern repeating every 7 residues has content at 3.5 as
    well, so it scored 0.368 in the helical band and looked like a helix face.
    The fundamental does dominate - a period-7 pattern peaks at 7.1, not 3.5 -
    so reporting which period wins separates a real helical face from a
    harmonic of something else.

    A high band power is therefore only interesting when the dominant period
    is *in* that band, and both numbers are reported so the distinction cannot
    be lost downstream.
    """
    series = np.asarray(series, dtype=float)
    low, high = float(band[0]), float(band[1])
    if not low < high:
        raise SequenceError(f"a band needs low below high, got {band}")
    periods = _grid(low, high, int(probes))
    powers = probe_powers(series, periods)
    best = int(np.argmax(powers))
    period = float(periods[best])
    return {
        "period": period,
        "power": float(powers[best]),
        "in_helix_band": float(HELIX_BAND[0] <= period <= HELIX_BAND[1]),
        "in_strand_band": float(STRAND_BAND[0] <= period <= STRAND_BAND[1]),
    }


def dominant_charge_period(sequence: str) -> Dict[str, float]:
    """The dominant period of a window's charge track."""
    return dominant_period(as_track(sequence, "charge"))


def helix_charge_periodicity(sequence: str) -> float:
    """Fraction of the charge variation in the helical band.

    This is the feature with the clearest structural reading: a helix whose
    face lines a pocket repeats its lining residues every 3.6 positions, so a
    basic face shows up as charge power at that period. It is also destroyed by
    shuffling, which is what makes it a candidate for a signature that
    composition cannot explain.
    """
    return band_periodicity(as_track(sequence, "charge"), HELIX_BAND)


def strand_charge_periodicity(sequence: str) -> float:
    """The same in the extended band, where residues alternate faces."""
    return band_periodicity(as_track(sequence, "charge"), STRAND_BAND)


def dominant_hydropathy_period(sequence: str) -> Dict[str, float]:
    """The dominant period of a window's hydropathy track."""
    return dominant_period(as_track(sequence, "hydropathy"))


def helix_hydropathy_periodicity(sequence: str) -> float:
    """Fraction of the hydropathy variation in the helical band.

    The charge statistics above are the right track for an inositol-phosphate
    site, which is basic. **Pocket formation in general is a hydropathy
    question**: an amphipathic helix puts one face against a cavity and the
    other against solvent, and that is a periodicity in Kyte-Doolittle
    hydropathy at the helical repeat rather than in charge.

    Declared in `docs/POCKET_SEQUENCE_PLAN.md` before any sequence was read for
    study Q9, and defined exactly as its charge counterpart on the hydropathy
    scale already in this module - no new scale, no new estimator.
    """
    return band_periodicity(as_track(sequence, "hydropathy"), HELIX_BAND)


def strand_hydropathy_periodicity(sequence: str) -> float:
    """The same in the extended band: alternating in and out along a strand."""
    return band_periodicity(as_track(sequence, "hydropathy"), STRAND_BAND)


def spacing_spectrum(positions: Sequence[int]) -> Dict[str, float]:
    """What the gaps between coordinating residues look like.

    Computable from data this project already has: the candidate records carry
    coordinating residues as identity plus position. A pocket lined by one
    helix face has gaps of 3 and 4; one lined by a strand has gaps of 2; a loop
    has neither.
    """
    ordered = sorted(int(p) for p in positions)
    if len(ordered) < 2:
        raise SequenceError("a spacing spectrum needs at least two positions")
    gaps = np.diff(np.asarray(ordered, dtype=float))
    local = gaps[gaps <= 12.0]
    out = {
        "n_positions": float(len(ordered)),
        "median_gap": float(np.median(gaps)),
        "span": float(ordered[-1] - ordered[0]),
        "fraction_local": float(len(local) / len(gaps)),
    }
    for label, wanted in (("helical", (3.0, 4.0)), ("strand", (2.0,)), ("adjacent", (1.0,))):
        out[f"fraction_{label}"] = float(np.mean([g in wanted for g in gaps]))
    return out


# ----------------------------------------------------- the composition null

#: The statistics the null is applied to, by name.
STATISTICS: Dict[str, Callable[[str], float]] = {
    "entropy": lambda s: shannon_entropy(s),
    "entropy_rate_k2": lambda s: entropy_rate(s, 2),
    "excess_entropy": lambda s: excess_entropy(s, 3),
    "lz_complexity": lambda s: float(lempel_ziv_complexity(s)),
    "lz_normalised": lz_normalised,
    "helix_charge_periodicity": helix_charge_periodicity,
    "strand_charge_periodicity": strand_charge_periodicity,
    # The dominant period's own power, which a shuffle should destroy whatever
    # the period turns out to be.
    "dominant_period_power": lambda s: dominant_charge_period(s)["power"],
    # The hydropathy track, declared in docs/POCKET_SEQUENCE_PLAN.md before any
    # sequence was read for study Q9. Charge is the track for a basic site;
    # hydropathy is the track for a pocket.
    "helix_hydropathy_periodicity": helix_hydropathy_periodicity,
    "strand_hydropathy_periodicity": strand_hydropathy_periodicity,
    "dominant_hydropathy_period_power": lambda s: dominant_hydropathy_period(s)["power"],
}


#: Names `composition_z` refuses, because they are functions of the multiset of
#: a window's letters alone and so are invariant under its null.
#:
#: Such a statistic has **zero variance** under an order-permuting shuffle, so
#: its z is 0/0. The arithmetic below returns 0.0 for a zero standard deviation,
#: which is right for a statistic that happens to be constant on one window and
#: badly wrong for one that is constant on every window: it reads as "no
#: ordering signal" when the question is undefined. None of this module's own
#: statistics is order-invariant, so nothing here ever hit it; study Q10's
#: `composition.py` adds nine that would, and `docs/COMPOSITION_PLAN.md`
#: section 5 fixes this guard as part of adding them.
#:
#: Spelled out here rather than imported so that this module keeps no
#: dependency on `composition`, which imports *it*.
#: `tests/test_detect_composition.py` fails if the two sets drift apart, and
#: also if any name in this set turns out to read order after all.
ORDER_INVARIANT_STATISTICS = frozenset({
    "basic_fraction",
    "acidic_fraction",
    "aromatic_fraction",
    "glycine_fraction",
    "proline_fraction",
    "mean_hydropathy",
    "net_charge_per_residue",
    "composition_divergence",
    "nonbasic_divergence",
})


@dataclass(frozen=True)
class NullResult:
    """A statistic beside the composition-preserving null it was compared with."""

    statistic: str
    observed: float
    null_mean: float
    null_sd: float
    z: float
    p_two_sided: float
    n_shuffles: int

    @property
    def is_composition_explained(self) -> bool:
        """Whether a shuffle of the same letters reproduces the observation.

        If it does, whatever the statistic measured is a property of the
        window's composition - which this project already knows about these
        sites - and not of the order the residues are in.
        """
        return abs(self.z) < 2.0

    def as_dict(self) -> Dict[str, object]:
        return {
            "statistic": self.statistic,
            "observed": float(self.observed),
            "null_mean": float(self.null_mean),
            "null_sd": float(self.null_sd),
            "z": float(self.z),
            "p_two_sided": float(self.p_two_sided),
            "n_shuffles": int(self.n_shuffles),
            "composition_explained": self.is_composition_explained,
        }


def composition_z(
    sequence: str,
    statistic: str,
    *,
    n_shuffles: int = DEFAULT_SHUFFLES,
    seed: int = 20261012,
) -> NullResult:
    """A statistic's z-score against shuffles of the same window.

    The null holds the window's length and its letter counts fixed and destroys
    only the order. Two things follow, and together they are the reason this is
    the headline statistic rather than the raw measure:

    * It is free of the small-sample entropy bias, because the bias depends on
      the length and the counts and both are identical in the null.
    * It answers the only interesting form of the question. These sites are
      known to be arginine-rich, so a raw difference in any composition-driven
      measure is already explained. A non-zero z says the *ordering* carries
      something extra.
    """
    if statistic in ORDER_INVARIANT_STATISTICS:
        raise SequenceError(
            f"{statistic!r} is a function of the window's letters alone, so it is "
            "constant under this null and its z is undefined; measure it against a "
            "null that moves what it reads (docs/COMPOSITION_PLAN.md)")
    function = STATISTICS.get(statistic)
    if function is None:
        raise SequenceError(f"unknown statistic {statistic!r}; known: {sorted(STATISTICS)}")
    text = str(sequence).upper()
    if len(text) < 4:
        raise SequenceError("too short for a shuffle null")
    observed = float(function(text))
    letters = np.frombuffer(text.encode("ascii", "ignore"), dtype=np.uint8).copy()
    rng = np.random.default_rng(seed)
    draws = np.empty(n_shuffles)
    for index in range(n_shuffles):
        rng.shuffle(letters)
        draws[index] = function(letters.tobytes().decode("ascii"))
    mean = float(draws.mean())
    sd = float(draws.std(ddof=1))
    z = 0.0 if sd == 0.0 else float((observed - mean) / sd)
    # A permutation p-value, with the usual +1 so it can never be zero.
    extreme = int(np.sum(np.abs(draws - mean) >= abs(observed - mean)))
    p = float((extreme + 1) / (n_shuffles + 1))
    return NullResult(
        statistic=statistic, observed=observed, null_mean=mean, null_sd=sd,
        z=z, p_two_sided=p, n_shuffles=n_shuffles,
    )


# ------------------------------------------------------------- the bundle

@dataclass
class SequenceFeatures:
    """Every descriptor of one window, in every alphabet asked for."""

    sequence: str
    length: int
    raw: Dict[str, float] = field(default_factory=dict)
    nulls: Dict[str, NullResult] = field(default_factory=dict)
    per_alphabet: Dict[str, Dict[str, float]] = field(default_factory=dict)
    saturated_blocks: List[int] = field(default_factory=list)

    def ordering_signals(self) -> Dict[str, float]:
        """The statistics a composition-preserving shuffle could not reproduce."""
        return {
            name: result.z for name, result in self.nulls.items()
            if not result.is_composition_explained
        }

    def as_dict(self) -> Dict[str, object]:
        return {
            "length": self.length,
            "raw": {k: float(v) for k, v in self.raw.items()},
            "nulls": {k: v.as_dict() for k, v in self.nulls.items()},
            "per_alphabet": {
                a: {k: float(v) for k, v in d.items()} for a, d in self.per_alphabet.items()
            },
            "saturated_blocks": list(self.saturated_blocks),
            "ordering_signals": self.ordering_signals(),
        }


def describe_sequence(
    sequence: str,
    *,
    alphabets: Sequence[str] = ("identity", "charge", "hydrophobicity"),
    statistics: Optional[Sequence[str]] = None,
    n_shuffles: int = 500,
    seed: int = 20261012,
) -> SequenceFeatures:
    """Describe one window: raw measures, composition nulls, reduced alphabets."""
    text = str(sequence).upper()
    cleaned = recode(text, "identity")
    if len(cleaned) < 4:
        raise SequenceError(
            f"a {len(cleaned)}-residue window is too short to describe; "
            "the estimators need at least 4"
        )
    out = SequenceFeatures(sequence=cleaned, length=len(cleaned))

    out.raw["entropy"] = shannon_entropy(cleaned)
    out.raw["entropy_plugin"] = shannon_entropy(cleaned, estimator="plugin")
    out.raw["entropy_bias_visible"] = out.raw["entropy"] - out.raw["entropy_plugin"]
    out.raw["lz_complexity"] = float(lempel_ziv_complexity(cleaned))
    out.raw["lz_normalised"] = lz_normalised(cleaned)
    out.raw["compressed_bits_zlib"] = compressed_bits(cleaned, "zlib")
    for k in (2, 3):
        if len(cleaned) > k:
            out.raw[f"entropy_rate_k{k}"] = entropy_rate(cleaned, k)
            if block_entropy_is_saturated(cleaned, k):
                out.saturated_blocks.append(k)
    if len(cleaned) >= 6:
        out.raw["excess_entropy"] = excess_entropy(cleaned, 3)
    out.raw["helix_charge_periodicity"] = helix_charge_periodicity(cleaned)
    out.raw["strand_charge_periodicity"] = strand_charge_periodicity(cleaned)
    dominant = dominant_charge_period(cleaned)
    out.raw["dominant_charge_period"] = dominant["period"]
    out.raw["dominant_period_power"] = dominant["power"]
    # A band power only means "helical face" when the band also holds the
    # dominant period; otherwise it is a harmonic of something else.
    out.raw["helix_face"] = float(
        dominant["in_helix_band"] and dominant["power"] > 0.3
    )
    out.raw["strand_face"] = float(
        dominant["in_strand_band"] and dominant["power"] > 0.3
    )

    for alphabet in alphabets:
        recoded = recode(cleaned, alphabet)
        if len(recoded) < 4:
            continue
        entry = {
            "entropy": shannon_entropy(recoded),
            "lz_complexity": float(lempel_ziv_complexity(recoded)),
            "lz_normalised": lz_normalised(recoded),
            "alphabet_size": float(len(set(recoded))),
        }
        if len(recoded) > 2:
            entry["entropy_rate_k2"] = entropy_rate(recoded, 2)
        out.per_alphabet[alphabet] = entry

    for name in (statistics if statistics is not None else sorted(STATISTICS)):
        try:
            out.nulls[name] = composition_z(
                cleaned, name, n_shuffles=n_shuffles, seed=seed)
        except SequenceError:
            continue
    return out


# ------------------------------------------------------- association testing

@dataclass
class FeatureComparison:
    """One feature compared between annotated sites and matched controls."""

    feature: str
    n_site_groups: int
    n_control_groups: int
    site_mean: float
    control_mean: float
    difference: float
    low: float
    high: float
    p_value: float
    minimum_detectable: float

    @property
    def is_evidence(self) -> bool:
        """Whether the comparison had the groups to mean anything.

        The project's standing floor: fewer than `MIN_GROUPS` homology groups
        on either side is not evidence whichever way it points.
        """
        return min(self.n_site_groups, self.n_control_groups) >= MIN_GROUPS

    @property
    def exceeds_power(self) -> bool:
        """Whether the observed difference is one this design could detect."""
        return abs(self.difference) >= self.minimum_detectable

    def as_dict(self) -> Dict[str, object]:
        return {
            "feature": self.feature,
            "n_site_groups": self.n_site_groups,
            "n_control_groups": self.n_control_groups,
            "site_mean": float(self.site_mean),
            "control_mean": float(self.control_mean),
            "difference": float(self.difference),
            "interval": [float(self.low), float(self.high)],
            "p_value": float(self.p_value),
            "minimum_detectable": float(self.minimum_detectable),
            "is_evidence": self.is_evidence,
            "exceeds_power": self.exceeds_power,
        }


def minimum_detectable_difference(
    site_values: Sequence[float],
    control_values: Sequence[float],
    *,
    power: float = 0.8,
) -> float:
    """The smallest difference in group means this design could detect.

    Computed from the group counts and the observed spread, and meant to be
    read *before* the difference itself. The project has been here before: a
    held-out AUC of 0.499 was taken for a modelling failure when at five
    homology groups nothing below 0.87 was detectable at all.
    """
    site = np.asarray(site_values, dtype=float)
    control = np.asarray(control_values, dtype=float)
    if site.size < 2 or control.size < 2:
        return float("inf")
    if power != 0.8:
        raise SequenceError("only power=0.8 is implemented")
    z_alpha, z_beta = 1.959963985, 0.8416212336
    pooled = float(np.sqrt(site.var(ddof=1) / site.size + control.var(ddof=1) / control.size))
    return float((z_alpha + z_beta) * pooled)


def compare_feature(
    feature: str,
    site_values: Mapping[object, float],
    control_values: Mapping[object, float],
    *,
    n_bootstrap: int = 2000,
    seed: int = 20261012,
) -> FeatureComparison:
    """Compare one feature between sites and controls, resampling whole groups.

    ``site_values`` and ``control_values`` map a homology group to that group's
    mean value. Resampling groups rather than windows is the project's standing
    rule and the reason it is not possible to report a confident difference
    from 134 copies of five proteins.
    """
    sites = list(site_values.values())
    controls = list(control_values.values())
    if not sites or not controls:
        raise SequenceError("both arms need at least one group")
    site_keys = list(site_values)
    control_keys = list(control_values)
    rng = np.random.default_rng(seed)
    observed = float(np.mean(sites) - np.mean(controls))
    draws = np.empty(n_bootstrap)
    for index in range(n_bootstrap):
        a = [site_values[site_keys[i]] for i in rng.integers(0, len(site_keys), len(site_keys))]
        b = [control_values[control_keys[i]]
             for i in rng.integers(0, len(control_keys), len(control_keys))]
        draws[index] = float(np.mean(a) - np.mean(b))
    low, high = (float(v) for v in np.percentile(draws, [2.5, 97.5]))
    # A two-sided permutation p-value over group labels.
    pooled = sites + controls
    labels = np.array([1] * len(sites) + [0] * len(controls))
    extreme = 0
    for _ in range(n_bootstrap):
        rng.shuffle(labels)
        a = [v for v, lab in zip(pooled, labels) if lab == 1]
        b = [v for v, lab in zip(pooled, labels) if lab == 0]
        if abs(float(np.mean(a) - np.mean(b))) >= abs(observed):
            extreme += 1
    return FeatureComparison(
        feature=feature,
        n_site_groups=len(sites), n_control_groups=len(controls),
        site_mean=float(np.mean(sites)), control_mean=float(np.mean(controls)),
        difference=observed, low=low, high=high,
        p_value=float((extreme + 1) / (n_bootstrap + 1)),
        minimum_detectable=minimum_detectable_difference(sites, controls),
    )


def holm(p_values: Mapping[str, float]) -> Dict[str, float]:
    """Holm-corrected p-values, for the features tried in parallel.

    Trying every feature in every alphabet and reporting the best is how a
    signature gets manufactured. The correction is across all of them.
    """
    items = sorted(p_values.items(), key=lambda kv: kv[1])
    n = len(items)
    out: Dict[str, float] = {}
    running = 0.0
    for index, (name, p) in enumerate(items):
        adjusted = min(1.0, (n - index) * float(p))
        running = max(running, adjusted)
        out[name] = running
    return out
