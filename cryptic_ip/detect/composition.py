"""The composition half of the sequence decomposition (study Q10).

A window of sequence factorises as **composition** (which letters) **⊕ order**
(in what arrangement). `sequence.py` holds the first fixed and asks about the
second: its `composition_z` permutes a window's own letters, so composition is
invariant under its null and invisible to every statistic it carries.

This module is the other half. Every statistic here is a function of the
**multiset** of letters alone and is therefore invariant to their order, which
has two consequences fixed in `docs/COMPOSITION_PLAN.md`:

* These statistics are sound on the **gap-closed** `lining` window, which
  `signature.POSITIONAL_STATISTICS` exists to keep periodicity statistics off.
  Closing the gaps reorders nothing that any of them reads.
* They may **never** be passed to `composition_z`. An order-invariant statistic
  has zero variance under an order-permuting null, so its z is 0/0 - and a
  naive implementation returns a silent zero that reads as "no ordering
  signal" when the truth is that the question is undefined. `sequence.py`
  refuses them by name, using the registry below.

The null these are measured against is not here but in the study driver: random
same-size residue subsets of the same chain, which holds the chain, the fold and
the lining's size fixed and randomises only which residues.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, Mapping, Sequence

import numpy as np

from .sequence import CHARGE_OF, HYDROPATHY, SequenceError

#: The twenty standard one-letter codes, in a fixed order so that a
#: distribution over them is comparable between windows.
ALPHABET = tuple("ACDEFGHIKLMNPQRSTVWY")

#: The residues this study calls basic. Histidine is left out on purpose: it is
#: only fractionally protonated at pH 7.4, `CHARGE_OF` gives it 0.1 for exactly
#: that reason, and counting it whole would put a judgement call inside the
#: positive control's statistic.
BASIC = frozenset("RK")

#: And the ones it calls acidic and aromatic.
ACIDIC = frozenset("DE")
AROMATIC = frozenset("FWY")

# No smoothing constant appears in this module, and that is deliberate. A
# twelve-residue lining cannot populate twenty letters, so a one-sided
# Kullback-Leibler divergence against a chain background would be infinite and
# would need a pseudocount to tame - a free parameter inside the statistic. The
# Jensen-Shannon divergence needs none: where one side is zero the mid-point is
# still positive, so every term is finite and the whole quantity is bounded by
# one bit. Dropping the pseudocount also restores the identity a divergence
# must have - a distribution compared with itself is exactly zero - which an
# inconsistently smoothed pair does not satisfy.


def letter_counts(sequence: str, alphabet: Sequence[str] = ALPHABET) -> np.ndarray:
    """Counts of each letter of `alphabet`, in its order."""
    text = str(sequence).upper()
    if not text:
        raise SequenceError("an empty window has no composition")
    index = {letter: position for position, letter in enumerate(alphabet)}
    counts = np.zeros(len(alphabet), dtype=float)
    for letter in text:
        position = index.get(letter)
        if position is not None:
            counts[position] += 1.0
    if counts.sum() == 0.0:
        raise SequenceError("no standard residues in the window")
    return counts


def distribution(sequence: str, alphabet: Sequence[str] = ALPHABET) -> np.ndarray:
    """The window's letter distribution, unsmoothed."""
    counts = letter_counts(sequence, alphabet)
    return counts / counts.sum()


def jensen_shannon(left: np.ndarray, right: np.ndarray) -> float:
    """Jensen-Shannon divergence in bits, which is bounded by 1 and symmetric.

    Symmetry matters: a lining is not privileged over the chain it sits in, and
    a one-sided KL divergence would blow up on whichever side happens to have a
    zero. Boundedness matters because the group bootstrap averages these.
    """
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape:
        raise SequenceError(f"distributions differ in length: {left.shape} vs {right.shape}")
    mid = 0.5 * (left + right)
    out = 0.0
    for part in (left, right):
        mask = part > 0.0
        out += 0.5 * float(np.sum(part[mask] * np.log2(part[mask] / mid[mask])))
    # Floating point can take a divergence of zero a hair below it.
    return float(min(1.0, max(0.0, out)))


def fraction_of(sequence: str, letters: frozenset) -> float:
    """The share of the window's standard residues drawn from `letters`."""
    counts = letter_counts(sequence)
    total = float(counts.sum())
    wanted = sum(counts[ALPHABET.index(letter)] for letter in letters
                 if letter in ALPHABET)
    return float(wanted / total)


def basic_fraction(sequence: str) -> float:
    """The fraction of arginine and lysine - the positive control's statistic."""
    return fraction_of(sequence, BASIC)


def acidic_fraction(sequence: str) -> float:
    return fraction_of(sequence, ACIDIC)


def aromatic_fraction(sequence: str) -> float:
    return fraction_of(sequence, AROMATIC)


def glycine_fraction(sequence: str) -> float:
    return fraction_of(sequence, frozenset("G"))


def proline_fraction(sequence: str) -> float:
    return fraction_of(sequence, frozenset("P"))


def mean_hydropathy(sequence: str) -> float:
    """Mean Kyte-Doolittle hydropathy, on the scale already in `sequence.py`."""
    counts = letter_counts(sequence)
    weights = np.array([HYDROPATHY[letter] for letter in ALPHABET], dtype=float)
    return float(np.dot(counts, weights) / counts.sum())


def net_charge_per_residue(sequence: str) -> float:
    """Mean formal charge at pH 7.4, on `sequence.CHARGE_OF`."""
    counts = letter_counts(sequence)
    weights = np.array([CHARGE_OF[letter] for letter in ALPHABET], dtype=float)
    return float(np.dot(counts, weights) / counts.sum())


def composition_divergence(sequence: str, background: Mapping[str, float]) -> float:
    """How far the window's letters sit from a background distribution, in bits.

    The omnibus statistic: it fires on any compositional difference in any
    direction, so it does not need a direction named in advance.
    """
    reference = np.array([float(background.get(letter, 0.0)) for letter in ALPHABET])
    if reference.sum() <= 0.0:
        raise SequenceError("the background distribution is empty")
    return jensen_shannon(distribution(sequence), reference / reference.sum())


def nonbasic_divergence(sequence: str, background: Mapping[str, float]) -> float:
    """The same divergence with arginine and lysine removed and the rest renormalised.

    The one primary statistic of study Q10 whose answer is not already known.
    That an inositol-phosphate shell is basic is textbook; whether anything
    *else* about its composition is unusual for its own chain is not.
    """
    keep = tuple(letter for letter in ALPHABET if letter not in BASIC)
    counts = letter_counts(sequence, keep)
    reference = np.array([float(background.get(letter, 0.0)) for letter in keep])
    if reference.sum() <= 0.0:
        raise SequenceError("the background distribution has no non-basic residues")
    return jensen_shannon(counts / counts.sum(), reference / reference.sum())


def background_of(sequence: str) -> Dict[str, float]:
    """A chain's own letter distribution, to be the background for its linings."""
    counts = letter_counts(sequence)
    total = float(counts.sum())
    return {letter: float(counts[position] / total)
            for position, letter in enumerate(ALPHABET)}


#: Statistics of a window alone, needing no background.
SIMPLE_STATISTICS: Dict[str, Callable[[str], float]] = {
    "basic_fraction": basic_fraction,
    "acidic_fraction": acidic_fraction,
    "aromatic_fraction": aromatic_fraction,
    "glycine_fraction": glycine_fraction,
    "proline_fraction": proline_fraction,
    "mean_hydropathy": mean_hydropathy,
    "net_charge_per_residue": net_charge_per_residue,
}

#: Statistics that compare a window with the chain it came from.
BACKGROUND_STATISTICS: Dict[str, Callable[[str, Mapping[str, float]], float]] = {
    "composition_divergence": composition_divergence,
    "nonbasic_divergence": nonbasic_divergence,
}

#: Every statistic this module carries. Each one is invariant to the order of a
#: window's letters, which is what `sequence.composition_z` checks against.
COMPOSITION_STATISTICS = frozenset(SIMPLE_STATISTICS) | frozenset(BACKGROUND_STATISTICS)


def evaluate(statistic: str, sequence: str,
             background: Mapping[str, float]) -> float:
    """One statistic of one window, whichever family it belongs to."""
    if statistic in SIMPLE_STATISTICS:
        return float(SIMPLE_STATISTICS[statistic](sequence))
    if statistic in BACKGROUND_STATISTICS:
        return float(BACKGROUND_STATISTICS[statistic](sequence, background))
    raise SequenceError(
        f"unknown composition statistic {statistic!r}; "
        f"known: {sorted(COMPOSITION_STATISTICS)}")


def is_order_invariant(statistic: str, sequence: str, *, draws: int = 32,
                       seed: int = 20261025,
                       background: Mapping[str, float] = None) -> bool:
    """Check empirically that a statistic ignores the order of its letters.

    Used by the tests rather than by the pipeline: the registry above is the
    claim, and this is how the claim is checked, so a statistic that reads
    order cannot be added to that registry and go unnoticed.
    """
    letters = list(str(sequence).upper())
    if background is None:
        background = background_of(sequence)
    observed = evaluate(statistic, "".join(letters), background)
    rng = np.random.default_rng(seed)
    for _ in range(draws):
        rng.shuffle(letters)
        if not math.isclose(evaluate(statistic, "".join(letters), background),
                            observed, rel_tol=1e-12, abs_tol=1e-12):
            return False
    return True
