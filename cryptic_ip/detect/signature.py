"""From a pocket record to a sequence window, with the plan's rules in code.

`sequence.py` measures a string. This module decides *which* string, which is
where the mistakes live. A pocket's lining is a set of residue numbers, and
there are two very different sequences you can build from them:

* the lining residues concatenated in sequence order, gaps closed, and
* the contiguous stretch they span.

They answer different questions and only the second one can carry a period.
Closing the gaps destroys exactly the positional information a periodicity
measures, so a helical-face score computed on a gap-closed lining window is
meaningless - and it is the easiest mistake to make with this data, because
both windows are strings and every function in `sequence.py` will happily
accept either. `docs/SEQUENCE_PLAN.md` writes the rule down;
`statistic_for_window` refuses it.

The second thing here is a refusal. `PRIMARY_TESTS` is the pre-registered test
family, transcribed from the plan so a test can check the two agree, and
`ESTABLISHED_SIGNATURES` is empty. Until study Q6 runs, the detector gets no
sequence term at all: `sequence_contribution` abstains and says why. A feature
that has not been shown to carry signal does not get to move a score, however
good the idea is - this project has a held-out AUC of 0.499 on record from the
last time an untested feature was trusted.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .sequence import (
    AMINO_ACIDS,
    SequenceError,
    composition_z,
    describe_sequence,
    spacing_spectrum,
)

#: The three window definitions, fixed by the plan.
WINDOWS = ("lining", "local", "spacing")

#: The longest contiguous stretch a local window may span. Past this a "period"
#: is being fitted to a whole domain rather than to a pocket's lining.
LOCAL_CAP = 40

#: The shortest window any descriptor is computed on. The estimators need 4;
#: the floor is 8 so a k=3 block entropy is not read off 6 points.
MIN_WINDOW = 8

#: A lining with more than this fraction of non-standard residues is dropped
#: rather than measured, because `recode` drops them and the window silently
#: stops being the lining.
MAX_NONSTANDARD_FRACTION = 0.10

#: Statistics that depend on where a residue sits, and so may not be computed
#: on a gap-closed window. This is the rule the plan states, enforced.
POSITIONAL_STATISTICS = frozenset({
    "helix_charge_periodicity",
    "strand_charge_periodicity",
    "dominant_period_power",
})

#: Statistics that are functions of the multiset of letters and their order
#: *within the window given*, and so are meaningful on either window. Spelled
#: out rather than taken as "everything else" so that a statistic added to
#: `sequence.STATISTICS` lands in neither set and fails
#: `tests/test_detect_signature.py` until someone classifies it. A new
#: positional statistic that silently defaulted to "safe on the lining window"
#: is exactly the failure this fence exists to prevent.
NON_POSITIONAL_STATISTICS = frozenset({
    "entropy",
    "entropy_rate_k2",
    "excess_entropy",
    "lz_complexity",
    "lz_normalised",
})


#: The share of a pocket's lining that must sit on one chain before a single
#: sequence window can stand for it. A pocket at an interface has no such
#: window - its lining is two stretches of two different sequences - and
#: concatenating them would invent a neighbourhood that does not exist.
MIN_CHAIN_FRACTION = 0.8


class WindowError(SequenceError):
    """A pocket from which no honest sequence window can be built."""


def dominant_chain(residues: Sequence[Tuple[str, int]]) -> Tuple[str, float, List[int]]:
    """Which chain holds the lining, what share of it, and those positions.

    `PocketAnalyzer.get_pocket_residues` returns residue *numbers* and says in
    its own docstring that they are chain-blind. On a multi-chain structure the
    same number appears on every chain, so a window built from the numbers
    alone silently mixes chains. Lining residues therefore arrive here as
    ``(chain_id, residue_number)`` pairs, and the share on the winning chain is
    returned so an interface pocket can be excluded rather than measured.
    """
    if not residues:
        raise WindowError("no lining residues")
    counts: Dict[str, List[int]] = {}
    for chain_id, number in residues:
        counts.setdefault(str(chain_id), []).append(int(number))
    # max() keeps the first maximal item, and the items are sorted, so a tie
    # breaks to the alphabetically first chain rather than to dict order.
    chain_id, positions = max(sorted(counts.items()), key=lambda kv: len(kv[1]))
    return chain_id, len(positions) / float(len(residues)), sorted(set(positions))


def spans_chains(residues: Sequence[Tuple[str, int]]) -> bool:
    """Whether the lining is too evenly split between chains to window at all."""
    _chain_id, fraction, _positions = dominant_chain(residues)
    return fraction < MIN_CHAIN_FRACTION


@dataclass(frozen=True)
class SiteWindows:
    """The three windows of one pocket, with what it cost to build them."""

    lining: str
    local: str
    positions: Tuple[int, ...]
    local_start: int
    local_end: int
    n_nonstandard: int
    truncated: bool

    @property
    def nonstandard_fraction(self) -> float:
        total = len(self.lining) + self.n_nonstandard
        return 0.0 if total == 0 else self.n_nonstandard / float(total)

    def as_dict(self) -> Dict[str, object]:
        return {
            "lining": self.lining,
            "local": self.local,
            "positions": list(self.positions),
            "local_span": [self.local_start, self.local_end],
            "lining_length": len(self.lining),
            "local_length": len(self.local),
            "n_nonstandard": self.n_nonstandard,
            "nonstandard_fraction": self.nonstandard_fraction,
            "truncated": self.truncated,
        }


def site_windows(
    sequence: str,
    positions: Sequence[int],
    *,
    offset: int = 1,
    local_cap: int = LOCAL_CAP,
) -> SiteWindows:
    """Build the lining and local windows of one pocket.

    ``sequence`` is the chain's sequence and ``offset`` the residue number of
    its first residue - 1 for a deposited chain that starts at residue 1, and
    the real first number otherwise. Getting this wrong shifts every window by
    a constant and would be invisible in the output, so an out-of-range
    position raises rather than being clamped.

    A local window longer than ``local_cap`` is truncated around the median
    lining position and the record says so, because a 300-residue "local"
    window is a domain and a period fitted to it means nothing.
    """
    chain = str(sequence).upper()
    if not chain:
        raise WindowError("no chain sequence")
    wanted = sorted({int(p) for p in positions})
    if len(wanted) < 2:
        raise WindowError(f"a pocket needs at least 2 lining positions, got {len(wanted)}")
    low, high = offset, offset + len(chain) - 1
    outside = [p for p in wanted if not low <= p <= high]
    if outside:
        raise WindowError(
            f"lining positions {outside[:4]} fall outside the chain's residues "
            f"{low}-{high}; the numbering offset is probably wrong"
        )

    letters = [chain[p - offset] for p in wanted]
    lining = "".join(c for c in letters if c in AMINO_ACIDS)
    n_nonstandard = len(letters) - len(lining)

    start, end = wanted[0], wanted[-1]
    truncated = False
    if end - start + 1 > local_cap:
        centre = wanted[len(wanted) // 2]
        half = local_cap // 2
        start = max(low, centre - half)
        end = min(high, start + local_cap - 1)
        start = max(low, end - local_cap + 1)
        truncated = True
    local = "".join(
        c for c in chain[start - offset:end - offset + 1] if c in AMINO_ACIDS
    )
    return SiteWindows(
        lining=lining, local=local, positions=tuple(wanted),
        local_start=start, local_end=end,
        n_nonstandard=n_nonstandard, truncated=truncated,
    )


def is_measurable(windows: SiteWindows) -> Tuple[bool, str]:
    """Whether the plan allows this pocket into the study, and why not if not.

    Returned as a reason rather than a bare False so an excluded pocket can be
    counted by cause in the results, which is the difference between "we
    excluded 40 pockets" and a usable exclusion table.
    """
    # The non-standard fraction is tested first because it is the cause of a
    # short window when both apply: a lining of ten residues half of which are
    # X comes out as a five-residue window, and reporting that as "too short"
    # would put it in the wrong column of the exclusion table.
    if windows.nonstandard_fraction > MAX_NONSTANDARD_FRACTION:
        return False, f"non-standard fraction {windows.nonstandard_fraction:.2f}"
    if len(windows.lining) < MIN_WINDOW:
        return False, f"lining window {len(windows.lining)} < {MIN_WINDOW}"
    if len(windows.local) < MIN_WINDOW:
        return False, f"local window {len(windows.local)} < {MIN_WINDOW}"
    return True, "measurable"


def statistic_for_window(windows: SiteWindows, window: str, statistic: str,
                         *, n_shuffles: int = 500, seed: int = 20261012):
    """One statistic's composition z on one named window, with the rule enforced.

    Refuses a positional statistic on the ``lining`` window. The gap-closed
    lining window has no positions left, so a period fitted to it is a period
    of the concatenation and not of the protein.
    """
    if window not in ("lining", "local"):
        raise WindowError(f"unknown window {window!r}; use 'lining' or 'local'")
    if window == "lining" and statistic in POSITIONAL_STATISTICS:
        raise WindowError(
            f"{statistic!r} is positional and the lining window is gap-closed; "
            "compute it on the 'local' window (docs/SEQUENCE_PLAN.md)"
        )
    text = windows.lining if window == "lining" else windows.local
    return composition_z(text, statistic, n_shuffles=n_shuffles, seed=seed)


def describe_site(
    windows: SiteWindows,
    *,
    n_shuffles: int = 500,
    seed: int = 20261012,
) -> Dict[str, object]:
    """Every descriptor of one pocket, by window, with the lining window fenced.

    The lining bundle carries no periodicity at all - not a zero, which would
    be read as "no helical face", but nothing, because the question cannot be
    asked of that window.
    """
    ok, reason = is_measurable(windows)
    out: Dict[str, object] = {
        "windows": windows.as_dict(),
        "measurable": ok,
        "exclusion_reason": None if ok else reason,
    }
    if not ok:
        return out
    local = describe_sequence(windows.local, n_shuffles=n_shuffles, seed=seed)
    lining = describe_sequence(
        windows.lining,
        n_shuffles=n_shuffles,
        seed=seed,
        statistics=[s for s in ("entropy", "entropy_rate_k2", "excess_entropy",
                                "lz_complexity", "lz_normalised")],
    )
    lining_raw = {
        k: v for k, v in lining.as_dict()["raw"].items()
        if not any(k.startswith(prefix) for prefix in
                   ("helix", "strand", "dominant"))
    }
    out["local"] = local.as_dict()
    out["lining"] = {**lining.as_dict(), "raw": lining_raw}
    out["spacing"] = spacing_spectrum(windows.positions)
    return out


# ------------------------------------------------- the pre-registered family

@dataclass(frozen=True)
class PrimaryTest:
    """One of the six tests `docs/SEQUENCE_PLAN.md` fixes in advance."""

    number: int
    window: str
    statistic: str
    hypothesis: str


#: Transcribed from the plan's table. `tests/test_detect_signature.py` parses
#: the document and checks the two have not drifted apart, which is the only
#: way a pre-registration stays honest once there is code beside it.
PRIMARY_TESTS: Tuple[PrimaryTest, ...] = (
    PrimaryTest(1, "local", "entropy_rate_k2",
                "the ordering is more predictable at a site than its letters imply"),
    PrimaryTest(2, "local", "helix_charge_periodicity",
                "a basic helix face lines the pocket"),
    PrimaryTest(3, "local", "strand_charge_periodicity",
                "a basic strand face lines it instead"),
    PrimaryTest(4, "local", "dominant_period_power",
                "some charge period dominates, whichever it is"),
    PrimaryTest(5, "lining", "lz_normalised",
                "the lining letters are more repetitive than chance"),
    PrimaryTest(6, "spacing", "fraction_helical",
                "the coordinating residues are 3-4 apart"),
)


def primary_statistics() -> List[Tuple[str, str]]:
    """The (window, statistic) pairs of the primary family, in plan order."""
    return [(t.window, t.statistic) for t in PRIMARY_TESTS]


# --------------------------------------------------------------- the refusal

#: The ordering-signal z worth full credit in a detector score. 4 is chosen,
#: not fitted: `composition_z` calls anything under 2 composition-explained, so
#: a feature has to clear twice the threshold it is judged by before it is
#: worth as much as a whole physical term. It is not tuned against any outcome.
Z_FOR_FULL_CREDIT = 4.0

#: Signatures study Q6 has established, as ``"<window>:<statistic>"``. Empty,
#: and it stays empty until the study reports: a name may be added here only
#: alongside the results file that cleared all five of the plan's decision
#: rules.
ESTABLISHED_SIGNATURES: Tuple[str, ...] = ()


def sequence_contribution(
    site: Dict[str, object],
    *,
    established: Sequence[str] = ESTABLISHED_SIGNATURES,
) -> Tuple[Optional[float], str]:
    """What the sequence descriptors contribute to a detector score.

    ``(None, reason)`` while nothing is established, which is the current state
    and the point of the function: the descriptors are computed, recorded and
    available, and they move no score. `monotone.py` treats a ``None`` as an
    abstention rather than a zero, so an absent sequence term cannot quietly
    push a candidate down.

    Once a signature clears the plan, it goes in `ESTABLISHED_SIGNATURES` as
    ``"<window>:<statistic>"`` - the window is part of the name because the
    plan pins each test to one window and the same statistic on the other
    window is a different test that was not pre-registered. This then returns
    the mean ordering-signal z over the established features, clipped to
    [-1, 1] after division by `Z_FOR_FULL_CREDIT`.
    """
    if not established:
        return None, (
            "no sequence signature has been established; study Q6 "
            "(docs/SEQUENCE_PLAN.md) has not reported"
        )
    signals: List[float] = []
    for name in established:
        if ":" not in name:
            raise WindowError(
                f"an established signature is named '<window>:<statistic>', got {name!r}"
            )
        window, statistic = name.split(":", 1)
        bundle = site.get(window)
        if not isinstance(bundle, dict):
            continue
        nulls = bundle.get("nulls")
        if isinstance(nulls, dict) and statistic in nulls:
            signals.append(float(nulls[statistic]["z"]))
    if not signals:
        return None, f"none of {list(established)} was measurable on this pocket"
    mean = sum(signals) / len(signals)
    scaled = mean / float(Z_FOR_FULL_CREDIT)
    return max(-1.0, min(1.0, scaled)), f"{len(signals)} established signal(s)"
