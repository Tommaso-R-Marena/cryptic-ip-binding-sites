#!/usr/bin/env python3
"""Study Q9: does the sequence predict pocket formation, type, or function?

Pre-registered in `docs/POCKET_SEQUENCE_PLAN.md`, dated before any sequence was
read from any structure for this study.

A window of sequence factorises the way a lining of positions does -
**composition** (which letters) **⊕ order** (in what arrangement) - and the
canonical null for order is permutation of the window's own letters. That is
the same quotient `architecture.py` takes on positions, and it is why every
primary test here is on `composition_z` rather than on a raw measure: a raw
difference is already explained by composition, since these sites are basic
and pocket linings are hydrophobic.

    pocket_sequence.py sequences --structures DIR --out sequences.csv
    pocket_sequence.py analyse --sequences sequences.csv --pockets pockets.csv ...

The `sequences` step reads structures and runs no pocket detection, so it is
minutes rather than hours; `analyse` is a pure function of committed tables.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.detect import sequence as sq  # noqa: E402

SEED = 20261020

#: Shuffles per window for the order null. The plan's number.
SHUFFLES = 2000

#: Longest `local` window, from `docs/SEQUENCE_PLAN.md`: the contiguous stretch
#: a lining spans, capped so one sprawling pocket cannot dominate.
MAX_WINDOW = 40

#: Shortest window worth a null at all, also from that plan.
MIN_WINDOW = 8

#: fpocket's own rank cut for the formation contrast. Geometric, reads no
#: sequence, fixed in `docs/POCKET_SEQUENCE_PLAN.md` before any data.
TOP_RANK = 5

#: Window lengths must match within this many residues, or the comparison is
#: reported as length-confounded.
LENGTH_TOLERANCE = 2

#: Amendment 1 section 4: the `dense` window's fixed length. The same cap as
#: `local`, used as a width rather than as a refusal, so every window of this
#: kind is exactly this long and the length confound is zero by construction.
DENSE_WINDOW = MAX_WINDOW

#: Amendment 1 section 3: at most this many matched pairs from one protein, so
#: a model with two hundred pockets cannot carry an arm on its own.
MAX_PAIRS_PER_PROTEIN = 12

#: And at most this many from an arm that is not split by protein.
MAX_PAIRS_PER_ARM = 400

#: Permutation draws per architecture null in the `classes` step. The
#: architecture module's own default, fixed here so the class table is
#: reproducible from the committed command line alone.
CLASS_DRAWS = 800

#: The two window definitions. `local` is the plan's; `dense` is amendment 1's
#: and is exploratory on every arm.
WINDOW_KINDS = ("local", "dense")

#: What the pre-registered arm reported in the one run amendment 1 section 1
#: discloses, so that this run can say whether it is looking at the same data.
#: If these do not reproduce, the pocket table has changed and the held-out
#: arm's claim to be held out from anything has to be re-established.
RECORDED_LOOK = {"arm": "formation", "statistic": "helix_hydropathy_periodicity",
                 "n_site": 10, "n_groups": 9, "difference": 1.317}

#: How many Monte-Carlo standard errors a difference may move and still count as
#: the same number. Amendment 2 section 3: a re-seeded null is a fresh draw on
#: the same windows, so the right unit is the spread the shuffle count implies
#: and not a fixed constant - which is what the first version of this check used,
#: and why it reported a gap of 0.41 standard errors as a failure to reproduce.
LOOK_SPREAD_BOUND = 3.0

#: Arms that may become a finding. The plan's two, plus the one confirmatory
#: arm amendment 1 adds; everything else is reported and never promoted.
CONFIRMATORY_ARMS = ("formation", "function", "type", "formation_heldout")

#: The primary family, fixed by the plan: two statistics on two arms.
PRIMARY = ("helix_charge_periodicity", "helix_hydropathy_periodicity")

#: Reported in full, never promoted to a finding.
SECONDARY = ("strand_charge_periodicity", "strand_hydropathy_periodicity",
             "dominant_period_power", "dominant_hydropathy_period_power",
             "entropy_rate_k2", "lz_normalised")

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


# ----------------------------------------------------------- the sequences

def chain_sequences(path: Path) -> Dict[str, Dict[int, str]]:
    """Every polymer chain's residue number to one-letter code, from a structure.

    Read straight from the file the pockets were detected in, so a residue
    number and its letter cannot disagree - which they can, and silently, when
    a sequence comes from a database entry and the numbering comes from a
    model.
    """
    out: Dict[str, Dict[int, str]] = {}
    seen: Dict[Tuple[str, int, str], bool] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            resname = line[17:20].strip().upper()
            letter = THREE_TO_ONE.get(resname)
            if letter is None:
                continue
            chain = line[21:22].strip() or " "
            try:
                number = int(line[22:26].strip())
            except ValueError:
                continue
            icode = line[26:27].strip()
            key = (chain, number, icode)
            if key in seen:
                continue
            seen[key] = True
            # An insertion code shares a residue number with its parent; the
            # first one wins, because the window is indexed by number.
            out.setdefault(chain, {}).setdefault(number, letter)
    return out


def cmd_sequences(args: argparse.Namespace) -> int:
    rows = []
    paths = []
    for directory in args.structures:
        if not directory.is_dir():
            print(f"no such directory: {directory}", flush=True)
            continue
        paths += sorted(path for path in directory.iterdir()
                        if path.suffix.lower() in {".pdb", ".ent"})
    for index, path in enumerate(paths, start=1):
        try:
            chains = chain_sequences(path)
        except Exception as error:  # noqa: BLE001 - one bad file is not the run
            print(f"[{index}/{len(paths)}] {path.name}: {type(error).__name__}", flush=True)
            continue
        for chain, residues in chains.items():
            if not residues:
                continue
            ordered = sorted(residues)
            rows.append({
                "structure": path.name,
                "stem": path.stem,
                "chain": chain,
                "first": ordered[0],
                "last": ordered[-1],
                "n_residues": len(ordered),
                "numbers": ",".join(str(n) for n in ordered),
                "letters": "".join(residues[n] for n in ordered),
            })
        print(f"[{index}/{len(paths)}] {path.name}: {len(chains)} chains", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(json.dumps({"structures": len(paths), "chains": len(rows)}, indent=2))
    return 0


# -------------------------------------------------------------- the windows

def local_window(numbers: Dict[int, str], positions: Sequence[int]) -> Optional[str]:
    """The contiguous stretch a lining spans, capped, as letters.

    Returns None when the stretch is too short, too broken, or runs past the
    cap - the plan's `local` window, which is the only one a periodicity
    statistic may use, because closing the gaps would destroy what a period
    means.
    """
    wanted = sorted({int(p) for p in positions})
    if not wanted:
        return None
    first, last = wanted[0], wanted[-1]
    if last - first + 1 > MAX_WINDOW:
        return None
    letters = [numbers.get(n) for n in range(first, last + 1)]
    if any(letter is None for letter in letters):
        return None
    text = "".join(letter for letter in letters if letter)
    if len(text) < MIN_WINDOW or not set(text) <= STANDARD:
        return None
    return text


def dense_starts(numbers: Dict[int, str]) -> Tuple[np.ndarray, Dict[int, str]]:
    """Every start of a complete, standard, `DENSE_WINDOW`-long stretch of a chain.

    Built once per chain, because a 2500-residue model has as many starts as it
    has residues and as many pockets as fpocket cares to find.
    """
    ordered = sorted(numbers)
    texts: Dict[int, str] = {}
    if len(ordered) < DENSE_WINDOW:
        return np.empty(0, dtype=int), texts
    for start in ordered:
        letters = [numbers.get(n) for n in range(start, start + DENSE_WINDOW)]
        if any(letter is None for letter in letters):
            continue
        text = "".join(letter for letter in letters if letter)
        if not set(text) <= STANDARD:
            continue
        texts[start] = text
    return np.asarray(sorted(texts), dtype=int), texts


def dense_window(starts: np.ndarray, texts: Dict[int, str],
                 positions: Sequence[int]) -> Optional[str]:
    """Amendment 1's window: the densest fixed-width stretch the lining sits in.

    The `DENSE_WINDOW`-residue contiguous stretch containing the most of the
    lining's positions, ties going to the lowest start. It reads positions only
    - never the letters it will be scored on - and every window it returns is
    exactly `DENSE_WINDOW` long, so the length confound the plan matches away
    is zero here by construction. Exploratory on every arm; see the amendment.
    """
    if starts.size == 0:
        return None
    wanted = np.asarray(sorted({int(p) for p in positions}), dtype=int)
    if wanted.size == 0:
        return None
    inside = (np.searchsorted(wanted, starts + DENSE_WINDOW, side="left")
              - np.searchsorted(wanted, starts, side="left"))
    best = int(np.argmax(inside))
    if int(inside[best]) <= 0:
        return None
    return texts[int(starts[best])]


def refusal(numbers: Dict[int, str], positions: Sequence[int]) -> Optional[str]:
    """Why `local_window` refused a lining, for the yield diagnostics.

    Amendment 1 section 6. Reporting only: it reproduces the refusal `local_window`
    already made rather than deciding anything.
    """
    wanted = sorted({int(p) for p in positions})
    if not wanted:
        return "no positions"
    first, last = wanted[0], wanted[-1]
    if last - first + 1 > MAX_WINDOW:
        return "span over the cap"
    letters = [numbers.get(n) for n in range(first, last + 1)]
    if any(letter is None for letter in letters):
        return "gap in the model"
    text = "".join(letter for letter in letters if letter)
    if len(text) < MIN_WINDOW:
        return "under the floor"
    if not set(text) <= STANDARD:
        return "non-standard letter"
    return None


def window_z(text: str, statistic: str, *, shuffles: int, seed: int) -> Optional[float]:
    """One window's order-null z, or None when the window cannot carry it."""
    try:
        return float(sq.composition_z(text, statistic, n_shuffles=shuffles, seed=seed).z)
    except sq.SequenceError:
        return None


def match_on_length(left: List[Dict[str, object]], right: List[Dict[str, object]],
                    rng: np.random.Generator) -> Tuple[List[Dict[str, object]],
                                                       List[Dict[str, object]]]:
    """Pair the two arms on window length, to within the tolerance.

    Greedy and seeded: every window of the smaller arm takes the closest
    unused partner from the other. Length is the covariate that would
    otherwise explain a difference in any length-sensitive statistic, and the
    plan fixes the tolerance rather than letting a run choose it.
    """
    short, long = (left, right) if len(left) <= len(right) else (right, left)
    pool = list(long)
    rng.shuffle(pool)
    paired_short, paired_long = [], []
    for item in short:
        best, distance = None, None
        for candidate in pool:
            gap = abs(int(candidate["length"]) - int(item["length"]))
            if distance is None or gap < distance:
                best, distance = candidate, gap
        if best is None or distance is None or distance > LENGTH_TOLERANCE:
            continue
        pool.remove(best)
        paired_short.append(item)
        paired_long.append(best)
    if short is left:
        return paired_short, paired_long
    return paired_long, paired_short


def contrast(name: str, site: List[Dict[str, object]], control: List[Dict[str, object]],
             statistic: str, *, n_bootstrap: int) -> Dict[str, object]:
    """One arm's comparison, grouped, with the plan's floors applied."""
    by_group_site: Dict[str, List[float]] = {}
    by_group_control: Dict[str, List[float]] = {}
    for rows, target in ((site, by_group_site), (control, by_group_control)):
        for row in rows:
            value = row.get(f"z_{statistic}")
            if value is None:
                continue
            target.setdefault(str(row["group"]), []).append(float(value))
    shared = sorted(set(by_group_site) & set(by_group_control))
    if not shared:
        return {"arm": name, "statistic": statistic, "runnable": False,
                "reason": "no group has both arms"}
    observed = {g: float(np.mean(by_group_site[g])) for g in shared}
    controls = {g: float(np.mean(by_group_control[g])) for g in shared}
    result = sq.compare_feature(f"{name}:{statistic}", observed, controls,
                                n_bootstrap=n_bootstrap, seed=SEED)
    entry = result.as_dict()
    lengths_site = float(np.mean([int(r["length"]) for r in site])) if site else float("nan")
    lengths_control = (float(np.mean([int(r["length"]) for r in control]))
                       if control else float("nan"))
    entry.update({
        "arm": name, "statistic": statistic, "runnable": True,
        "n_site": len(site), "n_control": len(control), "n_groups": len(shared),
        "mean_site_z": float(np.mean(list(observed.values()))),
        "mean_control_z": float(np.mean(list(controls.values()))),
        "mean_length_site": lengths_site,
        "mean_length_control": lengths_control,
        "length_confounded": bool(abs(lengths_site - lengths_control) > LENGTH_TOLERANCE),
    })
    return entry


def sequence_index(sequences: pd.DataFrame) -> Dict[Tuple[str, str], Dict[int, str]]:
    """Residue number to letter, per structure stem and chain."""
    out: Dict[Tuple[str, str], Dict[int, str]] = {}
    for row in sequences.itertuples():
        numbers = [int(n) for n in str(row.numbers).split(",") if n.strip()]
        letters = str(row.letters)
        if len(numbers) != len(letters):
            continue
        out[(str(row.stem), str(row.chain))] = dict(zip(numbers, letters))
    return out


def windows_for_pockets(pockets: pd.DataFrame,
                        index: Dict[Tuple[str, str], Dict[int, str]],
                        *, key_of, chain_of, group_of, kind: str = "local",
                        yields: Optional[Dict[str, int]] = None
                        ) -> List[Dict[str, object]]:
    """One row per pocket whose `local` window can carry an order null.

    Windows only: the order null costs 2000 shuffles of eight statistics per
    window, and a 2500-residue model has well over a hundred pockets of which
    this study scores ten. Scoring is deferred to `score_windows`, after the
    arms have chosen, which is the difference between minutes and hours.

    The seed travels with the row, so a window's null does not depend on which
    arm selected it or on how many pockets were scored before it.
    """
    rows: List[Dict[str, object]] = []
    dense_cache: Dict[Tuple[str, str], Tuple[np.ndarray, Dict[int, str]]] = {}
    for position, pocket in enumerate(pockets.itertuples()):
        key = (key_of(pocket), chain_of(pocket))
        numbers = index.get(key)
        if numbers is None:
            if yields is not None:
                yields["chain not in the sequence table"] = yields.get(
                    "chain not in the sequence table", 0) + 1
            continue
        residues = [int(r) for r in str(pocket.residues).split(",") if r.strip()]
        if kind == "local":
            text = local_window(numbers, residues)
            if text is None and yields is not None:
                reason = refusal(numbers, residues) or "refused"
                yields[reason] = yields.get(reason, 0) + 1
        else:
            if key not in dense_cache:
                dense_cache[key] = dense_starts(numbers)
            starts, texts = dense_cache[key]
            text = dense_window(starts, texts, residues)
            if text is None and yields is not None:
                yields["no dense stretch"] = yields.get("no dense stretch", 0) + 1
        if text is None:
            continue
        rows.append({
            "id": str(getattr(pocket, "uniprot_id", "")),
            "group": str(group_of(pocket)),
            "pocket_id": int(pocket.pocket_id),
            "rank": int(pocket.pocket_id),
            "length": len(text),
            "n_lining": len(residues),
            "window": text,
            # The seed is the window's, not the arm's: a window keyed by its
            # own letters gives the same null wherever it is selected, and the
            # same null for `local` and `dense` when they coincide.
            "seed": SEED + (hash_window(text) % 1_000_003),
        })
    return rows


def hash_window(text: str) -> int:
    """A stable integer for a window's letters, across processes and runs.

    `hash()` is salted per process, so it cannot seed anything that has to be
    reproducible from a committed command line.
    """
    return int(hashlib.sha256(text.encode("ascii")).hexdigest()[:12], 16)


def cap_pairs(left: List[Dict[str, object]], right: List[Dict[str, object]],
              limit: int, rng: np.random.Generator
              ) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Keep at most `limit` of the matched pairs, seeded (amendment 1 section 3)."""
    if limit <= 0 or len(left) <= limit:
        return left, right
    keep = sorted(rng.choice(len(left), size=limit, replace=False).tolist())
    return [left[i] for i in keep], [right[i] for i in keep]


def score_windows(rows: Sequence[Dict[str, object]], statistics: Sequence[str],
                  *, shuffles: int) -> None:
    """Fill each row's order-null z, in place, once and only once."""
    for row in rows:
        if row.get("scored"):
            continue
        for statistic in statistics:
            row[f"z_{statistic}"] = window_z(
                str(row["window"]), statistic, shuffles=shuffles,
                seed=int(row["seed"]))
        row["scored"] = True


def by_protein_of(rows: Sequence[Dict[str, object]]) -> Dict[str, List[Dict[str, object]]]:
    out: Dict[str, List[Dict[str, object]]] = {}
    for row in rows:
        out.setdefault(str(row["id"]), []).append(row)
    return out


def formation(rows: List[Dict[str, object]], *, n_bootstrap: int,
              statistics: Sequence[str], shuffles: int, arm: str = "formation"
              ) -> Tuple[List[Dict[str, object]], frozenset]:
    """Q9a as pre-registered: fpocket's top-ranked pockets against its bottom.

    Unchanged by amendment 1, which re-specifies the arm beside this one rather
    than in place of it. The second return value is the set of resampling groups
    this rule selected, which is what defines the amendment's held-out set - so
    that set is a function of the plan and not of me.
    """
    rng = np.random.default_rng(SEED)
    top: List[Dict[str, object]] = []
    bottom: List[Dict[str, object]] = []
    for _, found in sorted(by_protein_of(rows).items()):
        ordered = sorted(found, key=lambda r: int(r["rank"]))
        if len(ordered) < 2 * TOP_RANK:
            continue
        best = [r for r in ordered if int(r["rank"]) <= TOP_RANK]
        worst = ordered[-TOP_RANK:]
        if not best or not worst:
            continue
        matched_best, matched_worst = match_on_length(best, worst, rng)
        top.extend(matched_best)
        bottom.extend(matched_worst)
    score_windows(top + bottom, statistics, shuffles=shuffles)
    selected = frozenset(str(row["group"]) for row in top + bottom)
    return ([contrast(arm, top, bottom, statistic, n_bootstrap=n_bootstrap)
             for statistic in statistics], selected)


def formation_split(rows: List[Dict[str, object]], *, n_bootstrap: int,
                    statistics: Sequence[str], shuffles: int, arm: str,
                    skip_groups: frozenset = frozenset()) -> List[Dict[str, object]]:
    """Amendment 1 section 3: the rank split taken inside the surviving windows.

    Each protein's windows are ordered by fpocket rank and cut at the median,
    then length-matched across the cut. Rank is no longer collinear with length,
    because the comparison is made among the windows that exist rather than
    between a size class that survives the span cap and one that does not.

    `skip_groups` drops whole resampling groups, which is how the held-out arm
    is built: the groups the pre-registered rule selected are exactly the ones
    whose values have been looked at.
    """
    rng = np.random.default_rng(SEED + 3)
    better: List[Dict[str, object]] = []
    worse: List[Dict[str, object]] = []
    for _, found in sorted(by_protein_of(rows).items()):
        if any(str(row["group"]) in skip_groups for row in found):
            continue
        ordered = sorted(found, key=lambda r: int(r["rank"]))
        if len(ordered) < 2:
            continue
        half = len(ordered) // 2
        matched_better, matched_worse = match_on_length(
            ordered[:half], ordered[len(ordered) - half:], rng)
        matched_better, matched_worse = cap_pairs(
            matched_better, matched_worse, MAX_PAIRS_PER_PROTEIN, rng)
        better.extend(matched_better)
        worse.extend(matched_worse)
    better, worse = cap_pairs(better, worse, MAX_PAIRS_PER_ARM, rng)
    score_windows(better + worse, statistics, shuffles=shuffles)
    return [contrast(arm, better, worse, statistic, n_bootstrap=n_bootstrap)
            for statistic in statistics]


def pocket_type(rows: List[Dict[str, object]], *, n_bootstrap: int,
                statistics: Sequence[str], shuffles: int, arm: str = "type"
                ) -> List[Dict[str, object]]:
    """Q9b: do the letters anticipate the architecture class of the positions?

    The class is computed from positions alone and reads no letters, so the two
    sides of this comparison are independent by construction.
    """
    rng = np.random.default_rng(SEED + 1)
    faces: List[Dict[str, object]] = []
    others: List[Dict[str, object]] = []
    for row in rows:
        element = row.get("element")
        if element == "helix_face":
            faces.append(row)
        elif element == "none":
            others.append(row)
    matched_faces, matched_others = match_on_length(faces, others, rng)
    matched_faces, matched_others = cap_pairs(
        matched_faces, matched_others, MAX_PAIRS_PER_ARM, rng)
    score_windows(matched_faces + matched_others, statistics, shuffles=shuffles)
    return [contrast(arm, matched_faces, matched_others, statistic,
                     n_bootstrap=n_bootstrap)
            for statistic in statistics]


def sampled_windows(numbers: Dict[int, str], exclude: Sequence[int], length: int,
                    count: int, rng: np.random.Generator) -> List[str]:
    """Length-matched windows from the same chain, clear of the site."""
    ordered = sorted(numbers)
    banned = {int(p) for p in exclude}
    out: List[str] = []
    if len(ordered) <= length:
        return out
    for _ in range(count * 20):
        if len(out) >= count:
            break
        start = int(rng.integers(0, len(ordered) - length))
        chosen = ordered[start:start + length]
        if chosen[-1] - chosen[0] + 1 != length:
            continue
        if banned & set(chosen):
            continue
        text = "".join(numbers[n] for n in chosen)
        if len(text) < MIN_WINDOW or not set(text) <= STANDARD:
            continue
        out.append(text)
    return out


def function(sites: pd.DataFrame, index: Dict[Tuple[str, str], Dict[int, str]],
             *, statistics: Sequence[str], shuffles: int, per_site: int,
             kind: str = "local", yields: Optional[Dict[str, int]] = None
             ) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Q9c: an annotated contact shell against windows of its own chain.

    The controls are drawn from the shell's own chain, length-matched to it and
    clear of it, so the comparison is within one sequence and one organism.
    """
    rng = np.random.default_rng(SEED + 2)
    site_rows: List[Dict[str, object]] = []
    control_rows: List[Dict[str, object]] = []
    dense_cache: Dict[Tuple[str, str], Tuple[np.ndarray, Dict[int, str]]] = {}
    for row in sites.itertuples():
        stem = str(row.pdb_id).upper()
        key = (stem, str(row.chain))
        numbers = index.get(key)
        if numbers is None:
            key = (stem.lower(), str(row.chain))
            numbers = index.get(key)
        if numbers is None:
            if yields is not None:
                yields["chain not in the sequence table"] = yields.get(
                    "chain not in the sequence table", 0) + 1
            continue
        residues = [int(r) for r in str(row.residues).split(",") if r.strip()]
        if kind == "local":
            text = local_window(numbers, residues)
            if text is None and yields is not None:
                reason = refusal(numbers, residues) or "refused"
                yields[reason] = yields.get(reason, 0) + 1
        else:
            if key not in dense_cache:
                dense_cache[key] = dense_starts(numbers)
            starts, texts = dense_cache[key]
            text = dense_window(starts, texts, residues)
            if text is None and yields is not None:
                yields["no dense stretch"] = yields.get("no dense stretch", 0) + 1
        if text is None:
            continue
        group = str(getattr(row, "group", row.pdb_id))
        site_rows.append({"id": str(row.copy), "group": group,
                          "length": len(text), "window": text,
                          "n_lining": len(residues),
                          "seed": SEED + (hash_window(text) % 1_000_003)})
        for control in sampled_windows(numbers, residues, len(text), per_site, rng):
            control_rows.append({"id": f"{row.copy}:control", "group": group,
                                 "length": len(control), "window": control,
                                 "seed": SEED + (hash_window(control) % 1_000_003)})
    score_windows(site_rows + control_rows, statistics, shuffles=shuffles)
    return site_rows, control_rows


def holm_over(entries: List[Dict[str, object]], family: Sequence[str]) -> None:
    """Holm across one declared family, in place, with the plan's decision rule."""
    raw = {f"{e['arm']}:{e['statistic']}": float(e["p_value"])
           for e in entries
           if e.get("runnable") and f"{e['arm']}:{e['statistic']}" in family}
    corrected = sq.holm(raw) if raw else {}
    for entry in entries:
        if not entry.get("runnable"):
            continue
        key = f"{entry['arm']}:{entry['statistic']}"
        if key not in corrected:
            # Exploratory by declaration, not by outcome: an arm the amendment
            # named exploratory is in no family, so no Holm value exists for it
            # and `is_finding` cannot become true however the row reads.
            exploratory = str(entry["arm"]) not in CONFIRMATORY_ARMS
            entry["family"] = "exploratory" if exploratory else "secondary"
            entry["is_finding"] = False
            continue
        entry["family"] = "primary"
        entry["p_holm"] = corrected[key]
        entry["is_finding"] = bool(
            str(entry["arm"]) in CONFIRMATORY_ARMS
            and corrected[key] < 0.05 and entry["is_evidence"] and entry["exceeds_power"]
            and not (entry["interval"][0] <= 0.0 <= entry["interval"][1])
            and not entry["length_confounded"])


def null_spread(shuffles: int, n_windows: int, z: float = 1.0) -> float:
    """How far a fresh draw of the same windows' nulls moves a contrast.

    A ``composition_z`` against ``m`` shuffles carries Monte-Carlo error of about
    ``sqrt((1 + z**2 / 2) / m)`` from its null mean and standard deviation. An
    arm averages that over its windows, the contrast subtracts two such arms,
    and comparing two independent draws doubles the variance again - so two
    factors of root two on top of the per-window error.
    """
    if shuffles <= 0 or n_windows <= 0:
        return float("inf")
    per_window = math.sqrt((1.0 + z * z / 2.0) / float(shuffles))
    return 2.0 * per_window / math.sqrt(float(n_windows))


def look_reproduction(entries: Sequence[Dict[str, object]], *, shuffles: int
                      ) -> Dict[str, object]:
    """Whether this run is looking at the data amendment 1 section 1 disclosed.

    Selection and Monte-Carlo noise are judged separately, because only one of
    them can be exact. The selection rule reads ranks and lengths and never a
    seed, so the window and group counts must match with no tolerance at all;
    the difference is a mean of permutation z-scores, so it moves whenever the
    nulls are redrawn and is judged against the spread its own shuffle count
    implies. Both numbers are reported, so the arithmetic can be checked
    instead of the verdict being taken on trust.
    """
    found = [e for e in entries
             if e.get("runnable") and e.get("arm") == RECORDED_LOOK["arm"]
             and e.get("statistic") == RECORDED_LOOK["statistic"]]
    if not found:
        return {"recorded": RECORDED_LOOK, "verdict": "not reproduced",
                "reason": "the pre-registered arm did not run at all",
                "selection_reproduced": False}
    entry = found[0]
    n_site, n_groups = int(entry["n_site"]), int(entry["n_groups"])
    difference = float(entry["difference"])
    selection = (n_site == RECORDED_LOOK["n_site"]
                 and n_groups == RECORDED_LOOK["n_groups"])
    gap = abs(difference - float(RECORDED_LOOK["difference"]))
    spread = null_spread(shuffles, n_site)
    within = gap <= LOOK_SPREAD_BOUND * spread
    if selection and within:
        reason = ("the same windows in the same groups, and a difference within "
                  "the spread a redrawn null implies")
    elif not selection:
        reason = "the selection differs, which no re-seeding can cause"
    else:
        reason = "the selection holds but the difference moved further than a redrawn null explains"
    return {
        "recorded": RECORDED_LOOK, "n_site": n_site, "n_groups": n_groups,
        "difference": difference, "selection_reproduced": selection,
        "difference_gap": gap, "expected_spread": spread,
        "gap_in_spreads": (gap / spread) if spread > 0 else float("inf"),
        "spread_bound": LOOK_SPREAD_BOUND,
        "verdict": "reproduced" if (selection and within) else "not reproduced",
        "reason": reason,
    }


def run(sequences: pd.DataFrame, pockets: pd.DataFrame, sites: Optional[pd.DataFrame],
        classes: Optional[pd.DataFrame], *, shuffles: int, n_bootstrap: int,
        per_site: int) -> Dict[str, object]:
    """Every arm, on the statistics the plan fixed, in both window definitions.

    The plan's arms (`formation`, `type`, `function`) run on the plan's `local`
    window and are the only ones that may become a finding, together with the
    one confirmatory arm amendment 1 adds. Everything suffixed `_dense`, and
    `formation_matched`, are exploratory by declaration and can never be
    promoted however they come out.
    """
    statistics = list(PRIMARY) + list(SECONDARY)
    index = sequence_index(sequences)
    by_pocket_class: Dict[Tuple[str, int], str] = {}
    if classes is not None and not classes.empty:
        for row in classes.itertuples():
            by_pocket_class[(str(row.uniprot_id), int(row.pocket_id))] = str(row.element)

    yields: Dict[str, Dict[str, int]] = {}
    pocket_rows: Dict[str, List[Dict[str, object]]] = {}
    for kind in WINDOW_KINDS:
        yields[f"pockets/{kind}"] = {}
        pocket_rows[kind] = windows_for_pockets(
            pockets, index,  # windows only; the null runs after the arms choose
            key_of=lambda p: str(p.model).replace(".pdb", "").replace(".ent", ""),
            # An AlphaFold model is one chain, always A, so the lookup is exact
            # rather than a guess; the deposited structures take their chain
            # from the census row instead, in `function`.
            chain_of=lambda p: "A",
            group_of=lambda p: str(getattr(p, "cluster", p.uniprot_id)),
            kind=kind, yields=yields[f"pockets/{kind}"])
        for row in pocket_rows[kind]:
            row["element"] = by_pocket_class.get((str(row["id"]), int(row["pocket_id"])))

    entries: List[Dict[str, object]] = []
    shared = dict(n_bootstrap=n_bootstrap, statistics=statistics, shuffles=shuffles)

    # The pre-registered rule first, because its selection defines the held-out
    # set the amendment's confirmatory arm runs on.
    pre_entries, selected = formation(pocket_rows["local"], arm="formation", **shared)
    entries += pre_entries
    entries += formation(pocket_rows["dense"], arm="formation_dense", **shared)[0]
    entries += formation_split(pocket_rows["local"], arm="formation_matched", **shared)
    entries += formation_split(pocket_rows["dense"], arm="formation_matched_dense",
                               **shared)
    entries += formation_split(pocket_rows["local"], arm="formation_heldout",
                               skip_groups=selected, **shared)
    entries += pocket_type(pocket_rows["local"], arm="type", **shared)
    entries += pocket_type(pocket_rows["dense"], arm="type_dense", **shared)

    counts: Dict[str, int] = {}
    if sites is not None and not sites.empty:
        for kind in WINDOW_KINDS:
            arm = "function" if kind == "local" else "function_dense"
            yields[f"sites/{kind}"] = {}
            site_rows, control_rows = function(
                sites, index, statistics=statistics, shuffles=shuffles,
                per_site=per_site, kind=kind, yields=yields[f"sites/{kind}"])
            counts[f"n_site_windows_{kind}"] = len(site_rows)
            counts[f"n_control_windows_{kind}"] = len(control_rows)
            entries += [contrast(arm, site_rows, control_rows, statistic,
                                 n_bootstrap=n_bootstrap) for statistic in statistics]

    # The plan's family, untouched by the amendment.
    primary_family = [f"{arm}:{statistic}"
                      for arm in ("formation", "function") for statistic in PRIMARY]
    holm_over(entries, primary_family)
    # Q9b asks a different question of a different contrast, so its own family.
    type_family = [f"type:{statistic}" for statistic in PRIMARY]
    holm_over([e for e in entries if e.get("arm") == "type"], type_family)
    # And amendment 1's one confirmatory arm, in a third family of its own.
    heldout_family = [f"formation_heldout:{statistic}" for statistic in PRIMARY]
    holm_over([e for e in entries if e.get("arm") == "formation_heldout"], heldout_family)

    reproduction = look_reproduction(entries, shuffles=shuffles)

    report: Dict[str, object] = {
        "seed": SEED, "shuffles": shuffles,
        "look_reproduction": reproduction,
        "primary_family": primary_family,
        "type_family": type_family,
        "heldout_family": heldout_family,
        "confirmatory_arms": list(CONFIRMATORY_ARMS),
        "max_window": MAX_WINDOW, "min_window": MIN_WINDOW,
        "dense_window": DENSE_WINDOW,
        "top_rank": TOP_RANK, "length_tolerance": LENGTH_TOLERANCE,
        "max_pairs_per_protein": MAX_PAIRS_PER_PROTEIN,
        "max_pairs_per_arm": MAX_PAIRS_PER_ARM,
        "heldout_groups_excluded": sorted(selected),
        "window_yield": yields,
        "n_pocket_windows": len(pocket_rows["local"]),
        "n_dense_windows": len(pocket_rows["dense"]),
        "n_scored_windows": int(sum(1 for rows in pocket_rows.values()
                                    for r in rows if r.get("scored"))),
        "n_scored_local": int(sum(1 for r in pocket_rows["local"] if r.get("scored"))),
        "n_scored_dense": int(sum(1 for r in pocket_rows["dense"] if r.get("scored"))),
        "n_classified": int(sum(1 for r in pocket_rows["local"] if r.get("element"))),
        "n_site_windows": counts.get("n_site_windows_local", 0),
        "n_control_windows": counts.get("n_control_windows_local", 0),
        "tests": entries,
    }
    report.update(counts)
    return report


def markdown(report: Dict[str, object]) -> str:
    """The report, with the decomposition stated where it is read."""
    lines = [
        "# Does the sequence predict pocket formation, type, or function?",
        "",
        f"**Run.** `scripts/pocket_sequence.py analyse`, seed {report['seed']}, "
        f"{report['shuffles']} shuffles per window.",
        "",
        "**Pre-registered** in `docs/POCKET_SEQUENCE_PLAN.md`, dated before any sequence "
        "was read from any structure for this study.",
        "",
        "## What is being tested, and what is not",
        "",
        "A window of sequence factorises the way a lining of positions does — "
        "**composition** (which letters) **⊕ order** (in what arrangement) — and the null "
        "for order is permutation of the window's own letters. That is the same quotient "
        "`architecture.py` takes on positions.",
        "",
        "**Every primary test below is on `composition_z`, never on a raw measure.** A raw "
        "difference is already explained: these sites are basic and pocket linings are "
        "hydrophobic. A non-zero z says the *arrangement* carries something the letters "
        "alone do not, which is the only form of the question this study calls open.",
        "",
        f"{report['n_pocket_windows']} `local` pocket windows and "
        f"{report.get('n_dense_windows', 0)} `dense` ones, of which "
        f"{report['n_classified']} carry an architecture class; "
        f"{report['n_site_windows']} annotated-site windows against "
        f"{report['n_control_windows']} length-matched controls from the same chains. "
        f"{report.get('n_scored_windows', 0)} windows carried an order null.",
        "",
        "**Amendment 1** (`docs/POCKET_SEQUENCE_PLAN_AMENDMENT_1.md`, dated after the "
        "first run and before any of what it adds was read) re-specifies the formation "
        "arm, which the pre-registered rule could only run on ten windows out of two "
        "thousand, and adds a second window definition. Of what it adds, **only "
        "`formation_heldout` may become a finding** — the rule restricted to the "
        "resampling groups the pre-registered rule selected nothing from, so disjoint "
        "from everything already looked at. Every arm named `*_dense`, and "
        "`formation_matched`, is **exploratory by declaration**: no Holm family, and "
        "`is_finding` cannot become true for it however the row reads.",
        "",
        "The three arms: **formation** is fpocket's best-ranked pockets against its worst "
        f"within the same protein (rank ≤ {report['top_rank']}), which asks whether the "
        "sequence anticipates a geometric judgement that reads no sequence; **type** is "
        "pockets the architecture module calls a helical face against those it calls none, "
        "which asks whether the letters anticipate the positions; **function** is an "
        "annotated inositol-phosphate contact shell against length-matched windows of its "
        "own chain.",
        "",
        "## Primary",
        "",
        "| arm | statistic | windows | groups | mean z (site) | mean z (control) | "
        "difference | 95% interval | Holm | could have detected | finding |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]

    def row_for(entry: Dict[str, object]) -> str:
        low, high = entry["interval"]
        holm = entry.get("p_holm")
        verdict = "**yes**" if entry.get("is_finding") else (
            f"**no evidence** — {entry['n_groups']} groups" if not entry["is_evidence"]
            else "no")
        return (
            f"| {entry['arm']} | `{entry['statistic']}` | {entry['n_site']} | "
            f"{entry['n_groups']} | {entry['mean_site_z']:+.2f} | "
            f"{entry['mean_control_z']:+.2f} | {entry['difference']:+.3f} | "
            f"[{low:+.3f}, {high:+.3f}] | "
            f"{'—' if holm is None else format(holm, '.4f')} | "
            f"{entry['minimum_detectable']:.3f} | {verdict} |"
        )

    primary = [e for e in report["tests"]  # type: ignore[union-attr]
               if e.get("runnable") and e.get("family") == "primary"]
    for entry in primary:
        entry.setdefault("confirmatory", str(entry["arm"]) in CONFIRMATORY_ARMS)
    for entry in sorted(primary, key=lambda e: (str(e["arm"]), str(e["statistic"]))):
        lines.append(row_for(entry))
    if not primary:
        lines.append("| — | — | — | — | — | — | — | — | — | — | nothing runnable |")

    findings = [e for e in primary if e.get("is_finding")]
    lines += [
        "",
        f"**{len(findings)} of {len(primary)}** primary tests is a finding. A finding needs "
        "a Holm-corrected p below 0.05, an interval excluding zero, five resampling groups, "
        "an observed difference larger than the minimum detectable one, and matched window "
        "lengths — all five, fixed by the plan before any of this ran.",
        "",
        "## Secondary and exploratory, reported and never promoted",
        "",
        "| arm | statistic | windows | groups | mean z (site) | mean z (control) | "
        "difference | 95% interval | could have detected | family |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for entry in sorted((e for e in report["tests"]  # type: ignore[union-attr]
                         if e.get("runnable") and e.get("family") != "primary"),
                        key=lambda e: (str(e["arm"]), str(e["statistic"]))):
        low, high = entry["interval"]
        lines.append(
            f"| {entry['arm']} | `{entry['statistic']}` | {entry['n_site']} | "
            f"{entry['n_groups']} | {entry['mean_site_z']:+.2f} | "
            f"{entry['mean_control_z']:+.2f} | {entry['difference']:+.3f} | "
            f"[{low:+.3f}, {high:+.3f}] | {entry['minimum_detectable']:.3f} | "
            f"{entry.get('family')}"
            f"{' — length-confounded' if entry['length_confounded'] else ''} |"
        )
    skipped = [e for e in report["tests"] if not e.get("runnable")]  # type: ignore[union-attr]
    if skipped:
        lines += ["", "Not runnable:", ""]
        lines += [f"- {e['arm']}:{e['statistic']} — {e['reason']}" for e in skipped]
    yields = report.get("window_yield") or {}
    if yields:
        lines += ["", "## Why a lining produced no window", "",
                  "Amendment 1 section 6. A refused lining is not a pocket without "
                  "sequence signal; it is a pocket whose positions the window "
                  "definition cannot carry a period over.", "",
                  "| source | reason | linings |", "|---|---|---|"]
        for source in sorted(yields):
            for reason, count in sorted(yields[source].items(),
                                        key=lambda kv: (-kv[1], kv[0])):
                lines.append(f"| `{source}` | {reason} | {count} |")
    excluded = report.get("heldout_groups_excluded") or []
    reproduction = report.get("look_reproduction") or {}
    verdict = str(reproduction.get("verdict", "unknown"))
    recorded = reproduction.get("recorded") or {}
    arithmetic = (
        f"Amendment 1 section 1 records {recorded.get('n_site')} windows across "
        f"{recorded.get('n_groups')} groups at a difference of "
        f"{recorded.get('difference')}; this run selected "
        f"{reproduction.get('n_site')} across {reproduction.get('n_groups')} and "
        f"gives {float(reproduction.get('difference', float('nan'))):.6f}, a gap of "
        f"{float(reproduction.get('difference_gap', float('nan'))):.5f} against an "
        f"expected redraw spread of "
        f"{float(reproduction.get('expected_spread', float('nan'))):.5f} "
        f"({float(reproduction.get('gap_in_spreads', float('nan'))):.2f} of it, bound "
        f"{reproduction.get('spread_bound')}).")
    if verdict == "reproduced":
        audit = (
            "The pre-registered arm reproduced the disclosed look: the selection is "
            "identical, and the difference moved only as far as a redrawn null "
            f"explains. {arithmetic} So this is the same pocket table that was looked "
            "at, and the exclusion above is the set of groups that look covered."
        )
    elif reproduction.get("selection_reproduced"):
        audit = (
            "**The selection reproduced but the difference did not.** "
            f"{arithmetic} The windows and groups are the ones that were looked at, so "
            "the exclusion above is the right set; but something other than the null's "
            "seed has moved the contrast, and that **must be explained before the "
            "held-out arm is read as confirmatory**."
        )
    else:
        audit = (
            "**The pre-registered arm did not select what it selected before** "
            f"({reproduction.get('reason')}). {arithmetic} The selection rule reads "
            "ranks and lengths and never a seed, so the pocket table is not the one "
            "that was looked at, the exclusion above is not the set that look covered, "
            "and the held-out arm **must be read as exploratory** until the difference "
            "is explained."
        )
    lines += [
        "",
        "## The held-out arm",
        "",
        f"The pre-registered rule selected windows from {len(excluded)} resampling "
        "groups, and those are exactly the groups whose values were read before "
        "amendment 1 was written. `formation_heldout` runs the amendment's rule with "
        "all of them removed, so it is a confirmatory test on data no decision in the "
        "amendment could have been tuned to. Its Holm family is its own two primary "
        "statistics, and the plan's five decision rules apply to it unchanged.",
        "",
        "**The held-out claim is only as good as that reproduction.** " + audit,
        "",
        "## What this cannot conclude",
        "",
        "* fpocket's rank is a geometric judgement of a predicted structure. The formation "
        "arm predicts fpocket, not nature.",
        "* A crystal contact shell is one deposited conformation. The function arm predicts "
        "where an inositol phosphate sat in that crystal, not affinity and not specificity.",
        "* Every arm is within-protein, so none of this is a sequence-only predictor of "
        "anything, and none of it generalises to sequences with no structure.",
        "* The transfer ceiling in `results/architecture/ARCHITECTURE.md` stands over the "
        "periodicity statistics here as it does there.",
    ]
    return "\n".join(lines) + "\n"


def read_sequences(paths: Sequence[Path]) -> pd.DataFrame:
    """One sequence table from however many the data steps wrote.

    The predicted models and the deposited structures are read by different
    jobs, because only one of them needs the network, and the arms need both.
    """
    frames = [pd.read_csv(path) for path in paths if Path(path).exists()]
    if not frames:
        raise SystemExit(f"no sequence table exists among {[str(p) for p in paths]}")
    return pd.concat(frames, ignore_index=True)


def cmd_sites(args: argparse.Namespace) -> int:
    """The function arm's data step: contact shells, and nothing else.

    Q9c needs each annotated copy's shell and its chain's sequence. It does not
    need a pocket, so this runs no detector - which is the difference between
    minutes and the within-protein study's hours - and it reuses that study's
    census filter rather than writing a second one, so the two arms select the
    same copies for the same stated reason.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import within_protein as wp  # noqa: E402  (a sibling script, not a package)
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    copies = wp.select_copies(pd.read_csv(args.census), args.comp_id)
    if args.limit:
        copies = copies.head(int(args.limit))
    rows, excluded = [], {}
    for index, row in enumerate(copies.itertuples(), start=1):
        key = f"{row.pdb_id}_{row.chain}_{row.resseq}"
        stem = str(row.pdb_id).lower()
        found = sorted(args.structures.glob(f"{stem}.*")) + \
            sorted(args.structures.glob(f"{stem.upper()}.*"))
        if not found:
            excluded["structure not downloaded"] = excluded.get(
                "structure not downloaded", 0) + 1
            continue
        icode = "" if pd.isna(getattr(row, "icode", "")) else str(getattr(row, "icode", ""))
        try:
            arrays = load_structure_arrays(found[0])
            contacts = wp.contact_residues(arrays, str(row.chain), int(row.resseq),
                                           icode, comp_id=args.comp_id)
        except Exception as error:  # noqa: BLE001 - one bad entry is not the run
            reason = f"{type(error).__name__}"
            excluded[reason] = excluded.get(reason, 0) + 1
            print(f"[{index}/{len(copies)}] {key}: {reason}", flush=True)
            continue
        by_chain = contacts["by_chain"]
        if not contacts["found"] or str(row.chain) not in by_chain:
            reason = ("ligand copy not present in the structure" if not contacts["found"]
                      else "no contact with the ligand's own chain")
            excluded[reason] = excluded.get(reason, 0) + 1
            continue
        residues = list(by_chain[str(row.chain)])
        rows.append({"copy": key, "pdb_id": row.pdb_id, "chain": row.chain,
                     "group": getattr(row, "homology_group_strict", row.pdb_id),
                     "burial_class": getattr(row, "burial_class", ""),
                     "n_site_residues": len(residues),
                     "residues": ",".join(str(int(r)) for r in residues)})
        print(f"[{index}/{len(copies)}] {key}: {len(residues)} site residues", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    report = {"n_copies": int(len(copies)), "n_measured": len(rows),
              "n_groups": len({r["group"] for r in rows}), "excluded": excluded}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True, default=str))
    return 0


def cmd_classes(args: argparse.Namespace) -> int:
    """The type arm's data step: each pocket's architecture class, from positions.

    Only pockets that produced a window are classified, because only those can
    enter an arm; the class reads positions and never letters, which is what
    makes Q9b's two sides independent by construction.
    """
    from cryptic_ip.detect.architecture import ArchitectureError, describe_architecture

    sequences = read_sequences(args.sequences)
    pockets = pd.read_csv(args.pockets)
    index = sequence_index(sequences)
    wanted = set()
    for kind in WINDOW_KINDS:
        for row in windows_for_pockets(
                pockets, index,
                key_of=lambda p: str(p.model).replace(".pdb", "").replace(".ent", ""),
                chain_of=lambda p: "A",
                group_of=lambda p: str(getattr(p, "cluster", p.uniprot_id)),
                kind=kind):
            wanted.add((str(row["id"]), int(row["pocket_id"])))

    rows, failures = [], {}
    for pocket in pockets.itertuples():
        key = (str(getattr(pocket, "uniprot_id", "")), int(pocket.pocket_id))
        if key not in wanted:
            continue
        residues = sorted({int(r) for r in str(pocket.residues).split(",") if r.strip()})
        try:
            described = describe_architecture(residues, n_draws=args.draws, seed=SEED)
        except ArchitectureError as error:
            reason = str(error)[:80]
            failures[reason] = failures.get(reason, 0) + 1
            continue
        rows.append({"uniprot_id": key[0], "pocket_id": key[1],
                     "n_positions": len(residues), "element": described.element,
                     "architecture_class": described.architecture_class})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(args.out, index=False)
    counts = ({} if frame.empty
              else {str(k): int(v) for k, v in frame["element"].value_counts().items()})
    print(json.dumps({"n_windows": len(wanted), "n_classified": len(rows),
                      "elements": counts, "failures": failures},
                     indent=2, sort_keys=True))
    return 0


def cmd_analyse(args: argparse.Namespace) -> int:
    sequences = read_sequences(args.sequences)
    pockets = pd.read_csv(args.pockets)
    sites = pd.read_csv(args.sites) if args.sites and args.sites.exists() else None
    classes = pd.read_csv(args.classes) if args.classes and args.classes.exists() else None
    report = run(sequences, pockets, sites, classes, shuffles=args.shuffles,
                 n_bootstrap=args.bootstrap, per_site=args.controls_per_site)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")
    text = markdown(report)
    if args.markdown:
        args.markdown.write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    grab = sub.add_parser("sequences", help="chain sequences from structures, no detection")
    grab.add_argument("--structures", type=Path, nargs="+", required=True)
    grab.add_argument("--out", type=Path, required=True)

    shells = sub.add_parser("sites", help="annotated contact shells, no pocket detection")
    shells.add_argument("--census", type=Path, default=Path("results/redocking/census.csv"))
    shells.add_argument("--structures", type=Path, required=True)
    shells.add_argument("--comp-id", default="IHP")
    shells.add_argument("--out", type=Path, required=True)
    shells.add_argument("--report", type=Path, default=None)
    shells.add_argument("--limit", type=int, default=0)

    label = sub.add_parser("classes", help="each pocket's architecture class, from positions")
    label.add_argument("--sequences", type=Path, nargs="+", required=True)
    label.add_argument("--pockets", type=Path, required=True)
    label.add_argument("--out", type=Path, required=True)
    label.add_argument("--draws", type=int, default=CLASS_DRAWS)

    analyse = sub.add_parser("analyse", help="the three arms")
    analyse.add_argument("--sequences", type=Path, nargs="+", required=True)
    analyse.add_argument("--pockets", type=Path, required=True)
    analyse.add_argument("--sites", type=Path, default=None)
    analyse.add_argument("--classes", type=Path, default=None)
    analyse.add_argument("--out", type=Path, required=True)
    analyse.add_argument("--markdown", type=Path, default=None)
    analyse.add_argument("--shuffles", type=int, default=SHUFFLES)
    analyse.add_argument("--bootstrap", type=int, default=2000)
    analyse.add_argument("--controls-per-site", type=int, default=5)

    args = parser.parse_args(argv)
    return {"sequences": cmd_sequences, "sites": cmd_sites,
            "classes": cmd_classes, "analyse": cmd_analyse}[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
