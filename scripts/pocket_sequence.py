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
import json
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
    paths = sorted(p for p in args.structures.iterdir()
                   if p.suffix.lower() in {".pdb", ".ent"})
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
                        *, key_of, chain_of, group_of) -> List[Dict[str, object]]:
    """One row per pocket whose `local` window can carry an order null.

    Windows only: the order null costs 2000 shuffles of eight statistics per
    window, and a 2500-residue model has well over a hundred pockets of which
    this study scores ten. Scoring is deferred to `score_windows`, after the
    arms have chosen, which is the difference between minutes and hours.

    The seed travels with the row, so a window's null does not depend on which
    arm selected it or on how many pockets were scored before it.
    """
    rows: List[Dict[str, object]] = []
    for position, pocket in enumerate(pockets.itertuples()):
        numbers = index.get((key_of(pocket), chain_of(pocket)))
        if numbers is None:
            continue
        residues = [int(r) for r in str(pocket.residues).split(",") if r.strip()]
        text = local_window(numbers, residues)
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
            "seed": SEED + position,
        })
    return rows


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


def formation(rows: List[Dict[str, object]], *, n_bootstrap: int,
              statistics: Sequence[str], shuffles: int) -> List[Dict[str, object]]:
    """Q9a: fpocket's best-ranked pockets against its worst, within each protein."""
    rng = np.random.default_rng(SEED)
    top: List[Dict[str, object]] = []
    bottom: List[Dict[str, object]] = []
    by_protein: Dict[str, List[Dict[str, object]]] = {}
    for row in rows:
        by_protein.setdefault(str(row["id"]), []).append(row)
    for protein, found in by_protein.items():
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
    return [contrast("formation", top, bottom, statistic, n_bootstrap=n_bootstrap)
            for statistic in statistics]


def pocket_type(rows: List[Dict[str, object]], *, n_bootstrap: int,
                statistics: Sequence[str], shuffles: int) -> List[Dict[str, object]]:
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
    score_windows(matched_faces + matched_others, statistics, shuffles=shuffles)
    return [contrast("type", matched_faces, matched_others, statistic,
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
             *, n_bootstrap: int, statistics: Sequence[str], shuffles: int,
             per_site: int) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Q9c: an annotated contact shell against windows of its own chain."""
    rng = np.random.default_rng(SEED + 2)
    site_rows: List[Dict[str, object]] = []
    control_rows: List[Dict[str, object]] = []
    for position, row in enumerate(sites.itertuples()):
        stem = str(row.pdb_id).upper()
        numbers = index.get((stem, str(row.chain))) or index.get(
            (stem.lower(), str(row.chain)))
        if numbers is None:
            continue
        residues = [int(r) for r in str(row.residues).split(",") if r.strip()]
        text = local_window(numbers, residues)
        if text is None:
            continue
        group = str(getattr(row, "group", row.pdb_id))
        entry: Dict[str, object] = {"id": str(row.copy), "group": group,
                                    "length": len(text), "window": text,
                                    "n_lining": len(residues)}
        for statistic in statistics:
            entry[f"z_{statistic}"] = window_z(text, statistic, shuffles=shuffles,
                                               seed=SEED + position)
        site_rows.append(entry)
        for offset, control in enumerate(
                sampled_windows(numbers, residues, len(text), per_site, rng)):
            item: Dict[str, object] = {"id": f"{row.copy}:control{offset}",
                                       "group": group, "length": len(control),
                                       "window": control}
            for statistic in statistics:
                item[f"z_{statistic}"] = window_z(
                    control, statistic, shuffles=shuffles,
                    seed=SEED + 10_000 + position * per_site + offset)
            control_rows.append(item)
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
            entry["family"] = "secondary"
            entry["is_finding"] = False
            continue
        entry["family"] = "primary"
        entry["p_holm"] = corrected[key]
        entry["is_finding"] = bool(
            corrected[key] < 0.05 and entry["is_evidence"] and entry["exceeds_power"]
            and not (entry["interval"][0] <= 0.0 <= entry["interval"][1])
            and not entry["length_confounded"])


def run(sequences: pd.DataFrame, pockets: pd.DataFrame, sites: Optional[pd.DataFrame],
        classes: Optional[pd.DataFrame], *, shuffles: int, n_bootstrap: int,
        per_site: int) -> Dict[str, object]:
    """All three arms, on the statistics the plan fixed."""
    statistics = list(PRIMARY) + list(SECONDARY)
    index = sequence_index(sequences)
    by_pocket_class: Dict[Tuple[str, int], str] = {}
    if classes is not None and not classes.empty:
        for row in classes.itertuples():
            by_pocket_class[(str(row.uniprot_id), int(row.pocket_id))] = str(row.element)

    pocket_rows = windows_for_pockets(
        pockets, index,  # windows only; the null runs after the arms choose
        key_of=lambda p: str(p.model).replace(".pdb", "").replace(".ent", ""),
        # An AlphaFold model is one chain, always A, so the lookup is exact
        # rather than a guess; the deposited structures take their chain from
        # the census row instead, in `function`.
        chain_of=lambda p: "A",
        group_of=lambda p: str(getattr(p, "cluster", p.uniprot_id)))
    for row in pocket_rows:
        row["element"] = by_pocket_class.get((str(row["id"]), int(row["pocket_id"])))

    entries: List[Dict[str, object]] = []
    entries += formation(pocket_rows, n_bootstrap=n_bootstrap, statistics=statistics,
                         shuffles=shuffles)
    entries += pocket_type(pocket_rows, n_bootstrap=n_bootstrap, statistics=statistics,
                           shuffles=shuffles)
    site_rows: List[Dict[str, object]] = []
    control_rows: List[Dict[str, object]] = []
    if sites is not None and not sites.empty:
        site_rows, control_rows = function(
            sites, index, n_bootstrap=n_bootstrap, statistics=statistics,
            shuffles=shuffles, per_site=per_site)
        entries += [contrast("function", site_rows, control_rows, statistic,
                             n_bootstrap=n_bootstrap) for statistic in statistics]

    primary_family = [f"{arm}:{statistic}"
                      for arm in ("formation", "function") for statistic in PRIMARY]
    holm_over(entries, primary_family)
    type_family = [f"type:{statistic}" for statistic in PRIMARY]
    type_entries = [e for e in entries if e.get("arm") == "type"]
    holm_over(type_entries, type_family)

    return {
        "seed": SEED, "shuffles": shuffles,
        "primary_family": primary_family,
        "type_family": type_family,
        "max_window": MAX_WINDOW, "min_window": MIN_WINDOW,
        "top_rank": TOP_RANK, "length_tolerance": LENGTH_TOLERANCE,
        "n_pocket_windows": len(pocket_rows),
        "n_scored_windows": int(sum(1 for r in pocket_rows if r.get("scored"))),
        "n_site_windows": len(site_rows),
        "n_control_windows": len(control_rows),
        "n_classified": int(sum(1 for r in pocket_rows if r.get("element"))),
        "tests": entries,
    }


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
        f"{report['n_pocket_windows']} pocket windows, of which "
        f"{report['n_classified']} carry an architecture class; "
        f"{report['n_site_windows']} annotated-site windows against "
        f"{report['n_control_windows']} length-matched controls from the same chains.",
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
        "## Secondary, reported and never promoted",
        "",
        "| arm | statistic | windows | groups | difference | Holm family | length-confounded |",
        "|---|---|---|---|---|---|---|",
    ]
    for entry in sorted((e for e in report["tests"]  # type: ignore[union-attr]
                         if e.get("runnable") and e.get("family") != "primary"),
                        key=lambda e: (str(e["arm"]), str(e["statistic"]))):
        lines.append(
            f"| {entry['arm']} | `{entry['statistic']}` | {entry['n_site']} | "
            f"{entry['n_groups']} | {entry['difference']:+.3f} | {entry.get('family')} | "
            f"{'yes' if entry['length_confounded'] else 'no'} |"
        )
    skipped = [e for e in report["tests"] if not e.get("runnable")]  # type: ignore[union-attr]
    if skipped:
        lines += ["", "Not runnable:", ""]
        lines += [f"- {e['arm']}:{e['statistic']} — {e['reason']}" for e in skipped]
    lines += [
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


def cmd_analyse(args: argparse.Namespace) -> int:
    sequences = pd.read_csv(args.sequences)
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
    grab.add_argument("--structures", type=Path, required=True)
    grab.add_argument("--out", type=Path, required=True)

    analyse = sub.add_parser("analyse", help="the three arms")
    analyse.add_argument("--sequences", type=Path, required=True)
    analyse.add_argument("--pockets", type=Path, required=True)
    analyse.add_argument("--sites", type=Path, default=None)
    analyse.add_argument("--classes", type=Path, default=None)
    analyse.add_argument("--out", type=Path, required=True)
    analyse.add_argument("--markdown", type=Path, default=None)
    analyse.add_argument("--shuffles", type=int, default=SHUFFLES)
    analyse.add_argument("--bootstrap", type=int, default=2000)
    analyse.add_argument("--controls-per-site", type=int, default=5)

    args = parser.parse_args(argv)
    return cmd_sequences(args) if args.command == "sequences" else cmd_analyse(args)


if __name__ == "__main__":
    raise SystemExit(main())
