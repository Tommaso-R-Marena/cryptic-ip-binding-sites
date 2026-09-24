#!/usr/bin/env python3
"""Pre-registered exploration of zero-parameter descriptors (docs/EXPLORATION_PLAN.md).

Two stages, each appending one line to a ledger:

``explore``
    Development rows only. Every descriptor, signed by its a-priori
    direction, is scored on pooled AUROC and within-structure recovery, with
    group-bootstrap intervals over strict homology groups and
    Benjamini-Hochberg across the whole family. The confirmatory set is chosen
    by the plan's rule and written to the ledger.

``confirm``
    Holdout rows, for the confirmatory set read back from the ledger, once.

    python scripts/explore_descriptors.py explore --table table.csv.gz \\
        --ledger results/exploration/ledger.jsonl --out-dir results/exploration
    python scripts/explore_descriptors.py confirm --table table.csv.gz \\
        --ledger results/exploration/ledger.jsonl --out-dir results/exploration
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from math import comb
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.benchmark import protocol  # noqa: E402

SEED = 20260924
Q = 0.05
MAX_CONFIRMATORY = 3
TOP_K = 3
GROUP = "group_strict"

# +1: higher predicts binding; -1: lower does; 0: no prior (scored as +1, two-sided).
DIRECTIONS: Dict[str, int] = {
    **{name: +1 for name in (
        "n_basic_residues", "n_strong_basic_residues", "basic_fraction", "n_basic_residues_core",
        "n_basic_nitrogens", "net_formal_charge", "positive_charge_density", "charge_balance",
        "coulomb_potential_kt", "enclosure", "burial_depth", "buried_residue_fraction", "hull_depth",
        "n_hydroxyl_residues", "hydroxyl_fraction", "polar_fraction", "n_aromatic_residues",
        "aromatic_fraction", "alpha_sphere_density",
    )},
    **{name: -1 for name in (
        "n_acidic_residues", "acidic_fraction", "n_acidic_residues_core", "basic_nitrogen_min_distance",
        "basic_nitrogen_mean_distance", "basic_nitrogen_dispersion", "mean_relative_sasa", "sasa_mean",
        "sasa_median", "sasa_min", "sasa_max", "hydropathy_mean",
    )},
    **{name: 0 for name in (
        "pocket_volume", "hull_volume", "n_alpha_spheres", "radius_of_gyration", "asphericity",
        "max_extent", "sasa_total", "n_residues",
    )},
}
COMPOSITES: Dict[str, Sequence[str]] = {
    "electropositive_enclosure": ("coulomb_potential_kt", "enclosure"),
    "basic_cluster": ("n_basic_nitrogens", "basic_nitrogen_dispersion"),
}
# Evaluated on these kinds of data before this plan: reported, never confirmatory.
NOT_BLIND = frozenset({"rule_score", "burial_depth", "hull_depth"})
REFERENCE = "rule_score"
TASKS = {"ip_site": "inferential", "cryptic_ip_site": "descriptive"}


class LedgerError(RuntimeError):
    """The ledger forbids this run."""


# ------------------------------------------------------------------ scores
def signed(values: pd.Series, sign: int) -> np.ndarray:
    """``sign x value`` with missing values ranked least binding-like."""
    s = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float) * (1 if sign >= 0 else -1)
    finite = np.isfinite(s)
    floor = (s[finite].min() - 1.0) if finite.any() else 0.0
    return np.where(finite, s, floor)


def within_structure_rank(values: np.ndarray, structures: np.ndarray) -> np.ndarray:
    """Percentile rank of each pocket among its own structure's pockets (ties averaged)."""
    return pd.Series(values).groupby(pd.Series(structures)).rank(pct=True, method="average").to_numpy()


def descriptor_scores(rows: pd.DataFrame) -> Dict[str, np.ndarray]:
    """Every pre-registered score, signed so that higher means more binding-like."""
    missing = [name for name in DIRECTIONS if name not in rows.columns]
    if missing:
        raise SystemExit(f"table lacks descriptors {missing}")
    assert set(DIRECTIONS) == set(protocol.BENCHMARK_FEATURES), "directions must cover every descriptor"
    scores = {name: signed(rows[name], sign) for name, sign in DIRECTIONS.items()}
    structures = rows["structure_id"].to_numpy()
    for name, parts in COMPOSITES.items():
        scores[name] = sum(within_structure_rank(scores[p], structures) for p in parts)
    scores[REFERENCE] = signed(rows[REFERENCE], +1)
    return scores


def prior(name: str) -> int:
    if name in COMPOSITES:
        return +1
    if name == REFERENCE:
        return +1
    return DIRECTIONS[name]


# ---------------------------------------------------------------- recovery
def _hit_probability(ordered_labels: Sequence[np.ndarray], k: int) -> float:
    """P(any positive in the top ``k``) with random tie-breaking.

    ``ordered_labels`` holds each tie block's labels, highest score first.
    """
    taken = 0
    for block in ordered_labels:
        m, p = len(block), int(np.sum(block))
        if taken + m <= k:
            if p:
                return 1.0
            taken += m
            if taken == k:
                return 0.0
            continue
        r = k - taken  # draw r of this block's m at random
        return 1.0 - comb(m - p, r) / comb(m, r)
    return 0.0


def recovery_table(y: np.ndarray, score: np.ndarray, structures: np.ndarray, groups: np.ndarray,
                   k: int) -> pd.DataFrame:
    """Per structure with a positive: top-``k`` hit probability, and its chance rate."""
    frame = pd.DataFrame({"y": y, "s": score, "structure": structures, "group": groups})
    out = []
    for structure, part in frame.groupby("structure", sort=True):
        labels = part["y"].to_numpy()
        n, p = len(labels), int(labels.sum())
        if p == 0:
            continue
        blocks = [g["y"].to_numpy() for _, g in part.groupby("s", sort=True)][::-1]
        kk = min(k, n)
        chance = 1.0 - comb(n - p, kk) / comb(n, kk)
        out.append({
            "structure": structure, "group": part["group"].iloc[0], "pockets": n, "positives": p,
            "hit": _hit_probability(blocks, kk), "chance": chance,
        })
    return pd.DataFrame(out)


def group_mean_statistic(values: np.ndarray, groups: np.ndarray) -> Callable[[np.ndarray], float]:
    """Mean of ``values`` giving each group equal weight, as a function of bootstrap weights."""
    sizes = pd.Series(groups).map(pd.Series(groups).value_counts()).to_numpy(dtype=float)

    def statistic(w: np.ndarray) -> float:
        weight = w / sizes
        total = weight.sum()
        return float(np.sum(weight * values) / total) if total > 0 else float("nan")

    return statistic


# ------------------------------------------------------------- inference
def benjamini_hochberg(p_values: Mapping[str, float]) -> Dict[str, float]:
    """BH-adjusted p-values (step-up, monotone)."""
    items = sorted(p_values.items(), key=lambda kv: kv[1])
    m = len(items)
    adjusted, running = {}, 1.0
    for i in range(m - 1, -1, -1):
        name, p = items[i]
        running = min(running, p * m / (i + 1))
        adjusted[name] = min(1.0, running)
    return adjusted


def evaluate(rows: pd.DataFrame, task: str, n_bootstrap: int, inferential: bool) -> Dict[str, dict]:
    """Pooled AUROC and within-structure recovery for every score."""
    y = rows[f"label_{task}"].to_numpy(dtype=int)
    groups = rows[GROUP].astype(str).to_numpy()
    structures = rows["structure_id"].astype(str).to_numpy()
    results: Dict[str, dict] = {}
    for name, score in descriptor_scores(rows).items():
        ranked = protocol.Ranked(y, score)
        entry: Dict[str, object] = {"prior": prior(name)}
        rec1 = recovery_table(y, score, structures, groups, 1)
        rec3 = recovery_table(y, score, structures, groups, TOP_K)
        if inferential:
            entry["auroc"] = protocol.bootstrap_statistic(
                ranked.roc_auc, groups, null=0.5, n_bootstrap=n_bootstrap, seed=SEED).as_dict()
            for label, rec in (("top1", rec1), (f"top{TOP_K}", rec3)):
                excess = (rec["hit"] - rec["chance"]).to_numpy()
                stat = group_mean_statistic(excess, rec["group"].to_numpy())
                est = protocol.bootstrap_statistic(stat, rec["group"].to_numpy(), n_bootstrap=n_bootstrap, seed=SEED)
                entry[f"{label}_excess"] = est.as_dict()
                entry[f"{label}_hit"] = float(group_mean_statistic(rec["hit"].to_numpy(), rec["group"].to_numpy())(
                    np.ones(len(rec))))
                entry[f"{label}_chance"] = float(group_mean_statistic(
                    rec["chance"].to_numpy(), rec["group"].to_numpy())(np.ones(len(rec))))
        else:
            entry["auroc_point"] = ranked.roc_auc(np.ones(len(y)))
            per_group = rec1.assign(excess=rec1["hit"] - rec1["chance"]).groupby("group").agg(
                structures=("structure", "size"), top1_hit=("hit", "mean"), top1_chance=("chance", "mean"))
            entry["top1_by_group"] = {g: {k: float(v) for k, v in r.items()} for g, r in per_group.iterrows()}
        results[name] = entry
    return results


def select_confirmatory(results: Mapping[str, dict]) -> Dict[str, object]:
    """Apply BH across the family and the plan's rule for the confirmatory set."""
    family = {n: r for n, r in results.items() if n != REFERENCE}
    p = {}
    for name, r in family.items():
        p[f"{name}:auroc"] = r["auroc"]["p_value"]
        p[f"{name}:top1"] = r["top1_excess"]["p_value"]
    q = benjamini_hochberg(p)
    eligible = []
    for name, r in family.items():
        r["q_auroc"], r["q_top1"] = q[f"{name}:auroc"], q[f"{name}:top1"]
        above = r["auroc"]["point"] > 0.5
        r["verdict"] = (
            "no prior" if r["prior"] == 0 else
            ("matches prior" if above else "contradicts prior") if r["q_auroc"] < Q else "not significant"
        )
        if (r["q_auroc"] < Q and r["prior"] != 0 and above and name not in NOT_BLIND
                and r["top1_excess"]["point"] > 0):
            eligible.append(name)
    eligible.sort(key=lambda n: -family[n]["auroc"]["low"])
    return {"family_size": len(p), "q": Q, "eligible": eligible, "confirmatory_set": eligible[:MAX_CONFIRMATORY]}


def paired_vs_reference(rows: pd.DataFrame, task: str, names: Sequence[str], n_bootstrap: int) -> Dict[str, dict]:
    """AUROC(name) - AUROC(rule_score), paired over the same rows and resamples."""
    y = rows[f"label_{task}"].to_numpy(dtype=int)
    groups = rows[GROUP].astype(str).to_numpy()
    scores = descriptor_scores(rows)
    reference = protocol.Ranked(y, scores[REFERENCE])
    out = {}
    for name in names:
        ranked = protocol.Ranked(y, scores[name])
        out[name] = protocol.bootstrap_statistic(
            lambda w, a=ranked: a.roc_auc(w) - reference.roc_auc(w), groups,
            n_bootstrap=n_bootstrap, seed=SEED).as_dict()
    return out


# ----------------------------------------------------------------- ledger
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def read_ledger(path: Path) -> List[dict]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def append_ledger(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as fh:
        fh.write(json.dumps(record, sort_keys=True, default=float) + "\n")


def _record(stage: str, table: Path, **fields) -> dict:
    return {
        "stage": stage, "table_sha256": sha256(table), "commit": git_commit(),
        "time": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), **fields,
    }


def _task_rows(table: pd.DataFrame, task: str, holdout: bool) -> pd.DataFrame:
    rows = table[(table[f"label_{task}"] >= 0) & (table["holdout"].astype(bool) == holdout)]
    return rows.reset_index(drop=True)


# ---------------------------------------------------------------- stages
def cmd_explore(args: argparse.Namespace) -> int:
    table = pd.read_csv(args.table, low_memory=False)
    ledger = read_ledger(args.ledger)
    digest = sha256(args.table)
    if any(r["stage"] == "explore" and r["table_sha256"] == digest for r in ledger) and not args.rerun:
        raise LedgerError("this table was already explored; results are in the ledger (use --rerun to append again)")
    dev = table[~table["holdout"].astype(bool)].reset_index(drop=True)
    assert not dev["holdout"].astype(bool).any(), "exploration must never see holdout rows"
    tasks = {}
    for task, role in TASKS.items():
        rows = _task_rows(dev, task, holdout=False)
        if rows[f"label_{task}"].nunique() < 2:
            tasks[task] = {"role": role, "not_evaluable": "one class"}
            continue
        results = evaluate(rows, task, args.n_bootstrap, inferential=(role == "inferential"))
        entry = {"role": role, "rows": int(len(rows)), "positives": int(rows[f"label_{task}"].sum()),
                 "positive_groups": int(rows.loc[rows[f"label_{task}"] == 1, GROUP].nunique()),
                 "results": results}
        if role == "inferential":
            selection = select_confirmatory(results)
            eligible = list(selection["eligible"])
            selection["vs_rule_score"] = paired_vs_reference(rows, task, eligible, args.n_bootstrap)
            entry["selection"] = selection
        tasks[task] = entry
    record = _record("explore", args.table, n_bootstrap=args.n_bootstrap, seed=SEED, tasks=tasks)
    append_ledger(args.ledger, record)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "EXPLORATION.md").write_text(render_explore(record))
    print((args.out_dir / "EXPLORATION.md").read_text())
    return 0


def cmd_confirm(args: argparse.Namespace) -> int:
    table = pd.read_csv(args.table, low_memory=False)
    digest = sha256(args.table)
    ledger = read_ledger(args.ledger)
    explored = [r for r in ledger if r["stage"] == "explore" and r["table_sha256"] == digest]
    if not explored:
        raise LedgerError("no exploration of this table in the ledger; explore first")
    if any(r["stage"] == "confirm" and r["table_sha256"] == digest for r in ledger):
        raise LedgerError("confirmation already ran on this table; it runs once")
    selection = explored[0]["tasks"]["ip_site"].get("selection", {})
    names = list(selection.get("confirmatory_set", []))
    hold = _task_rows(table, "ip_site", holdout=True)
    result: Dict[str, object] = {"confirmatory_set": names, "rows": int(len(hold)),
                                 "positives": int(hold["label_ip_site"].sum()) if len(hold) else 0,
                                 "positive_groups": int(hold.loc[hold["label_ip_site"] == 1, GROUP].nunique())}
    if names and hold["label_ip_site"].nunique() == 2:
        y = hold["label_ip_site"].to_numpy(dtype=int)
        groups = hold[GROUP].astype(str).to_numpy()
        scores = descriptor_scores(hold)
        estimates = {n: protocol.bootstrap_statistic(protocol.Ranked(y, scores[n]).roc_auc, groups, null=0.5,
                                                     n_bootstrap=args.n_bootstrap, seed=SEED) for n in names}
        adjusted = protocol.holm({n: e.p_value for n, e in estimates.items()})
        result["results"] = {
            n: {**e.as_dict(), "holm_p": adjusted[n],
                "decision": "confirmed" if e.low > 0.5 and adjusted[n] < Q else "not confirmed (underpowered)"}
            for n, e in estimates.items()
        }
    record = _record("confirm", args.table, n_bootstrap=args.n_bootstrap, seed=SEED, confirm=result)
    append_ledger(args.ledger, record)
    print(json.dumps(record, indent=2, default=float))
    return 0


# ---------------------------------------------------------------- report
def _ci(est: Mapping[str, float], digits: int = 3) -> str:
    return f"{est['point']:.{digits}f} [{est['low']:.{digits}f}, {est['high']:.{digits}f}]"


def render_explore(record: Mapping[str, object]) -> str:
    lines = [f"## Descriptor exploration (table {record['table_sha256'][:12]}, commit {record['commit'][:10]})", ""]
    for task, entry in record["tasks"].items():
        if "not_evaluable" in entry:
            lines += [f"### {task}: not evaluable ({entry['not_evaluable']})", ""]
            continue
        lines += [f"### {task} ({entry['role']}): {entry['rows']} pockets, {entry['positives']} positives "
                  f"in {entry['positive_groups']} strict groups", ""]
        results = entry["results"]
        if entry["role"] == "inferential":
            lines += ["| score | prior | AUROC [95% CI] | q | top-1 hit / chance | top-1 excess [95% CI] | q "
                      "| verdict |",
                      "|---|---|---|---|---|---|---|---|"]
            order = sorted(results, key=lambda n: -results[n]["auroc"]["point"])
            for n in order:
                r = results[n]
                lines.append(
                    f"| {n} | {r['prior']:+d} | {_ci(r['auroc'])} | {r.get('q_auroc', float('nan')):.2g} | "
                    f"{r['top1_hit']:.2f} / {r['top1_chance']:.2f} | {_ci(r['top1_excess'])} | "
                    f"{r.get('q_top1', float('nan')):.2g} | {r.get('verdict', 'reference')} |")
            sel = entry["selection"]
            lines += ["", f"BH across {sel['family_size']} tests at q = {sel['q']}. "
                      f"Eligible: {', '.join(sel['eligible']) or 'none'}. "
                      f"Confirmatory set: {', '.join(sel['confirmatory_set']) or 'none'}.", ""]
            if sel["vs_rule_score"]:
                lines += ["| eligible score | AUROC - rule_score [95% CI] | p |", "|---|---|---|"]
                for n, est in sel["vs_rule_score"].items():
                    lines.append(f"| {n} | {_ci(est)} | {est['p_value']:.2g} |")
                lines.append("")
        else:
            lines += ["Descriptive only (plan): pooled AUROC and top-1 hits per strict group; no p-values.", "",
                      "| score | prior | AUROC | top-1 hits by group (hit / chance, structures) |", "|---|---|---|---|"]
            for n in sorted(results, key=lambda n: -results[n]["auroc_point"]):
                r = results[n]
                cells = "; ".join(f"{g}: {v['top1_hit']:.2f}/{v['top1_chance']:.2f} ({int(v['structures'])})"
                                  for g, v in r["top1_by_group"].items())
                lines.append(f"| {n} | {r['prior']:+d} | {r['auroc_point']:.3f} | {cells} |")
            lines.append("")
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("explore", "confirm"):
        p = sub.add_parser(name)
        p.add_argument("--table", type=Path, required=True)
        p.add_argument("--ledger", type=Path, default=ROOT / "results/exploration/ledger.jsonl")
        p.add_argument("--out-dir", type=Path, default=ROOT / "results/exploration")
        p.add_argument("--n-bootstrap", type=int, default=2000)
        if name == "explore":
            p.add_argument("--rerun", action="store_true", help="append another exploration of the same table")
    args = parser.parse_args(argv)
    return {"explore": cmd_explore, "confirm": cmd_confirm}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
