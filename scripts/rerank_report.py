#!/usr/bin/env python3
"""Study F, report (docs/RERANK_PLAN.md): F1-F4 from the docking shards of scripts/rerank.py.

    python scripts/rerank_report.py --census census.csv --arms-dir arms --out-dir out [--primary-dir p]

``--primary-dir`` holds the redocking primary arm's shards (``primary_*.jsonl``):
Vina's top-pose RMSD per seed there and here are compared as a reproducibility
check. Writes rerank.json, RERANK.md and report.html.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from cryptic_ip.rescoring import crossfit as cf  # noqa: E402

ARRESTIN_ENTRIES = {"1ZSH", "5TV1", "7F1W", "7F1X", "7JTB", "7JXA", "7MOR"}


def read_records(directory: Path, pattern: str) -> List[dict]:
    out = []
    for path in sorted(directory.rglob(pattern)):
        out += [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    return out


def build_runs(records: Sequence[dict], groups: Dict[str, str]) -> List[cf.Run]:
    runs = []
    for rec in records:
        for r in rec.get("runs", []):
            if r.get("seed", 0) < 1 or not r.get("vina"):
                continue
            runs.append(cf.Run(rec["copy_key"], groups[rec["copy_key"]], int(r["seed"]),
                               np.asarray(r["vina"], float), np.asarray(r["rmsd"], float),
                               np.asarray(r["eel"], float)))
    return runs


def reproducibility(records: Sequence[dict], primary: Sequence[dict]) -> Dict[str, object]:
    theirs = {}
    for rec in primary:
        for r in rec.get("runs", []):
            if str(r.get("arm", "")).startswith("vina_s") and r.get("top_rmsd") is not None:
                theirs[(rec["copy_key"], int(r["seed"]))] = float(r["top_rmsd"])
    pairs = []
    for rec in records:
        for r in rec.get("runs", []):
            key = (rec["copy_key"], int(r.get("seed", 0)))
            if key in theirs and r.get("rmsd"):
                pairs.append((theirs[key], float(r["rmsd"][0])))
    if not pairs:
        return {"pairs": 0}
    a, b = np.array(pairs).T
    return {"pairs": len(pairs), "same_success": float(np.mean((a <= cf.SUCCESS) == (b <= cf.SUCCESS))),
            "within_0.5A": float(np.mean(np.abs(a - b) <= 0.5))}


def report(census: pd.DataFrame, records: Sequence[dict], primary: Sequence[dict], n_bootstrap: int,
           n_permutations: int) -> Dict[str, object]:
    census = census.copy()
    census["copy_key"] = census["copy_key"].astype(str)
    groups = dict(zip(census["copy_key"], census["homology_group_strict"].astype(str)))
    burial = dict(zip(census["copy_key"], census["burial_class"].astype(str)))
    runs = build_runs(records, groups)
    result = cf.evaluate(runs, n_bootstrap=n_bootstrap, n_permutations=n_permutations)
    errors = [r for r in records if "error" in r]
    result["accounting"] = {"records": len(records), "errors": len(errors),
                            "not_reached": sum("not reached" in r["error"] for r in errors),
                            "error_examples": [r["error"][:120] for r in errors[:10]]}
    copies = sorted({r.copy_key for r in runs})
    if copies:
        codes = pd.factorize(pd.Series([groups[c] for c in copies]))[0]
        folds = cf.fold_of_groups(codes)
        chosen = result["chosen_w_by_fold"]
        chosen_w = {c: chosen[folds[i]] for i, c in enumerate(copies)}
        strata = {cls: {k for k in copies if burial.get(k) == cls} for cls in ("surface", "semi_cryptic", "cryptic")}
        strata["classic arrestins"] = {k for k in copies if k.split(":")[0].upper() in ARRESTIN_ENTRIES}
        result["F4"] = {name: cf.stratum(runs, keys, chosen_w) for name, keys in strata.items()}
        rho = []
        for r in runs:
            if len(r.eel) >= 3 and np.ptp(r.eel) > 0 and np.ptp(r.rmsd) > 0:
                from scipy.stats import spearmanr
                rho.append(spearmanr(r.eel, r.rmsd).correlation)
        result["spearman_eel_rmsd"] = {"runs": len(rho), "median": float(np.median(rho)) if rho else None}
    crystal = [(r["crystal_minimised_vina"], r["crystal_minimised_eel"]) for rec in records
               for r in rec.get("runs", []) if r.get("seed") == 0]
    if crystal:
        v, e = np.array(crystal).T
        result["crystal_minimised"] = {"copies": len(crystal), "median_vina": float(np.median(v)),
                                       "median_eel": float(np.median(e))}
    result["reproducibility_vs_redocking_primary"] = reproducibility(records, primary)
    result["plan"] = "docs/RERANK_PLAN.md"
    return result


def fmt(est: Optional[dict]) -> str:
    if not est or "point" not in est or est["point"] is None:
        return "–"
    return f"{est['point']:.3f} [{est['low']:.3f}, {est['high']:.3f}]"


def markdown(res: Dict[str, object]) -> str:
    f1 = res.get("F1", {})
    lines = ["## Electrostatic re-ranking of docked IP poses (docs/RERANK_PLAN.md)", "",
             f"{res.get('copies', 0)} copies in {res.get('groups', 0)} strict groups, {res.get('runs', 0)} seed runs.",
             "", f"**F1:** {f1.get('decision')}.", ""]
    if "group" in f1:
        lines += ["| estimand | Vina | re-ranked (out of fold) | difference |", "|---|---|---|---|"]
        for e in ("group", "copy"):
            lines.append(f"| {e} | {fmt(f1[e]['vina'])} | {fmt(f1[e]['reranked'])} | {fmt(f1[e]['difference'])} |")
        lines += ["", f"Chosen w by fold: {res.get('chosen_w_by_fold')}.", ""]
    if "F2" in res:
        lines.append(f"**F2:** per-run shares {json.dumps(res['F2']['per_run_shares'])}; sampling ceiling "
                     f"{fmt(res['F2']['sampling_ceiling_group'])}.")
    if "F3" in res:
        f3 = res["F3"]
        lines.append(f"**F3:** {f3['n']} permutations, mean {f3['mean']:.3f}, max {f3['max']:.3f}, "
                     f"fraction ≥ observed {f3['fraction_at_least_observed']:.2f}.")
    for name, s in (res.get("F4") or {}).items():
        note = "" if s.get("evidence") else " (fewer than 5 groups: not evidence)"
        if s.get("copies"):
            lines.append(f"- {name}: {s['copies']} copies, {s['groups']} groups; Vina {s['vina']:.3f}, "
                         f"re-ranked {s['reranked']:.3f}{note}")
    lines += ["", f"Accounting: {json.dumps(res.get('accounting'))}",
              f"Reproducibility vs the redocking primary arm: "
              f"{json.dumps(res.get('reproducibility_vs_redocking_primary'))}"]
    return "\n".join(lines) + "\n"


def page(res: Dict[str, object]) -> str:
    from study_report_page import decisions, forest
    from study_report_page import page as make_page

    f1 = res.get("F1", {})
    parts = [decisions([("F1 electrostatic re-ranking", str(f1.get("decision")),
                         f"group difference {fmt(f1.get('group', {}).get('difference'))}")])]
    if "group" in f1:
        parts.append(forest([("difference (group)", f1["group"]["difference"], ""),
                             ("difference (copy)", f1["copy"]["difference"], "")], null=0.0,
                            caption="re-ranked minus Vina top-pose success at 2 Å"))
        parts.append(forest([("Vina", f1["group"]["vina"], ""), ("re-ranked", f1["group"]["reranked"], ""),
                             ("sampling ceiling", res["F2"]["sampling_ceiling_group"], "")], null=0.5,
                            caption="top-pose success, group estimand"))
    return make_page("Electrostatic re-ranking", "docs/RERANK_PLAN.md", [("Decision", parts)])


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--arms-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--primary-dir", type=Path, default=None)
    parser.add_argument("--n-bootstrap", type=int, default=cf.N_BOOTSTRAP)
    parser.add_argument("--n-permutations", type=int, default=cf.N_PERMUTATIONS)
    args = parser.parse_args(argv)
    census = pd.read_csv(args.census, dtype={"icode": str})
    records = read_records(args.arms_dir, "rerank_*.jsonl")
    primary = read_records(args.primary_dir, "primary_*.jsonl") if args.primary_dir else []
    res = report(census, records, primary, args.n_bootstrap, args.n_permutations)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "rerank.json").write_text(json.dumps(res, indent=2, default=float))
    text = markdown(res)
    (args.out_dir / "RERANK.md").write_text(text)
    (args.out_dir / "report.html").write_text(page(res))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
