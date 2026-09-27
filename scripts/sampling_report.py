#!/usr/bin/env python3
"""Study G, report (docs/SAMPLING_PLAN.md): G1-G4 from the arms of scripts/sampling.py.

    python scripts/sampling_report.py --census census.csv --e32-dir rerank_arms \
        --arms-dir arms --out-dir results/sampling

``--e32-dir`` holds study F's pose lists (``rerank_*.jsonl``), which are this
study's exhaustiveness-32 arm: the same receptor, ligand, box, seeds and starting
poses, so they are reused rather than re-docked.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from cryptic_ip.benchmark import protocol  # noqa: E402
from cryptic_ip.docking import stats  # noqa: E402

SUCCESS = 2.0
#: docs/SAMPLING_PLAN.md freezes the re-ranking weight at study F's modal fold choice.
FROZEN_W = 0.1
SEED = 20260930


def read_records(directory: Optional[Path], pattern: str) -> List[dict]:
    if directory is None or not directory.exists():
        return []
    out = []
    for path in sorted(directory.rglob(pattern)):
        out += [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    return out


def run_metrics(run: dict, w: float = FROZEN_W) -> Optional[Dict[str, float]]:
    """Ceiling, Vina top-pose and re-ranked top-pose outcomes for one seed run."""
    vina, rmsd, eel = (np.asarray(run.get(k) or [], dtype=float) for k in ("vina", "rmsd", "eel"))
    if not len(vina) or not np.isfinite(rmsd).all():
        return None
    top = int(np.argmin(vina))
    reranked = int(np.argmin(vina + w * eel))
    return {"ceiling": float(np.min(rmsd) <= SUCCESS), "vina": float(rmsd[top] <= SUCCESS),
            "reranked": float(rmsd[reranked] <= SUCCESS), "poses": float(len(vina)),
            "seconds": float(run.get("seconds") or np.nan)}


def per_copy(records: Sequence[dict]) -> Dict[Tuple[str, int], Dict[str, float]]:
    """(copy_key, seed) -> metrics, for every seed run with a usable pose list."""
    out = {}
    for rec in records:
        for run in rec.get("runs", []) or []:
            m = run_metrics(run)
            if m is not None:
                out[(str(rec["copy_key"]), int(run["seed"]))] = m
    return out


def arm_table(runs: Dict[Tuple[str, int], Dict[str, float]], keys: Sequence[str]) -> pd.DataFrame:
    """Per copy, the mean over its seed runs of each metric (nan where the copy has none)."""
    rows = {}
    for key in keys:
        mine = [m for (k, _), m in runs.items() if k == key]
        rows[key] = ({name: float(np.mean([m[name] for m in mine])) for name in
                      ("ceiling", "vina", "reranked", "poses")} | {"seeds": len(mine),
                     "seconds": float(np.nansum([m["seconds"] for m in mine]))}) if mine else {}
    columns = ["ceiling", "vina", "reranked", "poses", "seeds", "seconds"]
    return pd.DataFrame.from_dict(rows, orient="index").reindex(index=keys, columns=columns)


def decide(est: Optional[dict], n_groups: int, holm_p: Optional[float], *, good: str, null_label: str,
           margin: float = 0.05) -> str:
    if n_groups < stats.MIN_GROUPS:
        return f"not evaluable: {n_groups} strict groups (fewer than {stats.MIN_GROUPS})"
    if not est or not np.isfinite(est.get("low", np.nan)):
        return "not evaluable"
    if est["low"] > 0 and (holm_p is not None and holm_p < 0.05):
        return good
    if est["high"] < margin:
        return null_label
    return "inconclusive"


def build(census: pd.DataFrame, e32: Sequence[dict], arms: Dict[str, Sequence[dict]],
          n_bootstrap: int) -> Dict[str, object]:
    census = census.copy()
    census["copy_key"] = census["copy_key"].astype(str)
    group_of = dict(zip(census["copy_key"], census["homology_group_strict"].astype(str)))
    burial = dict(zip(census["copy_key"], census["burial_class"].astype(str)))

    runs = {"e32": per_copy(e32), **{name: per_copy(recs) for name, recs in arms.items()}}
    shared = sorted({k for k, _ in runs["e32"]} & {k for k, _ in runs.get("e128", {})})
    tables = {name: arm_table(r, shared) for name, r in runs.items()}
    groups = np.array([group_of.get(k, "?") for k in shared])
    n_groups = int(pd.Series(groups).nunique())

    report: Dict[str, object] = {"plan": "docs/SAMPLING_PLAN.md", "frozen_w": FROZEN_W,
                                 "copies_compared": len(shared), "groups": n_groups,
                                 "arms": {name: {"copies": int(t["seeds"].notna().sum()),
                                                 "seed_runs": int(np.nansum(t["seeds"])),
                                                 "mean_poses": float(np.nanmean(t["poses"])),
                                                 "median_seconds_per_copy": float(np.nanmedian(t["seconds"]))}
                                          for name, t in tables.items()}}
    report["accounting"] = {name: {"records": len(recs), "errors": sum(1 for r in recs if "error" in r),
                                   "not_reached": sum(1 for r in recs
                                                      if "not reached" in str(r.get("error", "")))}
                            for name, recs in ({"e32": e32} | dict(arms)).items()}
    if not shared:
        report["decisions"] = {"G1": {"decision": "not evaluable: no copy in both arms"}}
        return report

    def mean(name: str, metric: str) -> dict:
        return stats.mean_estimates(tables[name][metric].to_numpy(), groups, null=0.5,
                                    n_bootstrap=n_bootstrap, seed=SEED)

    report["levels"] = {name: {m: mean(name, m) for m in ("ceiling", "vina", "reranked")}
                        for name in tables if name in ("e32", "e128")}

    pairs = {"G1": ("ceiling", "e128", "ceiling", "e32"), "G2": ("vina", "e128", "vina", "e32"),
             "G3": ("reranked", "e128", "vina", "e32")}
    questions = {"G1": "sampling ceiling, E128 - E32", "G2": "top-pose success, E128 - E32",
                 "G3": "re-ranked (w = 0.1) top-pose success at E128 - Vina top-pose at E32"}
    estimates, p_values = {}, {}
    for name, (m1, a1, m0, a0) in pairs.items():
        est = stats.paired_difference(tables[a1][m1].to_numpy(), tables[a0][m0].to_numpy(), groups,
                                      n_bootstrap=n_bootstrap, seed=SEED)
        estimates[name] = est
        if est.get("per_group") and n_groups >= stats.MIN_GROUPS:
            p_values[name] = est["per_group"]["p_value"]
    holm = protocol.holm(p_values) if p_values else {}
    labels = {"G1": ("budget-limited", "search-saturated"), "G2": ("improves", "no gain"),
              "G3": ("improves", "no gain")}
    report["decisions"] = {
        name: {"question": questions[name], "estimate": estimates[name].get("per_group"),
               "per_copy": estimates[name].get("per_copy"), "holm_p": holm.get(name),
               "decision": decide(estimates[name].get("per_group"), n_groups, holm.get(name),
                                  good=labels[name][0], null_label=labels[name][1])}
        for name in pairs}

    # G4: rescue rate among seed runs whose E32 list held no near-native pose (descriptive).
    missed = [(k, s) for (k, s), m in runs["e32"].items() if m["ceiling"] == 0.0 and k in set(shared)]
    rescue: Dict[str, object] = {"e32_seed_runs_without_a_near_native_pose": len(missed)}
    for name in ("e128", "e512"):
        got = [(k, s) for (k, s) in missed if (k, s) in runs.get(name, {})]
        if got:
            y = np.array([runs[name][ks]["ceiling"] for ks in got])
            g = np.array([group_of.get(k, "?") for k, _ in got])
            rescue[name] = {"covered_runs": len(got), "groups": int(pd.Series(g).nunique()),
                            "rescued": stats.mean_estimates(y, g, null=0.0, n_bootstrap=n_bootstrap, seed=SEED)}
    report["G4_rescue"] = rescue

    # The E512 ladder on the subset, seed 1 only, and the strata.
    sub = sorted({k for k, s in runs.get("e512", {}) if s == 1})
    if sub:
        ladder = {}
        for name in ("e32", "e128", "e512"):
            have = [k for k in sub if (k, 1) in runs.get(name, {})]
            g = np.array([group_of.get(k, "?") for k in have])
            ladder[name] = {m: stats.mean_estimates(np.array([runs[name][(k, 1)][m] for k in have]), g,
                                                    n_bootstrap=n_bootstrap, seed=SEED)
                            for m in ("ceiling", "vina")} | {"copies": len(have)}
        report["E512_ladder_seed1"] = ladder
    report["strata"] = {}
    for label_, mask in {cls: np.array([burial.get(k) == cls for k in shared])
                         for cls in ("surface", "semi_cryptic", "cryptic")}.items():
        if mask.sum():
            report["strata"][label_] = {
                "copies": int(mask.sum()),
                "ceiling_e32": stats.mean_estimates(tables["e32"]["ceiling"].to_numpy()[mask], groups[mask],
                                                    n_bootstrap=n_bootstrap, seed=SEED),
                "ceiling_e128": stats.mean_estimates(tables["e128"]["ceiling"].to_numpy()[mask], groups[mask],
                                                     n_bootstrap=n_bootstrap, seed=SEED),
                "vina_e128": stats.mean_estimates(tables["e128"]["vina"].to_numpy()[mask], groups[mask],
                                                  n_bootstrap=n_bootstrap, seed=SEED)}
    report["notes"] = [
        "G4 conditions on an E32 outcome, so it describes where a gain lands rather than evidencing it.",
        "Docking scores for a -9 polyanion from functions without explicit electrostatics are weak evidence "
        "whatever the search budget.",
    ]
    return report


def fmt(est: Optional[dict]) -> str:
    if not est or "point" not in est or est["point"] is None:
        return "-"
    return f"{est['point']:.3f} [{est['low']:.3f}, {est['high']:.3f}]"


def markdown(r: Dict[str, object]) -> str:
    lines = ["## The sampling ceiling (docs/SAMPLING_PLAN.md)", "",
             f"{r['copies_compared']} copies in {r['groups']} strict groups compared across arms; "
             f"the re-ranking weight is frozen at w = {r['frozen_w']}.", "",
             "| arm | copies | seed runs | mean poses | median s/copy |", "|---|---|---|---|---|"]
    for name, a in r["arms"].items():
        lines.append(f"| {name} | {a['copies']} | {a['seed_runs']} | {a['mean_poses']:.1f} | "
                     f"{a['median_seconds_per_copy']:.0f} |")
    lines += ["", "### Decisions (Holm across G1-G3)", "",
              "| | question | difference (group) | Holm p | decision |", "|---|---|---|---|---|"]
    for name, d in r["decisions"].items():
        hp = d.get("holm_p")
        lines.append(f"| {name} | {d.get('question', '')} | {fmt(d.get('estimate'))} | "
                     f"{'-' if hp is None else f'{hp:.3g}'} | **{d['decision']}** |")
    if "levels" in r:
        lines += ["", "### Levels (group estimand)", "", "| arm | ceiling | Vina top pose | re-ranked top pose |",
                  "|---|---|---|---|"]
        for name, m in r["levels"].items():
            lines.append(f"| {name} | {fmt(m['ceiling'].get('per_group'))} | {fmt(m['vina'].get('per_group'))} | "
                         f"{fmt(m['reranked'].get('per_group'))} |")
    g4 = r.get("G4_rescue", {})
    if g4:
        lines += ["", f"**G4 (descriptive).** {g4['e32_seed_runs_without_a_near_native_pose']} seed runs had no "
                  "near-native pose at E32."]
        for name in ("e128", "e512"):
            if name in g4:
                lines.append(f"- {name}: {fmt(g4[name]['rescued'].get('per_group'))} of "
                             f"{g4[name]['covered_runs']} such runs find one.")
    for name, lad in (r.get("E512_ladder_seed1") or {}).items():
        lines.append(f"- ladder {name} (seed 1, {lad['copies']} copies): ceiling "
                     f"{fmt(lad['ceiling'].get('per_group'))}, top pose {fmt(lad['vina'].get('per_group'))}")
    lines += ["", "### Strata (ceiling and success)", "",
              "| stratum | copies | ceiling E32 | ceiling E128 | top pose E128 |", "|---|---|---|---|---|"]
    for name, s in (r.get("strata") or {}).items():
        note = "" if s["ceiling_e128"].get("evidence") else " (not evidence)"
        lines.append(f"| {name}{note} | {s['copies']} | {fmt(s['ceiling_e32'].get('per_group'))} | "
                     f"{fmt(s['ceiling_e128'].get('per_group'))} | {fmt(s['vina_e128'].get('per_group'))} |")
    lines += ["", f"Accounting: {json.dumps(r.get('accounting'))}", "", *[f"- {n}" for n in r.get("notes", [])]]
    return "\n".join(lines) + "\n"


def page(r: Dict[str, object]) -> str:
    from study_report_page import decisions, forest
    from study_report_page import page as make_page

    parts = [decisions([(n, d["decision"], d.get("question", "")) for n, d in r["decisions"].items()])]
    if "levels" in r:
        parts.append(forest([(n, d.get("estimate"), "") for n, d in r["decisions"].items()], null=0.0,
                            caption="paired differences at 2 A (group estimand)"))
        rows = []
        for name, m in r["levels"].items():
            rows += [(f"{name} ceiling", m["ceiling"].get("per_group"), ""),
                     (f"{name} top pose", m["vina"].get("per_group"), "")]
        parts.append(forest(rows, null=0.5, caption="levels: sampling ceiling and top-pose success"))
    return make_page("The sampling ceiling", "docs/SAMPLING_PLAN.md", [("Decisions", parts)])


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--e32-dir", type=Path, required=True, help="study F's pose lists (rerank_*.jsonl)")
    parser.add_argument("--arms-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    args = parser.parse_args(argv)
    census = pd.read_csv(args.census, dtype={"icode": str})
    e32 = read_records(args.e32_dir, "rerank_*.jsonl")
    arms = {name: read_records(args.arms_dir, f"{name}_*.jsonl") for name in ("e128", "e512")}
    report = build(census, e32, {k: v for k, v in arms.items() if v}, args.n_bootstrap)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "sampling.json").write_text(json.dumps(report, indent=2, default=float))
    text = markdown(report)
    (args.out_dir / "SAMPLING.md").write_text(text)
    (args.out_dir / "report.html").write_text(page(report))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
