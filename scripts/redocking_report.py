#!/usr/bin/env python3
"""Merge the redocking arms into per-copy outcomes, estimands and decisions (docs/REDOCKING_PLAN.md).

    python scripts/redocking_report.py --census census.csv --results-dir arms/ --out-dir results/redocking
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
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from cryptic_ip.benchmark import protocol  # noqa: E402
from cryptic_ip.docking import stats  # noqa: E402

SUCCESS = 2.0
SEEDS = (1, 2, 3)
RESOLUTION_BINS = [(0.0, 2.0, "≤ 2.0 Å"), (2.0, 2.5, "2.0-2.5 Å"), (2.5, 3.0, "2.5-3.0 Å"), (3.0, 99.0, "> 3.0 Å")]


def _truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["true", "1", "1.0"])


def read_records(*results_dirs: Path) -> Dict[str, List[dict]]:
    """arm group -> records (one per copy; a docked record beats a failed or not-reached one)."""
    out: Dict[str, Dict[str, dict]] = {}
    paths = [p for d in results_dirs if d is not None and d.exists() for p in sorted(d.rglob("*.jsonl"))]
    for path in paths:
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            group = r.get("arm_group") or path.stem.rsplit("_", 1)[0]
            previous = out.setdefault(group, {}).get(r["copy_key"])
            # A completed record replaces a not-reached one from an earlier dispatch.
            if previous is None or "runs" in r or "runs" not in previous:
                out[group][r["copy_key"]] = r
    return {g: list(v.values()) for g, v in out.items()}


def _run(record: dict, arm: str) -> Optional[dict]:
    for run in record.get("runs", []) or []:
        if run.get("arm") == arm:
            return run
    return None


def _ok(run: Optional[dict], key: str = "top_rmsd", cutoff: float = SUCCESS) -> float:
    if not run or run.get("error") or run.get(key) is None:
        return float("nan")
    return float(float(run[key]) <= cutoff)


def per_copy_table(census: pd.DataFrame, records: Dict[str, List[dict]]) -> pd.DataFrame:
    table = census[_truthy(census["selected"])].copy().set_index("copy_key")
    table["primary_set"] = _truthy(table["primary_set"])
    by_arm = {g: {r["copy_key"]: r for r in rs} for g, rs in records.items()}
    cols: Dict[str, Dict[str, object]] = {}

    def put(key, name, value):
        cols.setdefault(name, {})[key] = value

    for key in table.index:
        prim = by_arm.get("primary", {}).get(key)
        put(key, "primary_error", (prim or {}).get("error", "not run" if prim is None else ""))
        if prim and "runs" in prim:
            runs = [_run(prim, f"vina_s{s}") for s in SEEDS]
            for s, run in zip(SEEDS, runs):
                put(key, f"top_rmsd_s{s}", (run or {}).get("top_rmsd"))
                put(key, f"top_score_s{s}", (run or {}).get("top_score"))
            succ = [_ok(r) for r in runs]
            put(key, "n_seeds", int(np.sum(np.isfinite(succ))))
            put(key, "success", float(np.nanmean(succ)) if np.isfinite(succ).any() else np.nan)
            for cutoff in (1.0, 3.0):
                vals = [_ok(r, cutoff=cutoff) for r in runs]
                put(key, f"success_{cutoff:g}A", float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan)
            best = [_ok(r, "best_rmsd") for r in runs]
            put(key, "best20_success", float(np.nanmean(best)) if np.isfinite(best).any() else np.nan)
            pp = [_ok(r, "top_rmsd_p") for r in runs]
            put(key, "p_success", float(np.nanmean(pp)) if np.isfinite(pp).any() else np.nan)
            rho = [r.get("spearman") for r in runs if r and r.get("spearman") is not None]
            rho = [x for x in rho if np.isfinite(x)]
            put(key, "spearman", float(np.mean(rho)) if rho else np.nan)
            scores = [r.get("top_score") for r in runs if r and r.get("top_score") is not None]
            put(key, "top_score_sd", float(np.std(scores, ddof=1)) if len(scores) > 1 else np.nan)
            put(key, "seed_agreement", float(len(set(s for s in succ if np.isfinite(s))) <= 1))
            ctrl = _run(prim, "crystal_control")
            if ctrl and runs[0] and runs[0].get("top_score") is not None:
                put(key, "crystal_minimised_score", ctrl.get("crystal_minimised_score"))
                put(key, "crystal_minimised_rmsd", ctrl.get("minimised_rmsd"))
                if runs[0].get("top_rmsd") is not None and runs[0]["top_rmsd"] > SUCCESS:
                    kind = ("sampling" if ctrl["crystal_minimised_score"] < runs[0]["top_score"] else "scoring")
                    put(key, "failure_kind", kind)
            put(key, "receptor_his", json.dumps((runs[0] or {}).get("receptor", {}).get("his", {})))
        sec = by_arm.get("secondary", {}).get(key)
        put(key, "secondary_error", (sec or {}).get("error", "not run" if sec is None else ""))
        if sec and "runs" in sec:
            for arm, name in (("vinardo_s1", "vinardo_success"), ("ad4_s1", "ad4_success"),
                              ("vina_deprotonated_s1", "deprotonated_success")):
                run = _run(sec, arm)
                put(key, name, _ok(run))
                if run and run.get("error"):
                    put(key, name.replace("success", "note"), run["error"])
            metal = [_ok(_run(sec, f"vina_metals_s{s}")) for s in SEEDS]
            if np.isfinite(metal).any():
                put(key, "metals_success", float(np.nanmean(metal)))
        pock = by_arm.get("pockets", {}).get(key)
        if pock and "runs" in pock:
            decoy = _run(pock, "decoy_s1")
            put(key, "decoy_top_score", (decoy or {}).get("top_score") if decoy and not decoy.get("error") else None)
            site = _run(pock, "site_finding")
            if site:
                put(key, "site_finding_lands", site.get("lands_in_true_site"))
                put(key, "site_finding_any_positive_pocket", site.get("any_positive_pocket"))
        elif pock:
            put(key, "pockets_error", pock.get("error"))
        af = by_arm.get("alphafold", {}).get(key)
        if af and "runs" in af:
            run = _run(af, "alphafold_s1") or {}
            put(key, "af_error", run.get("error", ""))
            put(key, "af_success", _ok(run))
            for k in ("top_rmsd", "top_score", "site_ca_rmsd", "accession", "identity"):
                put(key, f"af_{k}", run.get(k))
        elif af:
            put(key, "af_error", af.get("error"))
    for name, values in cols.items():
        table[name] = pd.Series(values)
    return table.reset_index()


def strata(frame: pd.DataFrame) -> Dict[str, Dict[str, pd.Series]]:
    res = pd.to_numeric(frame["resolution"], errors="coerce")
    method = frame["experimental_method"].astype(str).str.upper()
    species = frame["species"].where(frame["species"].isin(["InsP3", "InsP4", "InsP5", "InsP6"]), "other")
    out = {
        "burial class": {c: frame["burial_class"] == c for c in ("cryptic", "semi_cryptic", "surface")},
        "metal": {"metal within 3 Å": _truthy(frame["metal"]), "no metal": ~_truthy(frame["metal"])},
        "interface": {"interface": _truthy(frame["interface"]), "single chain": ~_truthy(frame["interface"])},
        "species": {s: species == s for s in ("InsP3", "InsP4", "InsP5", "InsP6", "other")},
        "resolution": {label: (res > lo) & (res <= hi) if lo else res <= hi for lo, hi, label in RESOLUTION_BINS},
        "method": {"X-ray": method.str.contains("X-RAY"), "cryo-EM": method.str.contains("ELECTRON")},
    }
    return out


def build_report(census: pd.DataFrame, records: Dict[str, List[dict]], n_bootstrap: int):
    table = per_copy_table(census, records)
    prim = table[table["primary_set"]].copy()
    groups = prim["homology_group_strict"].astype(str).to_numpy()
    report: Dict[str, object] = {"plan": "docs/REDOCKING_PLAN.md"}

    statuses = census["status"].fillna("unknown").value_counts().to_dict()
    report["accounting"] = {
        "entries_with_copies": int(census["pdb_id"].nunique()),
        "copies_found": int(len(census)),
        "copies_by_status": {str(k): int(v) for k, v in statuses.items()},
        "copies_eligible": int((census["status"] == "eligible").sum()),
        "copies_selected": int(_truthy(census["selected"]).sum()),
        "copies_unselected_eligible": int(((census["status"] == "eligible") & ~_truthy(census["selected"])).sum()),
        "primary_set": int(len(prim)),
        "incomplete_selected": int(len(table) - len(prim)),
        "arms": {},
    }
    for group, rs in records.items():
        errors = pd.Series([r.get("error") for r in rs if r.get("error")]).str.slice(0, 60).value_counts()
        report["accounting"]["arms"][group] = {"records": len(rs), "docked": sum(1 for r in rs if "runs" in r),
                                               "failed_or_not_reached": {str(k): int(v) for k, v in errors.items()}}

    def mean(col, rows=prim, null=0.5):
        g = rows["homology_group_strict"].astype(str).to_numpy()
        return stats.mean_estimates(pd.to_numeric(rows[col], errors="coerce").to_numpy(), g, null=null,
                                    n_bootstrap=n_bootstrap)

    outcomes = {c: mean(c) for c in ("success", "best20_success", "success_1A", "success_3A", "p_success")
                if c in prim}
    if "spearman" in prim:
        outcomes["spearman"] = mean("spearman", null=0.0)
    report["outcomes"] = outcomes
    report["strata"] = {}
    for name, masks in strata(prim).items():
        report["strata"][name] = {label: mean("success", prim[mask.to_numpy()]) for label, mask in masks.items()
                                  if mask.sum() > 0}

    # R1
    r1 = outcomes.get("success", {})
    decisions: Dict[str, dict] = {}
    p_values: Dict[str, float] = {}
    decisions["R1"] = {"question": "protocol reliability (group estimand)", "estimate": r1.get("per_group"),
                       "decision": stats.label(r1.get("per_group"), good="reliable", bad="unreliable")}
    if r1.get("per_group"):
        p_values["R1"] = r1["per_group"]["p_value"]
    # R2
    burial = prim["burial_class"]
    r2 = stats.difference_estimates(pd.to_numeric(prim["success"], errors="coerce"), burial == "cryptic",
                                    burial == "surface", groups, n_bootstrap=n_bootstrap)
    decisions["R2"] = {"question": "success(cryptic) - success(surface)", "estimate": r2.get("per_group"),
                       "groups": r2.get("groups")}
    if not r2.get("evidence"):
        decisions["R2"]["decision"] = f"not evaluable: {r2.get('groups')} strict groups (cryptic, surface); need 5"
    else:
        p_values["R2"] = r2["per_group"]["p_value"]
    # R3
    af = prim.iloc[0:0]
    if "af_success" in prim:
        af = prim[pd.to_numeric(prim["af_success"], errors="coerce").notna()]
    r3 = mean("af_success", af) if len(af) else {}
    decisions["R3"] = {"question": "AlphaFold cross-docking success (group estimand)",
                       "estimate": r3.get("per_group"), "copies": int(len(af)),
                       "decision": stats.label(r3.get("per_group"), good="trustworthy", bad="not trustworthy")}
    if r3.get("per_group"):
        p_values["R3"] = r3["per_group"]["p_value"]
    report["alphafold"] = {"overall": r3, "by_burial_class": {
        c: mean("af_success", af[af["burial_class"] == c]) for c in ("cryptic", "semi_cryptic", "surface")
        if (af["burial_class"] == c).any()}}
    if len(af):
        report["alphafold"]["paired_crystal_minus_af_seed1"] = stats.paired_difference(
            (pd.to_numeric(af["top_rmsd_s1"], errors="coerce") <= SUCCESS).astype(float),
            af["af_success"].astype(float), af["homology_group_strict"].astype(str))
    # R4
    r4 = {}
    if "decoy_top_score" in prim:
        both = prim[pd.to_numeric(prim["decoy_top_score"], errors="coerce").notna()
                    & pd.to_numeric(prim["top_score_s1"], errors="coerce").notna()]
        if len(both):
            labels = np.r_[np.ones(len(both)), np.zeros(len(both))]
            scores = -np.r_[pd.to_numeric(both["top_score_s1"]).to_numpy(),
                            pd.to_numeric(both["decoy_top_score"]).to_numpy()]
            g = np.r_[both["homology_group_strict"].astype(str), both["homology_group_strict"].astype(str)]
            r4 = stats.auc_estimate(labels, scores, g, n_bootstrap=n_bootstrap)
    decisions["R4"] = {"question": "Vina score separates true site from decoy (ROC-AUC)",
                       "estimate": r4.get("roc_auc"), "copies": r4.get("positives")}
    if r4.get("roc_auc"):
        p_values["R4"] = r4["roc_auc"]["p_value"]
    holm = protocol.holm(p_values) if p_values else {}
    for name, d in decisions.items():
        d["holm_p"] = holm.get(name)
    if "decision" not in decisions["R2"]:
        e, hp = r2["per_group"], holm.get("R2", 1.0)
        decisions["R2"]["decision"] = ("buried harder" if e["high"] < 0 and hp < 0.05 else
                                       "buried easier" if e["low"] > 0 and hp < 0.05 else "no difference detected")
    if r4.get("roc_auc"):
        e, hp = r4["roc_auc"], holm.get("R4", 1.0)
        decisions["R4"]["decision"] = ("discriminates" if e["low"] > 0.5 and hp < 0.05 else
                                       "does not discriminate" if e["high"] < 0.6 else "inconclusive")
    else:
        decisions["R4"]["decision"] = "not evaluable: no decoy scores"
    for name, source in (("R1", r1), ("R3", r3), ("R4", r4)):
        if source and source.get("evidence") is False:
            decisions[name]["decision"] = f"not evaluable: {source.get('groups')} strict groups; need 5"
    report["decisions"] = decisions

    # Controls and secondary arms (descriptive).
    controls: Dict[str, object] = {}
    if "failure_kind" in prim:
        kinds = prim["failure_kind"].dropna()
        controls["failure_decomposition"] = {"failed_seed1": int(len(kinds)),
                                             **{str(k): int(v) for k, v in kinds.value_counts().items()}}
    controls["seed_noise"] = {"top_score_sd": mean("top_score_sd", null=0.0) if "top_score_sd" in prim else {},
                              "seed_agreement": mean("seed_agreement") if "seed_agreement" in prim else {}}
    s1 = (pd.to_numeric(prim.get("top_rmsd_s1"), errors="coerce") <= SUCCESS).astype(float)
    s1[pd.to_numeric(prim.get("top_rmsd_s1"), errors="coerce").isna()] = np.nan
    for col in ("vinardo_success", "ad4_success", "deprotonated_success", "metals_success"):
        if col in prim:
            ok = pd.to_numeric(prim[col], errors="coerce")
            keep = ok.notna() & s1.notna()
            if keep.sum():
                controls[col] = {"arm": mean(col, prim[keep]),
                                 "arm_minus_primary_seed1": stats.paired_difference(
                                     ok[keep], s1[keep], prim.loc[keep, "homology_group_strict"].astype(str))}
    if "site_finding_lands" in prim:
        sf = prim[prim["site_finding_lands"].notna()]
        controls["site_finding"] = {
            "lands_in_true_site": mean("site_finding_lands", sf.assign(
                site_finding_lands=sf["site_finding_lands"].astype(float))) if len(sf) else {},
            "any_positive_pocket_fraction": float(sf["site_finding_any_positive_pocket"].astype(float).mean())
            if len(sf) else None}
    report["controls"] = controls
    report["notes"] = [
        "The cryptic class holds few strict groups; any stratum marked 'evidence: false' has fewer than 5 "
        "independent groups and its interval is not evidence.",
        "Docking scores for a -9 polyanion from a scoring function without electrostatics are weak evidence.",
    ]
    return report, table


def markdown(report: Dict[str, object]) -> str:
    from study_report_page import fmt

    lines = ["## Redocking benchmark (docs/REDOCKING_PLAN.md)", ""]
    acc = report["accounting"]
    lines += [f"Copies found: {acc['copies_found']} in {acc['entries_with_copies']} entries; eligible "
              f"{acc['copies_eligible']}; selected {acc['copies_selected']}; primary set {acc['primary_set']}; "
              f"incomplete (flagged stratum) {acc['incomplete_selected']}.", "",
              "| status | copies |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in acc["copies_by_status"].items()]
    lines += ["", "| arm group | records | docked | failed or not reached |", "|---|---|---|---|"]
    lines += [f"| {g} | {a['records']} | {a['docked']} | {a['failed_or_not_reached']} |"
              for g, a in acc["arms"].items()]
    lines += ["", "### Decisions (Holm across R1-R4)", "", "| | question | estimate (group) | Holm p | decision |",
              "|---|---|---|---|---|"]
    for name, d in report["decisions"].items():
        hp = d.get("holm_p")
        lines.append(f"| {name} | {d['question']} | {fmt(d.get('estimate'))} | "
                     f"{'–' if hp is None else f'{hp:.3g}'} | **{d['decision']}** |")
    lines += ["", "### Outcomes (primary set)", "", "| outcome | copies | groups | per copy | per group |",
              "|---|---|---|---|---|"]
    for name, e in report["outcomes"].items():
        lines.append(f"| {name} | {e.get('copies')} | {e.get('groups')} | {fmt(e.get('per_copy'))} | "
                     f"{fmt(e.get('per_group'))} |")
    lines += ["", "### Top-pose success by stratum", "",
              "| stratum | level | copies | groups | per copy | per group | |",
              "|---|---|---|---|---|---|---|"]
    for name, levels in report["strata"].items():
        for level, e in levels.items():
            note = "" if e.get("evidence") else "fewer than 5 groups: not evidence"
            lines.append(f"| {name} | {level} | {e.get('copies')} | {e.get('groups')} | {fmt(e.get('per_copy'))} | "
                         f"{fmt(e.get('per_group'))} | {note} |")
    lines += ["", "### Controls", "", "```json", json.dumps(report["controls"], indent=1, default=str)[:6000], "```"]
    lines += ["", *[f"- {n}" for n in report["notes"]]]
    return "\n".join(lines) + "\n"


def html_page(report: Dict[str, object]) -> str:
    from study_report_page import as_rows, decisions, forest, page

    dec = decisions([(n, d["decision"], d["question"]) for n, d in report["decisions"].items()])
    outcome_rows = as_rows(report["outcomes"])
    parts = [dec, forest([(n, d.get("estimate"), "") for n, d in report["decisions"].items()
                          if n in ("R1", "R3", "R4")], null=0.5, caption="R1, R3 (success) and R4 (ROC-AUC)")]
    strata_parts = [forest(as_rows(levels), null=0.5, caption=f"top-pose success at 2 Å by {name} (group estimand)")
                    for name, levels in report["strata"].items()]
    return page("Redocking inositol phosphates", "docs/REDOCKING_PLAN.md",
                [("Decisions", parts), ("Outcomes", [forest(outcome_rows, null=0.5, caption="group estimand")]),
                 ("Strata", strata_parts)], footer="Generated from redocking.json.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--earlier-dir", type=Path, default=None, help="arms of an earlier dispatch")
    args = parser.parse_args(argv)
    census = pd.read_csv(args.census, dtype={"icode": str})
    records = read_records(args.results_dir, args.earlier_dir)
    report, table = build_report(census, records, args.n_bootstrap)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "redocking.json").write_text(json.dumps(report, indent=2, default=float))
    table.to_csv(args.out_dir / "copies.csv", index=False)
    text = markdown(report)
    (args.out_dir / "REPORT.md").write_text(text)
    (args.out_dir / "report.html").write_text(html_page(report))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
