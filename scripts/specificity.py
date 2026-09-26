#!/usr/bin/env python3
"""Inositol phosphate versus other polyanion sites (docs/SPECIFICITY_PLAN.md).

``prepare``  one table per variant in the benchmark's format: IP pockets (label 1)
             and other-polyanion pockets (label 0) with joint homology groups and
             the temporal holdout; ``scripts/benchmark.py run`` evaluates them
             unchanged (the label column is ``label_ip_site``)
``decide``   each variant's decision from the comparison report and the
             ten-permutation null, and the gate for the re-ranking
``rerank``   (conditional) P(site) x P(IP | site) over the screened proteomes,
             L3a/L3b and candidates
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from cryptic_ip.benchmark import protocol  # noqa: E402

VARIANTS = ("all", "no_supergroup", "xray25")
MAX_SHARE = 0.40
MIN_HOLDOUT_GROUPS = 5
NOT_LEARNABLE_UPPER = 0.60
PERM_MEAN = 0.52
PERM_HIGH = 0.60
PERM_MAX_HIGH = 1
SEED = 20260930
EXPLAINED = re.compile(
    r"ATP-binding|GTP-binding|Nucleotide-binding|NAD|FAD|FMN|Coenzyme A|sulfotransferase|glycosyltransferase|"
    r"UDP-|diphosphate|pyrophosphate|bisphosphate|phosphoglycer|mitochondrial carrier|Solute carrier family 25|"
    r"ATPase|kinase|3'-phosphoadenos", re.IGNORECASE)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------- prepare
def build_table(ip: pd.DataFrame, other: pd.DataFrame, joint: pd.DataFrame) -> pd.DataFrame:
    """IP-site pockets (1) and other-polyanion-site pockets (0), with joint groups."""
    joint = joint.copy()
    joint["pdb_id"] = joint["pdb_id"].str.upper()
    groups = joint.drop_duplicates("pdb_id").set_index("pdb_id")
    a = ip[ip["label_ip_site"] == 1].copy()
    a["source"] = "inositol"
    b = other[other["label_ip_site"] == 1].copy()
    b["source"] = "transfer"
    shared = set(a["structure_id"].str.upper()) & set(b["structure_id"].str.upper())
    if shared:
        raise SystemExit(f"entries in both tables: {sorted(shared)[:5]}")
    table = pd.concat([a, b], ignore_index=True)
    table["structure_id"] = table["structure_id"].str.upper()
    missing = sorted(set(table["structure_id"]) - set(groups.index))
    if missing:
        raise SystemExit(f"{len(missing)} entries lack a joint group, e.g. {missing[:5]}")
    table["group_sequence"] = table["structure_id"].map(groups["homology_group"]).astype(str)
    table["group_strict"] = table["structure_id"].map(groups["homology_group_strict"]).astype(str)
    table["release_date"] = table["structure_id"].map(groups["release_date"])
    table["label_ip_site"] = (table["source"] == "inositol").astype(int)
    for col in ("label_cryptic_ip_site", "label_burial"):
        table[col] = -1  # not used here; kept so the benchmark's tools see their columns
    if table["row_id"].duplicated().any():
        raise SystemExit("duplicate pocket identifiers")
    entries = table.drop_duplicates("structure_id").rename(columns={"structure_id": "pdb_id"})
    held = protocol.temporal_holdout(entries[["pdb_id", "group_strict", "release_date"]],
                                     group_column="group_strict")
    table["holdout"] = table["structure_id"].isin(held)
    keep = ["row_id", "structure_id", "source", "release_date", "group_sequence", "group_strict", "holdout",
            "label_ip_site", "label_cryptic_ip_site", "label_burial", "rule_score", *protocol.BENCHMARK_FEATURES]
    return table[keep]


def supergroup(table: pd.DataFrame) -> str:
    dev = table[~table["holdout"] & (table["label_ip_site"] == 0)]
    return str(dev["group_strict"].value_counts().idxmax())


def xray25(table: pd.DataFrame, ip_entries: pd.DataFrame) -> pd.DataFrame:
    meta = ip_entries.assign(pdb_id=ip_entries["pdb_id"].str.upper()).drop_duplicates("pdb_id").set_index("pdb_id")
    res = pd.to_numeric(table["structure_id"].map(meta["resolution"]), errors="coerce")
    method = table["structure_id"].map(meta["experimental_method"]).fillna("").str.upper()
    ok_ip = method.str.contains("X-RAY") & (res <= 2.5)
    return table[(table["source"] != "inositol") | ok_ip].reset_index(drop=True)


def summary(table: pd.DataFrame) -> Dict[str, object]:
    out: Dict[str, object] = {"pockets": int(len(table)), "entries": int(table["structure_id"].nunique()),
                              "holdout_entries": int(table.loc[table["holdout"], "structure_id"].nunique())}
    for split, mask in (("development", ~table["holdout"]), ("holdout", table["holdout"])):
        part = table[mask]
        out[split] = {cls: {"pockets": int((part["label_ip_site"] == lab).sum()),
                            "groups": {g: int(part.loc[part["label_ip_site"] == lab, f"group_{g}"].nunique())
                                       for g in ("sequence", "strict")}}
                      for cls, lab in (("ip", 1), ("other", 0))}
    dev = table[~table["holdout"]]
    out["largest_group_share"] = {
        cls: {g: float(dev.loc[dev["label_ip_site"] == lab, f"group_{g}"].value_counts(normalize=True).max())
              if (dev["label_ip_site"] == lab).any() else float("nan") for g in ("sequence", "strict")}
        for cls, lab in (("ip", 1), ("other", 0))}
    # The benchmark's compare reads largest_group_positive_share per task.
    out["tasks"] = {"ip_site": {"largest_group_positive_share": out["largest_group_share"]["ip"]}}
    return out


def cmd_prepare(args: argparse.Namespace) -> int:
    ip = pd.read_csv(args.ip_table, low_memory=False)
    other = pd.read_csv(args.transfer_table, low_memory=False)
    joint = pd.read_csv(args.joint_entries, dtype=str)
    table = build_table(ip, other, joint)
    sg = supergroup(table)
    variants = {"all": table, "no_supergroup": table[table["group_strict"] != sg].reset_index(drop=True),
                "xray25": xray25(table, pd.read_csv(args.ip_entries, dtype=str))}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    meta: Dict[str, object] = {"plan": "docs/SPECIFICITY_PLAN.md", "ip_table_sha256": _sha256(args.ip_table),
                               "transfer_table_sha256": _sha256(args.transfer_table), "supergroup": sg,
                               "supergroup_entries": sorted(table.loc[table["group_strict"] == sg, "structure_id"]
                                                            .unique().tolist()),
                               "variants": {}}
    for name, frame in variants.items():
        frame.to_csv(args.out_dir / f"table_{name}.csv.gz", index=False)
        s = summary(frame)
        (args.out_dir / f"summary_{name}.json").write_text(json.dumps(s, indent=2))
        meta["variants"][name] = s
    (args.out_dir / "prepare.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps({k: v for k, v in meta.items() if k != "supergroup_entries"}, indent=1)[:5000])
    return 0


# ----------------------------------------------------------------- decide
def perm_verdict(null: Mapping[str, object]) -> Optional[str]:
    task = (null.get("tasks") or {}).get("ip_site")
    if not task:
        return None
    pooled = [p["pooled"] for p in task["per_permutation"]]
    if len(pooled) < 10:
        return f"incomplete ({len(pooled)} permutations)"
    v = np.asarray(pooled, dtype=float)
    return "leak" if (v.mean() > PERM_MEAN or int(np.sum(v > PERM_HIGH)) > PERM_MAX_HIGH) else "chance"


def decide_variant(report: Mapping[str, object], null: Mapping[str, object], summ: Mapping[str, object]
                   ) -> Dict[str, object]:
    ev = report.get("evaluations", {})
    seq = ev.get("ip_site/full/sequence/cv")
    strict = ev.get("ip_site/full/strict/cv")
    hold = ev.get("ip_site/full/sequence/holdout")
    verdict = perm_verdict(null)
    keys = ("roc_auc", "pr_auc", "pockets", "positives", "groups", "positive_groups")

    def pick(e):
        return {k: e[k] for k in keys} if e else None

    out: Dict[str, object] = {"permutation_null": verdict, "sequence_cv": pick(seq), "strict_cv": pick(strict),
                              "holdout": pick(hold),
                              "largest_group_share": summ.get("largest_group_share"),
                              "not_evaluable_runs": [k for k in report.get("not_evaluable", {})
                                                     if "__permuted" not in k]}
    share = (summ.get("largest_group_share") or {}).get("ip", {})
    over = [g for g, v in share.items() if v is not None and v > MAX_SHARE]
    hold_groups = (summ.get("holdout") or {})
    hold_ip = hold_groups.get("ip", {}).get("groups", {}).get("sequence", 0)
    hold_other = hold_groups.get("other", {}).get("groups", {}).get("sequence", 0)
    out["holdout_powered"] = bool(hold and hold_ip >= MIN_HOLDOUT_GROUPS and hold_other >= MIN_HOLDOUT_GROUPS)
    if verdict != "chance":
        out["decision"] = f"not evaluable: permutation control {verdict}"
    elif over:
        out["decision"] = f"not evaluable: one group holds > 40 % of IP pockets under {over}"
        out["reason"] = "forty_percent_rule"
    elif out["not_evaluable_runs"] or not seq or not strict:
        out["decision"] = "not evaluable: a required grouping cannot be split"
    else:
        s, t = seq["roc_auc"], strict["roc_auc"]
        hold_ok = (not out["holdout_powered"]) or hold["roc_auc"]["low"] > 0.5
        if s["low"] > 0.5 and t["low"] > 0.5 and hold_ok:
            out["decision"] = "learnable"
        elif s["high"] <= NOT_LEARNABLE_UPPER and t["high"] <= NOT_LEARNABLE_UPPER:
            out["decision"] = "not learnable"
        else:
            out["decision"] = "inconclusive"
    out["p_value_sequence"] = seq["roc_auc"]["p_value"] if seq else None
    return out


def gate(decisions: Mapping[str, Mapping[str, object]]) -> Dict[str, object]:
    s1, s1b, s1r = decisions.get("all", {}), decisions.get("no_supergroup", {}), decisions.get("xray25", {})
    source = None
    if s1.get("decision") == "learnable":
        source = "all"
    elif s1.get("reason") == "forty_percent_rule" and s1b.get("decision") == "learnable":
        source = "no_supergroup"
    r_point = ((s1r.get("sequence_cv") or {}).get("roc_auc") or {}).get("point")
    comparable = r_point is not None and r_point > 0.5
    opened = source is not None and comparable
    return {"open": opened, "variant": source if opened else None, "s1r_point": r_point,
            "reason": ("opened" if opened else "S1r point estimate not above 0.5" if source else
                       "neither S1 nor (S1 not evaluable by the 40 % rule and S1b) is learnable")}


def cmd_decide(args: argparse.Namespace) -> int:
    decisions = {}
    for v in VARIANTS:
        rep = args.reports_dir / v / "benchmark_report.json"
        nul = args.reports_dir / v / "permutation_null.json"
        summ = json.loads((args.prepare_dir / f"summary_{v}.json").read_text())
        if not rep.exists():
            decisions[v] = {"decision": "not evaluable: no comparison report"}
            continue
        decisions[v] = decide_variant(json.loads(rep.read_text()),
                                      json.loads(nul.read_text()) if nul.exists() else {}, summ)
    out = {"plan": "docs/SPECIFICITY_PLAN.md", "decisions": decisions, "gate": gate(decisions),
           "prepare": json.loads((args.prepare_dir / "prepare.json").read_text())}
    out["prepare"].pop("supergroup_entries", None)
    args.output.write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps({v: d.get("decision") for v, d in decisions.items()}), json.dumps(out["gate"]))
    return 0


# ----------------------------------------------------------------- rerank
def lock_specificity(table: pd.DataFrame, n_draws: int = 10, seed: int = 0) -> Dict[str, object]:
    features = list(protocol.ARMS["full"])
    X = table[features].astype(float)
    y = table["label_ip_site"].to_numpy(dtype=int)
    groups = table["group_strict"].astype(str).to_numpy()
    specs = protocol.benchmark_specs()
    candidates = protocol.draw_candidates(n_draws, 7, specs)
    selection = protocol.select_candidate(X, y, groups, candidates, n_inner=3, seed=seed * 1000 + 999, specs=specs)
    model = protocol.fit_candidate(selection.candidate, X, y, seed, specs)
    return {"model": model, "selection": selection, "features": features, "candidate": selection.candidate.as_dict()}


def explained(row: Mapping[str, object]) -> bool:
    text = " ".join(str(row.get(c) or "") for c in ("gene", "protein_name", "keywords", "function"))
    return bool(EXPLAINED.search(text))


def cmd_rerank(args: argparse.Namespace) -> int:
    import joblib

    import learned_screen as ls
    from cryptic_ip.analysis.proteome_stats import known_ip_annotation

    decisions = json.loads(args.decisions.read_text())
    g = decisions["gate"]
    if not g["open"]:
        out = {"ran": False, "reason": g["reason"]}
        args.output.write_text(json.dumps(out, indent=2))
        print(json.dumps(out))
        return 0
    table = pd.read_csv(args.prepare_dir / f"table_{g['variant']}.csv.gz", low_memory=False)
    spec = lock_specificity(table)
    site = joblib.load(args.site_model)
    pockets = ls._read_many(sorted(args.shards_dir.rglob("*_pockets_part*.csv.gz")), low_memory=False)
    catalog = ls._read_many(sorted(args.catalog_dir.rglob("*_catalog.csv")))
    annotations = ls._read_many(sorted(args.catalog_dir.rglob("*_uniprot.tsv")), sep="\t", dtype=str)
    pockets["p_site"] = ls.score_pockets(site, pockets)
    raw = protocol._scores(spec["model"], pockets[spec["features"]].astype(float))
    pockets["p_ip_given_site"] = protocol._calibrate(spec["selection"], raw)
    pockets["combined"] = pockets["p_site"] * pockets["p_ip_given_site"]
    from cryptic_ip.analysis.scorer import PocketScorer

    pockets["composite_score"] = PocketScorer().score_frame(pockets)
    conf = pockets[pockets["plddt_mean"].astype(float).fillna(0) >= ls.MIN_PLDDT]
    prot = conf.groupby("uniprot_id").agg(combined=("combined", "max"), learned_score=("p_site", "max")).reset_index()
    top = conf.sort_values("combined", ascending=False).drop_duplicates("uniprot_id")
    prot = prot.merge(top[["uniprot_id", "pocket_id", "pocket_residues", "p_ip_given_site", "hull_depth",
                           "plddt_mean"]].rename(columns={"pocket_id": "top_pocket", "pocket_residues":
                                                          "top_pocket_residues"}), on="uniprot_id", how="left")
    organism = catalog.drop_duplicates("uniprot_id").set_index("uniprot_id")["organism_key"]
    prot["organism_key"] = prot["uniprot_id"].map(organism)
    known = known_ip_annotation(annotations) if not annotations.empty else pd.Series(dtype=bool)
    prot["annotated"] = prot["uniprot_id"].map(known).fillna(False).astype(bool)
    seen = ls.seen_proteins(args.hits_benchmark) | ls.seen_proteins(args.hits_transfer)
    prot["seen"] = prot["uniprot_id"].isin(seen)
    clusters = ls.read_clusters(args.clusters)
    prot["cluster"] = prot["uniprot_id"].map(clusters).fillna(prot["uniprot_id"])
    unseen = prot[~prot["seen"]].reset_index(drop=True)
    y = unseen["annotated"].to_numpy(dtype=int)
    grp = unseen["cluster"].astype(str).to_numpy()
    comb = protocol.Ranked(y, unseen["combined"].to_numpy(dtype=float))
    l1 = protocol.Ranked(y, unseen["learned_score"].to_numpy(dtype=float))
    l3a = protocol.bootstrap_statistic(comb.roc_auc, grp, null=0.5, n_bootstrap=args.n_bootstrap, seed=SEED)
    l3b = protocol.bootstrap_statistic(lambda w: comb.roc_auc(w) - l1.roc_auc(w), grp, n_bootstrap=args.n_bootstrap,
                                       seed=SEED)
    l1e = protocol.bootstrap_statistic(l1.roc_auc, grp, null=0.5, n_bootstrap=args.n_bootstrap, seed=SEED)
    s1_p = decisions["decisions"][g["variant"]]["p_value_sequence"]
    holm = protocol.holm({"S1": s1_p, "L3a": l3a.p_value, "L3b": l3b.p_value})
    out: Dict[str, object] = {
        "ran": True, "variant": g["variant"], "specificity_model": spec["candidate"],
        "proteins_scored": int(len(prot)), "proteins_seen": int(prot["seen"].sum()),
        "unseen": int(len(unseen)), "annotated_unseen": int(y.sum()),
        "L1_same_proteins": l1e.as_dict(), "L3a_combined": l3a.as_dict(), "L3b_combined_minus_L1": l3b.as_dict(),
        "holm": holm,
        "decisions": {"L3a": "supported" if l3a.low > 0.5 and holm["L3a"] < 0.05 else "not supported",
                      "L3b": "supported" if l3b.low > 0 and holm["L3b"] < 0.05 else "not supported",
                      "S1_holm_p": holm["S1"]},
        "per_organism": {},
    }
    for org, part in unseen.groupby("organism_key"):
        yy = part["annotated"].to_numpy(dtype=int)
        if 0 < yy.sum() < len(yy):
            r = protocol.Ranked(yy, part["combined"].to_numpy(dtype=float))
            out["per_organism"][org] = {"binders": int(yy.sum()), "roc_auc": protocol.bootstrap_statistic(
                r.roc_auc, part["cluster"].astype(str).to_numpy(), null=0.5, n_bootstrap=args.n_bootstrap,
                seed=SEED).as_dict()}
    ann = annotations.drop_duplicates("uniprot_id").set_index("uniprot_id") if not annotations.empty else None
    cands: List[pd.DataFrame] = []
    for org, part in unseen.groupby("organism_key"):
        ranked = part.sort_values("combined", ascending=False).reset_index(drop=True)
        ranked["rank"] = np.arange(1, len(ranked) + 1)
        ranked["precision_at_rank"] = ranked["annotated"].cumsum() / ranked["rank"]
        cands.append(ranked[~ranked["annotated"]].head(25))
    cand = pd.concat(cands, ignore_index=True)
    if ann is not None:
        for col in ("gene", "protein_name", "keywords", "function"):
            if col in ann:
                cand[col] = cand["uniprot_id"].map(ann[col])
    cand["explained"] = cand.apply(explained, axis=1)
    out["candidates_explained"] = int(cand["explained"].sum())
    out["candidates_unexplained"] = int((~cand["explained"]).sum())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cand.to_csv(args.output.with_name("candidates.csv"), index=False)
    prot.to_csv(args.output.with_name("proteins_combined.csv.gz"), index=False)
    out["candidates"] = cand[[c for c in ("organism_key", "rank", "uniprot_id", "gene", "protein_name", "combined",
                                          "learned_score", "p_ip_given_site", "precision_at_rank", "explained")
                              if c in cand]].to_dict(orient="records")
    args.output.write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps({k: v for k, v in out.items() if k != "candidates"}, indent=1, default=float))
    return 0


def markdown(decisions: Mapping[str, object], rerank: Mapping[str, object]) -> str:
    from study_report_page import fmt

    lines = ["## IP versus other polyanion sites (docs/SPECIFICITY_PLAN.md)", "",
             f"Super-group removed in S1b: {decisions['prepare']['supergroup']}.", "",
             "| variant | 10-permutation control | sequence CV ROC-AUC | strict CV ROC-AUC | holdout ROC-AUC | "
             "largest IP group share (seq / strict) | decision |", "|---|---|---|---|---|---|---|"]
    names = {"all": "S1 (all)", "no_supergroup": "S1b (without super-group)", "xray25": "S1r (X-ray ≤ 2.5 Å)"}
    for v, d in decisions["decisions"].items():
        share = (d.get("largest_group_share") or {}).get("ip", {})
        lines.append(f"| {names.get(v, v)} | {d.get('permutation_null')} | "
                     f"{fmt((d.get('sequence_cv') or {}).get('roc_auc'))} | "
                     f"{fmt((d.get('strict_cv') or {}).get('roc_auc'))} | "
                     f"{fmt((d.get('holdout') or {}).get('roc_auc'))}"
                     f"{'' if d.get('holdout_powered') else ' (underpowered)'} | "
                     f"{share.get('sequence', float('nan')):.2f} / {share.get('strict', float('nan')):.2f} | "
                     f"**{d.get('decision')}** |")
    lines += ["", f"Gate for the re-ranking: {decisions['gate']}.", ""]
    if rerank.get("ran"):
        lines += [f"Re-ranking with {rerank['variant']}: {rerank['unseen']} unseen proteins, "
                  f"{rerank['annotated_unseen']} annotated binders.", "",
                  "| | ROC-AUC [95 %] | Holm p | decision |", "|---|---|---|---|",
                  f"| L1 (same proteins) | {fmt(rerank['L1_same_proteins'])} | – | – |",
                  f"| L3a combined | {fmt(rerank['L3a_combined'])} | {rerank['holm']['L3a']:.3g} | "
                  f"**{rerank['decisions']['L3a']}** |",
                  f"| L3b combined − L1 | {fmt(rerank['L3b_combined_minus_L1'])} | {rerank['holm']['L3b']:.3g} | "
                  f"**{rerank['decisions']['L3b']}** |", ""]
        cand = pd.DataFrame(rerank["candidates"])
        if len(cand):
            lines += ["### Candidates (hypotheses, not findings)", "",
                      cand.to_markdown(index=False, floatfmt=".3f"), ""]
    else:
        lines += [f"Re-ranking not run: {rerank.get('reason')}.", ""]
    return "\n".join(lines) + "\n"


def cmd_report(args: argparse.Namespace) -> int:
    decisions = json.loads(args.decisions.read_text())
    rerank = {"ran": False, "reason": "not run"}
    if args.rerank and args.rerank.exists():
        rerank = json.loads(args.rerank.read_text())
    report = {**decisions, "rerank": rerank}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "specificity.json").write_text(json.dumps(report, indent=2, default=float))
    text = markdown(decisions, rerank)
    (args.out_dir / "SPECIFICITY.md").write_text(text)
    from study_report_page import decisions as dec_table, forest, page

    rows = []
    for key, label in (("sequence_cv", "sequence CV"), ("strict_cv", "strict CV")):
        rows += [(f"{v} {label}", (d.get(key) or {}).get("roc_auc"), "") for v, d in decisions["decisions"].items()]
    parts = [dec_table([(v, d.get("decision", ""), "") for v, d in decisions["decisions"].items()]),
             forest(rows, null=0.5, caption="IP site versus other polyanion site: ROC-AUC")]
    if rerank.get("ran"):
        parts.append(forest([("L1 (same proteins)", rerank["L1_same_proteins"], ""),
                             ("L3a combined", rerank["L3a_combined"], "")], null=0.5, caption="protein ranking"))
    (args.out_dir / "report.html").write_text(page("IP versus other polyanion sites", "docs/SPECIFICITY_PLAN.md",
                                                   [("Decisions", parts)]))
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--ip-table", type=Path, required=True)
    p.add_argument("--transfer-table", type=Path, required=True)
    p.add_argument("--joint-entries", type=Path, required=True)
    p.add_argument("--ip-entries", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    d = sub.add_parser("decide")
    d.add_argument("--reports-dir", type=Path, required=True)
    d.add_argument("--prepare-dir", type=Path, required=True)
    d.add_argument("--output", type=Path, required=True)
    r = sub.add_parser("rerank")
    r.add_argument("--decisions", type=Path, required=True)
    r.add_argument("--prepare-dir", type=Path, required=True)
    r.add_argument("--site-model", type=Path, required=True)
    r.add_argument("--shards-dir", type=Path, required=True)
    r.add_argument("--catalog-dir", type=Path, required=True)
    r.add_argument("--hits-benchmark", type=Path, required=True)
    r.add_argument("--hits-transfer", type=Path, required=True)
    r.add_argument("--clusters", type=Path, required=True)
    r.add_argument("--output", type=Path, required=True)
    r.add_argument("--n-bootstrap", type=int, default=2000)
    o = sub.add_parser("report")
    o.add_argument("--decisions", type=Path, required=True)
    o.add_argument("--rerank", type=Path, default=None)
    o.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    return {"prepare": cmd_prepare, "decide": cmd_decide, "rerank": cmd_rerank, "report": cmd_report}[
        args.command](args)


if __name__ == "__main__":
    sys.exit(main())
