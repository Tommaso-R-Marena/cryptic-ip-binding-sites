#!/usr/bin/env python3
"""Study L (docs/COEVOLUTION_PLAN.md): does buried IP-site frequency track IP concentration?

Dictyostelium holds roughly ten times the IP6 of human or yeast, and far more IP7 and
IP8. If protein architecture co-evolved with inositol-phosphate availability, the
organism swimming in IPs should bury more of them. Every earlier study took candidates
as a rank cutoff per organism, which fixes the count by construction; a rate needs one
absolute threshold applied to all three proteomes.

    python scripts/coevolution.py table --proteins p.csv.gz --shards-dir screen/shards \
        --fasta seq/screened.fasta --out coev_table.csv.gz
    python scripts/coevolution.py report --table coev_table.csv.gz --out-dir results/coevolution

The primary test is the **matched** comparison (C2), not the raw one (C1): AlphaFold
confidence, protein length and amino-acid composition all differ between these proteomes
and any of them would produce a rate difference on its own.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

LOGGER = logging.getLogger("coevolution")

SEED = 20261005
N_BOOTSTRAP = 2000
MIN_PLDDT = 70.0                 # the founding document's own filter
PRIMARY_FPR = 0.01               # the anchor the decision rests on
FPR_CURVE = (0.001, 0.005, 0.01, 0.02, 0.05)
NO_DIFFERENCE = 0.002            # 0.2 percentage points on a hit rate
FOCUS = "dictyostelium"
BASIC = ("K", "R", "H")

#: Matching bin widths, fixed in the plan.
PLDDT_BIN = 5.0
LOG2_LENGTH_BIN = 0.5
BASIC_BIN = 0.02


# ----------------------------------------------------------------- the table
def basic_fraction(sequence: str) -> float:
    if not sequence:
        return float("nan")
    return sum(1 for c in sequence.upper() if c in BASIC) / len(sequence)


def rule_scores(shards_dir: Path) -> pd.Series:
    """Protein-level rule score: the best composite score over pLDDT-eligible pockets.

    The rule carries no training distribution, so it is the arm that cannot inherit the
    learned model's bias towards the proteomes its training data came from.
    """
    from cryptic_ip.analysis.scorer import PocketScorer
    import learned_screen as ls

    pockets = ls._read_many(sorted(shards_dir.rglob("*_pockets_part*.csv.gz")), low_memory=False)
    if pockets.empty:
        raise ValueError(f"no pocket shards under {shards_dir}; the rule-only arm (C4) needs them")
    pockets["rule_score"] = PocketScorer().score_frame(pockets)
    eligible = pockets[pd.to_numeric(pockets["plddt_mean"], errors="coerce").fillna(0) >= MIN_PLDDT]
    return eligible.groupby("uniprot_id")["rule_score"].max()


def build_table(proteins: pd.DataFrame, sequences: Dict[str, str], rules: pd.Series) -> pd.DataFrame:
    required = ("uniprot_id", "combined", "organism_key", "plddt_mean", "annotated", "seen", "cluster")
    missing = [c for c in required if c not in proteins.columns]
    if missing:
        raise ValueError(f"the proteins table lacks {missing}; study L needs study C's "
                         f"proteins_combined.csv.gz (the specificity-rerank artifact)")
    out = proteins.copy()
    out["sequence_length"] = out["uniprot_id"].map(lambda a: len(sequences.get(str(a), "")) or np.nan)
    out["basic_fraction"] = out["uniprot_id"].map(lambda a: basic_fraction(sequences.get(str(a), "")))
    out["rule_score"] = out["uniprot_id"].map(rules)
    for col in ("combined", "rule_score", "plddt_mean", "hull_depth"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    out["annotated"] = out["annotated"].astype(bool)
    out["seen"] = out["seen"].astype(bool)
    out["eligible"] = out["plddt_mean"].ge(MIN_PLDDT).fillna(False)
    return out


# ------------------------------------------------------------- the threshold
def threshold_at_fpr(frame: pd.DataFrame, score_col: str, fpr: float) -> float:
    """The score whose exceedance rate among pooled non-annotated, unseen proteins is ``fpr``.

    Pooling the anchor across organisms leaves their rates free to differ; it does not
    force them equal, which anchoring within organism would.
    """
    pool = frame[frame["eligible"] & ~frame["annotated"] & ~frame["seen"]][score_col].dropna()
    if pool.empty:
        raise ValueError(f"no pooled background protein carries {score_col!r}; no threshold can be set")
    return float(pool.quantile(1.0 - fpr))


def hits(frame: pd.DataFrame, score_col: str, threshold: float) -> pd.DataFrame:
    part = frame[frame["eligible"] & frame[score_col].notna()].copy()
    part["hit"] = (part[score_col] >= threshold).astype(float)
    return part


# --------------------------------------------------------------- the estimands
def organism_rates(part: pd.DataFrame, n_bootstrap: int) -> Dict[str, object]:
    from cryptic_ip.docking import stats

    out: Dict[str, object] = {}
    for org, sub in part.groupby("organism_key"):
        out[str(org)] = {
            "proteins": int(len(sub)), "hits": int(sub["hit"].sum()),
            "estimate": stats.estimate(lambda w, y=sub["hit"].to_numpy(): stats.per_copy_mean(y, w),
                                       sub["cluster"].astype(str).to_numpy(), null=0.0,
                                       n_bootstrap=n_bootstrap, seed=SEED),
            "clusters": int(sub["cluster"].astype(str).nunique()),
        }
    return out


def _bin_codes(part: pd.DataFrame) -> np.ndarray:
    """Matching cell: pLDDT, log2 length and basic fraction, at the plan's widths."""
    plddt = np.floor(part["plddt_mean"].to_numpy(dtype=float) / PLDDT_BIN)
    length = np.floor(np.log2(np.clip(part["sequence_length"].to_numpy(dtype=float), 1, None)) / LOG2_LENGTH_BIN)
    basic = np.floor(part["basic_fraction"].to_numpy(dtype=float) / BASIC_BIN)
    key = pd.Series([f"{a}|{b}|{c}" for a, b, c in zip(plddt, length, basic)])
    return pd.factorize(key)[0]


def matched_difference(part: pd.DataFrame, focus: str, n_bootstrap: int,
                       seed: int = SEED) -> Dict[str, object]:
    """Bin-size-weighted mean of within-bin (focus − comparator) hit-rate differences."""
    from cryptic_ip.docking import stats

    usable = part.dropna(subset=["plddt_mean", "sequence_length", "basic_fraction"])
    if usable.empty:
        return {"pairs": 0, "bins": 0, "clusters": 0,
                "decision": "not evaluable: no protein carries the matching covariates"}
    codes = _bin_codes(usable)
    n_bins = int(codes.max()) + 1 if len(codes) else 0
    is_focus = (usable["organism_key"].astype(str) == focus).to_numpy(dtype=float)
    is_comp = 1.0 - is_focus
    y = usable["hit"].to_numpy(dtype=float)

    def stat(w: np.ndarray) -> float:
        wf, wc = w * is_focus, w * is_comp
        sf = np.bincount(codes, wf, n_bins)
        sc = np.bincount(codes, wc, n_bins)
        yf = np.bincount(codes, wf * y, n_bins)
        yc = np.bincount(codes, wc * y, n_bins)
        ok = (sf > 0) & (sc > 0)
        if not ok.any():
            return float("nan")
        mass = sf[ok] + sc[ok]
        total = mass.sum()
        if total <= 0:
            return float("nan")
        return float((mass * (yf[ok] / sf[ok] - yc[ok] / sc[ok])).sum() / total)

    groups = usable["cluster"].astype(str).to_numpy()
    shared = _shared_bins(codes, n_bins, is_focus)
    contributing = usable.loc[np.isin(codes, shared)]
    n_clusters = int(contributing["cluster"].astype(str).nunique())
    est = stats.estimate(stat, groups, null=0.0, n_bootstrap=n_bootstrap, seed=seed)
    return {"estimate": est, "bins": int(len(shared)), "clusters": n_clusters,
            "focus_proteins": int(is_focus.sum()), "comparator_proteins": int(is_comp.sum()),
            "contributing_proteins": int(len(contributing)),
            "decision": decide(est, n_clusters)}


def _shared_bins(codes: np.ndarray, n_bins: int, is_focus: np.ndarray) -> np.ndarray:
    f = np.bincount(codes, is_focus, n_bins)
    c = np.bincount(codes, 1.0 - is_focus, n_bins)
    return np.flatnonzero((f > 0) & (c > 0))


def decide(est: Optional[Dict[str, float]], clusters: int) -> str:
    from cryptic_ip.docking import stats

    if not est or not np.isfinite(est.get("point", float("nan"))):
        return "not evaluable: the statistic is undefined (no bin holds both arms)"
    if clusters < stats.MIN_GROUPS:
        return f"not evaluable: {clusters} clusters (fewer than {stats.MIN_GROUPS})"
    if est["low"] > 0:
        return f"higher in {FOCUS}"
    if est["high"] < 0:
        return f"lower in {FOCUS}"
    if est["low"] >= -NO_DIFFERENCE and est["high"] <= NO_DIFFERENCE:
        return "no material difference"
    return "inconclusive"


def raw_difference(part: pd.DataFrame, focus: str, n_bootstrap: int) -> Dict[str, object]:
    """C1: focus minus the mean of the comparator organisms, unmatched."""
    from cryptic_ip.docking import stats

    org = part["organism_key"].astype(str).to_numpy()
    y = part["hit"].to_numpy(dtype=float)
    others = sorted({o for o in org if o != focus})

    def stat(w: np.ndarray) -> float:
        def rate(mask):
            ww = w * mask
            return float((ww * y).sum() / ww.sum()) if ww.sum() > 0 else float("nan")
        comp = [rate(org == o) for o in others]
        comp = [c for c in comp if np.isfinite(c)]
        f = rate(org == focus)
        return f - float(np.mean(comp)) if comp and np.isfinite(f) else float("nan")

    est = stats.estimate(stat, part["cluster"].astype(str).to_numpy(), null=0.0,
                         n_bootstrap=n_bootstrap, seed=SEED)
    return {"estimate": est, "comparators": others,
            "note": "unmatched and confounded by model confidence, length and composition; "
                    "reported so the reader can see how much the matching removes"}


def plddt_deciles(part: pd.DataFrame, focus: str, n_bootstrap: int) -> Dict[str, object]:
    """C3: a difference that lives only in the low-confidence deciles is an artefact."""
    out: Dict[str, object] = {}
    edges = part["plddt_mean"].quantile(np.linspace(0, 1, 11)).to_numpy()
    for i in range(10):
        lo, hi = edges[i], edges[i + 1]
        sub = part[(part["plddt_mean"] >= lo) & (part["plddt_mean"] <= hi)]
        if sub.empty:
            continue
        r = matched_difference(sub, focus, n_bootstrap)
        out[f"{lo:.1f}-{hi:.1f}"] = {"proteins": int(len(sub)), "clusters": r.get("clusters"),
                                     "estimate": r.get("estimate"), "decision": r.get("decision")}
    return out


# ------------------------------------------------------------------- report
def build(table: pd.DataFrame, n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    from cryptic_ip.benchmark import protocol

    report: Dict[str, object] = {"plan": "docs/COEVOLUTION_PLAN.md", "focus": FOCUS,
                                 "min_plddt": MIN_PLDDT, "primary_fpr": PRIMARY_FPR}
    report["eligibility"] = {
        str(org): {"proteins": int(len(sub)), "eligible": int(sub["eligible"].sum()),
                   "excluded_low_plddt": int((~sub["eligible"]).sum())}
        for org, sub in table.groupby("organism_key")}

    threshold = threshold_at_fpr(table, "combined", PRIMARY_FPR)
    part = hits(table, "combined", threshold)
    report["threshold"] = {"score": "combined", "fpr_anchor": PRIMARY_FPR, "value": threshold}
    report["rates"] = organism_rates(part, n_bootstrap)

    report["C1_raw"] = raw_difference(part, FOCUS, n_bootstrap)
    c2 = matched_difference(part, FOCUS, n_bootstrap)
    c2["question"] = "hit rate, dictyostelium - comparators, matched on pLDDT, length and basic fraction"
    report["C2_matched"] = c2
    report["C3_plddt_deciles"] = plddt_deciles(part, FOCUS, n_bootstrap)

    # C4: the same test on the rule score, which carries no training distribution.
    p_values = {}
    if table["rule_score"].notna().any():
        rule_threshold = threshold_at_fpr(table, "rule_score", PRIMARY_FPR)
        rule_part = hits(table, "rule_score", rule_threshold)
        c4 = matched_difference(rule_part, FOCUS, n_bootstrap)
        c4["threshold"] = rule_threshold
        report["C4_rule_only"] = c4
        report["robust"] = _agree(c2.get("estimate"), c4.get("estimate"))
        for key, res in (("C2", c2), ("C4", c4)):
            if res.get("estimate"):
                p_values[key] = res["estimate"]["p_value"]
    else:
        report["C4_rule_only"] = {"decision": "not evaluable: no rule score in the table"}
        report["robust"] = None
        if c2.get("estimate"):
            p_values["C2"] = c2["estimate"]["p_value"]
    report["holm"] = protocol.holm(p_values) if p_values else {}

    # C5: the pyrophosphate prediction - the difference should be larger where pockets are deepest.
    if "hull_depth" in part.columns and part["hull_depth"].notna().any():
        cut = part["hull_depth"].quantile(0.75)
        deep = matched_difference(part[part["hull_depth"] >= cut], FOCUS, n_bootstrap)
        deep["hull_depth_cutoff"] = float(cut)
        report["C5_deep_pockets"] = deep
    else:
        report["C5_deep_pockets"] = {"decision": "not evaluable: no hull depth in the table"}

    report["sensitivity_fpr"] = {}
    for fpr in FPR_CURVE:
        try:
            t = threshold_at_fpr(table, "combined", fpr)
        except ValueError:
            continue
        r = matched_difference(hits(table, "combined", t), FOCUS, n_bootstrap=max(200, n_bootstrap // 5))
        report["sensitivity_fpr"][str(fpr)] = {"threshold": t, "estimate": r.get("estimate"),
                                               "decision": r.get("decision")}

    report["notes"] = [
        "The primary test is C2, the matched comparison. C1 is reported only to show how much "
        "of the raw difference is model confidence, length and composition.",
        "Three organisms give two contrasts and no replication, so even a clean positive is a "
        "correlation across three points, not a demonstration of co-evolution.",
        "If C2 and C4 disagree in direction the primary is not robust to the learned model's "
        "training distribution and no co-evolution claim is made.",
    ]
    return report


def _agree(a: Optional[Dict[str, float]], b: Optional[Dict[str, float]]) -> Optional[bool]:
    if not a or not b:
        return None
    pa, pb = a.get("point"), b.get("point")
    if pa is None or pb is None or not np.isfinite(pa) or not np.isfinite(pb):
        return None
    return bool(np.sign(pa) == np.sign(pb) or pa == 0 or pb == 0)


def markdown(r: Dict[str, object]) -> str:
    def fmt(e):
        if not e or "point" not in e or not np.isfinite(e.get("point", float("nan"))):
            return "-"
        return f"{100 * e['point']:.3f} [{100 * e['low']:.3f}, {100 * e['high']:.3f}] pp"

    c2 = r.get("C2_matched", {})
    c4 = r.get("C4_rule_only", {})
    lines = ["## Does buried IP-site frequency track IP concentration? (docs/COEVOLUTION_PLAN.md)", "",
             f"Threshold: `combined` >= {r.get('threshold', {}).get('value', float('nan')):.6g} "
             f"(the {100 * r.get('primary_fpr', 0):.1f} % pooled false-positive anchor), "
             f"pLDDT floor {r.get('min_plddt')}.", "",
             f"**C2 (primary, matched):** {c2.get('decision')} - {fmt(c2.get('estimate'))} over "
             f"{c2.get('bins')} shared bins and {c2.get('clusters')} clusters.", "",
             "| organism | proteins | hits | hit rate |", "|---|---|---|---|"]
    for org, v in (r.get("rates") or {}).items():
        e = v.get("estimate") or {}
        rate = "-" if "point" not in e else f"{100 * e['point']:.3f} % [{100 * e['low']:.3f}, {100 * e['high']:.3f}]"
        lines.append(f"| {org} | {v.get('proteins')} | {v.get('hits')} | {rate} |")
    c1 = r.get("C1_raw", {})
    lines += ["", f"**C1 (raw, not the test):** {fmt(c1.get('estimate'))}.", "",
              f"**C4 (rule-only):** {c4.get('decision')} - {fmt(c4.get('estimate'))}. "
              f"Directions agree: {r.get('robust')}.", "",
              f"**C5 (deepest quartile):** {r.get('C5_deep_pockets', {}).get('decision')} - "
              f"{fmt(r.get('C5_deep_pockets', {}).get('estimate'))}.", "",
              "| pLDDT decile | proteins | matched difference | decision |", "|---|---|---|---|"]
    for band, v in (r.get("C3_plddt_deciles") or {}).items():
        lines.append(f"| {band} | {v.get('proteins')} | {fmt(v.get('estimate'))} | {v.get('decision')} |")
    lines += ["", *[f"- {n}" for n in r.get("notes", [])]]
    return "\n".join(lines) + "\n"


# -------------------------------------------------------------------- commands
def cmd_table(args: argparse.Namespace) -> int:
    from learned_screen import parse_fasta

    proteins = pd.read_csv(args.proteins)
    sequences = parse_fasta(args.fasta.read_text()) if args.fasta and args.fasta.exists() else {}
    rules = rule_scores(args.shards_dir) if args.shards_dir else pd.Series(dtype=float)
    table = build_table(proteins, sequences, rules)
    table.to_csv(args.out, index=False)
    LOGGER.info("%d proteins, %d with a sequence, %d with a rule score",
                len(table), int(table["sequence_length"].notna().sum()), int(table["rule_score"].notna().sum()))
    print(json.dumps({"proteins": len(table),
                      "organisms": table["organism_key"].value_counts().to_dict()}))
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    table = pd.read_csv(args.table)
    report = build(table, args.n_bootstrap)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "coevolution.json").write_text(json.dumps(report, indent=2, default=str))
    text = markdown(report)
    (args.out_dir / "COEVOLUTION.md").write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    t = sub.add_parser("table")
    t.add_argument("--proteins", type=Path, required=True)
    t.add_argument("--shards-dir", type=Path)
    t.add_argument("--fasta", type=Path)
    t.add_argument("--out", type=Path, required=True)

    r = sub.add_parser("report")
    r.add_argument("--table", type=Path, required=True)
    r.add_argument("--out-dir", type=Path, required=True)
    r.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)

    args = parser.parse_args(argv)
    return cmd_table(args) if args.command == "table" else cmd_report(args)


if __name__ == "__main__":
    sys.exit(main())
