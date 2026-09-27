#!/usr/bin/env python3
"""Keep, relax or remove the calibrated screen's hull-depth gate (docs/HULL_GATE_PLAN.md).

    python scripts/hull_gate.py --shards-dir screen/shards --proteins proteins.csv.gz --out-dir results/hull_gate
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from cryptic_ip.analysis.proteome_stats import CALIBRATED_CRITERIA, HitCriteria  # noqa: E402
from cryptic_ip.benchmark import protocol  # noqa: E402

ARMS: Dict[str, Optional[float]] = {"gate10": 10.0, "gate5": 5.0, "none": None}
MARGIN = 0.01
SEED = 20260929
NO_POCKET = -1.0


def arm_criteria(min_hull_depth: Optional[float]) -> HitCriteria:
    """The calibrated criteria with only the hull-depth gate changed (every other gate as calibrated)."""
    base = CALIBRATED_CRITERIA.replace(min_hull_depth=10.0)  # the plan's reference arm, whatever the current value
    return base.replace(min_hull_depth=min_hull_depth)


def protein_scores(pockets: pd.DataFrame, criteria: HitCriteria) -> pd.DataFrame:
    """Per protein: best composite score among pockets passing every non-score gate, and hit status."""
    gates = criteria.gates(pockets)
    non_score = gates.drop(columns=["score"]).all(axis=1)
    frame = pockets.assign(_eligible=non_score, _hit=gates.all(axis=1))
    eligible = frame[frame["_eligible"]]
    best = eligible.groupby("uniprot_id")["composite_score"].max()
    hits = frame.groupby("uniprot_id")["_hit"].any()
    return pd.DataFrame({"rank_score": best, "hit": hits}).reindex(hits.index)


def evaluate(proteins: pd.DataFrame, arms: Mapping[str, pd.DataFrame], n_bootstrap: int) -> Dict[str, object]:
    unseen = proteins[~proteins["seen"].astype(str).str.lower().isin(["true", "1"])].copy()
    unseen["annotated"] = unseen["annotated"].astype(str).str.lower().isin(["true", "1"])
    y = unseen["annotated"].to_numpy(dtype=int)
    groups = unseen["cluster"].fillna(unseen["uniprot_id"]).astype(str).to_numpy()
    scores, hits = {}, {}
    for name, table in arms.items():
        s = unseen["uniprot_id"].map(table["rank_score"]).fillna(NO_POCKET).to_numpy(dtype=float)
        scores[name] = s
        hits[name] = unseen["uniprot_id"].map(table["hit"]).fillna(False).to_numpy(dtype=bool)
    k = int(hits["gate10"].sum())
    top: Dict[str, np.ndarray] = {}
    order_key = unseen["uniprot_id"].to_numpy()
    for name, s in scores.items():
        order = np.lexsort((order_key, -s))  # descending score, ties by accession
        mask = np.zeros(len(s), dtype=bool)
        mask[order[:k]] = True
        top[name] = mask
    ranked = {name: protocol.Ranked(y, s) for name, s in scores.items()}

    def recall(name, w):
        denom = float(np.sum(w * y))
        return float(np.sum(w * y * top[name]) / denom) if denom > 0 else float("nan")

    def est(fn, null=0.0):
        return protocol.bootstrap_statistic(fn, groups, null=null, n_bootstrap=n_bootstrap, seed=SEED).as_dict()

    out: Dict[str, object] = {"unseen_proteins": int(len(unseen)), "annotated_unseen": int(y.sum()),
                              "clusters": int(pd.Series(groups).nunique()), "K": k, "arms": {}}
    for name in arms:
        n_hits = int(hits[name].sum())
        n_true = int((hits[name] & (y == 1)).sum())
        out["arms"][name] = {
            "roc_auc": est(ranked[name].roc_auc, 0.5),
            "recall_at_K": est(lambda w, n=name: recall(n, w)),
            "hits": n_hits, "annotated_hits": n_true,
            "precision_among_hits": n_true / n_hits if n_hits else None,
            "recall_among_hits": n_true / int(y.sum()) if y.sum() else None,
            "per_organism_hits": unseen.assign(_h=hits[name]).groupby("organism_key")["_h"].sum().astype(int).to_dict()
            if "organism_key" in unseen else {},
        }
    out["delta_none"] = est(lambda w: ranked["none"].roc_auc(w) - ranked["gate10"].roc_auc(w))
    out["delta_relax"] = est(lambda w: ranked["gate5"].roc_auc(w) - ranked["gate10"].roc_auc(w))
    out["recall_delta_none"] = est(lambda w: recall("none", w) - recall("gate10", w))
    out["recall_delta_relax"] = est(lambda w: recall("gate5", w) - recall("gate10", w))
    return out


def decide(result: Mapping[str, object], margin: float = MARGIN) -> Dict[str, str]:
    arms = result["arms"]
    r10 = arms["gate10"]["recall_at_K"]["point"]
    if result["delta_none"]["low"] >= -margin and arms["none"]["recall_at_K"]["point"] >= r10:
        return {"decision": "remove", "reason": "removing the gate is non-inferior (margin 0.01) and recall@K "
                "does not fall"}
    if result["delta_relax"]["low"] >= -margin and arms["gate5"]["recall_at_K"]["point"] >= r10:
        return {"decision": "relax", "reason": "a 5 A gate is non-inferior (margin 0.01) and recall@K does not fall"}
    return {"decision": "keep", "reason": "neither removing nor relaxing the gate is shown non-inferior"}


def markdown(result: Mapping[str, object], decision: Mapping[str, str]) -> str:
    from study_report_page import fmt

    lines = ["## The screen's hull-depth gate (docs/HULL_GATE_PLAN.md)", "",
             f"{result['unseen_proteins']} unseen proteins, {result['annotated_unseen']} annotated IP binders, "
             f"{result['clusters']} sequence clusters; K = {result['K']} (hits under the 10 Å gate).", "",
             "| arm | ROC-AUC [95 %] | recall@K [95 %] | hits | annotated hits |", "|---|---|---|---|---|"]
    for name, a in result["arms"].items():
        lines.append(f"| {name} | {fmt(a['roc_auc'])} | {fmt(a['recall_at_K'])} | {a['hits']} | "
                     f"{a['annotated_hits']} |")
    lines += ["", f"- Δ ROC-AUC, no gate − 10 Å: {fmt(result['delta_none'])}",
              f"- Δ ROC-AUC, 5 Å − 10 Å: {fmt(result['delta_relax'])}",
              f"- Δ recall@K, no gate − 10 Å: {fmt(result['recall_delta_none'])}",
              f"- Δ recall@K, 5 Å − 10 Å: {fmt(result['recall_delta_relax'])}", "",
              f"**Decision: {decision['decision']}** ({decision['reason']}).", ""]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards-dir", type=Path, required=True)
    parser.add_argument("--proteins", type=Path, required=True, help="learned-screen proteins.csv.gz")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    args = parser.parse_args(argv)
    from cryptic_ip.analysis.scorer import PocketScorer

    paths = sorted(args.shards_dir.rglob("*_pockets_part*.csv.gz"))
    pockets = pd.concat([pd.read_csv(p, low_memory=False) for p in paths], ignore_index=True)
    pockets["composite_score"] = PocketScorer().score_frame(pockets)
    proteins = pd.read_csv(args.proteins, low_memory=False)
    arms = {name: protein_scores(pockets, arm_criteria(h)) for name, h in ARMS.items()}
    result = evaluate(proteins, arms, args.n_bootstrap)
    decision = decide(result)
    report = {"plan": "docs/HULL_GATE_PLAN.md", "pockets": int(len(pockets)), "shards": len(paths),
              "criteria": {n: dict(arm_criteria(h).__dict__) for n, h in ARMS.items()}, **result, **decision}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "hull_gate.json").write_text(json.dumps(report, indent=2, default=float))
    text = markdown(result, decision)
    (args.out_dir / "HULL_GATE.md").write_text(text)
    from study_report_page import decisions as dec_table, forest, page

    parts = [dec_table([("hull-depth gate", decision["decision"], decision["reason"])]),
             forest([("no gate − 10 Å", result["delta_none"], ""), ("5 Å − 10 Å", result["delta_relax"], "")],
                    null=0.0, caption="paired ROC-AUC difference (non-inferiority margin −0.01)"),
             forest([(n, a["roc_auc"], "") for n, a in result["arms"].items()], null=0.5, caption="ROC-AUC")]
    (args.out_dir / "report.html").write_text(page("The screen's hull-depth gate", "docs/HULL_GATE_PLAN.md",
                                                   [("Decision", parts)]))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
