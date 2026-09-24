#!/usr/bin/env python3
"""The transfer plan's last two steps (docs/TRANSFER_PLAN.md).

``decide``
    T1 and T2 are decided by the benchmark's own rule on the transfer
    dataset's comparison report, except that the permutation control is judged
    on ten permutations (diagnostic D1's lesson) instead of one.

``external``
    T3: a model locked on the transfer dataset - minus every entry sharing a
    joint strict homology group with the inositol phosphate benchmark - scores
    every pocket of that benchmark. Descriptive only.

    python scripts/transfer.py decide --report report/benchmark_report.json \\
        --null null/permutation_null.json --output transfer_decisions.json
    python scripts/transfer.py external --transfer-table t.csv.gz --ip-table ip.csv.gz \\
        --joint-entries joint_grouped.csv --output external.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.benchmark import protocol  # noqa: E402

MAX_GROUP_POSITIVE_SHARE = 0.40  # as scripts/benchmark.py
MIN_HOLDOUT_POSITIVE_GROUPS = 10
HYPOTHESES = {"T1": "cryptic_ip_site", "T2": "burial"}
SEED = 20260925


def _estimate(d: Mapping[str, float]) -> protocol.Estimate:
    return protocol.Estimate(float(d["point"]), float(d["low"]), float(d["high"]), float(d["p_value"]), 0)


def decide_transfer(report: Mapping[str, object], null: Mapping[str, object]) -> Dict[str, dict]:
    """Re-decide each hypothesis with the ten-permutation control."""
    hypotheses = report.get("hypotheses", {})
    summary = report.get("data", {})
    not_evaluable = report.get("not_evaluable", {})
    out: Dict[str, dict] = {}
    for name, task in HYPOTHESES.items():
        source = next((h for h in hypotheses.values() if h.get("task") == task), None)
        verdict = null.get("tasks", {}).get(task, {}).get("verdict")
        entry: Dict[str, object] = {"task": task, "permutation_null": verdict,
                                    "benchmark_decision": source.get("decision") if source else None}
        largest = summary.get("tasks", {}).get(task, {}).get("largest_group_positive_share", {})
        unsupported = [g for g, share in largest.items() if share > MAX_GROUP_POSITIVE_SHARE]
        missing = [s for s in not_evaluable if s.startswith(f"{task}__") and "__permuted" not in s]
        if source is None:
            entry["decision"] = "not evaluable: no comparison was produced"
        elif verdict is None:
            entry["decision"] = "not evaluable: permutation control missing"
        elif verdict != "chance":
            entry["decision"] = "not evaluable: permutation control failed (10 permutations)"
        elif missing:
            entry["decision"] = "not evaluable: a required grouping cannot be split"
        elif unsupported:
            entry["decision"] = (f"not evaluable: one group holds > {MAX_GROUP_POSITIVE_SHARE:.0%} "
                                 f"of positives under {unsupported}")
        else:
            dev = {g: _estimate(v["roc_auc"]) for g, v in source.get("development", {}).items()}
            hold = source.get("holdout")
            entry["decision"] = protocol.decide(
                dev, float(source.get("holm_p", float("nan"))),
                _estimate(hold["roc_auc"]) if hold else None,
                int(hold["positive_groups"]) if hold else 0,
                min_holdout_groups=MIN_HOLDOUT_POSITIVE_GROUPS,
            )
        out[name] = entry
    return out


def cmd_decide(args: argparse.Namespace) -> int:
    report = json.loads(args.report.read_text())
    null = json.loads(args.null.read_text())
    decisions = decide_transfer(report, null)
    args.output.write_text(json.dumps(decisions, indent=2))
    lines = ["## Transfer hypotheses (docs/TRANSFER_PLAN.md)", "",
             "| | task | 10-permutation control | benchmark rule (1 permutation) | decision |", "|---|---|---|---|---|"]
    for name, d in decisions.items():
        lines.append(f"| {name} | {d['task']} | {d['permutation_null']} | {d['benchmark_decision']} | "
                     f"**{d['decision']}** |")
    text = "\n".join(lines) + "\n"
    args.output.with_suffix(".md").write_text(text)
    print(text)
    return 0


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cmd_external(args: argparse.Namespace) -> int:
    transfer = pd.read_csv(args.transfer_table, low_memory=False)
    ip = pd.read_csv(args.ip_table, low_memory=False)
    joint = pd.read_csv(args.joint_entries)
    joint_group = dict(zip(joint["pdb_id"].str.upper(), joint["homology_group_strict"].astype(str)))
    ip_ids = set(ip["structure_id"].str.upper())
    ip_groups = {joint_group[i] for i in ip_ids if i in joint_group}
    if not ip_groups:
        raise SystemExit("no inositol phosphate entry in the joint grouping")

    dev = transfer[~transfer["holdout"].astype(bool) & (transfer["label_cryptic_ip_site"] >= 0)].copy()
    dev["joint_group"] = dev["structure_id"].str.upper().map(joint_group)
    if dev["joint_group"].isna().any():
        raise SystemExit(f"{int(dev['joint_group'].isna().sum())} transfer rows lack a joint group")
    shared = dev["joint_group"].isin(ip_groups)
    train = dev[~shared].reset_index(drop=True)
    features = list(protocol.ARMS["full"])
    y = train["label_cryptic_ip_site"].to_numpy(dtype=int)

    ip = ip.copy()
    ip["joint_group"] = ip["structure_id"].str.upper().map(joint_group).fillna("IP:" + ip["structure_id"])
    protocol.assert_disjoint(train["joint_group"], ip["joint_group"], "transfer training vs inositol table")
    preds, locked = protocol.run_locked(
        train[features].astype(float), y, train["joint_group"].to_numpy(), ip[features].astype(float),
        ip["joint_group"].to_numpy(), train["rule_score"].to_numpy(dtype=float),
        n_inner=3, n_draws=args.n_draws,
    )
    result: Dict[str, object] = {
        "plan": "docs/TRANSFER_PLAN.md (T3, descriptive)",
        "ip_table_sha256": _sha256(args.ip_table), "transfer_table_sha256": _sha256(args.transfer_table),
        "training": {"rows": int(len(train)), "positives": int(y.sum()),
                     "positive_groups": int(train.loc[y == 1, "joint_group"].nunique()),
                     "removed_rows_sharing_a_group_with_ip": int(shared.sum())},
        "locked_model": locked, "tasks": {},
    }
    groups_ip = ip["group_strict"].astype(str).to_numpy()
    for task in ("cryptic_ip_site", "ip_site"):
        keep = (ip[f"label_{task}"] >= 0).to_numpy()
        yt = ip.loc[keep, f"label_{task}"].to_numpy(dtype=int)
        entry: Dict[str, object] = {"pockets": int(keep.sum()), "positives": int(yt.sum()),
                                    "positive_groups": int(pd.Series(groups_ip[keep][yt == 1]).nunique())}
        for label, scores in (("transfer_model", preds["score"].to_numpy()[keep]),
                              ("rule_score", ip.loc[keep, "rule_score"].to_numpy(dtype=float))):
            ranked = protocol.Ranked(yt, scores)
            for metric, fn in protocol.METRICS.items():
                entry[f"{label}_{metric}"] = protocol.bootstrap_statistic(
                    lambda w, r=ranked, f=fn: f(r, w), groups_ip[keep], null=0.5 if metric == "roc_auc" else 0.0,
                    n_bootstrap=args.n_bootstrap, seed=SEED).as_dict()
        result["tasks"][task] = entry
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, default=float))
    lines = ["## T3: a transfer-trained model on the inositol phosphate benchmark (descriptive)", "",
             f"Training: {result['training']}", "",
             "| task | pockets (positives, groups) | transfer model ROC-AUC | PR-AUC | rule ROC-AUC | rule PR-AUC |",
             "|---|---|---|---|---|---|"]
    for task, e in result["tasks"].items():
        f = lambda d: f"{d['point']:.3f} [{d['low']:.3f}, {d['high']:.3f}]"  # noqa: E731
        lines.append(f"| {task} | {e['pockets']} ({e['positives']}, {e['positive_groups']}) | "
                     f"{f(e['transfer_model_roc_auc'])} | {f(e['transfer_model_pr_auc'])} | "
                     f"{f(e['rule_score_roc_auc'])} | {f(e['rule_score_pr_auc'])} |")
    text = "\n".join(lines) + "\n"
    args.output.with_suffix(".md").write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("decide")
    d.add_argument("--report", type=Path, required=True)
    d.add_argument("--null", type=Path, required=True)
    d.add_argument("--output", type=Path, required=True)
    e = sub.add_parser("external")
    e.add_argument("--transfer-table", type=Path, required=True)
    e.add_argument("--ip-table", type=Path, required=True)
    e.add_argument("--joint-entries", type=Path, required=True)
    e.add_argument("--output", type=Path, required=True)
    e.add_argument("--n-draws", type=int, default=10)
    e.add_argument("--n-bootstrap", type=int, default=2000)
    args = parser.parse_args(argv)
    return {"decide": cmd_decide, "external": cmd_external}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
