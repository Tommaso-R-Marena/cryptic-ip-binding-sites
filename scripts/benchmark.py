#!/usr/bin/env python3
"""Run the pre-registered benchmark (docs/ANALYSIS_PLAN.md).

prepare   one table: every pocket with its descriptors, each task's label, both
          homology groupings, the holdout flag and the rule-based score
run       one evaluation - a task, a descriptor arm, a grouping, a repeat - as
          out-of-fold predictions (plus the locked model's holdout predictions)
compare   every run together: metrics with group-bootstrap intervals, the
          paired hull-depth ablations, the permutation control, and the
          plan's decision for each hypothesis

``run`` is a separate step so each evaluation can be its own CI job; paired
runs share their seeds, so they see identical folds and candidates.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.analysis.labeling import TASKS, task_labels  # noqa: E402
from cryptic_ip.analysis.ml_classifier import expected_calibration_error  # noqa: E402
from cryptic_ip.analysis.scorer import PocketScorer  # noqa: E402
from cryptic_ip.benchmark import protocol  # noqa: E402
from cryptic_ip.utils.json_io import write_json_strict  # noqa: E402

LOGGER = logging.getLogger("benchmark")

GROUPINGS = {"sequence": "homology_group", "strict": "homology_group_strict"}
#: Primary hypotheses: the task each tests (section 6).
HYPOTHESES = {"H1": "cryptic_ip_site", "H2": "burial"}
#: Share of a task's positive pockets one group may hold before five-fold
#: evaluation under that grouping is reported as unsupported (section 3).
MAX_GROUP_POSITIVE_SHARE = 0.40
MIN_HOLDOUT_POSITIVE_GROUPS = 10


# ------------------------------------------------------------------ prepare
def cmd_prepare(args: argparse.Namespace) -> int:
    pockets = pd.read_csv(args.pockets_csv, low_memory=False)
    entries = pd.read_csv(args.entry_csv, dtype={"pdb_id": str})
    entries["pdb_id"] = entries["pdb_id"].str.upper()
    for column in list(GROUPINGS.values()) + ["release_date"]:
        if column not in entries.columns:
            raise SystemExit(f"{args.entry_csv} lacks {column!r}: run homology_groups.py on the builder's output")

    # An entry without a release date cannot be placed in time, so it cannot be
    # assigned to development or holdout: it is excluded, and counted.
    undated = entries["release_date"].isna() | (entries["release_date"].astype(str).str.strip() == "")
    excluded_undated = sorted(entries.loc[undated, "pdb_id"])
    entries = entries[~undated].copy()

    pockets["structure_id"] = pockets["structure_id"].astype(str).str.upper()
    pockets = pockets[~pockets["structure_id"].isin(excluded_undated)].reset_index(drop=True)
    unknown = sorted(set(pockets["structure_id"]) - set(entries["pdb_id"]))
    if unknown:
        # A pocket without an entry has no group; admitting it would leak.
        raise SystemExit(f"{len(unknown)} structures have no entry row, e.g. {unknown[:5]}")
    meta = entries.set_index("pdb_id")
    table = pd.DataFrame(
        {
            "row_id": pockets["structure_id"] + ":" + pockets["pocket_id"].astype(int).astype(str),
            "structure_id": pockets["structure_id"],
            "release_date": pockets["structure_id"].map(meta["release_date"]),
        }
    )
    for name, column in GROUPINGS.items():
        table[f"group_{name}"] = pockets["structure_id"].map(meta[column])
    holdout = protocol.temporal_holdout(entries, group_column=GROUPINGS["strict"])
    table["holdout"] = table["structure_id"].isin(holdout)
    for task in TASKS:
        table[f"label_{task}"] = task_labels(task, pockets)
    table["rule_score"] = PocketScorer().score_frame(pockets)
    for feature in protocol.BENCHMARK_FEATURES:
        table[feature] = pockets[feature] if feature in pockets.columns else np.nan
    if table["row_id"].duplicated().any():
        raise SystemExit("duplicate pocket identifiers")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)

    summary: Dict[str, object] = {
        "n_entries": int(table["structure_id"].nunique()),
        "n_pockets": int(len(table)),
        "holdout_entries": int(len(holdout)),
        "excluded_without_release_date": excluded_undated,
        "groups": {name: int(table[f"group_{name}"].nunique()) for name in GROUPINGS},
        "tasks": {},
    }
    for task in TASKS:
        label = table[f"label_{task}"]
        info: Dict[str, object] = {}
        for split, mask in (("development", ~table["holdout"]), ("holdout", table["holdout"])):
            rows = table[mask & (label >= 0)]
            positives = rows[rows[f"label_{task}"] == 1]
            info[split] = {
                "pockets": int(len(rows)),
                "positives": int(len(positives)),
                "positive_groups": {
                    name: int(positives[f"group_{name}"].nunique()) for name in GROUPINGS
                },
            }
        dev_pos = table[~table["holdout"] & (label == 1)]
        info["largest_group_positive_share"] = {
            name: float(dev_pos[f"group_{name}"].value_counts(normalize=True).max()) if len(dev_pos) else float("nan")
            for name in GROUPINGS
        }
        summary["tasks"][task] = info
    write_json_strict(args.summary_json, summary, indent=2)
    print(json.dumps(summary, indent=2, default=str))
    return 0


# ---------------------------------------------------------------------- run
def _task_rows(table: pd.DataFrame, task: str, holdout: bool) -> pd.DataFrame:
    rows = table[(table[f"label_{task}"] >= 0) & (table["holdout"] == holdout)]
    return rows.reset_index(drop=True)


def cmd_run(args: argparse.Namespace) -> int:
    table = pd.read_csv(args.table, low_memory=False)
    features = list(protocol.ARMS[args.arm])
    group_column = f"group_{args.grouping}"
    dev = _task_rows(table, args.task, holdout=False)
    y = dev[f"label_{args.task}"].to_numpy(dtype=int)
    if args.permute:
        y = np.random.default_rng(12345 + args.repeat).permutation(y)
    groups = dev[group_column].astype(str).to_numpy()
    X = dev[features].astype(float)
    rule = dev["rule_score"].to_numpy(dtype=float)
    LOGGER.info(
        "%s / %s / %s / repeat %d%s: %d pockets, %d positive, %d groups",
        args.task, args.arm, args.grouping, args.repeat, " (permuted)" if args.permute else "",
        len(y), int(y.sum()), len(np.unique(groups)),
    )
    meta: Dict[str, object] = {
        "task": args.task, "arm": args.arm, "grouping": args.grouping, "repeat": args.repeat,
        "permuted": bool(args.permute), "features": features,
        "n_outer": args.n_outer, "n_inner": args.n_inner, "n_draws": args.n_draws,
    }
    stem = f"{args.task}__{args.arm}__{args.grouping}__r{args.repeat}{'__permuted' if args.permute else ''}"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    try:
        oof = protocol.run_cv(
            X, y, groups, rule,
            n_outer=args.n_outer, n_inner=args.n_inner, n_draws=args.n_draws,
            fold_seed=args.repeat, model_seed=args.repeat, progress=LOGGER.info,
        )
    except protocol.SplitError as exc:
        # The grouping cannot support cross-validation for this task: recorded
        # and reported as not evaluable, never silently dropped.
        meta["not_evaluable"] = str(exc)
        write_json_strict(args.out_dir / f"{stem}.json", meta, indent=2)
        print(f"::warning::{stem}: not evaluable - {exc}")
        return 0
    frames = [_prediction_frame(dev, y, groups, rule, oof, "cv", args)]
    if args.holdout:
        hold = _task_rows(table, args.task, holdout=True)
        if len(hold) and len(np.unique(hold[f"label_{args.task}"])) == 2:
            preds, locked = protocol.run_locked(
                X, y, groups, hold[features].astype(float),
                hold[group_column].astype(str).to_numpy(), rule,
                n_inner=args.n_inner, n_draws=args.n_draws,
                fold_seed=args.repeat, model_seed=args.repeat,
            )
            frames.append(
                _prediction_frame(
                    hold, hold[f"label_{args.task}"].to_numpy(dtype=int),
                    hold[group_column].astype(str).to_numpy(),
                    hold["rule_score"].to_numpy(dtype=float), preds, "holdout", args,
                )
            )
            meta["locked_model"] = locked
        else:
            meta["locked_model"] = "holdout lacks both classes for this task"

    pd.concat(frames, ignore_index=True).to_csv(args.out_dir / f"{stem}.csv.gz", index=False)
    write_json_strict(args.out_dir / f"{stem}.json", meta, indent=2)
    return 0


def _prediction_frame(rows, y, groups, rule, preds, split, args) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "row_id": rows["row_id"].to_numpy(),
            "group": groups,
            "label": y,
            "split": split,
            "task": args.task,
            "arm": args.arm,
            "grouping": args.grouping,
            "repeat": args.repeat,
            "permuted": bool(args.permute),
            "fold": preds["fold"].to_numpy(),
            "score": preds["score"].to_numpy(),
            "calibrated": preds["calibrated"].to_numpy(),
            "threshold": preds["threshold"].to_numpy(),
            "family": preds["family"].to_numpy(),
            "rule_score": rule,
            "rule_threshold": preds["rule_threshold"].to_numpy(),
        }
    )


# ------------------------------------------------------------------ compare
def _mcc(y, s, t) -> float:
    from sklearn.metrics import matthews_corrcoef

    return float(matthews_corrcoef(y, (s >= t).astype(int))) if len(np.unique(y)) == 2 else float("nan")


def _metric_panel(frame: pd.DataFrame, n_bootstrap: int) -> Dict[str, object]:
    """Pooled out-of-fold metrics, averaged over repeats, with group-bootstrap intervals."""
    repeats = sorted(frame["repeat"].unique())
    wide = {r: frame[frame["repeat"] == r].set_index("row_id") for r in repeats}
    base = wide[repeats[0]]
    for r in repeats[1:]:
        if not wide[r].index.equals(base.index) or not (wide[r]["label"] == base["label"]).all():
            raise ValueError("repeats are not aligned on the same rows and labels")
    y = base["label"].to_numpy(dtype=int)
    groups = base["group"].to_numpy()
    out: Dict[str, object] = {
        "pockets": int(len(y)), "positives": int(y.sum()), "groups": int(len(np.unique(groups))),
        "positive_groups": int(len(np.unique(groups[y == 1]))), "repeats": len(repeats),
    }
    rule_ranked = protocol.Ranked(y, base["rule_score"].to_numpy())
    for name, fn in protocol.METRICS.items():
        ranked = [protocol.Ranked(y, wide[r]["score"].to_numpy()) for r in repeats]
        null = 0.5 if name == "roc_auc" else float(y.mean())

        def ml(w, fn=fn, ranked=ranked):
            return float(np.mean([fn(rk, w) for rk in ranked]))

        out[name] = protocol.bootstrap_statistic(ml, groups, null=null, n_bootstrap=n_bootstrap).as_dict()
        out[f"rule_{name}"] = protocol.bootstrap_statistic(
            lambda w, fn=fn: fn(rule_ranked, w), groups, null=null, n_bootstrap=n_bootstrap
        ).as_dict()
        out[f"ml_minus_rule_{name}"] = protocol.bootstrap_statistic(
            lambda w, fn=fn, ml=ml: ml(w) - fn(rule_ranked, w), groups, n_bootstrap=n_bootstrap
        ).as_dict()
    out["mcc"] = float(np.mean([_mcc(y, wide[r]["score"].to_numpy(), wide[r]["threshold"].to_numpy()) for r in repeats]))
    out["rule_mcc"] = _mcc(y, base["rule_score"].to_numpy(), base["rule_threshold"].to_numpy())
    calibrated = [wide[r]["calibrated"].to_numpy() for r in repeats]
    if all(np.isfinite(c).all() for c in calibrated):
        out["brier"] = float(np.mean([np.mean((c - y) ** 2) for c in calibrated]))
        out["ece"] = float(np.mean([expected_calibration_error(y, c) for c in calibrated]))
    out["families_chosen"] = frame.drop_duplicates(["repeat", "fold"])["family"].value_counts().to_dict()
    return out


def _paired_difference(full: pd.DataFrame, reduced: pd.DataFrame, metric: str, n_bootstrap: int) -> protocol.Estimate:
    """Mean over repeats of metric(full) - metric(reduced), resampling groups for both arms together."""
    repeats = sorted(set(full["repeat"]) & set(reduced["repeat"]))
    if not repeats:
        raise ValueError("no repeat present in both arms")
    a = {r: full[full["repeat"] == r].set_index("row_id") for r in repeats}
    b = {r: reduced[reduced["repeat"] == r].set_index("row_id").reindex(a[r].index) for r in repeats}
    base = a[repeats[0]]
    for r in repeats:
        if b[r]["score"].isna().any() or not (a[r]["label"] == b[r]["label"]).all():
            raise ValueError("paired arms are not aligned on the same rows and labels")
        if (a[r]["fold"] >= 0).all() and not (a[r]["fold"] == b[r]["fold"]).all():
            raise ValueError("paired arms used different folds")
    y = base["label"].to_numpy(dtype=int)
    fn = protocol.METRICS[metric]
    pairs = [
        (protocol.Ranked(y, a[r]["score"].to_numpy()), protocol.Ranked(y, b[r]["score"].to_numpy()))
        for r in repeats
    ]
    return protocol.bootstrap_statistic(
        lambda w: float(np.mean([fn(ra, w) - fn(rb, w) for ra, rb in pairs])),
        base["group"].to_numpy(), n_bootstrap=n_bootstrap,
    )


def cmd_compare(args: argparse.Namespace) -> int:
    frames = [pd.read_csv(path, low_memory=False) for path in sorted(args.predictions_dir.rglob("*.csv.gz"))]
    not_evaluable = {}
    for path in sorted(args.predictions_dir.rglob("*.json")):
        meta = json.loads(path.read_text())
        if meta.get("not_evaluable"):
            not_evaluable[path.stem] = meta["not_evaluable"]
    if not frames:
        raise SystemExit(f"no predictions under {args.predictions_dir}; not evaluable: {not_evaluable}")
    preds = pd.concat(frames, ignore_index=True)
    summary = json.loads(args.prepare_summary.read_text()) if args.prepare_summary else {}
    report: Dict[str, object] = {
        "data": summary, "evaluations": {}, "hypotheses": {}, "permutation": {},
        "not_evaluable": not_evaluable,
    }
    nb = args.n_bootstrap

    # A smoke run permutes every label, so it exercises the whole analysis
    # without revealing anything about the hypotheses.
    primary_permuted = bool(args.smoke)
    report["smoke"] = primary_permuted

    def select(task, arm, grouping, split="cv", permuted=primary_permuted):
        return preds[
            (preds["task"] == task) & (preds["arm"] == arm) & (preds["grouping"] == grouping)
            & (preds["split"] == split) & (preds["permuted"] == permuted)
        ]

    for (task, arm, grouping, split, permuted), frame in preds.groupby(["task", "arm", "grouping", "split", "permuted"]):
        key = f"{task}/{arm}/{grouping}/{split}{'/permuted' if permuted else ''}"
        LOGGER.info("metrics: %s", key)
        report["evaluations"][key] = _metric_panel(frame, nb)

    def perm_none():
        return preds.iloc[0:0]

    for task in TASKS:
        perm = select(task, "full", "sequence", permuted=True) if not primary_permuted else perm_none()
        if len(perm):
            roc = report["evaluations"][f"{task}/full/sequence/cv/permuted"]["roc_auc"]
            report["permutation"][task] = {
                "roc_auc": roc, "passes": bool(roc["low"] <= 0.5 <= roc["high"]),
            }

    raw_p: Dict[str, float] = {}
    for name, task in HYPOTHESES.items():
        result: Dict[str, object] = {"task": task, "development": {}}
        dev_estimates: Dict[str, protocol.Estimate] = {}
        for grouping in GROUPINGS:
            full, reduced = select(task, "full", grouping), select(task, "no_hull_depth", grouping)
            if full.empty or reduced.empty:
                continue
            roc = _paired_difference(full, reduced, "roc_auc", nb)
            ap = _paired_difference(full, reduced, "pr_auc", nb)
            dev_estimates[grouping] = roc
            result["development"][grouping] = {"roc_auc": roc.as_dict(), "pr_auc": ap.as_dict()}
        if "sequence" in dev_estimates:
            raw_p[name] = dev_estimates["sequence"].p_value
        hold_full = select(task, "full", "sequence", split="holdout")
        hold_reduced = select(task, "no_hull_depth", "sequence", split="holdout")
        holdout_estimate = None
        positive_groups = 0
        if len(hold_full) and len(hold_reduced):
            holdout_estimate = _paired_difference(hold_full, hold_reduced, "roc_auc", nb)
            positive_groups = int(hold_full[hold_full["label"] == 1]["group"].nunique())
            result["holdout"] = {"roc_auc": holdout_estimate.as_dict(), "positive_groups": positive_groups}
        result["_estimates"] = (dev_estimates, holdout_estimate, positive_groups)
        report["hypotheses"][name] = result

    adjusted = protocol.holm(raw_p) if raw_p else {}
    for name, result in report["hypotheses"].items():
        dev_estimates, holdout_estimate, positive_groups = result.pop("_estimates")
        permutation_ok = report["permutation"].get(result["task"], {}).get("passes", True)
        result["holm_p"] = adjusted.get(name, float("nan"))
        largest = summary.get("tasks", {}).get(result["task"], {}).get("largest_group_positive_share", {})
        unsupported = [g for g, share in largest.items() if share > MAX_GROUP_POSITIVE_SHARE]
        if unsupported:
            result["unsupported_groupings"] = unsupported
        missing = [
            stem for stem in not_evaluable
            if stem.startswith(f"{result['task']}__") and "__permuted" not in stem
        ]
        if missing:
            result["not_evaluable_runs"] = missing
        if not permutation_ok:
            result["decision"] = "not evaluable: permutation control failed"
        elif missing:
            result["decision"] = "not evaluable: a required grouping cannot be split (see not_evaluable_runs)"
        elif unsupported:
            result["decision"] = f"not evaluable: one group holds > {MAX_GROUP_POSITIVE_SHARE:.0%} of positives under {unsupported}"
        else:
            result["decision"] = protocol.decide(
                dev_estimates, result["holm_p"], holdout_estimate, positive_groups,
                min_holdout_groups=MIN_HOLDOUT_POSITIVE_GROUPS,
            )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_json_strict(args.out_dir / "benchmark_report.json", report, indent=2)
    (args.out_dir / "REPORT.md").write_text(_markdown(report), encoding="utf-8")
    print(_markdown(report))
    return 0


def _fmt(e: Dict[str, float], digits: int = 3) -> str:
    return f"{e['point']:.{digits}f} [{e['low']:.{digits}f}, {e['high']:.{digits}f}]"


def _markdown(report: Dict[str, object]) -> str:
    lines = ["## Benchmark (docs/ANALYSIS_PLAN.md)\n"]
    if report.get("smoke"):
        lines.append("**Smoke run: every label permuted. Nothing here is a result.**\n")
    data = report.get("data") or {}
    if data:
        lines.append(
            f"{data.get('n_entries')} entries, {data.get('n_pockets')} pockets; "
            f"groups: {data.get('groups')}; holdout entries: {data.get('holdout_entries')}.\n"
        )
    lines.append("### Primary hypotheses: adding hull depth (paired, ROC-AUC difference)\n")
    lines.append("| | task | sequence grouping | strict grouping | holdout | Holm p | decision |")
    lines.append("|---|---|---|---|---|---|---|")
    for name, r in report["hypotheses"].items():
        dev = r["development"]
        hold = r.get("holdout", {})
        hold_text = f"{_fmt(hold['roc_auc'])} ({hold['positive_groups']} positive groups)" if hold else "n/a"
        lines.append(
            f"| {name} | {r['task']} | "
            f"{_fmt(dev['sequence']['roc_auc']) if 'sequence' in dev else 'n/a'} | "
            f"{_fmt(dev['strict']['roc_auc']) if 'strict' in dev else 'n/a'} | "
            f"{hold_text} | "
            f"{r['holm_p']:.3g} | **{r['decision']}** |"
        )
    lines.append("\n### Every evaluation (ROC-AUC and PR-AUC, 95 % group-bootstrap intervals)\n")
    lines.append("| evaluation | pockets (pos) | groups (pos) | ROC-AUC | PR-AUC | rule ROC-AUC | rule PR-AUC | ML − rule ROC-AUC | MCC / rule | families |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for key, e in report["evaluations"].items():
        lines.append(
            f"| {key} | {e['pockets']} ({e['positives']}) | {e['groups']} ({e['positive_groups']}) | "
            f"{_fmt(e['roc_auc'])} | {_fmt(e['pr_auc'])} | {_fmt(e['rule_roc_auc'])} | {_fmt(e['rule_pr_auc'])} | "
            f"{_fmt(e['ml_minus_rule_roc_auc'])} | {e['mcc']:.3f} / {e['rule_mcc']:.3f} | {e['families_chosen']} |"
        )
    if report.get("not_evaluable"):
        lines.append("\n### Not evaluable\n")
        for stem, reason in report["not_evaluable"].items():
            lines.append(f"- {stem}: {reason}")
    if report["permutation"]:
        lines.append("\n### Permutation control (ROC-AUC must include 0.5)\n")
        for task, p in report["permutation"].items():
            lines.append(f"- {task}: {_fmt(p['roc_auc'])} - {'pass' if p['passes'] else '**FAIL**'}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------- cli
def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("prepare")
    p.add_argument("--pockets-csv", type=Path, required=True)
    p.add_argument("--entry-csv", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--summary-json", type=Path, required=True)
    p.set_defaults(func=cmd_prepare)

    p = sub.add_parser("run")
    p.add_argument("--table", type=Path, required=True)
    p.add_argument("--task", choices=TASKS, required=True)
    p.add_argument("--arm", choices=sorted(protocol.ARMS), required=True)
    p.add_argument("--grouping", choices=sorted(GROUPINGS), required=True)
    p.add_argument("--repeat", type=int, default=0)
    p.add_argument("--permute", action="store_true")
    p.add_argument("--holdout", action="store_true", help="also fit the locked model and score the holdout")
    p.add_argument("--n-outer", type=int, default=5)
    p.add_argument("--n-inner", type=int, default=3)
    p.add_argument("--n-draws", type=int, default=10)
    p.add_argument("--out-dir", type=Path, required=True)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("compare")
    p.add_argument("--predictions-dir", type=Path, required=True)
    p.add_argument("--prepare-summary", type=Path, default=None)
    p.add_argument("--n-bootstrap", type=int, default=2000)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--smoke", action="store_true", help="every run was permuted: analyse them as primary")
    p.set_defaults(func=cmd_compare)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    args = parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
