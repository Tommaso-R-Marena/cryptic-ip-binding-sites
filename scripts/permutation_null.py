#!/usr/bin/env python3
"""Null distribution of the permutation control (docs/ANALYSIS_PLAN.md, diagnostic D1).

Reads the out-of-fold predictions of many ``benchmark.py run --permute``
repeats and reports, per task, the pooled ROC-AUC and the mean within-fold
ROC-AUC of each permutation, with D1's leak/chance verdict.

    python scripts/permutation_null.py --predictions-dir preds --output null.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

LEAK_MEAN = 0.52
HIGH = 0.60
MAX_HIGH_FRACTION = 3 / 30


def per_permutation(frame: pd.DataFrame) -> Dict[str, float]:
    cv = frame[frame["split"] == "cv"]
    folds = [roc_auc_score(f["label"], f["score"]) for _, f in cv.groupby("fold") if f["label"].nunique() == 2]
    return {
        "pooled": float(roc_auc_score(cv["label"], cv["score"])),
        "within_fold_mean": float(np.mean(folds)) if folds else float("nan"),
        "folds": len(folds),
    }


def summarise(values: Sequence[float]) -> Dict[str, float]:
    v = np.asarray(values, dtype=float)
    return {"n": int(v.size), "mean": float(v.mean()), "sd": float(v.std(ddof=1)) if v.size > 1 else float("nan"),
            "p2.5": float(np.percentile(v, 2.5)), "p97.5": float(np.percentile(v, 97.5)),
            "above_0.60": int(np.sum(v > HIGH))}


def verdict(pooled: Sequence[float]) -> str:
    v = np.asarray(pooled, dtype=float)
    if v.mean() > LEAK_MEAN or np.mean(v > HIGH) > MAX_HIGH_FRACTION:
        return "leak"
    return "chance"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    rows: List[dict] = []
    for path in sorted(args.predictions_dir.glob("*__permuted.csv.gz")):
        frame = pd.read_csv(path)
        rows.append({"file": path.name, "task": frame["task"].iloc[0], "repeat": int(frame["repeat"].iloc[0]),
                     **per_permutation(frame)})
    skipped = []
    for meta_path in sorted(args.predictions_dir.glob("*__permuted.json")):
        meta = json.loads(meta_path.read_text())
        if "not_evaluable" in meta:
            skipped.append({"task": meta["task"], "repeat": meta["repeat"], "reason": meta["not_evaluable"]})
    if not rows:
        raise SystemExit("no permuted predictions found")
    table = pd.DataFrame(rows).sort_values(["task", "repeat"])
    report = {"rule": {"leak_mean": LEAK_MEAN, "high": HIGH, "max_high_fraction": MAX_HIGH_FRACTION},
              "tasks": {}, "not_evaluable": skipped}
    lines = ["## Permutation-control null (diagnostic D1)", "",
             "| task | permutations | pooled ROC-AUC mean ± sd [2.5, 97.5 pct] | > 0.60 "
             "| within-fold mean ± sd | verdict |",
             "|---|---|---|---|---|---|"]
    for task, part in table.groupby("task"):
        pooled, within = summarise(part["pooled"]), summarise(part["within_fold_mean"].dropna())
        report["tasks"][task] = {"pooled": pooled, "within_fold_mean": within, "verdict": verdict(part["pooled"]),
                                 "per_permutation": part.drop(columns=["task"]).to_dict(orient="records")}
        lines.append(
            f"| {task} | {pooled['n']} | {pooled['mean']:.3f} ± {pooled['sd']:.3f} "
            f"[{pooled['p2.5']:.3f}, {pooled['p97.5']:.3f}] | {pooled['above_0.60']} | "
            f"{within['mean']:.3f} ± {within['sd']:.3f} | **{report['tasks'][task]['verdict']}** |")
    if skipped:
        lines += ["", f"Not evaluable: {len(skipped)} permutation(s): "
                  + "; ".join(f"{s['task']} r{s['repeat']}" for s in skipped)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))
    text = "\n".join(lines) + "\n"
    args.output.with_suffix(".md").write_text(text)
    print(text)
    print(json.dumps({t: {k: v for k, v in r.items() if k != "per_permutation"} for t, r in report["tasks"].items()}))
    for task, r in report["tasks"].items():
        print(task, [round(p["pooled"], 3) for p in r["per_permutation"]])
    return 0


if __name__ == "__main__":
    sys.exit(main())
