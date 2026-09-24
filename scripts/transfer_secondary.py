#!/usr/bin/env python3
"""Transfer plan Amendment 1 (docs/TRANSFER_PLAN_AMENDMENT_1.md): drop the dominant group.

Removes every row of the single strict group holding the most
``cryptic_ip_site`` development positives, recomputes the prepare summary on
what remains (exactly as ``benchmark.py prepare`` does), and says whether
T1b/T2b are evaluable under the 40 % rule.

    python scripts/transfer_secondary.py --table table.csv.gz --summary prepare_summary.json \\
        --entries entries_grouped.csv --out-dir secondary
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence

import pandas as pd

MAX_SHARE = 0.40
TASKS = ("ip_site", "cryptic_ip_site", "burial")
GROUPINGS = ("sequence", "strict")


def summarise(table: pd.DataFrame) -> Dict[str, object]:
    """The prepare summary's task block, as scripts/benchmark.py computes it."""
    holdout = table["holdout"].astype(bool)
    tasks: Dict[str, object] = {}
    for task in TASKS:
        label = table[f"label_{task}"]
        info: Dict[str, object] = {}
        for split, mask in (("development", ~holdout), ("holdout", holdout)):
            rows = table[mask & (label >= 0)]
            positives = rows[rows[f"label_{task}"] == 1]
            info[split] = {"pockets": int(len(rows)), "positives": int(len(positives)),
                           "positive_groups": {g: int(positives[f"group_{g}"].nunique()) for g in GROUPINGS}}
        dev_pos = table[~holdout & (label == 1)]
        info["largest_group_positive_share"] = {
            g: float(dev_pos[f"group_{g}"].value_counts(normalize=True).max()) if len(dev_pos) else float("nan")
            for g in GROUPINGS
        }
        tasks[task] = info
    return tasks


def dominant_group(table: pd.DataFrame) -> str:
    dev = table[~table["holdout"].astype(bool) & (table["label_cryptic_ip_site"] == 1)]
    counts = dev["group_strict"].value_counts()
    return str(counts.index[0])


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--entries", type=Path, default=None, help="entries_grouped.csv, to name the removed entries")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    table = pd.read_csv(args.table, low_memory=False)
    summary = json.loads(args.summary.read_text())
    group = dominant_group(table)
    removed = table["group_strict"].astype(str) == group
    kept = table[~removed].reset_index(drop=True)
    entries = sorted(table.loc[removed, "structure_id"].unique())

    new_summary = dict(summary)
    new_summary["n_entries"] = int(kept["structure_id"].nunique())
    new_summary["n_pockets"] = int(len(kept))
    new_summary["groups"] = {g: int(kept[f"group_{g}"].nunique()) for g in GROUPINGS}
    new_summary["tasks"] = summarise(kept)
    shares = {t: new_summary["tasks"][t]["largest_group_positive_share"] for t in ("cryptic_ip_site", "burial")}
    evaluable = all(v <= MAX_SHARE for s in shares.values() for v in s.values())
    new_summary["amendment_1"] = {
        "removed_strict_group": group, "removed_entries": len(entries), "removed_rows": int(removed.sum()),
        "removed_cryptic_dev_positives": int((removed & ~table["holdout"].astype(bool)
                                              & (table["label_cryptic_ip_site"] == 1)).sum()),
        "shares_after": shares, "evaluable": evaluable,
    }
    if args.entries is not None and args.entries.exists():
        meta = pd.read_csv(args.entries, dtype=str)
        cols = [c for c in ("pdb_id", "title", "struct_title", "uniprot_ids", "organism") if c in meta.columns]
        named = meta[meta["pdb_id"].str.upper().isin(entries)][cols]
        args.out_dir.mkdir(parents=True, exist_ok=True)
        named.to_csv(args.out_dir / "removed_entries.csv", index=False)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    kept.to_csv(args.out_dir / "table.csv.gz", index=False)
    (args.out_dir / "prepare_summary.json").write_text(json.dumps(new_summary, indent=2))
    print(json.dumps(new_summary["amendment_1"], indent=2))
    return 0 if evaluable else 3


if __name__ == "__main__":
    sys.exit(main())
