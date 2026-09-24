"""Transfer plan Amendment 1: removing the dominant strict group."""

from __future__ import annotations

import json

import pandas as pd

from scripts import transfer_secondary as ts


def _table():
    rows = []
    # Group D holds 6 of 10 buried development positives; groups A-D otherwise.
    for i, (group, seq_group, label, holdout) in enumerate(
        [("D", "s1", 1, False)] * 6 + [("A", "s2", 1, False)] * 2 + [("B", "s3", 1, False)] * 2
        + [("C", "s4", 0, False)] * 5 + [("D", "s1", 0, True), ("E", "s5", 1, True)]
    ):
        rows.append({"row_id": i, "structure_id": f"{group}{i:03d}", "group_strict": group, "group_sequence": seq_group,
                     "holdout": holdout, "label_ip_site": label, "label_cryptic_ip_site": label,
                     "label_burial": label if label else -1})
    return pd.DataFrame(rows)


def test_dominant_group_is_removed_everywhere_and_summary_recomputed(tmp_path):
    table = _table()
    table_path, summary_path = tmp_path / "t.csv.gz", tmp_path / "s.json"
    table.to_csv(table_path, index=False)
    summary_path.write_text(json.dumps({"n_entries": 17, "tasks": {}}))
    code = ts.main(["--table", str(table_path), "--summary", str(summary_path), "--out-dir", str(tmp_path / "o")])
    kept = pd.read_csv(tmp_path / "o" / "table.csv.gz")
    assert "D" not in set(kept["group_strict"])  # development and holdout rows alike
    summary = json.loads((tmp_path / "o" / "prepare_summary.json").read_text())
    info = summary["amendment_1"]
    assert info["removed_strict_group"] == "D" and info["removed_cryptic_dev_positives"] == 6
    share = summary["tasks"]["cryptic_ip_site"]["largest_group_positive_share"]["strict"]
    assert share == 0.5 and info["evaluable"] is False and code == 3


def test_summary_matches_benchmark_prepare_definition():
    tasks = ts.summarise(_table())
    dev = tasks["cryptic_ip_site"]["development"]
    assert dev["positives"] == 10 and dev["positive_groups"] == {"sequence": 3, "strict": 3}
    assert tasks["cryptic_ip_site"]["largest_group_positive_share"]["strict"] == 0.6
