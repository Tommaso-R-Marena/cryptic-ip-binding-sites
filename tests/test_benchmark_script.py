"""End to end: prepare, run the paired arms and the permutation control, compare."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.benchmark import protocol
from scripts import benchmark


def _inputs(tmp_path, n_entries=48, pockets_per_entry=12, seed=0):
    rng = np.random.default_rng(seed)
    rows, entries = [], []
    for e in range(n_entries):
        pdb = f"{e + 1}X{e:02d}"[:4].upper()
        group = f"G:{e // 2:03d}"  # two entries per sequence group
        strict = f"G:{e // 4:03d}"  # four per strict group
        entries.append({
            "pdb_id": pdb, "homology_group": group, "homology_group_strict": strict,
            "release_date": f"{2000 + e // 2}-06-01",
        })
        for p in range(pockets_per_entry):
            site = p < 2
            cryptic = site and (e + p) % 2 == 0
            row = {name: rng.normal() for name in protocol.BENCHMARK_FEATURES}
            row.update({
                "structure_id": pdb, "pocket_id": p + 1,
                "label": 1 if site else 0,
                "matched_burial_class": ("cryptic" if cryptic else "surface") if site else "",
                "plddt_mean": 99.0 if site else 10.0,  # a B-factor leak the protocol must ignore
            })
            row["hull_depth"] = rng.normal() + (3.0 if cryptic else 0.0)
            row["enclosure"] = rng.normal() + (2.0 if site else 0.0)
            rows.append(row)
    pockets, entry_csv = tmp_path / "pockets.csv", tmp_path / "entries.csv"
    pd.DataFrame(rows).to_csv(pockets, index=False)
    pd.DataFrame(entries).to_csv(entry_csv, index=False)
    return pockets, entry_csv


@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("bench")
    pockets, entries = _inputs(tmp)
    table, summary = tmp / "table.csv.gz", tmp / "summary.json"
    assert benchmark.main([
        "prepare", "--pockets-csv", str(pockets), "--entry-csv", str(entries),
        "--output", str(table), "--summary-json", str(summary),
    ]) == 0
    return tmp, table, summary


def test_prepare_derives_labels_groups_and_holdout(prepared):
    tmp, table_path, summary_path = prepared
    table = pd.read_csv(table_path)
    summary = json.loads(summary_path.read_text())
    assert "plddt_mean" not in table.columns
    # Holdout: latest 20 % of the 12 strict groups -> 3 groups, 12 entries.
    assert summary["holdout_entries"] == 12
    held = table[table["holdout"]]
    assert set(held["group_strict"]).isdisjoint(set(table.loc[~table["holdout"], "group_strict"]))
    assert set(table["label_burial"]) <= {-1, 0, 1}
    assert (table.loc[table["label_ip_site"] == 0, "label_burial"] == -1).all()


def test_prepare_refuses_pockets_without_an_entry(tmp_path):
    pockets, entries = _inputs(tmp_path, n_entries=4)
    table = pd.read_csv(entries).iloc[1:]
    table.to_csv(entries, index=False)
    with pytest.raises(SystemExit):
        benchmark.main([
            "prepare", "--pockets-csv", str(pockets), "--entry-csv", str(entries),
            "--output", str(tmp_path / "t.csv"), "--summary-json", str(tmp_path / "s.json"),
        ])


def test_run_and_compare(prepared):
    tmp, table, summary = prepared
    preds = tmp / "preds"
    common = ["--table", str(table), "--n-outer", "3", "--n-inner", "2", "--n-draws", "1", "--out-dir", str(preds)]
    for task in ("cryptic_ip_site", "burial"):
        for arm in ("full", "no_hull_depth"):
            for grouping in ("sequence", "strict"):
                extra = ["--holdout"] if grouping == "sequence" else []
                assert benchmark.main(["run", "--task", task, "--arm", arm, "--grouping", grouping, *common, *extra]) == 0
        assert benchmark.main(["run", "--task", task, "--arm", "full", "--grouping", "sequence", "--permute", *common]) == 0

    out = tmp / "report"
    assert benchmark.main([
        "compare", "--predictions-dir", str(preds), "--prepare-summary", str(summary),
        "--n-bootstrap", "200", "--out-dir", str(out),
    ]) == 0
    report = json.loads((out / "benchmark_report.json").read_text())
    assert set(report["hypotheses"]) == {"H1", "H2"}
    h2 = report["hypotheses"]["H2"]
    # Hull depth carries the burial signal in this synthetic set.
    assert h2["development"]["sequence"]["roc_auc"]["point"] > 0
    assert {"sequence", "strict"} <= set(h2["development"])
    assert "holdout" in h2 and "decision" in h2
    for task in ("cryptic_ip_site", "burial"):
        assert report["permutation"][task]["passes"] in (True, False)
    assert (out / "REPORT.md").read_text().startswith("## Benchmark")


def test_a_task_that_cannot_be_split_is_reported_not_evaluable(prepared, monkeypatch):
    tmp, table, summary = prepared

    def refuse(*_a, **_k):
        raise benchmark.protocol.SplitError("no grouped split gives every training fold both classes")

    monkeypatch.setattr(benchmark.protocol, "run_cv", refuse)
    preds = tmp / "preds_unsplittable"
    common = ["--table", str(table), "--n-outer", "3", "--n-inner", "2", "--n-draws", "1", "--out-dir", str(preds)]
    assert benchmark.main(["run", "--task", "burial", "--arm", "full", "--grouping", "strict", *common]) == 0
    meta = json.loads((preds / "burial__full__strict__r0.json").read_text())
    assert "not_evaluable" in meta
    assert not list(preds.glob("*.csv.gz"))
