"""scripts/specificity.py (docs/SPECIFICITY_PLAN.md) on synthetic tables."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.benchmark import protocol

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def spec():
    sys.path.insert(0, str(ROOT / "scripts"))
    s = importlib.util.spec_from_file_location("specificity", ROOT / "scripts" / "specificity.py")
    module = importlib.util.module_from_spec(s)
    sys.modules["specificity"] = module
    s.loader.exec_module(module)
    return module


def _table(prefix: str, n_entries: int, label_frac: float, rng, start_year: int) -> pd.DataFrame:
    rows = []
    for e in range(n_entries):
        sid = f"{prefix}{e:03d}"
        for p in range(6):
            row = {"row_id": f"{sid}:{p}", "structure_id": sid, "label_ip_site": int(rng.random() < label_frac),
                   "rule_score": rng.random(), "group_sequence": "local", "group_strict": "local"}
            row.update({f: rng.normal() for f in protocol.BENCHMARK_FEATURES})
            rows.append(row)
    return pd.DataFrame(rows)


def _joint(ip: pd.DataFrame, other: pd.DataFrame, rng) -> pd.DataFrame:
    ids = sorted(set(ip["structure_id"]) | set(other["structure_id"]))
    return pd.DataFrame({"pdb_id": ids,
                         "homology_group": [f"S{i % 23}" for i in range(len(ids))],
                         "homology_group_strict": [f"T{i % 11}" for i in range(len(ids))],
                         "release_date": [f"{2000 + i % 20}-01-01" for i in range(len(ids))]})


def test_table_labels_joint_groups_and_holdout(spec):
    rng = np.random.default_rng(0)
    ip, other = _table("I", 30, 0.3, rng, 2000), _table("T", 40, 0.3, rng, 2000)
    joint = _joint(ip, other, rng)
    table = spec.build_table(ip, other, joint)
    assert set(table["label_ip_site"]) == {0, 1}
    assert (table.loc[table["source"] == "inositol", "label_ip_site"] == 1).all()
    assert (table.loc[table["source"] == "transfer", "label_ip_site"] == 0).all()
    # only polyanion-site pockets: every row came from a label-1 pocket of its own table
    assert len(table) == int((ip["label_ip_site"] == 1).sum() + (other["label_ip_site"] == 1).sum())
    # joint groups replace each table's own groups
    assert not (table["group_strict"] == "local").any()
    protocol.assert_disjoint(table.loc[table["holdout"], "group_strict"], table.loc[~table["holdout"], "group_strict"])


def test_shared_entries_and_missing_groups_stop_the_run(spec):
    rng = np.random.default_rng(1)
    ip, other = _table("I", 5, 0.5, rng, 2000), _table("I", 5, 0.5, rng, 2000)
    with pytest.raises(SystemExit):
        spec.build_table(ip, other, _joint(ip, other, rng))
    ip, other = _table("I", 5, 0.5, rng, 2000), _table("T", 5, 0.5, rng, 2000)
    joint = _joint(ip, other, rng).iloc[1:]
    with pytest.raises(SystemExit):
        spec.build_table(ip, other, joint)


def test_supergroup_and_xray_variant(spec):
    rng = np.random.default_rng(2)
    ip, other = _table("I", 20, 0.5, rng, 2000), _table("T", 30, 0.5, rng, 2000)
    table = spec.build_table(ip, other, _joint(ip, other, rng))
    sg = spec.supergroup(table)
    dev = table[~table["holdout"] & (table["label_ip_site"] == 0)]
    assert sg == dev["group_strict"].value_counts().idxmax()
    entries = pd.DataFrame({"pdb_id": [f"I{i:03d}" for i in range(20)],
                            "resolution": ["2.0" if i % 2 else "3.1" for i in range(20)],
                            "experimental_method": ["X-RAY DIFFRACTION" if i % 4 else "ELECTRON MICROSCOPY"
                                                    for i in range(20)]})
    x = spec.xray25(table, entries)
    kept_ip = set(x.loc[x["source"] == "inositol", "structure_id"])
    assert all(int(s[1:]) % 2 == 1 and int(s[1:]) % 4 != 0 for s in kept_ip)
    assert (x["source"] == "transfer").sum() == (table["source"] == "transfer").sum()


def _est(point, low, high, p=0.01):
    return {"point": point, "low": low, "high": high, "p_value": p}


def _report(seq, strict, hold=None, not_evaluable=None):
    ev = {"ip_site/full/sequence/cv": {"roc_auc": seq, "pr_auc": seq, "pockets": 10, "positives": 5, "groups": 9,
                                       "positive_groups": 6},
          "ip_site/full/strict/cv": {"roc_auc": strict, "pr_auc": strict, "pockets": 10, "positives": 5,
                                     "groups": 8, "positive_groups": 5}}
    if hold:
        ev["ip_site/full/sequence/holdout"] = {"roc_auc": hold, "pr_auc": hold, "pockets": 5, "positives": 2,
                                               "groups": 6, "positive_groups": 5}
    return {"evaluations": ev, "not_evaluable": not_evaluable or {}}


def _null(values):
    return {"tasks": {"ip_site": {"per_permutation": [{"pooled": v} for v in values]}}}


def _summary(ip_share=0.2, hold_ip=6, hold_other=6):
    return {"largest_group_share": {"ip": {"sequence": ip_share, "strict": ip_share}, "other": {}},
            "holdout": {"ip": {"groups": {"sequence": hold_ip}}, "other": {"groups": {"sequence": hold_other}}}}


def test_decision_rules(spec):
    chance = _null([0.5] * 10)
    d = spec.decide_variant(_report(_est(0.8, 0.7, 0.9), _est(0.75, 0.6, 0.85), _est(0.7, 0.55, 0.8)), chance,
                            _summary())
    assert d["decision"] == "learnable"
    d = spec.decide_variant(_report(_est(0.8, 0.7, 0.9), _est(0.75, 0.6, 0.85), _est(0.6, 0.45, 0.7)), chance,
                            _summary())
    assert d["decision"] == "inconclusive"  # powered holdout straddles 0.5
    d = spec.decide_variant(_report(_est(0.8, 0.7, 0.9), _est(0.75, 0.6, 0.85), _est(0.6, 0.45, 0.7)), chance,
                            _summary(hold_ip=3))
    assert d["decision"] == "learnable" and not d["holdout_powered"]
    d = spec.decide_variant(_report(_est(0.52, 0.45, 0.58), _est(0.5, 0.42, 0.6)), chance, _summary())
    assert d["decision"] == "not learnable"
    d = spec.decide_variant(_report(_est(0.8, 0.7, 0.9), _est(0.75, 0.6, 0.85)), chance, _summary(ip_share=0.5))
    assert d["decision"].startswith("not evaluable") and d["reason"] == "forty_percent_rule"
    d = spec.decide_variant(_report(_est(0.8, 0.7, 0.9), _est(0.75, 0.6, 0.85)), _null([0.5] * 8 + [0.65, 0.66]),
                            _summary())
    assert d["decision"] == "not evaluable: permutation control leak"
    d = spec.decide_variant(_report(_est(0.8, 0.7, 0.9), _est(0.75, 0.6, 0.85)), _null([0.5] * 3), _summary())
    assert d["decision"].startswith("not evaluable: permutation control incomplete")


def test_gate(spec):
    learn = {"decision": "learnable", "sequence_cv": {"roc_auc": _est(0.8, 0.7, 0.9)}}
    r_ok = {"sequence_cv": {"roc_auc": _est(0.7, 0.4, 0.8)}}
    r_bad = {"sequence_cv": {"roc_auc": _est(0.45, 0.3, 0.6)}}
    assert spec.gate({"all": learn, "no_supergroup": {}, "xray25": r_ok}) == {
        "open": True, "variant": "all", "s1r_point": 0.7, "reason": "opened"}
    assert not spec.gate({"all": learn, "xray25": r_bad})["open"]
    forty = {"decision": "not evaluable: ...", "reason": "forty_percent_rule"}
    assert spec.gate({"all": forty, "no_supergroup": learn, "xray25": r_ok})["variant"] == "no_supergroup"
    assert not spec.gate({"all": {"decision": "inconclusive"}, "no_supergroup": learn, "xray25": r_ok})["open"]


def test_explained_rule(spec):
    assert spec.explained({"protein_name": "Carbohydrate sulfotransferase 1"})
    assert spec.explained({"keywords": "ATP-binding;Kinase"})
    assert spec.explained({"protein_name": "Solute carrier family 25 member 16"})
    assert not spec.explained({"protein_name": "Arrestin domain-containing protein 2", "keywords": "Membrane"})


def test_benchmark_runner_reads_a_specificity_table(spec, tmp_path):
    """The benchmark's own run step evaluates the table unchanged (smallest settings)."""
    rng = np.random.default_rng(4)
    ip, other = _table("I", 30, 0.6, rng, 2000), _table("T", 40, 0.6, rng, 2000)
    table = spec.build_table(ip, other, _joint(ip, other, rng))
    path = tmp_path / "t.csv.gz"
    table.to_csv(path, index=False)
    s = importlib.util.spec_from_file_location("benchmark_script", ROOT / "scripts" / "benchmark.py")
    bench = importlib.util.module_from_spec(s)
    s.loader.exec_module(bench)
    assert bench.main(["run", "--table", str(path), "--task", "ip_site", "--arm", "full", "--grouping", "strict",
                       "--repeat", "0", "--n-draws", "1", "--n-outer", "3", "--n-inner", "2",
                       "--out-dir", str(tmp_path / "pred")]) == 0
    out = pd.read_csv(tmp_path / "pred" / "ip_site__full__strict__r0.csv.gz")
    assert set(out["label"]) == {0, 1} and len(out) == int((~table["holdout"]).sum())
    json.loads((tmp_path / "pred" / "ip_site__full__strict__r0.json").read_text())
