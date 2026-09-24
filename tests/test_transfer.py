"""The transfer plan (docs/TRANSFER_PLAN.md): ligand rule, decisions, external test."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from cryptic_ip.database import polyanion_ligands as pl
from scripts import benchmark, build_polyanion_dataset, transfer
from tests.test_benchmark_script import _inputs


def _payload(name, formula):
    return {"chem_comp": {"name": name, "formula": formula}}


@pytest.mark.parametrize("comp_id,name,formula,accepted", [
    ("ATP", "ADENOSINE-5'-TRIPHOSPHATE", "C10 H16 N5 O13 P3", True),
    ("ADP", "ADENOSINE-5'-DIPHOSPHATE", "C10 H15 N5 O10 P2", True),  # 2/27 = 0.074
    ("PPV", "PYROPHOSPHATE", "H4 O7 P2", True),
    ("NAP", "NADP NICOTINAMIDE-ADENINE-DINUCLEOTIDE PHOSPHATE", "C21 H28 N7 O17 P3", False),  # 3/48
    ("IHP", "INOSITOL HEXAKISPHOSPHATE", "C6 H18 O24 P6", False),
    ("AMP", "ADENOSINE MONOPHOSPHATE", "C10 H14 N5 O7 P", False),
    ("XYZ", "SOMETHING", "", False),
])
def test_class_rule(comp_id, name, formula, accepted):
    ligand, reason = pl.classify(comp_id, _payload(name, formula))
    assert (ligand is not None) == accepted, reason


def test_resolution_records_every_decision():
    known = {"ATP": _payload("ATP", "C10 H16 N5 O13 P3"), "NAP": _payload("NADP", "C21 H28 N7 O17 P3")}

    class Client:
        def fetch_chemcomp(self, comp_id):
            return known.get(comp_id)

    ligands, provenance = pl.resolve_polyanion_ligands(Client(), ["atp", "NAP", "ZZZ", "ATP"])
    assert [lig.comp_id for lig in ligands] == ["ATP"]
    assert provenance["decisions"]["ZZZ"] == "unknown identifier"
    assert provenance["decisions"]["NAP"].startswith("phosphorus per heavy atom")


def test_sample_ignores_input_order():
    ids = [f"{i}ABC" for i in range(50)]
    assert build_polyanion_dataset.sample(ids, 10, 1) == build_polyanion_dataset.sample(ids[::-1], 10, 1)
    assert len(build_polyanion_dataset.sample(ids, 100, 1)) == 50


def _est(point, low, high, p=0.001):
    return {"point": point, "low": low, "high": high, "p_value": p}


def _report(share=0.2, dev=(0.02, 0.01, 0.04)):
    hypothesis = lambda task: {  # noqa: E731
        "task": task, "holm_p": 0.004, "decision": "not evaluable: permutation control failed",
        "development": {g: {"roc_auc": _est(*dev)} for g in ("sequence", "strict")},
    }
    return {
        "data": {"tasks": {t: {"largest_group_positive_share": {"sequence": share, "strict": share}}
                           for t in ("cryptic_ip_site", "burial")}},
        "hypotheses": {"H1": hypothesis("cryptic_ip_site"), "H2": hypothesis("burial")},
        "not_evaluable": {},
    }


def test_ten_permutation_control_replaces_the_single_one():
    null = {"tasks": {"cryptic_ip_site": {"verdict": "chance"}, "burial": {"verdict": "leak"}}}
    decisions = transfer.decide_transfer(_report(), null)
    assert decisions["T1"]["decision"] == "supported"
    assert decisions["T2"]["decision"] == "not evaluable: permutation control failed (10 permutations)"


def test_concentrated_positives_stay_not_evaluable():
    null = {"tasks": {"cryptic_ip_site": {"verdict": "chance"}, "burial": {"verdict": "chance"}}}
    decisions = transfer.decide_transfer(_report(share=0.8), null)
    assert all(d["decision"].startswith("not evaluable: one group holds") for d in decisions.values())
    refuted = transfer.decide_transfer(_report(dev=(0.0, -0.004, 0.005), share=0.1), null)
    assert refuted["T1"]["decision"] in ("refuted", "inconclusive")


def test_external_training_excludes_groups_shared_with_the_inositol_table(tmp_path):
    pockets, entries = _inputs(tmp_path, n_entries=40)
    table = tmp_path / "t.csv.gz"
    assert benchmark.main(["prepare", "--pockets-csv", str(pockets), "--entry-csv", str(entries),
                           "--output", str(table), "--summary-json", str(tmp_path / "s.json")]) == 0
    t = pd.read_csv(table)
    ip = t.copy()
    ip["structure_id"] = "9" + ip["structure_id"].str[1:]
    ip_path = tmp_path / "ip.csv.gz"
    ip.to_csv(ip_path, index=False)
    ids = sorted(t["structure_id"].unique())
    joint = pd.DataFrame(
        [{"pdb_id": s, "homology_group_strict": f"J{i // 4}"} for i, s in enumerate(ids)]
        + [{"pdb_id": s, "homology_group_strict": f"K{i}"} for i, s in enumerate(sorted(ip["structure_id"].unique()))]
    )
    joint.loc[len(joint) - 1, "homology_group_strict"] = "J0"  # one inositol entry shares a transfer group
    joint_path = tmp_path / "joint.csv"
    joint.to_csv(joint_path, index=False)
    out = tmp_path / "external.json"
    assert transfer.main(["external", "--transfer-table", str(table), "--ip-table", str(ip_path),
                          "--joint-entries", str(joint_path), "--output", str(out),
                          "--n-draws", "1", "--n-bootstrap", "50"]) == 0
    result = json.loads(out.read_text())
    assert result["training"]["removed_rows_sharing_a_group_with_ip"] > 0
    assert {"cryptic_ip_site", "ip_site"} == set(result["tasks"])
    auc = result["tasks"]["cryptic_ip_site"]["transfer_model_roc_auc"]["point"]
    assert 0.0 <= auc <= 1.0
