"""scripts/triage.py (docs/TRIAGE_PLAN.md) on synthetic data."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def triage():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("triage", ROOT / "scripts" / "triage.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["triage"] = module
    spec.loader.exec_module(module)
    return module


def test_pocket_positions_parse(triage):
    assert triage.parse_positions("12,13,14") == [12, 13, 14]
    assert triage.parse_positions("12;13 ,14") == [12, 13, 14]
    assert triage.parse_positions(float("nan")) == [] and triage.parse_positions("") == []


def _candidates(n=8):
    return pd.DataFrame({"uniprot_id": [f"C{i}" for i in range(n)], "organism_key": "human",
                         "plddt_mean": 90.0, "hull_depth": 10.0, "combined": 0.7,
                         "top_pocket_residues": "10,20,30", "cluster": [f"g{i}" for i in range(n)],
                         "rank": range(1, n + 1), "gene": [f"G{i}" for i in range(n)],
                         "protein_name": "p", "explained": [i % 4 == 0 for i in range(n)]})


def _proteins(n=40):
    return pd.DataFrame({"uniprot_id": [f"P{i}" for i in range(n)], "organism_key": "human",
                         "plddt_mean": [88.0 + (i % 10) for i in range(n)],
                         "hull_depth": [9.0 + (i % 5) for i in range(n)],
                         "combined": [0.1 + 0.02 * i for i in range(n)],
                         "top_pocket_residues": "10,20,30", "cluster": [f"h{i}" for i in range(n)],
                         "annotated": [i < 3 for i in range(n)], "seen": False})


# ------------------------------------------------------------- matching
def test_controls_are_matched_low_ranking_and_distinct(triage):
    cands, prots = _candidates(), _proteins()
    ctrl = triage.matched_controls(cands, prots)
    assert len(ctrl) == len(ctrl["uniprot_id"].unique())          # never reuses a control
    assert set(ctrl["uniprot_id"]) & set(cands["uniprot_id"]) == set()
    cutoff = prots["combined"].quantile(0.5)
    picked = prots.set_index("uniprot_id").loc[ctrl["uniprot_id"]]
    assert (picked["combined"] <= cutoff).all()                   # below the median score
    assert (picked["plddt_mean"] - 90.0).abs().le(triage.PLDDT_TOLERANCE).all()
    assert (picked["hull_depth"] - 10.0).abs().le(triage.DEPTH_TOLERANCE).all()


def test_matching_is_deterministic(triage):
    cands, prots = _candidates(), _proteins()
    a = triage.matched_controls(cands, prots.sample(frac=1, random_state=1))
    b = triage.matched_controls(cands, prots.sample(frac=1, random_state=2))
    assert a.sort_values("matched_to")["uniprot_id"].tolist() == b.sort_values("matched_to")["uniprot_id"].tolist()


# ------------------------------------------------------------ decisions
def _conservation(cand_conserved, ctrl_conserved, positives_conserved, n=8):
    out = {}
    for i in range(n):
        out[f"C{i}"] = {"role": "candidate", "cluster": f"g{i}", "matched_to": None, "homologues": 50,
                        "basic_positions": [10, 20, 30], "conserved_positions": [10, 20, 30],
                        "basic_fraction": {}, "organism_key": "human",
                        "decision": "conserved" if i < cand_conserved else "not conserved"}
        out[f"X{i}"] = {"role": "control", "cluster": f"g{i}", "matched_to": f"C{i}", "homologues": 50,
                        "basic_positions": [], "conserved_positions": [], "basic_fraction": {},
                        "organism_key": "human",
                        "decision": "conserved" if i < ctrl_conserved else "not conserved"}
    for i in range(4):
        out[f"A{i}"] = {"role": "positive", "cluster": f"a{i}", "matched_to": None, "homologues": 50,
                        "basic_positions": [10], "conserved_positions": [], "basic_fraction": {},
                        "organism_key": "human",
                        "decision": "conserved" if i < positives_conserved else "not conserved"}
    return out


def test_enrichment_is_detected_and_paired(triage):
    r, _ = triage.build(_candidates(), _conservation(8, 0, 4), n_bootstrap=300)
    assert r["H1"]["decision"] == "enriched"
    assert r["H1"]["estimate"]["point"] == pytest.approx(1.0)
    assert r["H1"]["pairs"] == 8 and r["holm"]["H1"] <= 1.0


def test_no_enrichment_when_controls_match(triage):
    r, _ = triage.build(_candidates(), _conservation(4, 4, 4), n_bootstrap=300)
    assert r["H1"]["decision"] in ("not enriched", "inconclusive")
    assert r["H1"]["estimate"]["point"] == pytest.approx(0.0)


def test_few_clusters_are_not_evaluable(triage):
    r, _ = triage.build(_candidates(3), _conservation(3, 0, 4, n=3), n_bootstrap=100)
    assert r["H1"]["decision"].startswith("not evaluable")


def test_filter_that_rejects_known_binders_is_declared_uninformative(triage):
    """The guard that matters: a filter failing its positive controls must not triage anything."""
    r, _ = triage.build(_candidates(), _conservation(8, 0, 0), n_bootstrap=200)
    assert r["H2"]["informative"] is False
    assert "not informative" in r["H2"]["note"]
    r2, _ = triage.build(_candidates(), _conservation(8, 0, 4), n_bootstrap=200)
    assert r2["H2"]["informative"] is True


def test_explained_candidates_never_reach_the_shortlist(triage):
    r, merged = triage.build(_candidates(), _conservation(8, 0, 4), n_bootstrap=200)
    assert set(merged.loc[merged["explained"], "triage"]) == {"explained"}
    shortlisted = {c["uniprot_id"] for c in r["H3_shortlist"]}
    assert shortlisted == set(merged.loc[~merged["explained"], "uniprot_id"])
    assert r["H3_counts"]["conserved basic pocket"] == 6


def test_not_evaluable_candidates_are_labelled_not_shortlisted(triage):
    cons = _conservation(8, 0, 4)
    cons["C1"]["decision"] = "not evaluable: 4 homologues (fewer than 10)"
    r, merged = triage.build(_candidates(), cons, n_bootstrap=200)
    assert merged.set_index("uniprot_id").loc["C1", "triage"].startswith("not evaluable")
    assert "C1" not in {c["uniprot_id"] for c in r["H3_shortlist"]}
    assert r["H1"]["pairs"] == 7  # the pair with an unevaluable side drops out


def test_report_writes_its_outputs(triage, tmp_path):
    cands = _candidates()
    cands.to_csv(tmp_path / "c.csv", index=False)
    (tmp_path / "cons.json").write_text(json.dumps(_conservation(8, 0, 4)))
    assert triage.main(["report", "--candidates", str(tmp_path / "c.csv"), "--conservation",
                        str(tmp_path / "cons.json"), "--out-dir", str(tmp_path / "out"),
                        "--n-bootstrap", "100"]) == 0
    text = (tmp_path / "out" / "TRIAGE.md").read_text()
    assert "**H1:** enriched" in text and "H2 (filter calibration)" in text
    assert (tmp_path / "out" / "triaged_candidates.csv").exists()


def test_a_table_without_the_matching_columns_fails_loudly(triage):
    """The run that produced a silent empty control arm must now stop with a named cause."""
    thin = _proteins().drop(columns=["combined", "hull_depth"])
    with pytest.raises(ValueError) as err:
        triage.matched_controls(_candidates(), thin)
    assert "combined" in str(err.value) and "proteins_combined.csv.gz" in str(err.value)


def test_conserve_refuses_to_run_without_a_control_arm(triage, tmp_path, monkeypatch):
    monkeypatch.setattr(triage, "matched_controls", lambda c, p: pd.DataFrame())
    cands = _candidates()
    cands.to_csv(tmp_path / "c.csv", index=False)
    _proteins().to_csv(tmp_path / "p.csv", index=False)
    args = type("A", (), {"candidates": tmp_path / "c.csv", "proteins": tmp_path / "p.csv",
                          "out_dir": tmp_path / "w", "out": tmp_path / "o.json"})()
    with pytest.raises(RuntimeError, match="no matched control"):
        triage.cmd_conserve(args)
