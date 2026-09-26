"""scripts/orthogonal.py (docs/ORTHOGONAL_PLAN.md) on synthetic data."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def ortho():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("orthogonal", ROOT / "scripts" / "orthogonal.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["orthogonal"] = module
    spec.loader.exec_module(module)
    return module


def _candidates(n=8):
    return pd.DataFrame({"uniprot_id": [f"C{i}" for i in range(n)], "organism_key": "human",
                         "plddt_mean": 90.0, "hull_depth": 10.0, "combined": 0.7,
                         "top_pocket_residues": "10,20,30", "cluster": [f"g{i}" for i in range(n)],
                         "rank": range(1, n + 1), "gene": [f"G{i}" for i in range(n)],
                         "protein_name": "p", "explained": False})


def _proteins(n=40):
    return pd.DataFrame({"uniprot_id": [f"P{i}" for i in range(n)], "organism_key": "human",
                         "plddt_mean": [88.0 + (i % 10) for i in range(n)],
                         "hull_depth": [9.0 + (i % 5) for i in range(n)],
                         "combined": [0.1 + 0.02 * i for i in range(n)],
                         "top_pocket_residues": "10,20,30", "cluster": [f"h{i}" for i in range(n)],
                         "annotated": [i < 4 for i in range(n)], "seen": False})


def _disorder(cand, ctrl, positive=0.1, n=8):
    """Records shaped like the disorder step's output, with a known candidate-control gap."""
    out = {}
    for i in range(n):
        out[f"C{i}"] = {"role": "candidate", "cluster": f"g{i}", "matched_to": None,
                        "pocket_disorder": cand, "disordered": cand > 0.5}
        out[f"X{i}"] = {"role": "control", "cluster": f"g{i}", "matched_to": f"C{i}",
                        "pocket_disorder": ctrl, "disordered": ctrl > 0.5}
    for i in range(6):
        out[f"A{i}"] = {"role": "positive", "cluster": f"a{i}", "matched_to": None,
                        "pocket_disorder": positive, "disordered": positive > 0.5}
    return out


# ------------------------------------------------------------------- arms
def test_the_arms_are_study_h_s_arms(ortho):
    frame = ortho.roles_frame(_candidates(), _proteins())
    assert set(frame["role"]) == {"candidate", "control", "positive"}
    assert (frame[frame["role"] == "control"]["matched_to"].notna()).all()
    assert frame["uniprot_id"].is_unique       # a binder drawn as a control counts once, as a control


def test_a_thin_protein_table_fails_loudly(ortho):
    with pytest.raises(ValueError, match="proteins_combined"):
        ortho.roles_frame(_candidates(), _proteins().drop(columns=["combined"]))


def test_no_control_arm_refuses_to_run(ortho, monkeypatch):
    monkeypatch.setattr(ortho, "matched_controls", lambda c, p: pd.DataFrame(), raising=False)
    import triage
    monkeypatch.setattr(triage, "matched_controls", lambda c, p: pd.DataFrame())
    with pytest.raises(RuntimeError, match="no matched control"):
        ortho.roles_frame(_candidates(), _proteins())


# --------------------------------------------------------------- disorder
def test_only_basic_pocket_positions_are_kept(ortho):
    seq = "AKRDHEG"          # K2 R3 H5 are basic; D4 and the rest are not
    assert ortho.basic_pocket_positions(seq, [1, 2, 3, 4, 5, 99]) == [2, 3, 5]


def test_pocket_disorder_averages_only_the_pocket(ortho):
    pytest.importorskip("metapredict")
    ordered = "MKVLAAGIVGLNLGGSLAKELVKRGHEVTVYDVNQEAVDHLVAQGATAVASPAEAAKDADLVILAVPAEAVEAV"
    out = ortho.pocket_disorder(ordered, [5, 6, 7])
    assert 0.0 <= out["pocket_disorder"] <= 1.0 and out["residues_scored"] == 3
    assert out["length"] == len(ordered)
    assert out["disordered"] is (out["pocket_disorder"] > ortho.DISORDER_CUT)


def test_a_pocket_outside_the_sequence_is_an_error_not_a_score(ortho):
    pytest.importorskip("metapredict")
    assert "error" in ortho.pocket_disorder("MKVLAAG", [900])


# ----------------------------------------------------------- motif support
def test_motif_support_counts_positions_inside_a_match(ortho, tmp_path):
    table = tmp_path / "all_occurrences.tsv"
    table.write_text("sequenceID\treference_kmer\tmatch\tstart\tend\n"
                     "Q1\tKRK\tKRK\t10\t12\n"
                     "Q1\tKGR\tKGR\t40\t42\n"
                     "OTHER\tKRK\tKRK\t20\t22\n")
    out = ortho.motif_support(table, "Q1", [10, 20, 41])
    assert out["motif_support"] == pytest.approx(2 / 3)     # 10 and 41 are covered, 20 is not
    assert out["covered_positions"] == [10, 41] and out["motifs"] == 2


def test_motif_support_ignores_other_sequences(ortho, tmp_path):
    table = tmp_path / "all_occurrences.tsv"
    table.write_text("sequenceID\treference_kmer\tmatch\tstart\tend\nOTHER\tKRK\tKRK\t10\t12\n")
    assert ortho.motif_support(table, "Q1", [10])["motif_support"] == 0.0


# ---------------------------------------------------------------- decisions
def test_more_ordered_candidates_are_detected(ortho):
    r, _ = ortho.build(_candidates(), _disorder(0.10, 0.60), {}, {}, n_bootstrap=300)
    assert r["I1"]["decision"] == "candidates more ordered"
    assert r["I1"]["estimate"]["point"] == pytest.approx(-0.50)
    assert r["I1"]["pairs"] == 8 and r["holm"]["I1"] <= 1.0


def test_an_identical_arm_is_called_no_difference(ortho):
    r, _ = ortho.build(_candidates(), _disorder(0.30, 0.30), {}, {}, n_bootstrap=300)
    assert r["I1"]["decision"] == "no difference"
    assert r["I1"]["estimate"]["point"] == pytest.approx(0.0)


def test_few_clusters_are_not_evaluable(ortho):
    r, _ = ortho.build(_candidates(3), _disorder(0.1, 0.6, n=3), {}, {}, n_bootstrap=100)
    assert r["I1"]["decision"].startswith("not evaluable")


def test_no_pair_at_all_is_not_evaluable(ortho):
    r, _ = ortho.build(_candidates(), {}, {}, {}, n_bootstrap=100)
    assert r["I1"]["decision"].startswith("not evaluable")


def test_a_filter_that_flags_known_binders_demotes_nobody(ortho):
    """The guard that matters: a QC failing its positive controls must not triage anything."""
    r, _ = ortho.build(_candidates(), _disorder(0.9, 0.2, positive=0.9), {}, {}, n_bootstrap=200)
    assert r["I2"]["informative"] is False and r["I2"]["demoted"] == []
    assert "not informative" in r["I2"]["note"]


def test_an_informative_filter_demotes_the_disordered_candidates(ortho):
    r, _ = ortho.build(_candidates(), _disorder(0.9, 0.2, positive=0.1), {}, {}, n_bootstrap=200)
    assert r["I2"]["informative"] is True
    assert r["I2"]["demoted"] == [f"C{i}" for i in range(8)]


def test_i4_counts_only_the_unevaluable_and_does_not_re_run_conservation(ortho):
    dive = {"C1": {"remote_homologues": 7}, "C2": {"remote_homologues": 0},
            "C3": {"error": "RuntimeError: shark-dive produced no table"}}
    r, merged = ortho.build(_candidates(), _disorder(0.2, 0.2), {}, dive, n_bootstrap=100)
    assert r["I4"]["queries"] == 3 and r["I4"]["with_remote_homologues"] == 1
    assert "not re-run" in r["I4"]["note"]
    assert merged.set_index("uniprot_id").loc["C1", "remote_homologues"] == 7


def test_unevaluable_candidates_are_the_dive_queries(ortho):
    cons = {"C1": {"role": "candidate", "decision": "not evaluable: 4 homologues (fewer than 10)"},
            "C2": {"role": "candidate", "decision": "conserved"},
            "X1": {"role": "control", "decision": "not evaluable: 0 homologues (fewer than 10)"}}
    assert ortho.unevaluable_candidates(cons) == ["C1"]


# ------------------------------------------------------------------ outputs
def test_merge_joins_the_shards(ortho, tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "s0.json").write_text(json.dumps({"C0": {"motif_support": 1.0}}))
    (tmp_path / "a" / "s1.json").write_text(json.dumps({"C1": {"motif_support": 0.5}}))
    assert ortho.main(["merge", "--inputs", str(tmp_path / "a"),
                       "--out", str(tmp_path / "m.json")]) == 0
    assert set(json.loads((tmp_path / "m.json").read_text())) == {"C0", "C1"}


def test_report_writes_its_outputs(ortho, tmp_path):
    cands = _candidates()
    cands.to_csv(tmp_path / "c.csv", index=False)
    (tmp_path / "d.json").write_text(json.dumps(_disorder(0.1, 0.6)))
    (tmp_path / "m.json").write_text(json.dumps({}))
    (tmp_path / "v.json").write_text(json.dumps({}))
    assert ortho.main(["report", "--candidates", str(tmp_path / "c.csv"),
                       "--disorder", str(tmp_path / "d.json"), "--motifs", str(tmp_path / "m.json"),
                       "--dive", str(tmp_path / "v.json"), "--out-dir", str(tmp_path / "out"),
                       "--n-bootstrap", "100"]) == 0
    text = (tmp_path / "out" / "ORTHOGONAL.md").read_text()
    assert "**I1:** candidates more ordered" in text and "I2 (disorder QC)" in text
    assert (tmp_path / "out" / "orthogonal_candidates.csv").exists()
