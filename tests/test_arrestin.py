"""scripts/arrestin.py (docs/ARRESTIN_PLAN.md) on synthetic inputs."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def arrestin():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("arrestin", ROOT / "scripts" / "arrestin.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["arrestin"] = module
    spec.loader.exec_module(module)
    return module


USALIGN = """
Name of Structure_1: target.pdb:A (to be superimposed onto Structure_2)
Name of Structure_2: ref.pdb:A
Aligned length=    5, RMSD=   1.20, Seq_ID=n_identical/n_aligned= 0.400
TM-score= 0.61000 (normalized by length of Structure_1: L=8, d0=0.50)
TM-score= 0.72000 (normalized by length of Structure_2: L=7, d0=0.50)
(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
MK-TAYIAK
::.:: :.
MKRTA-IA-
"""


def test_usalign_parsing_and_mapping(arrestin):
    parsed = arrestin.parse_usalign(USALIGN)
    assert parsed["tm_norm_structure1"] == 0.61 and parsed["tm_norm_structure2"] == 0.72
    assert parsed["aligned_length"] == 5 and parsed["rmsd"] == 1.2
    target_numbers = [10, 11, 12, 13, 14, 15, 16, 17]  # MKTAYIAK
    ref_numbers = [100, 101, 102, 103, 104, 105, 106]    # MKRTAIA
    mapping = arrestin.mapping_from_alignment(parsed["seq1"], parsed["markers"], parsed["seq2"],
                                              target_numbers, ref_numbers)
    # columns: M/M ':', K/K ':', -/R '.', T/T ':', A/A ':', Y/- ' ', I/I ':', A/A '.', K/- (no marker)
    assert mapping == {100: 10, 101: 11, 103: 12, 104: 13, 105: 15}


def test_residue_lists_and_jaccard(arrestin):
    assert arrestin.parse_residue_list("12,13,14") == {12, 13, 14}
    assert arrestin.parse_residue_list("A:12;A:13 B:14A") == {12, 13, 14}
    assert arrestin.parse_residue_list(float("nan")) == set()
    assert arrestin.jaccard({1, 2, 3}, {2, 3, 4}) == pytest.approx(0.5)
    assert arrestin.jaccard(set(), set()) == 0.0


def test_family_definitions(arrestin):
    frame = pd.DataFrame({
        "Entry": ["P49407", "Q8TBH0", "P53244", "Q9H3M7", "P0XXXX"],
        "Gene Names (primary)": ["ARRB1", "ARRDC2", "ART5", "TXNIP", "VPS26A"],
        "Protein names": ["Beta-arrestin-1 (Arrestin beta-1)", "Arrestin domain-containing protein 2",
                          "Arrestin-related trafficking adapter 5", "Thioredoxin-interacting protein (VDUP1)",
                          "Vacuolar protein sorting-associated protein 26A"],
        "Pfam": ["PF02752;PF00339;", "PF02752;PF00339;", "PF02752;", "PF02752;PF00339;", "PF03643;"],
        "Organism (ID)": ["9606"] * 5, "Length": ["418"] * 5})
    c = arrestin.classify(frame).set_index("uniprot_id")
    assert c.loc["P49407", "classic"] and not c.loc["P49407", "alpha"]
    assert c.loc["Q8TBH0", "alpha"] and c.loc["P53244", "alpha"] and c.loc["Q9H3M7", "alpha"]
    assert not c.loc["P0XXXX", "arrestin_fold"]


def _proteins(n=400, members=(), scores=None, clusters=None, seen=()):
    rng = np.random.default_rng(0)
    ids = [f"X{i:05d}" for i in range(n)]
    s = rng.random(n) if scores is None else scores
    return pd.DataFrame({"uniprot_id": ids, "learned_score": s, "seen": [i in seen for i in ids],
                         "cluster": clusters if clusters is not None else ids, "organism_key": "human"})


def test_family_test_supported_when_members_rank_high(arrestin):
    scores = np.linspace(0, 1, 400)
    frame = _proteins(scores=scores)
    members = set(frame["uniprot_id"].iloc[-8:])  # the top 8, each its own cluster
    r = arrestin.family_test(frame, members, n_bootstrap=300)
    assert r["decision"] == "supported" and r["roc_auc"]["point"] > 0.95


def test_family_test_not_evaluable_with_few_clusters(arrestin):
    scores = np.linspace(0, 1, 400)
    clusters = [f"c{i}" for i in range(400)]
    for i in range(392, 400):
        clusters[i] = "paralogues" if i % 2 else "paralogues2"
    frame = _proteins(scores=scores, clusters=clusters)
    members = set(frame["uniprot_id"].iloc[-8:])
    r = arrestin.family_test(frame, members, n_bootstrap=100)
    assert r["member_clusters"] == 2 and r["decision"].startswith("not evaluable")


def test_family_test_excludes_seen_proteins(arrestin):
    frame = _proteins(scores=np.linspace(0, 1, 400))
    members = set(frame["uniprot_id"].iloc[-8:])
    frame["seen"] = frame["uniprot_id"].isin(members)
    r = arrestin.family_test(frame, members, n_bootstrap=50)
    assert r["members_unseen"] == 0 and r["decision"].startswith("not evaluable")


def test_family_test_not_supported_at_random(arrestin):
    frame = _proteins(scores=np.random.default_rng(3).random(400))
    members = set(frame["uniprot_id"].iloc[::50])
    assert arrestin.family_test(frame, members, n_bootstrap=300)["decision"] == "not supported"


def test_conservation_scoring(arrestin):
    alignment = {"T": "MK-RA", "h1": "MKQRA", "h2": "MR-KA", "h3": "M--EA"}
    scores = arrestin.conservation(alignment, "T", [2, 3])  # target K2 (column 1), R3 (column 3)
    assert scores[2] == pytest.approx(2 / 3)  # K, R, gap
    assert scores[3] == pytest.approx(2 / 3)  # R, K, E
    d = arrestin.conservation_decision("MKRA", [2, 3], {2: 0.9, 3: 0.9}, n_homologues=50)
    assert d["decision"] == "not conserved"  # only 2 basic positions: need 3
    d = arrestin.conservation_decision("MKRKA", [2, 3, 4], {2: 0.9, 3: 0.85, 4: 0.8}, n_homologues=50)
    assert d["decision"] == "conserved"
    d = arrestin.conservation_decision("MKRKA", [2, 3, 4], {2: 0.9, 3: 0.85, 4: 0.8}, n_homologues=4)
    assert d["decision"].startswith("not evaluable")


def _tasks():
    t = [
        {"task_id": 0, "protein": "Q8TBH0", "kind": "lead_site", "ligand": "IHP", "site": "s1"},
        {"task_id": 1, "protein": "Q8TBH0", "kind": "negative", "ligand": "IHP", "site": "p1"},
        {"task_id": 2, "protein": "Q8TBH0", "kind": "negative", "ligand": "IHP", "site": "p2"},
        {"task_id": 3, "protein": "P49407", "kind": "positive_crystal", "ligand": "IHP", "site": "x"},
        {"task_id": 4, "protein": "P49407", "kind": "positive_af", "ligand": "IHP", "site": "y"},
    ]
    return t


def test_decision_all_four_criteria(arrestin):
    tasks = _tasks()
    results = {0: {"best_score": -9.0, "convergence": {"converged": True, "fraction": 0.6}},
               1: {"best_score": -6.0}, 2: {"best_score": -7.0},
               3: {"best_score": -8.0, "success": 1.0}, 4: {"best_score": -7.5}}
    mapping = {"Q8TBH0": {"overlap": True, "has_mapped_site": True, "lead_jaccard": 0.3}}
    conserve = {"Q8TBH0": {"decision": "conserved"}}
    d = arrestin.decide(tasks, results, mapping, conserve)["proteins"]["Q8TBH0"]
    assert d["verdict"] == "supported for experiment" and d["failing"] == []
    # a negative pocket scoring better than the lead site fails criterion 4
    results[2]["best_score"] = -9.5
    d = arrestin.decide(tasks, results, mapping, conserve)["proteins"]["Q8TBH0"]
    assert d["verdict"] == "not supported" and d["failing"] == ["4_scores"]
    # weaker than the weakest positive control fails criterion 4 as well
    results[2]["best_score"] = -7.0
    results[0]["best_score"] = -7.2
    d = arrestin.decide(tasks, results, mapping, conserve)["proteins"]["Q8TBH0"]
    assert "4_scores" in d["failing"]


def test_invalid_protocol_makes_docking_criteria_not_evaluable(arrestin):
    tasks = _tasks()
    results = {0: {"best_score": -9.0, "convergence": {"converged": True}}, 1: {"best_score": -6.0},
               2: {"best_score": -7.0}, 3: {"best_score": -8.0, "success": 0.0}, 4: {"best_score": -7.5}}
    mapping = {"Q8TBH0": {"overlap": True, "has_mapped_site": True}}
    out = arrestin.decide(tasks, results, mapping, {"Q8TBH0": {"decision": "conserved"}})
    d = out["proteins"]["Q8TBH0"]
    assert not out["protocol_validity"]["valid"]
    assert d["criteria"]["3_convergence"]["decision"].startswith("not evaluable")
    assert d["verdict"] == "not supported" and set(d["failing"]) == {"3_convergence", "4_scores"}
    # ART5 has no tasks at all: every criterion fails, none silently passes
    assert arrestin.decide(tasks, results, mapping, {})["proteins"]["P53244"]["verdict"] == "not supported"


def test_a_copy_counts_once_for_protocol_validity(arrestin):
    tasks = _tasks() + [
        {"task_id": 5, "protein": "P49407", "kind": "positive_crystal", "ligand": "IHP", "site": "x"},
        {"task_id": 6, "protein": "P49407", "kind": "positive_crystal", "ligand": "IHP", "site": "z"},
    ]
    results = {0: {"best_score": -9.0, "convergence": {"converged": True}}, 1: {"best_score": -6.0},
               2: {"best_score": -7.0}, 3: {"best_score": -8.0, "success": 1.0}, 4: {"best_score": -7.5},
               5: {"best_score": -8.0, "success": 1.0}, 6: {"best_score": -5.0, "success": 0.0}}
    out = arrestin.decide(tasks, results, {}, {})
    # sites x (twice, success 1) and z (success 0): the duplicate of x must not tip the mean to 2/3
    assert out["protocol_validity"]["crystal_sites"] == 2
    assert out["protocol_validity"]["mean_top_pose_success"] == pytest.approx(0.5)
