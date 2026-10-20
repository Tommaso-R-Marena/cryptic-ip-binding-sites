"""Tests for study Q9: sequence against pocket formation, type and function.

The study's whole claim rests on one thing — that a difference it reports is a
difference in **arrangement** and not in composition — so the two tests that
matter most are the planted-signal and matched-noise pair at the bottom: the
pipeline must fire when an amphipathic stretch is planted in the pockets it
scores, and must stay quiet when the same construction carries no signal.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "pocket_sequence", ROOT / "scripts" / "pocket_sequence.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["pocket_sequence"] = script
SPEC.loader.exec_module(script)

ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
AMPHIPATHIC = "LKKLLKKLLKKLLKKLLKKL"


# ------------------------------------------------------- reading structures

PDB = """\
ATOM      1  N   MET A   1      11.104  13.207  10.000  1.00 90.00           N
ATOM      2  CA  MET A   1      12.104  13.207  10.000  1.00 90.00           C
ATOM      3  CA  LYS A   2      13.104  13.207  10.000  1.00 90.00           C
ATOM      4  CA  MSE A   3      14.104  13.207  10.000  1.00 90.00           C
ATOM      5  CA  LEU B   1      15.104  13.207  10.000  1.00 90.00           C
HETATM    6  O   HOH A  99      16.104  13.207  10.000  1.00 90.00           O
ATOM      7  CA  UNK A   4      17.104  13.207  10.000  1.00 90.00           C
"""


def test_sequences_are_read_per_chain_from_the_same_file_as_the_pockets(tmp_path):
    path = tmp_path / "model.pdb"
    path.write_text(PDB)
    out = script.chain_sequences(path)
    assert out["A"] == {1: "M", 2: "K", 3: "M"}      # MSE reads as methionine
    assert out["B"] == {1: "L"}


def test_waters_and_unknown_residues_are_not_sequence(tmp_path):
    path = tmp_path / "model.pdb"
    path.write_text(PDB)
    out = script.chain_sequences(path)
    assert 99 not in out["A"]
    assert 4 not in out["A"]


# --------------------------------------------------------------- windows

def test_the_local_window_is_the_contiguous_stretch_a_lining_spans():
    numbers = {n: ALPHABET[n % 20] for n in range(1, 60)}
    text = script.local_window(numbers, [10, 13, 17, 20])
    assert text == "".join(numbers[n] for n in range(10, 21))


def test_a_window_longer_than_the_cap_is_refused():
    numbers = {n: "A" for n in range(1, 200)}
    assert script.local_window(numbers, [1, 1 + script.MAX_WINDOW]) is None


def test_a_window_with_a_gap_in_the_structure_is_refused():
    """A missing residue is a hole in the window, not a shorter window."""
    numbers = {n: "A" for n in range(1, 40)}
    del numbers[15]
    assert script.local_window(numbers, [10, 20]) is None


def test_a_window_below_the_floor_is_refused():
    numbers = {n: "A" for n in range(1, 40)}
    assert script.local_window(numbers, [10, 12]) is None


def test_a_window_with_a_non_standard_letter_is_refused():
    numbers = {n: "A" for n in range(1, 40)}
    numbers[12] = "X"
    assert script.local_window(numbers, [10, 20]) is None


# ------------------------------------------------------- length matching

def test_matching_pairs_arms_on_window_length():
    left = [{"length": 10, "id": "a"}, {"length": 20, "id": "b"}]
    right = [{"length": 21, "id": "c"}, {"length": 11, "id": "d"},
             {"length": 40, "id": "e"}]
    matched_left, matched_right = script.match_on_length(
        left, right, np.random.default_rng(0))
    assert len(matched_left) == len(matched_right) == 2
    for one, other in zip(matched_left, matched_right):
        assert abs(one["length"] - other["length"]) <= script.LENGTH_TOLERANCE


def test_a_partner_outside_the_tolerance_is_not_paired():
    left = [{"length": 10, "id": "a"}]
    right = [{"length": 40, "id": "b"}]
    matched_left, matched_right = script.match_on_length(
        left, right, np.random.default_rng(0))
    assert matched_left == [] and matched_right == []


def test_matching_is_seeded():
    left = [{"length": 10, "id": f"a{i}"} for i in range(5)]
    right = [{"length": 10, "id": f"b{i}"} for i in range(20)]
    first = script.match_on_length(left, right, np.random.default_rng(3))[1]
    again = script.match_on_length(left, right, np.random.default_rng(3))[1]
    assert [r["id"] for r in first] == [r["id"] for r in again]


# ------------------------------------------------------------ the fences

def test_only_composition_z_enters_a_primary_test():
    """A raw measure is never promoted; the plan forbids it."""
    source = (ROOT / "scripts" / "pocket_sequence.py").read_text()
    assert "composition_z" in source
    # every statistic the contrast reads is a z column
    assert 'row.get(f"z_{statistic}")' in source


def test_the_primary_family_is_the_one_the_plan_fixed():
    assert script.PRIMARY == ("helix_charge_periodicity", "helix_hydropathy_periodicity")
    assert script.TOP_RANK == 5
    assert script.LENGTH_TOLERANCE == 2


def test_a_length_confounded_contrast_cannot_be_a_finding():
    entry = {"arm": "formation", "statistic": "helix_charge_periodicity",
             "runnable": True, "p_value": 0.0001, "is_evidence": True,
             "exceeds_power": True, "interval": [0.5, 1.5], "length_confounded": True,
             "n_groups": 9}
    script.holm_over([entry], ["formation:helix_charge_periodicity"])
    assert entry["is_finding"] is False


def test_an_underpowered_contrast_cannot_be_a_finding():
    entry = {"arm": "formation", "statistic": "helix_charge_periodicity",
             "runnable": True, "p_value": 0.0001, "is_evidence": True,
             "exceeds_power": False, "interval": [0.5, 1.5], "length_confounded": False,
             "n_groups": 9}
    script.holm_over([entry], ["formation:helix_charge_periodicity"])
    assert entry["is_finding"] is False


# --------------------------------------------- planted signal, matched noise

def _synthetic(planted: bool, n_proteins: int = 8, seed: int = 0):
    rng = np.random.default_rng(seed)
    sequences, pockets = [], []
    for index in range(n_proteins):
        stem = f"AF-P{index}-F1-model_v4"
        letters = list(rng.choice(ALPHABET, size=400))
        top = list(range(20, 220, 40))
        bottom = list(range(240, 400, 32))
        if planted:
            for start in top:
                for offset in range(20):
                    letters[start + offset - 1] = AMPHIPATHIC[offset]
        sequences.append({"structure": stem + ".pdb", "stem": stem, "chain": "A",
                          "first": 1, "last": 400, "n_residues": 400,
                          "numbers": ",".join(str(n) for n in range(1, 401)),
                          "letters": "".join(letters)})
        for rank, start in enumerate(top + bottom, start=1):
            pockets.append({"uniprot_id": f"P{index}", "model": stem + ".pdb",
                            "pocket_id": rank, "cluster": f"c{index}",
                            "residues": ",".join(str(start + i) for i in range(0, 20, 3))})
    return pd.DataFrame(sequences), pd.DataFrame(pockets)


def test_the_pipeline_finds_a_planted_arrangement():
    sequences, pockets = _synthetic(planted=True)
    report = script.run(sequences, pockets, None, None, shuffles=200,
                        n_bootstrap=400, per_site=3)
    primary = [t for t in report["tests"]
               if t.get("runnable") and t["statistic"] in script.PRIMARY]
    assert primary
    assert all(t["mean_site_z"] > 2.0 for t in primary)
    assert any(t["is_finding"] for t in primary)


def test_the_pipeline_stays_quiet_on_matched_noise():
    """Same construction, no planted arrangement: nothing may fire."""
    sequences, pockets = _synthetic(planted=False, seed=1)
    report = script.run(sequences, pockets, None, None, shuffles=200,
                        n_bootstrap=400, per_site=3)
    primary = [t for t in report["tests"]
               if t.get("runnable") and t["statistic"] in script.PRIMARY]
    assert primary
    assert not any(t["is_finding"] for t in primary)


def test_the_report_states_the_decomposition_and_its_limits():
    sequences, pockets = _synthetic(planted=False, seed=2)
    text = script.markdown(script.run(sequences, pockets, None, None, shuffles=100,
                                      n_bootstrap=200, per_site=3))
    assert "composition" in text and "order" in text
    assert "never on a raw measure" in text
    assert "predicts fpocket, not nature" in text
