"""Tests for study Q9: sequence against pocket formation, type and function.

The study's whole claim rests on one thing — that a difference it reports is a
difference in **arrangement** and not in composition — so the two tests that
matter most are the planted-signal and matched-noise pair at the bottom: the
pipeline must fire when an amphipathic stretch is planted in the pockets it
scores, and must stay quiet when the same construction carries no signal.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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

def _synthetic(planted: bool, n_proteins: int = 8, seed: int = 0, middle: int = 0,
               singletons: int = 0):
    """Eight proteins, five top-ranked pockets, five bottom-ranked, and
    optionally `middle` more plus `singletons` one-pocket proteins that no arm
    can select from - which is what makes the deferred scoring worth having."""
    rng = np.random.default_rng(seed)
    sequences, pockets = [], []
    for index in range(n_proteins + singletons):
        stem = f"AF-P{index}-F1-model_v4"
        letters = list(rng.choice(ALPHABET, size=400))
        top = list(range(20, 220, 40))
        bottom = list(range(240, 400, 32))
        between = [225 + 4 * i for i in range(middle)]
        if planted:
            for start in top:
                for offset in range(20):
                    letters[start + offset - 1] = AMPHIPATHIC[offset]
        sequences.append({"structure": stem + ".pdb", "stem": stem, "chain": "A",
                          "first": 1, "last": 400, "n_residues": 400,
                          "numbers": ",".join(str(n) for n in range(1, 401)),
                          "letters": "".join(letters)})
        starts = top + between + bottom if index < n_proteins else top[:1]
        for rank, start in enumerate(starts, start=1):
            pockets.append({"uniprot_id": f"P{index}", "model": stem + ".pdb",
                            "pocket_id": rank, "cluster": f"c{index}",
                            "residues": ",".join(str(start + i) for i in range(0, 20, 3))})
    return pd.DataFrame(sequences), pd.DataFrame(pockets)


@pytest.mark.slow
def test_the_pipeline_finds_a_planted_arrangement():
    sequences, pockets = _synthetic(planted=True)
    report = script.run(sequences, pockets, None, None, shuffles=200,
                        n_bootstrap=400, per_site=3)
    primary = [t for t in report["tests"]
               if t.get("runnable") and t["statistic"] in script.PRIMARY]
    assert primary
    assert all(t["mean_site_z"] > 2.0 for t in primary)
    assert any(t["is_finding"] for t in primary)


@pytest.mark.slow
def test_the_pipeline_stays_quiet_on_matched_noise():
    """Same construction, no planted arrangement: nothing may fire."""
    sequences, pockets = _synthetic(planted=False, seed=1)
    report = script.run(sequences, pockets, None, None, shuffles=200,
                        n_bootstrap=400, per_site=3)
    primary = [t for t in report["tests"]
               if t.get("runnable") and t["statistic"] in script.PRIMARY]
    assert primary
    assert not any(t["is_finding"] for t in primary)


@pytest.mark.slow
def test_the_report_states_the_decomposition_and_its_limits():
    sequences, pockets = _synthetic(planted=False, seed=2)
    text = script.markdown(script.run(sequences, pockets, None, None, shuffles=100,
                                      n_bootstrap=200, per_site=3))
    assert "composition" in text and "order" in text
    assert "never on a raw measure" in text
    assert "predicts fpocket, not nature" in text


@pytest.mark.slow
def test_only_the_windows_an_arm_selected_are_scored():
    """A window no arm can reach must not cost 2000 shuffles of eight statistics.

    The one-pocket proteins are unreachable by construction: the pre-registered
    rule needs ten windows in a protein and the amendment's needs two.
    """
    sequences, pockets = _synthetic(planted=False, seed=4, middle=6, singletons=3)
    report = script.run(sequences, pockets, None, None, shuffles=50,
                        n_bootstrap=100, per_site=2)
    assert report["n_scored_local"] < report["n_pocket_windows"]
    assert report["n_scored_dense"] < report["n_dense_windows"]
    assert report["n_scored_windows"] > 0


def test_a_window_s_null_does_not_depend_on_what_selected_it():
    """The seed travels with the row, so the z is a property of the window."""
    numbers = {n: ALPHABET[n % 20] for n in range(1, 60)}
    text = script.local_window(numbers, [10, 30])
    assert text is not None
    row = {"window": text, "seed": 7}
    other = {"window": text, "seed": 7}
    script.score_windows([row], ["helix_charge_periodicity"], shuffles=80)
    script.score_windows([other], ["helix_charge_periodicity"], shuffles=80)
    assert row["z_helix_charge_periodicity"] == other["z_helix_charge_periodicity"]


def test_scoring_twice_does_not_rescore():
    row = {"window": "LKKLLKKLLKKLLKKLLKKL", "seed": 1}
    script.score_windows([row], ["helix_charge_periodicity"], shuffles=50)
    first = row["z_helix_charge_periodicity"]
    row["window"] = "AAAAAAAAAAAAAAAAAAAA"      # would score differently if rescored
    script.score_windows([row], ["helix_charge_periodicity"], shuffles=50)
    assert row["z_helix_charge_periodicity"] == first


# ------------------------------------------------- amendment 1: the windows

def test_the_dense_window_is_always_the_same_length():
    """Amendment 1 section 4: a fixed width makes the length confound zero."""
    numbers = {n: ALPHABET[n % 20] for n in range(1, 200)}
    starts, texts = script.dense_starts(numbers)
    for lining in ([10, 90], [5, 7, 9], [1, 150], list(range(40, 120, 7))):
        text = script.dense_window(starts, texts, lining)
        assert text is not None
        assert len(text) == script.DENSE_WINDOW


def test_the_dense_window_is_the_stretch_holding_most_of_the_lining():
    numbers = {n: ALPHABET[n % 20] for n in range(1, 200)}
    starts, texts = script.dense_starts(numbers)
    lining = [100, 101, 102, 103, 104, 180]
    text = script.dense_window(starts, texts, lining)
    assert text is not None
    # The five clustered positions outvote the lone one, so the window holds them.
    windows = {start: body for start, body in texts.items()}
    chosen = [start for start, body in windows.items() if body == text]
    assert any(start <= 100 and 104 < start + script.DENSE_WINDOW for start in chosen)


def test_the_dense_window_reads_positions_and_never_letters():
    """Two chains with the same numbering pick the same stretch, whatever the letters."""
    first = {n: ALPHABET[n % 20] for n in range(1, 120)}
    second = {n: "A" if n % 2 else "K" for n in range(1, 120)}
    lining = [50, 52, 54, 56]
    one = script.dense_window(*script.dense_starts(first), lining)
    two = script.dense_window(*script.dense_starts(second), lining)
    assert one is not None and two is not None
    assert one != two                                   # different letters
    assert first[50] in one and second[50] in two       # same stretch


def test_a_chain_shorter_than_the_dense_window_yields_nothing():
    numbers = {n: "A" for n in range(1, script.DENSE_WINDOW)}
    starts, texts = script.dense_starts(numbers)
    assert starts.size == 0
    assert script.dense_window(starts, texts, [3, 5]) is None


def test_the_dense_window_refuses_a_gap_in_the_model():
    numbers = {n: "A" for n in range(1, 200)}
    for missing in range(1, 200, 3):
        numbers.pop(missing, None)
    starts, _ = script.dense_starts(numbers)
    assert starts.size == 0


def test_every_refusal_has_a_reason_and_it_agrees_with_the_window():
    """Amendment 1 section 6: the diagnostics reproduce the refusal, never decide it."""
    long = {n: "A" for n in range(1, 200)}
    gapped = {n: "A" for n in range(1, 40)}
    del gapped[15]
    short = {n: "A" for n in range(1, 40)}
    cases = {
        "span over the cap": (long, [1, 1 + script.MAX_WINDOW]),
        "gap in the model": (gapped, [10, 20]),
        "under the floor": (short, [10, 12]),
        "no positions": (short, []),
    }
    for expected, (numbers, lining) in cases.items():
        assert script.refusal(numbers, lining) == expected
        assert script.local_window(numbers, lining) is None
    accepted = {n: ALPHABET[n % 20] for n in range(1, 40)}
    assert script.refusal(accepted, [10, 20]) is None
    assert script.local_window(accepted, [10, 20]) is not None


def test_a_window_s_seed_is_a_function_of_its_letters():
    """Across processes too: `hash()` is salted and cannot seed a committed run."""
    assert script.hash_window("LKKLLKKL") == script.hash_window("LKKLLKKL")
    assert script.hash_window("LKKLLKKL") != script.hash_window("LKKLLKKK")
    # Pinned, so a change of hash cannot silently change every committed null.
    assert script.hash_window("LKKLLKKL") == 0xA42640709CC9


# --------------------------------------- amendment 1: the re-specified arm

def _collinear(n_proteins: int = 10, seed: int = 0):
    """Rank collinear with window length, which is what broke the plan's arm.

    Each protein's pocket 1 spans far more than the cap and its later pockets
    span less and less - exactly fpocket's volume ordering - so the plan's
    top-five-against-bottom-five contrast has almost nothing to match.
    """
    rng = np.random.default_rng(seed)
    sequences, pockets = [], []
    for index in range(n_proteins):
        stem = f"AF-Q{index}-F1-model_v4"
        letters = list(rng.choice(ALPHABET, size=1700))
        sequences.append({"structure": stem + ".pdb", "stem": stem, "chain": "A",
                          "first": 1, "last": 1700, "n_residues": 1700,
                          "numbers": ",".join(str(n) for n in range(1, 1701)),
                          "letters": "".join(letters)})
        for rank in range(1, 31):
            # Span falls with rank, with the scatter a real model has, so a
            # length class holds several ranks and the split rule can match
            # within one - which is the whole point of taking the split inside
            # the surviving population.
            span = max(8, int(round(100 - 3 * rank + rng.normal(0.0, 9.0))))
            start = 10 + 50 * rank
            pockets.append({"uniprot_id": f"Q{index}", "model": stem + ".pdb",
                            "pocket_id": rank, "cluster": f"g{index}",
                            "residues": f"{start},{start + span}"})
    return pd.DataFrame(sequences), pd.DataFrame(pockets)


def _windows(sequences, pockets):
    return script.windows_for_pockets(
        pockets, script.sequence_index(sequences),
        key_of=lambda p: str(p.model).replace(".pdb", ""),
        chain_of=lambda p: "A", group_of=lambda p: str(p.cluster))


def _selected(rows, rule, **kwargs):
    """Which windows a rule actually paid an order null for."""
    rule(rows, n_bootstrap=20, statistics=["lz_normalised"], shuffles=10, **kwargs)
    return [row for row in rows if row.get("scored")]


def test_the_plan_s_rule_starves_when_rank_is_collinear_with_length():
    """Amendment 1 section 2, as a test rather than as a claim.

    Every pocket's span falls with its rank, so the span cap keeps the
    top-ranked ones out entirely and the plan's top-against-bottom contrast has
    nothing at all to match. The split rule compares among the windows that do
    exist, and reaches them.
    """
    sequences, pockets = _collinear(seed=3)
    rows = _windows(sequences, pockets)
    assert rows                                      # windows did exist
    assert _selected(_windows(sequences, pockets), script.formation) == []
    split = _selected(_windows(sequences, pockets), script.formation_split,
                      arm="formation_matched")
    assert len(split) >= 2 * 10                      # both halves, every protein


def test_the_split_rule_reaches_more_windows_where_the_plan_s_rule_also_runs():
    sequences, pockets = _synthetic(planted=False, seed=5, middle=6)
    by_plan = _selected(_windows(sequences, pockets), script.formation)
    by_split = _selected(_windows(sequences, pockets), script.formation_split,
                         arm="formation_matched")
    assert by_plan
    assert len(by_split) > len(by_plan)


def test_the_heldout_arm_scores_no_window_of_a_group_the_plan_s_rule_saw():
    """The look is disclosed by construction: the seen groups are dropped.

    The two halves of the fixture are the two halves of the real dataset: one
    where rank and length are collinear and the plan's rule selects nothing,
    one where they are not and it selects. The held-out arm must reach only the
    first.
    """
    reachable, starved = _synthetic(planted=False, seed=7), _collinear(seed=8)
    sequences = pd.concat([reachable[0], starved[0]], ignore_index=True)
    pockets = pd.concat([reachable[1], starved[1]], ignore_index=True)
    rows = _windows(sequences, pockets)
    _, selected = script.formation(rows, n_bootstrap=20,
                                   statistics=["lz_normalised"], shuffles=10)
    assert selected                                  # the plan's rule did select some
    held = _selected(_windows(sequences, pockets), script.formation_split,
                     arm="formation_heldout", skip_groups=selected)
    assert held                                      # and the exclusion is not total
    assert all(str(row["group"]) not in selected for row in held)
    assert {str(row["group"]) for row in held} == {f"g{i}" for i in range(10)}


def test_an_exploratory_arm_can_never_be_a_finding():
    """Amendment 1 declares them exploratory, so no outcome may promote them."""
    for arm in ("formation_matched", "formation_dense", "type_dense",
                "function_dense", "formation_matched_dense"):
        entry = {"arm": arm, "statistic": "helix_charge_periodicity",
                 "runnable": True, "p_value": 1e-9, "is_evidence": True,
                 "exceeds_power": True, "interval": [2.0, 3.0],
                 "length_confounded": False, "n_groups": 40}
        # Even if someone puts it in a family by mistake, the arm fence holds.
        script.holm_over([entry], [f"{arm}:helix_charge_periodicity"])
        assert entry["is_finding"] is False
    assert "formation_matched" not in script.CONFIRMATORY_ARMS
    assert "formation_heldout" in script.CONFIRMATORY_ARMS


def test_a_cap_keeps_a_pair_count_and_its_pairing():
    rng = np.random.default_rng(0)
    left = [{"length": i} for i in range(20)]
    right = [{"length": i + 1} for i in range(20)]
    kept_left, kept_right = script.cap_pairs(left, right, 5, rng)
    assert len(kept_left) == len(kept_right) == 5
    for one, two in zip(kept_left, kept_right):
        assert int(two["length"]) - int(one["length"]) == 1
    assert script.cap_pairs(left, right, 50, rng)[0] is left


def test_several_sequence_tables_read_as_one(tmp_path):
    first, second = tmp_path / "a.csv", tmp_path / "b.csv"
    pd.DataFrame([{"stem": "x", "chain": "A", "numbers": "1,2", "letters": "AK"}]
                 ).to_csv(first, index=False)
    pd.DataFrame([{"stem": "y", "chain": "B", "numbers": "3,4", "letters": "LM"}]
                 ).to_csv(second, index=False)
    index = script.sequence_index(script.read_sequences([first, second, tmp_path / "no"]))
    assert index[("x", "A")] == {1: "A", 2: "K"}
    assert index[("y", "B")] == {3: "L", 4: "M"}


def test_the_amendment_is_committed_and_the_report_points_at_it():
    text = (ROOT / "docs" / "POCKET_SEQUENCE_PLAN_AMENDMENT_1.md").read_text()
    assert "+1.317" in text and "0.037" in text        # the look, disclosed
    assert "Nothing about the type arm or" in text
    source = (ROOT / "scripts" / "pocket_sequence.py").read_text()
    assert "POCKET_SEQUENCE_PLAN_AMENDMENT_1.md" in source


def test_the_classes_step_reads_positions_and_writes_an_element(tmp_path):
    sequences, pockets = _synthetic(planted=False, seed=11)
    sequences_path, pockets_path = tmp_path / "s.csv", tmp_path / "p.csv"
    sequences.to_csv(sequences_path, index=False)
    pockets.to_csv(pockets_path, index=False)
    out = tmp_path / "classes.csv"
    args = script.main(["classes", "--sequences", str(sequences_path),
                        "--pockets", str(pockets_path), "--out", str(out),
                        "--draws", "60"])
    assert args == 0
    frame = pd.read_csv(out)
    assert set(frame.columns) >= {"uniprot_id", "pocket_id", "element"}
    assert not frame.empty
    assert set(frame["element"]) <= {"none", "helix_face", "strand_face",
                                     "periodic_other"}


def _look(n_site, n_groups, difference):
    return [{"arm": script.RECORDED_LOOK["arm"],
             "statistic": script.RECORDED_LOOK["statistic"], "runnable": True,
             "n_site": n_site, "n_groups": n_groups, "difference": difference}]


def test_a_redrawn_null_still_counts_as_the_same_look():
    """Amendment 2 section 2: the real gap was 0.41 of one standard error."""
    recorded = script.RECORDED_LOOK
    check = script.look_reproduction(
        _look(recorded["n_site"], recorded["n_groups"], 1.309855329), shuffles=2000)
    assert check["selection_reproduced"] is True
    assert check["verdict"] == "reproduced"
    assert check["gap_in_spreads"] < 1.0
    # And the spread is the quantity the shuffle count implies, not a constant.
    assert check["expected_spread"] == pytest.approx(
        script.null_spread(2000, recorded["n_site"]))


def test_a_difference_a_redraw_cannot_explain_is_not_the_same_look():
    recorded = script.RECORDED_LOOK
    check = script.look_reproduction(
        _look(recorded["n_site"], recorded["n_groups"], 0.4), shuffles=2000)
    assert check["selection_reproduced"] is True
    assert check["verdict"] == "not reproduced"
    assert check["gap_in_spreads"] > script.LOOK_SPREAD_BOUND


def test_a_changed_selection_is_never_excused_by_noise():
    """No tolerance at all on the selection: no re-seeding can move it."""
    recorded = script.RECORDED_LOOK
    for n_site, n_groups in ((recorded["n_site"] + 1, recorded["n_groups"]),
                             (recorded["n_site"], recorded["n_groups"] - 1)):
        check = script.look_reproduction(
            _look(n_site, n_groups, recorded["difference"]), shuffles=2000)
        assert check["selection_reproduced"] is False
        assert check["verdict"] == "not reproduced"


def test_an_arm_that_did_not_run_cannot_reproduce_anything():
    check = script.look_reproduction([], shuffles=2000)
    assert check["verdict"] == "not reproduced"
    assert check["selection_reproduced"] is False
    assert "did not run" in str(check["reason"])


def test_the_redraw_spread_shrinks_with_shuffles_and_windows():
    assert script.null_spread(8000, 10) < script.null_spread(2000, 10)
    assert script.null_spread(2000, 40) < script.null_spread(2000, 10)
    assert script.null_spread(0, 10) == float("inf")
    assert script.null_spread(2000, 0) == float("inf")


def test_the_report_says_whether_the_disclosed_look_reproduced():
    """The held-out arm is only held out if the data underneath it is unchanged."""
    base = {"seed": 1, "shuffles": 2, "n_pocket_windows": 0, "n_dense_windows": 0,
            "n_classified": 0, "n_site_windows": 0, "n_control_windows": 0,
            "n_scored_windows": 0, "top_rank": 5, "heldout_groups_excluded": ["a"],
            "window_yield": {}, "tests": []}
    good = dict(base, look_reproduction=script.look_reproduction(
        _look(10, 9, 1.309855329), shuffles=2000))
    assert "reproduced the disclosed look" in script.markdown(good)

    # The selection held and the contrast moved further than a redraw explains:
    # the exclusion is right, but something unexplained moved the number.
    drifted = dict(base, look_reproduction=script.look_reproduction(
        _look(10, 9, 0.4), shuffles=2000))
    text = script.markdown(drifted)
    assert "selection reproduced but the difference did not" in text
    assert "before the held-out arm is read as confirmatory" in text

    # A changed selection is the serious case, and reads as exploratory.
    changed = dict(base, look_reproduction=script.look_reproduction(
        _look(4, 3, 0.2), shuffles=2000))
    text = script.markdown(changed)
    assert "did not select what it selected before" in text
    assert "must be read as exploratory" in text


def test_the_recorded_look_is_the_one_the_amendment_discloses():
    text = (ROOT / "docs" / "POCKET_SEQUENCE_PLAN_AMENDMENT_1.md").read_text()
    assert str(script.RECORDED_LOOK["n_site"]) in text
    assert str(script.RECORDED_LOOK["n_groups"]) in text
    assert f"+{script.RECORDED_LOOK['difference']}" in text


# ------------------------------------------- what the committed run must keep saying

def _committed():
    return json.loads(
        (ROOT / "results" / "pocket_sequence" / "pocket_sequence.json").read_text())


def test_the_committed_run_found_nothing_and_says_so():
    report = _committed()
    primary = [t for t in report["tests"]
               if t.get("runnable") and t.get("family") == "primary"]
    assert primary
    assert not any(t["is_finding"] for t in primary), \
        "a Q9 primary test became a finding; the reports and amendments say none did"


def test_the_confirmatory_heldout_arm_is_powered_and_null():
    """Amendment 2 section 1: this is the arm the whole re-specification was for."""
    report = _committed()
    rows = [t for t in report["tests"]
            if t.get("runnable") and t["arm"] == "formation_heldout"
            and t["statistic"] in script.PRIMARY]
    assert len(rows) == len(script.PRIMARY)
    for row in rows:
        assert row["n_groups"] >= 5, "below the plan's evidence floor"
        assert row["exceeds_power"] is False, \
            "the arm now detects its own difference, so it is no longer a null result"
        assert abs(row["difference"]) < row["minimum_detectable"]
        assert row["is_finding"] is False


def test_the_underpowered_look_did_not_survive_power():
    """The +1.3 of amendment 1 section 1 against every better-powered form of it."""
    report = _committed()
    by_arm = {t["arm"]: t for t in report["tests"]
              if t.get("runnable") and t["statistic"] == "helix_hydropathy_periodicity"}
    assert by_arm["formation"]["difference"] > 1.0          # the look, reproduced
    assert by_arm["formation"]["n_site"] == 10
    for arm in ("formation_matched", "formation_dense", "formation_matched_dense",
                "formation_heldout"):
        row = by_arm[arm]
        assert row["n_site"] >= 100, f"{arm} is not the powered form of the contrast"
        assert abs(row["difference"]) < 0.2, \
            f"{arm} no longer collapses the look; the amendments say it does"


def test_every_exploratory_arm_in_the_committed_run_is_marked_exploratory():
    report = _committed()
    for test in report["tests"]:
        if not test.get("runnable"):
            continue
        if test["arm"] in script.CONFIRMATORY_ARMS:
            continue
        assert test.get("family") == "exploratory", test["arm"]
        assert test.get("p_holm") is None
        assert test["is_finding"] is False


def test_the_committed_report_is_what_the_committed_numbers_render_to():
    """Nothing in the report was written by hand, and nothing has drifted from it."""
    report = _committed()
    committed = (ROOT / "results" / "pocket_sequence" / "POCKET_SEQUENCE.md").read_text()
    assert script.markdown(report) == committed


def test_the_committed_run_reproduced_the_disclosed_look():
    """Amendment 2 section 2: the held-out arm is held out from what it claims."""
    check = _committed()["look_reproduction"]
    assert check["selection_reproduced"] is True
    assert check["verdict"] == "reproduced"
    assert check["gap_in_spreads"] < 1.0
    assert check["n_site"] == script.RECORDED_LOOK["n_site"]
    assert check["n_groups"] == script.RECORDED_LOOK["n_groups"]
