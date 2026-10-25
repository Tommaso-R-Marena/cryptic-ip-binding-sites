"""Tests for study Q10: the composition half, and the gate.

The study's whole architecture rests on one claim — that a positive control
which fails makes every other row unreadable rather than null — so the tests
that matter most are the gate fence and the planted/quiet pair that calibrate
it. The rest attack the null's correctness and the window's legitimacy.
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
    "pocket_composition", ROOT / "scripts" / "pocket_composition.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["pocket_composition"] = script
SPEC.loader.exec_module(script)

cp = script.cp
ps = script.ps
ALPHABET = list(cp.ALPHABET)


# --------------------------------------------------------- the lining window

def test_the_window_is_the_lining_s_own_residues_gap_closed():
    numbers = {n: ALPHABET[n % 20] for n in range(1, 200)}
    positions = [10, 20, 30, 40, 50, 60, 70, 80]
    assert script.lining_letters(numbers, positions) == \
        "".join(numbers[n] for n in positions)


def test_a_lining_under_the_floor_is_refused():
    numbers = {n: "A" for n in range(1, 200)}
    assert script.lining_letters(numbers, list(range(1, script.MIN_LINING))) is None
    assert script.lining_letters(numbers, list(range(1, script.MIN_LINING + 1))) is not None


def test_a_position_missing_from_the_model_refuses_the_lining():
    """A hole in the lining is not a shorter lining."""
    numbers = {n: "A" for n in range(1, 200)}
    del numbers[50]
    assert script.lining_letters(numbers, [10, 20, 30, 40, 50, 60, 70, 80]) is None


def test_the_span_cap_that_crippled_study_q9_does_not_apply_here():
    """The whole reason the window changed: composition cannot see a gap."""
    numbers = {n: ALPHABET[n % 20] for n in range(1, 1000)}
    spread = [1, 200, 400, 500, 600, 700, 800, 900]
    assert span_of(spread) > ps.MAX_WINDOW
    assert ps.local_window(numbers, spread) is None       # Q9 refuses it
    assert script.lining_letters(numbers, spread) is not None   # Q10 does not


def span_of(positions):
    return max(positions) - min(positions) + 1


def test_a_lining_larger_than_half_its_chain_is_refused():
    """Past a half the rest of the chain is not a background to sit against."""
    chain = script.Chain({n: ALPHABET[n % 20] for n in range(1, 21)})
    assert script.refusal(chain, chain.numbers, list(range(1, 15))) == \
        "more than half its chain"
    assert script.refusal(chain, chain.numbers, list(range(1, 10))) is None


def test_every_refusal_reason_is_reachable_and_agrees_with_the_window():
    long_chain = script.Chain({n: "A" for n in range(1, 400)})
    gapped = {n: "A" for n in range(1, 400)}
    del gapped[50]
    cases = {
        "no positions": (long_chain, long_chain.numbers, []),
        "under the floor": (long_chain, long_chain.numbers, [1, 2, 3]),
        "position not in the model": (long_chain, gapped, [10, 20, 30, 40, 50, 60, 70, 80]),
        "chain not in the sequence table": (None, None, [1, 2]),
    }
    for expected, (chain, numbers, positions) in cases.items():
        assert script.refusal(chain, numbers, positions) == expected


# ------------------------------------------------------------------ the null

def test_the_subset_null_is_centred_on_the_chain_s_own_composition():
    """A random same-size subset has the chain's composition in expectation."""
    letters = ("RK" * 20) + ("A" * 160)              # exactly 20% basic
    chain = script.Chain({n: letters[n - 1] for n in range(1, len(letters) + 1)})
    null = chain.null(20, ["basic_fraction"], draws=4000, seed=1)
    mean, deviation = null["basic_fraction"]
    assert mean == pytest.approx(0.2, abs=0.01)
    assert deviation > 0.0


def test_the_null_is_cached_per_size_and_not_per_lining():
    """What makes eight thousand linings cost seconds rather than hours."""
    chain = script.Chain({n: ALPHABET[n % 20] for n in range(1, 400)})
    first = chain.null(20, ["basic_fraction"], draws=200, seed=1)
    second = chain.null(20, ["basic_fraction"], draws=200, seed=99)
    assert first is second                       # the seed could not even be reread


def test_a_basic_lining_scores_positively_against_its_own_chain():
    letters = "A" * 300
    numbers = {n: letters[n - 1] for n in range(1, 301)}
    for n in range(100, 120):
        numbers[n] = "R"
    chain = script.Chain(numbers)
    scored = script.score("R" * 10, chain, ["basic_fraction"], draws=2000, seed=1)
    assert scored["basic_fraction"] > 3.0


def test_the_vectorised_null_matches_the_scalar_statistics():
    """The null evaluates statistics in array form; it must be the same function."""
    rng = np.random.default_rng(0)
    letters = "".join(rng.choice(ALPHABET, size=300))
    chain = script.Chain({n: letters[n - 1] for n in range(1, 301)})
    drawn = np.array([cp.letter_counts(letters[i:i + 20]) for i in range(0, 200, 20)])
    statistics = sorted(cp.COMPOSITION_STATISTICS)
    vectorised = script.counts_statistics(drawn, statistics, chain.background)
    for statistic in statistics:
        scalar = []
        for row in drawn:
            window = "".join(letter * int(count)
                             for letter, count in zip(cp.ALPHABET, row))
            scalar.append(cp.evaluate(statistic, window, chain.background))
        assert vectorised[statistic][0] == pytest.approx(float(np.mean(scalar)), rel=1e-9)


# -------------------------------------------------------------------- arms

def test_the_rank_split_selects_what_study_q9_s_rule_selects():
    """The two halves of the decomposition must mean the same thing by 'formation'."""
    rows = [{"id": f"P{p}", "group": f"c{p}", "rank": r, "length": 10 + (r % 5),
             "window": "A" * (10 + (r % 5))}
            for p in range(6) for r in range(1, 13)]
    mine = script.rank_split([dict(r) for r in rows], np.random.default_rng(script.SEED))
    theirs = [dict(r) for r in rows]
    ps.formation_split(theirs, n_bootstrap=10, statistics=[], shuffles=5,
                       arm="formation")
    chosen_by_them = {(r["id"], r["rank"]) for r in theirs if r.get("scored")}
    chosen_by_me = {(r["id"], r["rank"]) for r in mine[0] + mine[1]}
    assert chosen_by_me == chosen_by_them


# -------------------------------------------------------------------- the gate

def _entry(difference, low, high, p, **kwargs):
    base = {"arm": "function", "statistic": "basic_fraction", "runnable": True,
            "difference": difference, "interval": [low, high], "p_value": p,
            "is_evidence": True, "exceeds_power": True, "length_confounded": False,
            "n_groups": 16}
    base.update(kwargs)
    return base


def test_the_control_fires_only_in_the_pre_stated_direction():
    assert script.fired(_entry(+2.0, 1.0, 3.0, 0.001))
    assert not script.fired(_entry(-2.0, -3.0, -1.0, 0.001))


def test_the_control_needs_every_floor_and_not_just_a_p_value():
    assert not script.fired(_entry(+2.0, 1.0, 3.0, 0.001, is_evidence=False))
    assert not script.fired(_entry(+2.0, 1.0, 3.0, 0.001, exceeds_power=False))
    assert not script.fired(_entry(+2.0, -1.0, 3.0, 0.001))      # spans zero
    assert not script.fired(_entry(+2.0, 1.0, 3.0, 0.9))         # p
    assert not script.fired({"runnable": False})


def test_a_shut_gate_forbids_every_finding_however_strong_the_row():
    """The plan's central commitment, as a fence rather than a promise."""
    overwhelming = _entry(+9.0, 8.0, 10.0, 1e-12)
    script.holm_over([overwhelming], ["function:basic_fraction"], gate=False)
    assert overwhelming["is_finding"] is False
    assert overwhelming["p_holm"] < 0.05           # it would otherwise have fired

    again = _entry(+9.0, 8.0, 10.0, 1e-12)
    script.holm_over([again], ["function:basic_fraction"], gate=True)
    assert again["is_finding"] is True


def test_the_control_arm_is_in_no_family_and_can_never_be_a_finding():
    assert script.CONTROL_ARM not in script.CONFIRMATORY_ARMS
    assert script.CONTROL_STATISTIC in script.PRIMARY


def test_the_report_says_the_gate_is_shut_before_anything_else():
    report = {"seed": 1, "subset_draws": 10, "control_fired": False,
              "control_arm": script.CONTROL_ARM,
              "control_statistic": script.CONTROL_STATISTIC,
              "n_pocket_linings": 0, "n_shell_linings": 0,
              "n_deposited_pocket_linings": 0, "window_yield": {},
              "tests": [_entry(-0.1, -0.5, 0.3, 0.8, arm=script.CONTROL_ARM,
                               family="control", is_finding=False,
                               mean_site_z=-0.1, mean_control_z=0.0,
                               n_site=9, minimum_detectable=0.9)]}
    text = script.markdown(report)
    assert "did NOT fire" in text
    assert "nothing below is evidence of absence" in text
    assert "cannot be read as null results either" in text
    # And the gate paragraph comes before the primary table.
    assert text.index("## The gate") < text.index("## Primary")


# -------------------------------------------------- the committed validation

def test_the_committed_validation_passes_every_check():
    path = ROOT / "results" / "composition" / "validation.json"
    if not path.exists():
        pytest.skip("the validation has not been generated")
    report = json.loads(path.read_text())
    assert report["all_checks_pass"] is True
    assert report["n_noise_pools_fired"] == 0
    assert report["planted"]["fired"] is True
    assert abs(report["noise_mean_z"]) < script.NOISE_TOLERANCE


def test_the_committed_validation_page_is_what_its_numbers_render_to():
    page = ROOT / "results" / "composition" / "VALIDATION.md"
    data = ROOT / "results" / "composition" / "validation.json"
    if not (page.exists() and data.exists()):
        pytest.skip("the validation has not been generated")
    assert script.validation_markdown(json.loads(data.read_text())) == page.read_text()


@pytest.mark.slow
def test_the_gate_opens_on_a_planted_signal_and_shuts_on_matched_noise():
    report = script.validate(draws=400, n_bootstrap=400)
    assert report["all_checks_pass"] is True, report["checks"]
    assert report["planted"]["mean_z"] > 3.0
    assert report["n_noise_pools_fired"] == 0


# ------------------------------------------------------ end to end, synthetic

def _synthetic(planted: bool, seed: int = 0, proteins: int = 18):
    rng = np.random.default_rng(seed)
    sequences, pockets, classes, sites, siblings = [], [], [], [], []
    for index in range(proteins):
        stem = f"AF-P{index}-F1-model_v4"
        letters = list(rng.choice(ALPHABET, size=600))
        sequences.append({"structure": stem + ".pdb", "stem": stem, "chain": "A",
                          "first": 1, "last": 600, "n_residues": 600,
                          "numbers": ",".join(str(n) for n in range(1, 601)),
                          "letters": "".join(letters)})
        for rank in range(1, 13):
            start = 20 + 45 * rank
            residues = list(range(start, start + 14))
            pockets.append({"uniprot_id": f"P{index}", "model": stem + ".pdb",
                            "pocket_id": rank, "cluster": f"c{index}",
                            "residues": ",".join(str(n) for n in residues)})
            classes.append({"uniprot_id": f"P{index}", "pocket_id": rank,
                            "element": "helix_face" if rank % 3 == 0 else "none"})
        name = f"{index:02d}XY"
        deposited = list(rng.choice(ALPHABET, size=400))
        shell = list(range(100, 114))
        if planted:
            for position in range(8):
                deposited[shell[position] - 1] = "RK"[position % 2]
        sequences.append({"structure": name + ".pdb", "stem": name, "chain": "B",
                          "first": 1, "last": 400, "n_residues": 400,
                          "numbers": ",".join(str(n) for n in range(1, 401)),
                          "letters": "".join(deposited)})
        copy = f"{name}_B_900"
        sites.append({"copy": copy, "pdb_id": name, "chain": "B",
                      "group": f"G:{index // 2}",
                      "residues": ",".join(str(n) for n in shell)})
        for rank in range(1, 9):
            start = 180 + 20 * rank
            siblings.append({"uniprot_id": copy, "model": name + ".pdb",
                             "pocket_id": rank,
                             "residues": ",".join(str(n) for n in
                                                  range(start, start + 14))})
    return (pd.DataFrame(sequences), pd.DataFrame(pockets), pd.DataFrame(classes),
            pd.DataFrame(sites), pd.DataFrame(siblings))


@pytest.mark.slow
def test_the_pipeline_detects_a_planted_compositional_signal():
    report = script.run(*_synthetic(True, seed=3), draws=400, n_bootstrap=400)
    assert report["control_fired"] is True
    findings = [t for t in report["tests"] if t.get("is_finding")]
    assert findings, "the gate opened and a planted signal was not detected"
    assert all(t["arm"] in script.CONFIRMATORY_ARMS for t in findings)


@pytest.mark.slow
def test_the_pipeline_stays_quiet_on_matched_noise():
    report = script.run(*_synthetic(False, seed=7), draws=400, n_bootstrap=400)
    assert report["control_fired"] is False
    assert not [t for t in report["tests"] if t.get("is_finding")]


@pytest.mark.slow
def test_every_arm_the_plan_names_runs_or_says_why_not():
    report = script.run(*_synthetic(True, seed=11), draws=200, n_bootstrap=200)
    arms = {t["arm"] for t in report["tests"]}
    assert {"formation", "type", "function", script.CONTROL_ARM} <= arms
    for test in report["tests"]:
        assert test.get("runnable") or test.get("reason")


def test_the_plan_is_committed_and_the_driver_points_at_it():
    text = (ROOT / "docs" / "COMPOSITION_PLAN.md").read_text()
    assert "positive control" in text
    assert "instrument not demonstrated" in text
    source = (ROOT / "scripts" / "pocket_composition.py").read_text()
    assert "docs/COMPOSITION_PLAN.md" in source
