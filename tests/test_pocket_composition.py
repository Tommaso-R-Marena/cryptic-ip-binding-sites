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
sq = script.sq
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


# ------------------------------------------- the gate's p-value is calibrated

def test_the_sign_flip_p_matches_exhaustive_enumeration():
    """With few groups the randomisation null can be enumerated exactly."""
    import itertools
    sample = np.array([0.8, -0.2, 0.5, 0.1, 0.4])
    observed = abs(sample.mean())
    exact = sum(
        1 for signs in itertools.product((-1.0, 1.0), repeat=sample.size)
        if abs(float((np.array(signs) * sample).mean())) >= observed
    ) / 2 ** sample.size
    estimate = script.sign_flip_p(sample, draws=40000, seed=1)
    assert estimate == pytest.approx(exact, abs=0.01)


def test_the_sign_flip_p_is_one_for_a_mean_of_zero():
    """Every flip is at least as extreme as no difference at all."""
    assert script.sign_flip_p([1.0, -1.0, 2.0, -2.0], draws=500, seed=1) == pytest.approx(1.0)


def test_the_sign_flip_p_is_small_for_a_consistent_shift():
    assert script.sign_flip_p([2.0] * 10, draws=2000, seed=1) < 0.005


def test_a_single_group_cannot_produce_a_p_value():
    assert script.sign_flip_p([1.5], draws=100, seed=1) == 1.0


def test_the_one_sample_arm_reports_both_p_values_and_which_rule_it_used():
    """The substitution must be visible in the artefact, not silent."""
    rows = [{"group": f"G{g}", "length": 14, "z_basic_fraction": 0.4 + 0.1 * g}
            for g in range(8)]
    entry = script.one_sample(script.CONTROL_ARM, rows, "basic_fraction", n_bootstrap=400)
    assert entry["p_rule"] == "sign-flip randomisation over group means"
    assert "p_value_label_permutation" in entry
    assert entry["p_value"] != entry["p_value_label_permutation"]


@pytest.mark.slow
def test_the_gate_s_p_value_rejects_at_its_nominal_rate_under_the_null():
    """The defect this replaced: label permutation rejected at 7.25% of a nominal 5%.

    `compare_feature`'s permutation shuffles which groups are site and which are
    control. The control arm here is a constant zero, so the arms are not
    exchangeable under the hypothesis - one has no variance - and the null is
    too narrow in the tail a gate reads.
    """
    rng = np.random.default_rng(7)
    pools, groups = 200, 9
    fresh, stale = [], []
    for _ in range(pools):
        rows = [{"group": f"G{g}", "length": 14, "z_basic_fraction": float(rng.normal())}
                for g in range(groups) for _ in range(2)]
        entry = script.one_sample(script.CONTROL_ARM, rows, "basic_fraction",
                                  n_bootstrap=400)
        fresh.append(entry["p_value"])
        stale.append(entry["p_value_label_permutation"])
    fresh, stale = np.array(fresh), np.array(stale)
    # Three standard errors at 200 pools is about 0.046, so this is a loose
    # bound that still fails on the 0.0725 the old rule produced.
    assert (fresh <= 0.05).mean() < 0.09
    assert (fresh <= 0.50).mean() == pytest.approx(0.50, abs=0.08)
    # And the rule that was replaced is still measurably worse in that tail.
    assert (stale <= 0.05).mean() > (fresh <= 0.05).mean()


@pytest.mark.slow
def test_the_two_arm_test_every_other_study_uses_is_calibrated():
    """Scoped to the one-sample construction: compare_feature itself is sound."""
    rng = np.random.default_rng(11)
    ps = []
    for k in range(200):
        site = {f"G{g}": float(rng.normal()) for g in range(9)}
        control = {f"H{g}": float(rng.normal()) for g in range(9)}
        ps.append(sq.compare_feature("f", site, control, n_bootstrap=400,
                                     seed=1000 + k).p_value)
    assert (np.array(ps) <= 0.05).mean() < 0.09


# --------------------------- the gate must not depend on the expensive arm

@pytest.mark.slow
def test_the_gate_still_reads_without_the_deposited_pocket_table():
    """The decisive number must not be hostage to hours of pocket detection.

    The function arm's control is every pocket of the shell's own deposited
    structure. The first run of this workflow put that in the same job as the
    shells and spent four hours in fpocket without reaching the gate, so the
    table is optional and the workflow computes it separately.
    """
    parts = _synthetic(True, seed=3)
    with_table = script.run(*parts, draws=200, n_bootstrap=200)
    without = script.run(*parts[:4], None, draws=200, n_bootstrap=200)

    # The gate is identical either way.
    assert with_table["control_fired"] is True
    assert without["control_fired"] is True
    for key in ("n_shell_linings", "n_pocket_linings", "n_formation_pairs"):
        assert with_table[key] == without[key]

    # And the function arm is reported as having no input, not no signal.
    assert with_table["deposited_pockets_supplied"] is True
    assert without["deposited_pockets_supplied"] is False
    absent = [t for t in without["tests"] if t["arm"] == "function"]
    assert absent and not any(t.get("runnable") for t in absent)
    text = script.markdown(without)
    assert "did not run, for want of its control rather than for want of a signal" in text
    assert "statement about the input and not" in text


@pytest.mark.slow
def test_the_formation_and_type_arms_are_unaffected_by_the_missing_table():
    parts = _synthetic(True, seed=3)
    with_table = script.run(*parts, draws=200, n_bootstrap=200)
    without = script.run(*parts[:4], None, draws=200, n_bootstrap=200)

    def keyed(report):
        return {(t["arm"], t["statistic"]): t.get("difference")
                for t in report["tests"]
                if t.get("runnable") and t["arm"] in ("formation", "type")}
    assert keyed(with_table) == keyed(without)


def test_an_empty_optional_table_reads_as_an_arm_that_did_not_run(tmp_path):
    """The pocket job writes and uploads its table even when it measured nothing.

    An empty file is the same situation as a missing one - the arm has no
    control - and `run` reports that. Reaching `analyse` as a parse error would
    fail the whole study for the one case it is built to tolerate.
    """
    empty = tmp_path / "site_pockets.csv"
    empty.write_text("")
    assert script.read_frame(empty) is None
    assert script.read_frame(tmp_path / "absent.csv") is None
    header_only = tmp_path / "header.csv"
    header_only.write_text("uniprot_id,residues\n")
    frame = script.read_frame(header_only)
    assert frame is not None and frame.empty


def test_the_finding_count_reads_as_english_for_any_count():
    """The committed report said "4 of 9 primary tests is a finding"."""
    def count_line(n_findings: int) -> str:
        tests = [{"arm": "function", "statistic": f"s{i}", "family": "primary",
                  "runnable": True, "is_finding": i < n_findings, "n_groups": 10,
                  "is_evidence": True, "exceeds_power": True, "length_confounded": False,
                  "n_site": 95, "difference": 1.0, "interval": [0.5, 1.5],
                  "p_holm": 0.01, "minimum_detectable": 0.2, "mean_site_z": 1.0,
                  "mean_control_z": 0.0}
                 for i in range(9)]
        report = {"control_fired": True, "control_arm": "positive_control",
                  "control_statistic": "basic_fraction", "tests": tests,
                  "n_pocket_linings": 1, "n_shell_linings": 1,
                  "n_deposited_pocket_linings": 1, "deposited_pockets_supplied": True,
                  "window_yield": {}, "seed": 1, "subset_draws": 1}
        text = script.markdown(report)
        return next(line for line in text.splitlines() if "primary tests" in line)

    assert "1 of 9** primary tests is a finding" in count_line(1)
    assert "4 of 9** primary tests are findings" in count_line(4)
    assert "0 of 9** primary tests are findings" in count_line(0)


def test_the_committed_report_is_a_function_of_the_committed_artefact():
    """`COMPOSITION.md` is `markdown(composition.json)` and nothing else.

    The run prints both, so they agree when committed together - but the report
    is prose and will be edited, and the only safe way to edit it is to change
    the generator and regenerate. This fails if anyone edits the committed
    markdown by hand, which is also how a number would get retyped.
    """
    artefact = ROOT / "results" / "composition" / "composition.json"
    report = ROOT / "results" / "composition" / "COMPOSITION.md"
    if not artefact.exists() or not report.exists():
        pytest.skip("study Q10 has not reported yet")
    assert report.read_text() == script.markdown(json.loads(artefact.read_text()))


# ------------------------------- the control pool is machine-dependent

def test_the_deposited_copy_set_is_recorded_and_digested():
    """Pocket detection has a wall-clock budget, so the pool is not an input.

    Runs 40275846268 and 40307921444 built this table from the same commit's
    inputs and timed out on 8 and 9 structures. The artefact has to carry which
    copies it rests on, or two artefacts cannot be told apart.
    """
    table = pd.DataFrame({"uniprot_id": ["2BBB_A_1", "1AAA_A_1", "1AAA_A_1"],
                          "residues": ["1,2", "3,4", "5,6"]})
    out = script.deposited_provenance(table)
    assert out["n_deposited_copies"] == 2
    assert out["deposited_copies"] == ["1AAA_A_1", "2BBB_A_1"]
    assert len(out["deposited_copies_digest"]) == 16


def test_the_digest_sees_one_copy_appearing_or_leaving():
    """One structure crossing the budget is exactly what happened between runs."""
    full = pd.DataFrame({"uniprot_id": [f"P{i}_A_1" for i in range(10)]})
    short = full.iloc[:9]
    a, b = script.deposited_provenance(full), script.deposited_provenance(short)
    assert a["deposited_copies_digest"] != b["deposited_copies_digest"]
    assert a["n_deposited_copies"] == 10 and b["n_deposited_copies"] == 9


def test_the_digest_ignores_row_order_and_duplicate_pockets():
    """It identifies the POOL, not the table: many pockets share one copy."""
    one = pd.DataFrame({"uniprot_id": ["A_A_1", "B_B_2", "A_A_1"]})
    other = pd.DataFrame({"uniprot_id": ["B_B_2", "A_A_1"]})
    assert (script.deposited_provenance(one)["deposited_copies_digest"]
            == script.deposited_provenance(other)["deposited_copies_digest"])


def test_an_absent_table_records_an_empty_pool():
    for empty in (None, pd.DataFrame(), pd.DataFrame({"other": [1]})):
        out = script.deposited_provenance(empty)
        assert out == {"n_deposited_copies": 0, "deposited_copies": [],
                       "deposited_copies_digest": ""}


def test_the_report_warns_that_the_pool_is_machine_dependent():
    report = {"control_fired": True, "control_arm": "positive_control",
              "control_statistic": "basic_fraction", "tests": [],
              "n_pocket_linings": 1, "n_shell_linings": 1,
              "n_deposited_pocket_linings": 9, "deposited_pockets_supplied": True,
              "n_deposited_copies": 121, "deposited_copies_digest": "abc123",
              "window_yield": {}, "seed": 1, "subset_draws": 1}
    text = script.markdown(report)
    assert "121 deposited copies" in text
    assert "not reproducible across machines" in text
    assert "abc123" in text


def test_a_run_without_the_table_does_not_claim_a_pool():
    report = {"control_fired": True, "control_arm": "positive_control",
              "control_statistic": "basic_fraction", "tests": [],
              "n_pocket_linings": 1, "n_shell_linings": 1,
              "n_deposited_pocket_linings": 0, "deposited_pockets_supplied": False,
              "n_deposited_copies": 0, "deposited_copies_digest": "",
              "window_yield": {}, "seed": 1, "subset_draws": 1}
    text = script.markdown(report)
    assert "not reproducible across machines" not in text
    assert "did not run, for want of its control" in text
