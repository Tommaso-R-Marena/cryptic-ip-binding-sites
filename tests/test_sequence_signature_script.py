"""Tests for the study Q6 driver.

The report's two gates are what matter here: the power block has to be written
before the differences, and a test only becomes a signature when all five of
the plan's conditions hold. Both are checked by building reports whose
conditions fail one at a time.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "sequence_signature", ROOT / "scripts" / "sequence_signature.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(script)


# ------------------------------------------------------------- the pools

def test_the_planted_pool_has_a_basic_residue_on_a_helical_spacing():
    pool = script.planted_pool(6, 24, 3.6, seed=1)
    assert len(pool) == 6
    for window in pool:
        positions = [i for i, c in enumerate(window) if c == "R"]
        gaps = np.diff(positions)
        assert set(gaps) <= {3, 4}, window


def test_the_planted_pool_puts_no_basic_residue_in_its_filler():
    """So the composition of a window is a function of length and period only."""
    pool = script.planted_pool(8, 30, 3.6, seed=2)
    counts = {sum(c in "RK" for c in window) for window in pool}
    assert len(counts) <= 2, counts


def test_the_shuffle_control_keeps_every_letter():
    pool = script.planted_pool(5, 24, 3.6, seed=3)
    for before, after in zip(pool, script.shuffled_pool(pool, seed=4)):
        assert sorted(before) == sorted(after)


def test_the_composition_pool_reproduces_each_window_count_exactly():
    """A rounded mean fraction leaves a residual the pools can be told apart by."""
    counts = [5, 6, 7, 5]
    pool = script.composition_pool(counts, 20, seed=5)
    assert [sum(c in "RK" for c in window) for window in pool] == counts
    assert all(len(window) == 20 for window in pool)


def test_the_composition_pool_refuses_a_count_that_does_not_fit():
    with pytest.raises(ValueError):
        script.composition_pool([25], 20, seed=5)


def test_the_composition_control_matches_the_planted_pool_letter_for_letter():
    planted = script.planted_pool(6, 24, 3.6, seed=7)
    counts = [sum(c in "RK" for c in w) for w in planted]
    control = script.composition_pool(counts, 24, seed=8)
    assert [sum(c in "RK" for c in w) for w in control] == counts


def test_a_pool_summary_reports_what_it_was_matched_on():
    pool = script.planted_pool(4, 24, 3.6, seed=6)
    summary = script.pool_summary(
        "p", pool, "helix_charge_periodicity", n_shuffles=50)
    assert summary["n_windows"] == 4
    assert summary["length"] == 24
    assert 0.0 <= summary["fraction_firing"] <= 1.0


# ---------------------------------------------------------- the validation

def test_the_validation_finds_a_planted_signature_and_stays_quiet_otherwise():
    """The whole point of the module, in one measurement."""
    report = script.run_validate(n_windows=12, length=24, n_shuffles=150)
    planted, shuffle, basic = report["pools"]
    assert planted["mean_z"] > 2.0
    assert planted["fraction_firing"] > 0.5
    assert abs(shuffle["mean_z"]) < 2.0
    assert abs(basic["mean_z"]) < 2.0
    assert report["all_checks_pass"] is True


def test_the_shuffle_control_is_matched_on_composition_exactly():
    report = script.run_validate(n_windows=10, length=24, n_shuffles=100)
    planted, shuffle, _basic = report["pools"]
    assert planted["mean_entropy"] == pytest.approx(shuffle["mean_entropy"])
    assert planted["basic_fraction"] == pytest.approx(shuffle["basic_fraction"])


def test_the_validation_markdown_names_every_check():
    report = script.run_validate(n_windows=8, length=24, n_shuffles=80)
    text = script.validate_markdown(report)
    for name in report["checks"]:
        assert name.replace("_", " ") in text


# -------------------------------------------------------------- the inputs

def write_inputs(tmp_path, records):
    path = tmp_path / "inputs.jsonl"
    path.write_text("\n".join(json.dumps(r) for r in records) + "\n")
    return path


def a_record(**overrides):
    chain = "MKVARLSTGGRKLAAVRRTGGKLSAARVRKLGGSARVLKTRGASRVLKTGGASRVL"
    record = {
        "id": "p1", "group": "G1", "arm": "site", "sequence": chain,
        "residues": [["A", p] for p in (5, 9, 12, 16, 19, 23, 26, 30, 33, 37)],
        "offset": 1,
    }
    record.update(overrides)
    return record


def test_inputs_load(tmp_path):
    path = write_inputs(tmp_path, [a_record()])
    assert len(script.load_inputs(path)) == 1


def test_comments_and_blank_lines_are_skipped(tmp_path):
    path = tmp_path / "inputs.jsonl"
    path.write_text("# a comment\n\n" + json.dumps(a_record()) + "\n")
    assert len(script.load_inputs(path)) == 1


def test_a_missing_field_stops_the_run(tmp_path):
    record = a_record()
    del record["group"]
    with pytest.raises(SystemExit) as caught:
        script.load_inputs(write_inputs(tmp_path, [record]))
    assert "group" in str(caught.value)


def test_an_unknown_arm_stops_the_run(tmp_path):
    with pytest.raises(SystemExit):
        script.load_inputs(write_inputs(tmp_path, [a_record(arm="decoy")]))


def test_unparseable_json_names_the_line(tmp_path):
    path = tmp_path / "inputs.jsonl"
    path.write_text(json.dumps(a_record()) + "\n{not json\n")
    with pytest.raises(SystemExit) as caught:
        script.load_inputs(path)
    assert ":2:" in str(caught.value)


def test_an_empty_inputs_file_stops_the_run(tmp_path):
    path = tmp_path / "inputs.jsonl"
    path.write_text("\n")
    with pytest.raises(SystemExit):
        script.load_inputs(path)


# ------------------------------------------------------------ describe

def test_a_pocket_is_described_on_its_dominant_chain():
    out = script.describe_one(a_record(), n_shuffles=30)
    assert out["measurable"] is True
    assert out["chain"] == "A"
    assert out["chain_fraction"] == 1.0
    assert "local" in out and "lining" in out and "spacing" in out


def test_an_interface_pocket_is_excluded_rather_than_measured():
    """Its lining is two stretches of two sequences; there is no window."""
    record = a_record(residues=[["A", 5], ["A", 9], ["B", 5], ["B", 9]])
    out = script.describe_one(record, n_shuffles=10)
    assert out["measurable"] is False
    assert "spans chains" in out["exclusion_reason"]


def test_a_chain_keyed_sequence_is_looked_up_by_chain():
    chain_a = "MKVARLSTGGRKLAAVRRTGGKLSAARVRKLGGSARVLKTRGASRVLKTGGASRVL"
    record = a_record(sequence={"A": chain_a, "B": "G" * 60})
    out = script.describe_one(record, n_shuffles=30)
    assert out["measurable"] is True
    assert out["windows"]["lining"] == "".join(
        chain_a[p - 1] for p in (5, 9, 12, 16, 19, 23, 26, 30, 33, 37))


def test_a_bad_offset_is_reported_as_an_exclusion_not_a_crash():
    out = script.describe_one(a_record(offset=500), n_shuffles=10)
    assert out["measurable"] is False
    assert "offset" in out["exclusion_reason"]


def test_a_pocket_with_no_residues_is_excluded():
    out = script.describe_one(a_record(residues=[]), n_shuffles=10)
    assert out["measurable"] is False


# -------------------------------------------------- pulling out the numbers

def site_record(arm, group, *, z=1.0, helical=0.5, lining_length=12,
                local_length=30, measurable=True, truncated=False):
    return {
        "id": f"{arm}-{group}-{z}", "arm": arm, "group": group,
        "measurable": measurable,
        "windows": {"lining_length": lining_length, "local_length": local_length,
                    "truncated": truncated},
        "local": {"nulls": {name: {"z": z} for name in (
            "entropy_rate_k2", "helix_charge_periodicity",
            "strand_charge_periodicity", "dominant_period_power")}},
        "lining": {"nulls": {"lz_normalised": {"z": z}}},
        "spacing": {"fraction_helical": helical},
    }


def test_a_primary_value_comes_from_the_right_place_per_window():
    record = site_record("site", "G1", z=3.0, helical=0.75)
    assert script.value_for(record, "local", "entropy_rate_k2") == 3.0
    assert script.value_for(record, "lining", "lz_normalised") == 3.0
    assert script.value_for(record, "spacing", "fraction_helical") == 0.75


def test_an_unmeasurable_pocket_contributes_nothing():
    record = site_record("site", "G1", measurable=False)
    assert script.value_for(record, "local", "entropy_rate_k2") is None


def test_a_missing_statistic_contributes_nothing():
    assert script.value_for(site_record("site", "G1"), "local", "nope") is None


def test_copies_in_one_group_count_once():
    """Groups, not copies: the standing rule, in the aggregation."""
    sites = [site_record("site", "G1", z=1.0), site_record("site", "G1", z=3.0),
             site_record("site", "G2", z=5.0)]
    means = script.group_means(sites, "local", "entropy_rate_k2", "site")
    assert means == {"G1": 2.0, "G2": 5.0}


def test_the_other_arm_is_not_mixed_in():
    sites = [site_record("site", "G1", z=1.0), site_record("control", "H1", z=9.0)]
    assert script.group_means(sites, "local", "entropy_rate_k2", "site") == {"G1": 1.0}


# ---------------------------------------------------------- the power block

def balanced(n_groups=6, *, site_z=3.0, control_z=0.0, **kwargs):
    return (
        [site_record("site", f"G{i}", z=site_z, **kwargs) for i in range(n_groups)]
        + [site_record("control", f"H{i}", z=control_z, **kwargs) for i in range(n_groups)]
    )


def test_the_power_block_flags_the_evidence_floor():
    power = script.power_block(balanced(3))
    assert all(not row["clears_evidence_floor"] for row in power["tests"])
    power = script.power_block(balanced(6))
    assert all(row["clears_evidence_floor"] for row in power["tests"])


def test_the_power_block_counts_exclusions_by_cause():
    sites = balanced(5) + [
        {"id": "x", "arm": "site", "group": "G9", "measurable": False,
         "exclusion_reason": "lining spans chains (0.50 on A); an interface pocket"},
        {"id": "y", "arm": "site", "group": "G8", "measurable": False,
         "exclusion_reason": "lining spans chains (0.40 on B); an interface pocket"},
    ]
    power = script.power_block(sites)
    assert power["exclusions"] == {"lining spans chains": 2}
    assert power["n_pockets"] == 12
    assert power["n_measurable"] == 10


def test_the_power_block_reports_the_truncation_share():
    """A capped local window hides a length difference from the confound flag."""
    sites = (
        [site_record("site", f"G{i}", truncated=True, local_length=40) for i in range(5)]
        + [site_record("control", f"H{i}", truncated=False, local_length=25)
           for i in range(5)]
    )
    power = script.power_block(sites)
    assert power["windows"]["site"]["truncated_fraction"] == 1.0
    assert power["windows"]["control"]["truncated_fraction"] == 0.0


def test_a_minimum_detectable_difference_is_computed_for_every_test():
    power = script.power_block(balanced(6, site_z=3.0, control_z=0.0))
    assert len(power["tests"]) == 6
    for row in power["tests"]:
        assert row["minimum_detectable"] >= 0.0


# -------------------------------------------------------- length confounding

def test_arms_that_differ_in_window_length_are_flagged():
    sites = (
        [site_record("site", f"G{i}", lining_length=20) for i in range(5)]
        + [site_record("control", f"H{i}", lining_length=12) for i in range(5)]
    )
    assert script.length_confounded(sites, "lining") is True
    assert script.length_confounded(sites, "local") is False


def test_the_spacing_window_has_no_length_to_confound():
    assert script.length_confounded(balanced(5), "spacing") is None


def test_one_empty_arm_cannot_be_checked_for_confounding():
    sites = [site_record("site", f"G{i}") for i in range(5)]
    assert script.length_confounded(sites, "lining") is None


# -------------------------------------------------------------- the report

def jittered(n_groups, arm, mean, spread, **kwargs):
    rng = np.random.default_rng(4)
    return [
        site_record(arm, f"{arm[0].upper()}{i}",
                    z=float(mean + spread * rng.normal()), **kwargs)
        for i in range(n_groups)
    ]


def test_a_clear_difference_with_everything_in_order_is_a_signature():
    sites = jittered(8, "site", 4.0, 0.3) + jittered(8, "control", 0.0, 0.3)
    report = script.run_report(sites, n_bootstrap=300)
    primary = {e["test"]: e for e in report["primary"]}
    assert primary[1]["is_signature"] is True
    assert primary[1]["p_holm"] >= primary[1]["p_value"]


def test_a_difference_below_the_minimum_detectable_effect_is_not_a_signature():
    """Significant, interval excluding zero, and still not a finding."""
    sites = jittered(8, "site", 0.5, 2.0) + jittered(8, "control", 0.0, 2.0)
    report = script.run_report(sites, n_bootstrap=300)
    for entry in report["primary"]:
        if entry.get("runnable") and not entry["exceeds_power"]:
            assert entry["is_signature"] is False


def test_a_stratum_below_the_floor_is_never_a_signature():
    sites = jittered(3, "site", 8.0, 0.1) + jittered(3, "control", 0.0, 0.1)
    report = script.run_report(sites, n_bootstrap=200)
    for entry in report["primary"]:
        if entry.get("runnable"):
            assert entry["is_evidence"] is False
            assert entry["is_signature"] is False
    assert report["n_signatures"] == 0


def test_length_confounded_arms_are_never_a_signature():
    sites = (
        jittered(8, "site", 6.0, 0.2, lining_length=24)
        + jittered(8, "control", 0.0, 0.2, lining_length=10)
    )
    report = script.run_report(sites, n_bootstrap=300)
    lining = [e for e in report["primary"] if e.get("feature") == "lining:lz_normalised"]
    assert lining and lining[0]["length_confounded"] is True
    assert lining[0]["is_signature"] is False


def test_holm_covers_the_whole_primary_family():
    sites = jittered(8, "site", 4.0, 0.3) + jittered(8, "control", 0.0, 0.3)
    report = script.run_report(sites, n_bootstrap=200)
    runnable = [e for e in report["primary"] if e.get("runnable")]
    assert report["holm_family_size"] == len(runnable)
    for entry in runnable:
        assert entry["p_holm"] >= entry["p_value"]


def test_a_test_with_an_empty_arm_is_reported_as_not_runnable():
    sites = jittered(6, "site", 3.0, 0.2)
    report = script.run_report(sites, n_bootstrap=100)
    assert all(e["runnable"] is False for e in report["primary"])
    assert report["n_signatures"] == 0


def test_the_markdown_puts_power_before_the_differences():
    """The plan says the MDE is read first; the artefact's order is how."""
    sites = jittered(6, "site", 3.0, 0.3) + jittered(6, "control", 0.0, 0.3)
    text = script.report_markdown(script.run_report(sites, n_bootstrap=100))
    assert text.index("## Power, before any difference") < text.index("## Primary tests")


def test_the_markdown_names_the_plan_and_the_five_conditions():
    sites = jittered(6, "site", 3.0, 0.3) + jittered(6, "control", 0.0, 0.3)
    text = script.report_markdown(script.run_report(sites, n_bootstrap=100))
    assert "docs/SEQUENCE_PLAN.md" in text
    assert "All five" in text


def test_an_infinite_minimum_detectable_effect_prints_as_inf():
    sites = [site_record("site", "G1"), site_record("control", "H1")]
    text = script.report_markdown(script.run_report(sites, n_bootstrap=50))
    assert "inf" in text


# ------------------------------------------------------------------ cli

def test_the_validate_subcommand_writes_both_files(tmp_path):
    out = tmp_path / "validation.json"
    markdown = tmp_path / "VALIDATION.md"
    code = script.main([
        "validate", "--windows", "8", "--length", "24", "--shuffles", "80",
        "--out", str(out), "--markdown", str(markdown),
    ])
    assert code == 0
    assert json.loads(out.read_text())["all_checks_pass"] is True
    assert "Study Q6 method validation" in markdown.read_text()


def test_the_describe_and_report_subcommands_round_trip(tmp_path):
    records = [
        a_record(id=f"s{i}", group=f"G{i // 2}", arm="site") for i in range(10)
    ] + [
        a_record(id=f"c{i}", group=f"H{i // 2}", arm="control") for i in range(10)
    ]
    inputs = write_inputs(tmp_path, records)
    sites = tmp_path / "sites.json"
    report = tmp_path / "report.json"
    assert script.main([
        "describe", "--inputs", str(inputs), "--shuffles", "40", "--out", str(sites),
    ]) == 0
    assert script.main([
        "report", "--sites", str(sites), "--bootstrap", "100", "--out", str(report),
    ]) == 0
    payload = json.loads(report.read_text())
    assert payload["plan"] == "docs/SEQUENCE_PLAN.md"
    assert len(payload["primary"]) == 6


def test_the_committed_validation_passed_every_check():
    """The artefact in the repo, re-read rather than retyped."""
    path = ROOT / "results" / "sequence" / "validation.json"
    payload = json.loads(path.read_text())
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())


def test_the_committed_provenance_is_what_the_generator_produces():
    """A provenance file with hand-typed numbers proves nothing.

    The commit is read back out of the file rather than from git, because the
    file records the tree it was generated from, not the tree it is tested in.
    """
    import re

    results = ROOT / "results" / "sequence"
    payload = json.loads((results / "validation.json").read_text())
    ci = results / "validation_ci.json"
    confirmation = json.loads(ci.read_text()) if ci.exists() else None
    committed = (results / "PROVENANCE.md").read_text()
    commit = re.search(r"at commit `([^`]+)`", committed)
    assert commit, "the provenance file must record the commit it came from"
    assert script.provenance_markdown(
        payload, commit.group(1), confirmation) == committed


def test_the_ci_confirmation_agrees_with_the_committed_run():
    """CI re-runs the validation with more shuffles on a different machine.

    Its artefact is extracted from the run log by `extract_log_block.py`, whose
    gzip CRC verifies the copy, so this compares two independent runs rather
    than a file with itself.
    """
    results = ROOT / "results" / "sequence"
    ci = results / "validation_ci.json"
    if not ci.exists():
        pytest.skip("no CI artefact committed yet")
    mine = json.loads((results / "validation.json").read_text())
    theirs = json.loads(ci.read_text())
    assert theirs["all_checks_pass"] is True
    assert theirs["n_shuffles"] > mine["n_shuffles"]
    by_name = {p["pool"]: p for p in theirs["pools"]}
    for pool in mine["pools"]:
        other = by_name[pool["pool"]]
        assert other["mean_z"] == pytest.approx(pool["mean_z"], abs=0.5)
        assert other["basic_fraction"] == pytest.approx(pool["basic_fraction"])


def test_the_provenance_subcommand_writes_the_same_file(tmp_path):
    validation = tmp_path / "validation.json"
    code = script.main([
        "validate", "--windows", "6", "--length", "24", "--shuffles", "60",
        "--out", str(validation),
    ])
    assert code == 0
    out = tmp_path / "PROVENANCE.md"
    assert script.main([
        "provenance", "--validation", str(validation), "--commit", "abc1234",
        "--out", str(out),
    ]) == 0
    assert out.read_text() == script.provenance_markdown(
        json.loads(validation.read_text()), "abc1234")


def test_the_provenance_states_what_it_may_not_be_used_for():
    text = (ROOT / "results" / "sequence" / "PROVENANCE.md").read_text()
    assert "may not be used for" in text
    assert "ESTABLISHED_SIGNATURES" in text
    assert "has not run" in text


# ------------------------------------------------------- Q6b, by pocket type

def typed(arm, group, pocket_type, **kwargs):
    record = site_record(arm, group, **kwargs)
    record["pocket_type"] = pocket_type
    return record


def test_q6b_does_not_run_without_a_recorded_pocket_type():
    """The split must come from the inputs, not be computed in the analysis."""
    strata = script.run_strata(balanced(6), n_bootstrap=100)
    assert strata["runnable"] is False
    assert "pocket_type" in strata["reason"]


def test_q6b_lists_the_types_it_found():
    sites = (
        [typed("site", f"G{i}", "cryptic", z=3.0) for i in range(5)]
        + [typed("site", f"H{i}", "open", z=0.0) for i in range(5)]
        + [typed("site", f"I{i}", "groove", z=0.0) for i in range(5)]
    )
    strata = script.run_strata(sites, n_bootstrap=100)
    assert strata["types"] == ["cryptic", "groove", "open"]


def test_q6b_compares_a_type_against_the_rest_of_the_sites_only():
    """Controls do not enter Q6b: it is a question about sites."""
    sites = (
        [typed("site", f"G{i}", "cryptic", z=4.0) for i in range(5)]
        + [typed("site", f"H{i}", "open", z=0.0) for i in range(5)]
        + [typed("site", f"J{i}", "groove", z=0.0) for i in range(5)]
        + [site_record("control", f"C{i}", z=99.0) for i in range(5)]
    )
    inside = script.group_means_by_type(
        sites, "local", "entropy_rate_k2", "cryptic", inside=True)
    rest = script.group_means_by_type(
        sites, "local", "entropy_rate_k2", "cryptic", inside=False)
    assert set(inside) == {f"G{i}" for i in range(5)}
    assert set(rest) == {f"H{i}" for i in range(5)} | {f"J{i}" for i in range(5)}
    assert 99.0 not in rest.values()


def test_two_types_are_tested_once_because_the_second_is_a_mirror():
    """One against rest on two types is the same partition with a flipped sign."""
    sites = (
        [typed("site", f"G{i}", "cryptic", z=3.0) for i in range(6)]
        + [typed("site", f"H{i}", "open", z=0.0) for i in range(6)]
    )
    strata = script.run_strata(sites, n_bootstrap=200)
    assert strata["mirrored"] is True
    assert strata["types_tested"] == ["cryptic"]
    assert strata["holm_family_size"] == len(script.sg.PRIMARY_TESTS)


def test_three_types_are_all_tested():
    sites = (
        [typed("site", f"G{i}", "cryptic", z=3.0) for i in range(6)]
        + [typed("site", f"H{i}", "open", z=0.0) for i in range(6)]
        + [typed("site", f"J{i}", "groove", z=1.0) for i in range(6)]
    )
    strata = script.run_strata(sites, n_bootstrap=200)
    assert strata["mirrored"] is False
    assert strata["types_tested"] == ["cryptic", "groove", "open"]
    assert strata["holm_family_size"] == 3 * len(script.sg.PRIMARY_TESTS)


def test_a_stratum_below_the_floor_gives_no_evidence():
    sites = (
        [typed("site", f"G{i}", "cryptic", z=9.0) for i in range(3)]
        + [typed("site", f"H{i}", "open", z=0.0) for i in range(8)]
        + [typed("site", f"J{i}", "groove", z=0.0) for i in range(8)]
    )
    strata = script.run_strata(sites, n_bootstrap=200)
    cryptic = [r for r in strata["power"] if r["pocket_type"] == "cryptic"]
    assert cryptic and all(not r["clears_evidence_floor"] for r in cryptic)
    for entry in strata["tests"]:
        if entry.get("runnable") and entry["pocket_type"] == "cryptic":
            assert entry["is_signature"] is False


def test_q6b_power_is_reported_for_every_stratum_and_feature():
    sites = (
        [typed("site", f"G{i}", "cryptic", z=3.0) for i in range(6)]
        + [typed("site", f"H{i}", "open", z=0.0) for i in range(6)]
        + [typed("site", f"J{i}", "groove", z=1.0) for i in range(6)]
    )
    strata = script.run_strata(sites, n_bootstrap=100)
    assert len(strata["power"]) == 3 * len(script.sg.PRIMARY_TESTS)
    for row in strata["power"]:
        assert row["minimum_detectable"] >= 0.0


def test_the_q6b_markdown_puts_power_before_the_differences():
    sites = (
        [typed("site", f"G{i}", "cryptic", z=3.0) for i in range(6)]
        + [typed("site", f"H{i}", "open", z=0.0) for i in range(6)]
    )
    text = "\n".join(script.strata_markdown(script.run_strata(sites, n_bootstrap=100)))
    assert text.index("### Power, before any difference") < text.index("### Differences")


def test_the_q6b_markdown_explains_an_unrunnable_study():
    text = "\n".join(script.strata_markdown(
        script.run_strata(balanced(6), n_bootstrap=50)))
    assert "Not run" in text


def test_the_pocket_type_survives_into_the_site_record():
    record = a_record()
    record["pocket_type"] = "cryptic"
    out = script.describe_one(record, n_shuffles=20)
    assert out["pocket_type"] == "cryptic"


def test_a_report_carries_its_q6b_section():
    sites = (
        [typed("site", f"G{i}", "cryptic", z=3.0) for i in range(6)]
        + [typed("site", f"H{i}", "open", z=0.0) for i in range(6)]
        + [site_record("control", f"C{i}", z=0.0) for i in range(6)]
    )
    report = script.run_report(sites, n_bootstrap=100)
    assert report["strata"]["runnable"] is True
    assert "Q6b: by pocket type" in script.report_markdown(report)


# -------------------------------- the committed validation must stay committed

def test_the_committed_validation_still_reproduces_exactly(tmp_path):
    """The one real fence on a change to the periodicity machinery.

    `results/sequence/validation.json` is deterministic, offline and committed,
    so any change to `periodicity_power` and the band scan that moved a number
    would move these. It is the check the vectorised periodogram was accepted
    against: every committed figure reproduced to 3.6e-15, which is summation
    order and nothing else.
    """
    committed = json.loads(
        (ROOT / "results" / "sequence" / "validation.json").read_text())
    out = tmp_path / "validation.json"
    assert script.main(["validate", "--windows", "40", "--length", "24",
                        "--shuffles", str(committed["n_shuffles"]),
                        "--out", str(out)]) == 0
    fresh = json.loads(out.read_text())

    assert fresh["checks"] == committed["checks"]
    assert fresh["all_checks_pass"] is True
    assert fresh["statistic"] == committed["statistic"]
    assert [pool["pool"] for pool in fresh["pools"]] == \
        [pool["pool"] for pool in committed["pools"]]
    for was, is_now in zip(committed["pools"], fresh["pools"]):
        for key in ("mean_raw", "mean_z", "mean_entropy", "basic_fraction",
                    "fraction_firing", "distinct_letters", "n_windows"):
            assert float(is_now[key]) == pytest.approx(float(was[key]), abs=1e-9), \
                f"{was['pool']}: {key} moved"
