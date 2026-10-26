"""Tests for the generated project-level findings page.

The page's claim is that no number on it was retyped and that a claim its
artefact stops supporting fails the render rather than standing. Both ends are
attacked: a missing or renamed key must raise, and every fence that guards an
interpretive sentence must be shown to fire when its artefact inverts.
"""

from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "project_findings", ROOT / "scripts" / "project_findings.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["project_findings"] = script
SPEC.loader.exec_module(script)

PAGE = ROOT / "results" / "PROJECT_FINDINGS.md"


@pytest.fixture(scope="module")
def artefacts():
    missing = [str(p) for p in script.SOURCES.values() if not p.exists()]
    if missing:
        pytest.skip(f"artefacts not generated: {', '.join(missing)}")
    return script.load()


# ------------------------------------------------------------ loud about gaps

def test_a_missing_artefact_names_its_path(tmp_path):
    with pytest.raises(FileNotFoundError, match="has not been generated"):
        script.load({"nowhere": tmp_path / "absent.json"}, {})


def test_pick_names_the_whole_path_when_a_key_is_gone():
    with pytest.raises(script.Contradiction, match="rerank -> holm"):
        script.pick({"rerank": {}}, "rerank", "holm", "L3b")


def test_a_renamed_key_fails_the_render_rather_than_blanking_it(artefacts):
    data = copy.deepcopy(artefacts)
    del data["specificity"]["rerank"]["L3b_combined_minus_L1"]
    with pytest.raises(script.Contradiction):
        script.render(data)


def test_a_dataset_without_a_burial_column_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    data["dataset"] = pd.DataFrame({"pdb_id": ["1ZY7"]})
    with pytest.raises(script.Contradiction, match="burial classification"):
        script.render(data)


# -------------------------------------------- the fences guarding the verdict

def test_more_cryptic_examples_force_the_page_to_be_rewritten(artefacts):
    """The limitation the whole page is built around must not change silently."""
    data = copy.deepcopy(artefacts)
    column = ("classification" if "classification" in data["dataset"].columns
              else "burial_class")
    data["dataset"] = pd.DataFrame({column: ["Cryptic"] * 40})
    with pytest.raises(script.Contradiction, match="cryptic structures"):
        script.render(data)


def test_a_full_positive_panel_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    data["controls"] = data["controls"].assign(passed=True)
    with pytest.raises(script.Contradiction, match="every positive control"):
        script.render(data)


def test_the_combination_beating_its_baseline_fails_the_render(artefacts):
    """Section 1's third row is the page's sharpest claim, so it is fenced."""
    data = copy.deepcopy(artefacts)
    data["specificity"]["rerank"]["decisions"]["L3b"] = "supported"
    with pytest.raises(script.Contradiction, match="same-protein baseline"):
        script.render(data)


def test_a_within_protein_finding_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    data["within"]["tests"][0]["is_finding"] = True
    with pytest.raises(script.Contradiction, match="within-protein test"):
        script.render(data)


def test_a_sequence_finding_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    for test in data["order"]["tests"]:
        if test.get("family") == "primary":
            test["is_finding"] = True
            break
    with pytest.raises(script.Contradiction, match="Q9 test"):
        script.render(data)


def test_a_look_that_stopped_reproducing_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    data["order"]["look_reproduction"]["verdict"] = "not reproduced"
    with pytest.raises(script.Contradiction, match="held out"):
        script.render(data)


def test_a_failed_gate_calibration_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    data["gate"]["all_checks_pass"] = False
    data["gate"]["checks"]["planted_pool_opens_the_gate"] = False
    with pytest.raises(script.Contradiction, match="gate calibration"):
        script.render(data)


# ------------------------------------- the study that has not reported yet

def test_an_unreported_study_is_named_rather_than_left_out(artefacts):
    data = copy.deepcopy(artefacts)
    data["composition"] = None
    text = script.render(data)
    assert "has not reported yet" in text
    assert "not yet reported" in text          # and in the artefact table


def test_a_shut_gate_on_real_shells_unreads_the_sequence_nulls(artefacts):
    data = copy.deepcopy(artefacts)
    data["composition"] = {
        "control_fired": False, "control_arm": "positive_control",
        "tests": [{"arm": "positive_control", "statistic": "basic_fraction",
                   "runnable": True, "mean_site_z": -0.1, "p_value": 0.8,
                   "n_groups": 16, "family": "control"}]}
    text = script.render(data)
    assert "did not fire" in text
    assert "are not null results" in text
    assert "instrument not demonstrated" in text


def test_an_open_gate_on_real_shells_firms_up_the_sequence_nulls(artefacts):
    data = copy.deepcopy(artefacts)
    data["composition"] = {
        "control_fired": True, "control_arm": "positive_control",
        "tests": [{"arm": "positive_control", "statistic": "basic_fraction",
                   "runnable": True, "mean_site_z": 4.2, "p_value": 0.0005,
                   "n_groups": 16, "family": "control"},
                  {"arm": "function", "statistic": "basic_fraction", "runnable": True,
                   "family": "primary", "is_finding": False}]}
    text = script.render(data)
    assert "positive control fired" in text
    assert "statements about order and not about" in text


# ------------------------------------- what a finding of that study may claim


def _reported(**overrides) -> dict:
    """A Q10 artefact with the gate open, shaped like the committed one."""
    report = {
        "control_fired": True, "control_arm": "positive_control",
        "tests": [{"arm": "positive_control", "statistic": "basic_fraction",
                   "runnable": True, "mean_site_z": 4.28, "p_value": 0.001,
                   "n_groups": 13, "family": "control"},
                  {"arm": "type", "statistic": "composition_divergence", "runnable": True,
                   "family": "primary", "is_finding": True, "n_groups": 54,
                   "difference": 0.317, "interval": [0.106, 0.534], "p_holm": 0.0135},
                  {"arm": "formation", "statistic": "basic_fraction", "runnable": True,
                   "family": "primary", "is_finding": False},
                  {"arm": "function", "statistic": "basic_fraction", "runnable": False,
                   "family": "primary", "is_finding": False}]}
    report.update(overrides)
    return report


def test_a_finding_is_described_by_the_contrast_that_produced_it(artefacts):
    """The arms hold different things fixed, and the page must not blur them."""
    data = copy.deepcopy(artefacts)
    data["composition"] = _reported()
    text = script.render(data)
    assert "pooled across proteins" in text
    assert "the protein is not held fixed" in text
    assert "`composition_divergence`" in text
    assert "0.0135" in text or "0.01349" in text


def test_the_findings_table_does_not_claim_a_within_protein_control(artefacts):
    """`type` pairs across proteins; only `formation` and `function` hold one fixed."""
    data = copy.deepcopy(artefacts)
    data["composition"] = _reported()
    text = script.render(data)
    section = text.split("## 5.")[1].split("## 6.")[0]
    assert "holds the protein fixed" not in section


def test_an_arm_the_page_cannot_describe_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    report = _reported()
    report["tests"].append({"arm": "solvent", "statistic": "basic_fraction",
                            "runnable": True, "family": "primary", "is_finding": True})
    data["composition"] = report
    with pytest.raises(script.Contradiction, match="cannot describe"):
        script.render(data)


def test_an_arm_that_did_not_run_is_said_to_be_missing_its_control(artefacts):
    data = copy.deepcopy(artefacts)
    data["composition"] = _reported()
    text = script.render(data)
    assert "did not run in the committed artefact" in text
    assert "for want of the control table rather than for want of a signal" in text


def test_every_arm_running_leaves_no_missing_arm_paragraph(artefacts):
    data = copy.deepcopy(artefacts)
    report = _reported()
    for test in report["tests"]:
        test["runnable"] = True
    data["composition"] = report
    text = script.render(data)
    assert "did not run in the committed artefact" not in text


def test_the_feature_table_row_stops_saying_the_study_is_in_flight(artefacts):
    data = copy.deepcopy(artefacts)
    data["composition"] = _reported()
    text = script.render(data)
    assert "study in flight" not in text
    assert "the instrument check **fired** on real shells" in text

    unreported = copy.deepcopy(artefacts)
    unreported["composition"] = None
    assert "study in flight" in script.render(unreported)


def test_the_feature_table_row_reports_a_shut_gate_as_unreadable(artefacts):
    data = copy.deepcopy(artefacts)
    data["composition"] = _reported(control_fired=False)
    text = script.render(data)
    assert "instrument **not** demonstrated" in text


def test_the_finding_count_agrees_between_the_section_and_the_table(artefacts):
    data = copy.deepcopy(artefacts)
    data["composition"] = _reported()
    text = script.render(data)
    assert "**1** is a finding" in text
    assert "1 of 2 primary tests is a finding" in text


# ------------------------------------------------------------- the page itself

def test_the_page_is_not_stale(artefacts):
    """Regenerating from the committed artefacts must reproduce the committed page."""
    if not PAGE.exists():
        pytest.skip("PROJECT_FINDINGS.md has not been generated")
    assert PAGE.read_text() == script.render(artefacts)


def test_the_page_lists_every_artefact_it_reads(artefacts):
    text = script.render(artefacts)
    for path in list(script.SOURCES.values()) + list(script.OPTIONAL.values()):
        assert str(path) in text


def test_the_commands_cover_every_source():
    assert set(script.COMMANDS) == set(script.SOURCES) | set(script.OPTIONAL)


def test_the_page_says_what_is_not_established(artefacts):
    text = script.render(artefacts)
    assert "What is not established" in text
    assert "binds an inositol phosphate" in text
    assert "three untested feature families" in text


def test_both_burial_tables_are_reported_and_not_conflated(artefacts):
    """The error this page was corrected for: 1 of 136 is the benchmark's bound,
    and the census's cryptic copies are the ground-truth arms' bound."""
    text = script.render(artefacts)
    census = script.cryptic_census(artefacts["census"])
    assert census["cryptic"] > 1, "fixture no longer exercises the distinction"
    assert f"{census['cryptic']} cryptic annotated" in text
    assert f"{census['cryptic_groups']} homology groups" in text
    assert f"{census['semi_cryptic']} semi-cryptic" in text
    assert "corrects that to a per-copy count" in text
    # And neither overstatement may come back.
    assert "has one experimentally established cryptic site" not in text
    assert "which is what bounds the benchmark" not in text


def test_the_superseded_artefacts_are_named_before_any_result(artefacts):
    """The publication outputs predate corrections the project already made."""
    text = script.render(artefacts)
    assert text.index("## 0. Which artefacts are superseded") < text.index("## The short answer")
    assert "yeast_pilot_hits.csv" in text
    assert "withdrawn" in text
    # The note is quoted from the draft, not paraphrased into this page.
    assert script.revision_note(artefacts["revision"]) in text


def test_the_screening_filter_rejecting_adar2_is_stated(artefacts):
    """The sharpest thing in the whole project, and it was not on the page before."""
    text = script.render(artefacts)
    assert "rejects ADAR2's site on both" in text
    assert "ligand bound and is not the" in text


def test_a_draft_without_its_revision_note_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    data["revision"] = "# Results\n\nNo revision note here.\n"
    with pytest.raises(script.Contradiction, match="revision note"):
        script.render(data)


def test_a_revision_note_that_stops_withdrawing_the_candidate_fails(artefacts):
    data = copy.deepcopy(artefacts)
    data["revision"] = "> **Revision note.** Everything was fine all along.\n"
    with pytest.raises(script.Contradiction, match="withdraws the yeast candidate"):
        script.render(data)


def test_a_census_without_a_burial_column_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    data["census"] = pd.DataFrame({"pdb_id": ["1ZY7"]})
    with pytest.raises(script.Contradiction, match="burial class"):
        script.render(data)


def test_the_feature_table_says_which_families_were_actually_tried(artefacts):
    """The other error: coevolution and a learned score had been tried."""
    text = script.render(artefacts)
    estimate = artefacts["coevolution"]["C2_matched"]["estimate"]
    assert f"{estimate['point'] * 100:+.3f} pp" in text
    assert "cross-organism IP-concentration coevolution" in text
    assert "**not tested**" in text
    for tried in ("a learned pocket score", "orthogroup-level conservation"):
        assert tried in text


def test_the_coevolution_test_turning_supported_fails_the_render(artefacts):
    data = copy.deepcopy(artefacts)
    data["coevolution"]["C2_matched"]["decision"] = "supported"
    with pytest.raises(script.Contradiction, match="coevolution test"):
        script.render(data)


def test_the_page_points_at_the_older_summary_and_says_which_is_current(artefacts):
    text = script.render(artefacts)
    assert "results/publication/RESULTS_SUMMARY.md" in text
    assert "this page is the current one" in text.lower()


def test_the_headline_numbers_are_the_artefact_s_own(artefacts):
    """Spot-check the spine: nothing in the short answer is typed by hand."""
    text = script.render(artefacts)
    rerank = artefacts["specificity"]["rerank"]
    assert str(rerank["unseen"]) in text
    assert str(rerank["annotated_unseen"]) in text
    assert f"{rerank['L3a_combined']['point']:.3f}" in text
    assert f"{rerank['holm']['L3b']:.3f}" in text
    assert str(len(artefacts["candidates"])) in text


def test_the_older_summary_no_longer_overclaims():
    """The generator was fixed, so a regeneration cannot bring these back."""
    source = (ROOT / "scripts" / "run_publication_package.py").read_text()
    for gone in ("reliably separate",
                 "supports proceeding to proteome-scale",
                 "within the 0.2-0.8% range anticipated",
                 "exceeds the 0.30 discrimination minimum"):
        assert gone not in source, gone
    assert "_burial_caveat" in source
    assert "PROJECT_FINDINGS.md" in source


def test_the_decoy_contrast_is_reported_and_not_called_a_finding(artefacts):
    """Section 2's strongest number, with the reason it is weaker than it looks."""
    text = script.render(artefacts)
    rates = artefacts["architecture"]["rates"]
    assert f"{rates['real_either']:.1%}" in text
    assert f"{rates['decoy_either']:.1%}" in text
    assert str(rates["n_decoys"]) in text
    assert "a random position set is not a protein" in text
    assert "close to tautological" in text


def test_a_per_lining_verdict_is_not_presented_as_a_contrast(artefacts):
    """The page must point at the within-protein control as the test."""
    text = script.render(artefacts)
    two = text.index("## 2.")
    three = text.index("## 3.")
    section = text[two:three]
    assert "why section 3 exists" in section
    assert "null" in section


def test_a_thin_finding_is_flagged_against_the_floor(artefacts):
    """The arm with the best control has the fewest groups; the table hides that."""
    data = copy.deepcopy(artefacts)
    report = _reported()
    for test in report["tests"]:
        test["runnable"] = True
        if test["arm"] == "function":
            test.update({"is_finding": True, "n_groups": 10, "difference": 3.994,
                         "interval": [3.239, 4.747], "p_holm": 0.003})
    data["composition"] = report
    text = script.render(data)
    assert "on 10 groups** against a floor of 5" in text
    assert "`type`, on 54" not in text, "the thinnest arm is the one to name"


def test_a_comfortable_finding_is_not_flagged(artefacts):
    """The caveat has to be a fact about the artefact, not boilerplate."""
    data = copy.deepcopy(artefacts)
    report = _reported()
    for test in report["tests"]:
        test["runnable"] = True
        test["n_groups"] = 54
    data["composition"] = report
    text = script.render(data)
    assert "against a floor of" not in text


def test_the_floor_is_read_from_the_code_that_enforces_it():
    """Not written out again in the page, where it could drift."""
    from cryptic_ip.detect import pu
    from cryptic_ip.detect import sequence as sq
    assert script.MIN_GROUPS == sq.MIN_GROUPS == pu.MIN_GROUPS
    assert sq.MIN_GROUPS == 5
