"""Tests for the generated architecture findings page.

The page's whole claim is that no number on it was retyped: every figure comes
out of a committed artefact. These tests attack that claim from both ends —
a missing or renamed key must raise rather than print a blank, and a figure
that *is* on the page must be traceable to the artefact it came from.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "architecture_findings", ROOT / "scripts" / "architecture_findings.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["architecture_findings"] = script
SPEC.loader.exec_module(script)

PAGE = ROOT / "results" / "architecture" / "FINDINGS.md"


@pytest.fixture(scope="module")
def artefacts():
    missing = [str(p) for p in script.SOURCES.values() if not p.exists()]
    if missing:
        pytest.skip(f"artefacts not generated: {', '.join(missing)}")
    return script.load()


# ------------------------------------------------------- loud about gaps

def test_a_missing_artefact_names_its_path():
    with pytest.raises(FileNotFoundError) as excinfo:
        script.load({"probes": Path("results/architecture/does_not_exist.json")})
    assert "does_not_exist.json" in str(excinfo.value)


def test_pick_names_the_whole_path_when_a_key_is_gone():
    with pytest.raises(script.MissingFigure) as excinfo:
        script.pick({"power": {"shapes": []}}, "power", "transfer", "ideal_rate")
    assert "transfer" in str(excinfo.value)


def test_a_renamed_key_fails_the_render_rather_than_blanking_it(artefacts):
    damaged = copy.deepcopy(artefacts)
    damaged["probes"]["transfer"].pop("realistic_rate")
    with pytest.raises(script.MissingFigure):
        script.render(damaged)


def test_a_dropped_probe_row_fails_the_render(artefacts):
    damaged = copy.deepcopy(artefacts)
    damaged["homology"]["retrieval"] = [
        r for r in damaged["homology"]["retrieval"] if r["metric"] != "architecture"]
    with pytest.raises(script.MissingFigure):
        script.render(damaged)


def test_row_for_refuses_an_ambiguous_match():
    rows = [{"metric": "register"}, {"metric": "register"}]
    with pytest.raises(script.MissingFigure):
        script.row_for(rows, "metric", "register")


def test_the_barcode_identity_failing_is_not_quietly_rendered(artefacts):
    damaged = copy.deepcopy(artefacts)
    damaged["probes"]["identities"]["barcode_identity_holds"] = False
    with pytest.raises(script.MissingFigure):
        script.render(damaged)


def test_a_resolved_strand_face_invalidates_the_blind_spot_bullet(artefacts):
    """The bullet asserts register is unresolved there. If it stops being so, fail."""
    damaged = copy.deepcopy(artefacts)
    for row in damaged["probes"]["panel"]["rows"]:
        if row["lining"] == "strand face, every 2nd":
            row["register_p"] = 0.01
    with pytest.raises(script.MissingFigure):
        script.render(damaged)


# ------------------------------------------------------- the numbers tie back

def test_the_transfer_ceiling_on_the_page_is_the_measured_one(artefacts):
    text = script.render(artefacts)
    transfer = artefacts["probes"]["transfer"]
    for key in ("ideal_rate", "realistic_rate", "worst_rate"):
        assert f"{100.0 * transfer[key]:.0f}%" in text


def test_every_retrieval_row_reaches_the_page(artefacts):
    text = script.render(artefacts)
    for row in artefacts["homology"]["retrieval"]:
        assert f"`{row['metric']}`" in text
        assert f"{row['permutation_p']:.4f}" in text


def test_the_decoy_arm_is_on_the_page_beside_the_shapes(artefacts):
    text = script.render(artefacts)
    assert "decoy (no planted face)" in text
    for row in artefacts["probes"]["power"]["shapes"]:
        assert row["shape"] in text


def test_the_real_and_decoy_rates_are_both_reported(artefacts):
    text = script.render(artefacts)
    rates = artefacts["real"]["rates"]
    for key in ("real_either", "decoy_either", "real_topology", "decoy_topology"):
        assert f"{100.0 * rates[key]:.0f}%" in text


def test_the_paired_table_carries_every_runnable_test(artefacts):
    text = script.render(artefacts)
    runnable = [t for t in artefacts["paired"]["tests"] if t.get("runnable")]
    assert runnable, "the paired artefact has no runnable test"
    for entry in runnable:
        assert f"{entry['p_holm']:.4f}" in text
        assert f"{entry['mean_quantile']:.3f}" in text


def test_the_finding_count_matches_the_artefact(artefacts):
    text = script.render(artefacts)
    runnable = [t for t in artefacts["paired"]["tests"] if t.get("runnable")]
    hits = [t for t in runnable if t["is_finding"]]
    assert f"**{len(hits)} of {len(runnable)}**" in text


def test_the_dilution_edge_is_the_last_row_above_the_threshold(artefacts):
    rows = artefacts["probes"]["power"]["dilution"]
    edge, beyond = script.dilution_edge(rows, "random")
    assert edge["detected_either"] >= 0.9 > beyond["detected_either"]
    assert edge["fraction_of_face"] > beyond["fraction_of_face"]


def test_a_dilution_curve_that_never_crosses_raises(artefacts):
    damaged = copy.deepcopy(artefacts)
    for row in damaged["probes"]["power"]["dilution"]:
        row["detected_either"] = 1.0
    with pytest.raises(script.MissingFigure):
        script.render(damaged)


# ------------------------------------------------------- the committed page

def test_the_page_is_not_stale(artefacts):
    """Regenerating from the committed artefacts must reproduce the committed page."""
    if not PAGE.exists():
        pytest.skip("FINDINGS.md has not been generated")
    assert PAGE.read_text() == script.render(artefacts)


def test_the_page_states_what_it_cannot_conclude(artefacts):
    text = script.render(artefacts)
    assert "Nothing on this page is evidence that a pocket binds" in text
    assert "not evolutionarily related" in text or "evolutionarily related" in text


def test_the_page_lists_every_artefact_it_reads(artefacts):
    text = script.render(artefacts)
    for path in script.SOURCES.values():
        assert str(path) in text


def test_the_commands_cover_every_source():
    assert set(script.COMMANDS) == set(script.SOURCES)


def test_the_page_is_valid_markdown_tables(artefacts):
    """Every table row has the same cell count as its header."""
    text = script.render(artefacts)
    header = None
    for line in text.splitlines():
        if not line.startswith("|"):
            header = None
            continue
        cells = line.count("|")
        if header is None:
            header = cells
        else:
            assert cells == header, f"ragged table row: {line}"


def test_the_artefacts_are_json_as_committed():
    for name, path in script.SOURCES.items():
        if path.exists():
            assert isinstance(json.loads(path.read_text()), dict), name


# ------------------------------------------- the within-protein control

def test_the_page_no_longer_claims_the_control_does_not_exist(artefacts):
    """It was built; a page saying otherwise would be wrong, not merely stale."""
    text = script.render(artefacts)
    assert "no committed table carries residue lists" not in text
    assert "The control that did not exist, built" in text


def test_both_within_protein_arms_reach_the_page(artefacts):
    text = script.render(artefacts)
    for name in ("within", "sites"):
        for test in artefacts[name]["tests"]:
            if test.get("runnable"):
                assert f"{test['mean_quantile']:.3f}" in text
                assert f"{test['minimum_detectable']:.3f}" in text


def test_a_within_protein_finding_would_break_the_section(artefacts):
    """The section asserts nothing fired; if that changes it must be rewritten."""
    damaged = copy.deepcopy(artefacts)
    damaged["sites"]["tests"][0]["is_finding"] = True
    damaged["sites"]["tests"][0]["runnable"] = True
    with pytest.raises(script.MissingFigure):
        script.render(damaged)


def test_the_ground_truth_arm_is_named_as_what_decides_the_other(artefacts):
    text = script.render(artefacts)
    assert "decides how to read the other" in text
    assert "exploratory" in text
