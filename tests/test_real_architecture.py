"""Tests for the real-lining analysis.

The analysis reads committed real data, so the checks here are that it reads
the right file, that its decoy arm really is matched to each real lining, and
that the caveats which make the result readable are present in the output
rather than only in a commit message.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "real_architecture", ROOT / "scripts" / "real_architecture.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["real_architecture"] = script
SPEC.loader.exec_module(script)

INPUTS = ROOT / "results" / "architecture" / "real" / "architecture_inputs.jsonl"
REPORT = ROOT / "results" / "architecture" / "real" / "real.json"


def a_record(positions=(10, 14, 17, 21, 24, 28, 31, 35, 38, 42)):
    return {"id": "P1", "group": "C1", "label": "False",
            "residues": [["A", int(p)] for p in positions]}


def test_positions_come_back_sorted_and_unique():
    record = a_record((42, 10, 14, 10))
    assert script.positions_of(record) == [10, 14, 42]


def test_a_decoy_matches_its_lining_on_size_and_span():
    """The whole point of the arm: it controls both by construction."""
    record = a_record()
    positions = script.positions_of(record)
    rng = np.random.default_rng(0)
    for decoy in script.decoys_for(record, n_draws=60, rng=rng, count=4):
        assert decoy["measurable"]
        assert decoy["n_positions"] == len(positions)
        assert decoy["span"] == positions[-1] - positions[0]


def test_a_measurement_reports_both_halves_and_their_better():
    out = script.measure(script.positions_of(a_record()), n_draws=80, seed=0)
    assert out["measurable"]
    assert out["either_p"] <= out["topology_p"] + 1e-9
    assert out["either_p"] <= out["register_p"] + 1e-9


def test_an_unmeasurable_lining_says_why_rather_than_scoring_zero():
    out = script.measure([5], n_draws=20, seed=0)
    assert out["measurable"] is False
    assert out["reason"]


def test_the_rate_ignores_unmeasurable_rows():
    rows = [{"measurable": True, "either_p": 0.01},
            {"measurable": False, "reason": "x"},
            {"measurable": True, "either_p": 0.9}]
    assert script.rate(rows) == pytest.approx(0.5)


# ------------------------------------------------------- the committed run

def test_the_committed_inputs_are_what_the_analysis_reads():
    assert INPUTS.exists()
    records = [json.loads(line) for line in INPUTS.read_text().splitlines() if line.strip()]
    assert len(records) == 75
    for record in records:
        assert record["chain_is_implicit"] is True, "AlphaFold models are single-chain"


def test_the_committed_report_separates_real_from_matched_decoys():
    if not REPORT.exists():
        pytest.skip("report not built")
    payload = json.loads(REPORT.read_text())
    rates = payload["rates"]
    for key in ("real_either", "decoy_either", "real_topology", "decoy_topology"):
        assert key in rates
    assert rates["n_decoys"] == payload["n_linings"] * script.DECOYS_PER_LINING


def test_the_committed_report_keeps_unresolved_out_of_none():
    if not REPORT.exists():
        pytest.skip("report not built")
    payload = json.loads(REPORT.read_text())
    assert "verdicts" in payload
    assert set(payload["verdicts"]) <= {"clear", "weak", "none", "unresolved", "unmeasurable"}


def test_the_committed_report_runs_exactly_one_pre_registered_split():
    if not REPORT.exists():
        pytest.skip("report not built")
    payload = json.loads(REPORT.read_text())
    features = payload["stratified"]["features"]
    assert len(features) == 4
    for entry in features:
        if entry.get("runnable"):
            assert entry["arms"] == ["False", "True"]


def test_the_markdown_states_the_control_it_does_not_have():
    """A nearly tautological comparison must say so where it is read."""
    if not REPORT.exists():
        pytest.skip("report not built")
    text = script.markdown(json.loads(REPORT.read_text()))
    assert "tautological" in text
    assert "other pockets from the same proteins" in text
    assert "candidates, not annotated sites" in text


def test_the_amendment_was_written_before_the_read():
    amendment = ROOT / "docs" / "ARCHITECTURE_PLAN_AMENDMENT_1.md"
    # Markdown wraps, so the claim is checked against normalised whitespace.
    text = " ".join(amendment.read_text().split())
    assert "before any real lining has been passed to `describe_architecture`" in text
    assert "candidates, not annotated sites" in text
    assert "2026-10-12" in text
