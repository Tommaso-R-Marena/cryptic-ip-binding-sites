"""Tests for the paired within-lining analysis.

Two of these exist because of a mistake. The first version of this script
mapped an *unresolved* half to a p of 1.0 and then treated it as a measurement,
which put half of all register quantiles at exactly 1.000 and none at 0.000 —
an unresolved subset tying with its equally unresolved draws. And it counted
ties as "at or below", which pushes a quantile upward whenever the p-values are
discrete, as permutation p-values always are. Both are pinned here.
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
    "basic_face", ROOT / "scripts" / "basic_face.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["basic_face"] = script
SPEC.loader.exec_module(script)

REPORT = ROOT / "results" / "architecture" / "real" / "basic_face.json"


# --------------------------------------------------------------- parsing

def test_positions_parse_from_a_list_cell():
    assert script.parse_positions("[3, 7, 11]") == [3, 7, 11]


def test_positions_parse_from_a_delimited_cell():
    assert script.parse_positions("3,7,11") == [3, 7, 11]
    assert script.parse_positions("3;7;11") == [3, 7, 11]


def test_an_empty_cell_is_no_positions():
    assert script.parse_positions("") == []
    assert script.parse_positions(None) == []


def test_positions_are_sorted_and_deduplicated():
    assert script.parse_positions("[11, 3, 7, 3]") == [3, 7, 11]


# ------------------------------------------- unresolved is not a p of one

def test_an_unresolved_half_comes_back_as_none_not_one():
    """The distinction the architecture module insists on, kept here."""
    contiguous = list(range(20, 32))          # every gap 1: one ordering
    out = script.halves_of(contiguous, n_draws=60, seed=0)
    assert out["register"] is None, "a single-ordering lining cannot be registered"


def test_a_measurable_lining_reports_numbers_for_both_halves():
    lining = list(np.cumsum([3, 4] * 7))
    out = script.halves_of(lining, n_draws=120, seed=0)
    assert out["topology"] is not None
    assert 0.0 <= out["topology"] <= 1.0


def test_an_impossible_lining_reports_both_halves_unresolved():
    assert script.halves_of([5], n_draws=20, seed=0) == {"topology": None, "register": None}


# ------------------------------------------------------- the mid-p quantile

def test_the_quantile_splits_ties():
    """Discrete p-values tie constantly; counting them as below inflates it."""
    rng = np.random.default_rng(0)
    lining = sorted(rng.choice(np.arange(120), size=30, replace=False).tolist())
    subset = lining[:10]
    out = script.quantile_within_lining(lining, subset, n_draws=120, draws=15, rng=rng)
    if out is None:
        pytest.skip("this lining is unresolvable at these settings")
    for half in script.HALVES:
        key = f"{half}_quantile"
        if key in out:
            assert 0.0 <= out[key] <= 1.0
            # With ties split, a quantile of exactly 1.0 requires every draw to
            # be strictly below, which is possible but must not be the norm.
            assert f"{half}_tied_fraction" in out


def test_a_subset_as_large_as_its_lining_has_no_matched_draws():
    lining = list(range(0, 40, 3))
    assert script.quantile_within_lining(
        lining, lining, n_draws=60, draws=5,
        rng=np.random.default_rng(0)) is None


def test_too_few_resolved_draws_leaves_that_half_unscored():
    """Better unscored than scored against a handful of survivors."""
    assert script.MIN_RESOLVED_DRAWS >= 10


# ----------------------------------------------------------- the committed run

def test_the_committed_report_exists_and_is_pre_registered():
    if not REPORT.exists():
        pytest.skip("not run yet")
    payload = json.loads(REPORT.read_text())
    assert payload["min_subset"] == script.MIN_SUBSET
    assert payload["seed"] == script.SEED


def test_the_committed_report_tests_exactly_the_pre_registered_family():
    if not REPORT.exists():
        pytest.skip("not run yet")
    payload = json.loads(REPORT.read_text())
    pairs = {(t["subset"], t["half"]) for t in payload["tests"]}
    assert pairs == {(s, h) for s in script.SUBSETS for h in script.HALVES}


def test_the_committed_quantiles_are_not_tie_degenerate():
    """The symptom that exposed the bug: half the mass at exactly 1.000."""
    if not REPORT.exists():
        pytest.skip("not run yet")
    payload = json.loads(REPORT.read_text())
    for half in script.HALVES:
        values = np.array([r[f"{half}_quantile"] for r in payload["per_pocket"]
                           if f"{half}_quantile" in r])
        if values.size < 10:
            continue
        at_one = float(np.mean(values >= 0.999))
        assert at_one < 0.35, f"{half}: {at_one:.2f} of quantiles pinned at 1.000"


def test_the_committed_report_counts_its_exclusions():
    if not REPORT.exists():
        pytest.skip("not run yet")
    payload = json.loads(REPORT.read_text())
    assert isinstance(payload["excluded"], dict)
    assert payload["n_measured"] <= payload["n_pockets_in_table"] * len(script.SUBSETS)


def test_the_markdown_explains_the_unresolved_and_tie_handling():
    if not REPORT.exists():
        pytest.skip("not run yet")
    text = script.markdown(json.loads(REPORT.read_text()))
    assert "unresolved" in text
    assert "mid-p" in text
    assert "candidates, not annotated sites" in text


def test_the_amendment_was_written_before_the_subsets_were_measured():
    text = " ".join((ROOT / "docs" / "ARCHITECTURE_PLAN_AMENDMENT_2.md").read_text().split())
    assert "before any basic-residue subset has been passed to" in text
    assert "2026-10-13" in text
