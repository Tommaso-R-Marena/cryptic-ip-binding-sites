"""The scorer must refuse an input that cannot be the quantity its name claims.

This project lost a candidate to a silent mislabelling: the yeast pilot's single
hit scored 0.955 on a burial-depth term that was reading fpocket's mean local
hydrophobic density as a distance. The withdrawal is recorded in the revision
note of `manuscript/RESULTS_DRAFT.md`, and the committed artefact still carries
the row, with the misread value beside the true one. These tests replay it.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from cryptic_ip.analysis.scorer import (
    MAX_PLAUSIBLE_BURIAL_DEPTH,
    PocketScorer,
    ScoringInputError,
    check_burial_depth,
)

ROOT = Path(__file__).resolve().parents[1]
WITHDRAWN = ROOT / "results" / "yeast_pilot" / "yeast_pilot_hits.csv"


# ----------------------------------------------------- the incident, replayed

def test_the_withdrawn_hit_still_carries_both_values():
    """If the artefact is ever corrected, this test should be the thing that notices."""
    if not WITHDRAWN.exists():
        pytest.skip("the yeast pilot artefact is no longer committed")
    row = pd.read_csv(WITHDRAWN).iloc[0]
    assert float(row["depth"]) > MAX_PLAUSIBLE_BURIAL_DEPTH, (
        "the misread value is gone from the artefact; if the pilot was re-run, "
        "this test and the supersession note in results/PROJECT_FINDINGS.md "
        "should be revisited together")
    assert float(row["burial_depth"]) < MAX_PLAUSIBLE_BURIAL_DEPTH


def test_the_scorer_now_refuses_the_value_that_produced_the_withdrawn_hit():
    if not WITHDRAWN.exists():
        pytest.skip("the yeast pilot artefact is no longer committed")
    misread = float(pd.read_csv(WITHDRAWN).iloc[0]["depth"])
    with pytest.raises(ScoringInputError, match="not a distance"):
        PocketScorer().score_depth(misread)


def test_the_true_burial_depth_of_that_pocket_scores_without_complaint():
    """The guard must reject the mislabelled input and accept the real one."""
    if not WITHDRAWN.exists():
        pytest.skip("the yeast pilot artefact is no longer committed")
    true_depth = float(pd.read_csv(WITHDRAWN).iloc[0]["burial_depth"])
    score = PocketScorer().score_depth(true_depth)
    assert 0.0 <= score <= 1.0
    # And it is a low score, which is the point: the pocket is not buried.
    assert score < 0.2


def test_the_ambiguous_depth_alias_is_gone_from_the_fallback_chain():
    """`depth` is fpocket's own descriptor name, not a burial distance."""
    source = (ROOT / "cryptic_ip" / "analysis" / "scorer.py").read_text()
    assert '"burial_depth", "pocket_depth", "depth"' not in source
    assert '"burial_depth", "pocket_depth"' in source


def test_a_frame_with_only_the_ambiguous_column_scores_it_as_missing():
    """Better a missing component than a confident wrong one."""
    frame = pd.DataFrame([{"volume": 500.0, "depth": 32.109, "sasa": 7.3,
                           "basic_residues": 7, "enclosure": 0.9}])
    scores = PocketScorer().score_frame(frame)
    assert len(scores) == 1
    assert 0.0 <= float(scores[0]) <= 1.0


# ------------------------------------------------------------- the two guards

def test_a_depth_past_the_ceiling_is_refused():
    with pytest.raises(ScoringInputError, match="not a distance"):
        check_burial_depth(MAX_PLAUSIBLE_BURIAL_DEPTH + 0.1)


def test_a_depth_at_the_ceiling_is_accepted():
    check_burial_depth(MAX_PLAUSIBLE_BURIAL_DEPTH)


def test_no_relation_between_burial_depth_and_hull_depth_is_enforced():
    """The guard that looked right and was not.

    "The nearest solvent-exposed atom lies inside the convex hull, so the
    distance to it cannot exceed the distance to the hull" is false. Distance to
    the hull is a minimum over directions and the hull surface is dense, so it
    is small; distance to the nearest exposed atom is a minimum over a sparse
    set of atoms, so it can be larger. This repository's own fixtures carry both
    of the pairs below, and an earlier version of this module rejected them.
    """
    for depth, hull in ((2.390, 0.205), (1.119, -1.032), (10.0, 0.0)):
        check_burial_depth(depth)
    source = (ROOT / "cryptic_ip" / "analysis" / "scorer.py").read_text()
    assert "HULL_DEPTH_TOLERANCE" not in source
    assert "Neither bounds the other" in source


def test_plausible_depths_score_exactly_as_before():
    """The guards are a fence, not a change of scale: no committed score moves."""
    scorer = PocketScorer()
    assert scorer.score_depth(2.84) == pytest.approx(0.0160, abs=5e-4)
    assert scorer.score_depth(10.14) == pytest.approx(0.3022, abs=5e-4)
