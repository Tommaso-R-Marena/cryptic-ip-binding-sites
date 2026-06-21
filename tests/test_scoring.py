"""Tests for scoring module."""

import pytest

from cryptic_ip.analysis.scorer import PocketScorer


class TestPocketScorer:
    """Tests for threshold-based pocket scoring."""

    def test_cryptic_pocket_scores_high(self):
        scorer = PocketScorer()
        score = scorer.calculate_composite_score(
            volume=650.0,
            depth=18.0,
            sasa=2.0,
            basic_count=6,
            potential=7.5,
        )
        assert score > 0.55

    def test_surface_pocket_scores_lower(self):
        scorer = PocketScorer()
        score = scorer.calculate_composite_score(
            volume=400.0,
            depth=6.0,
            sasa=65.0,
            basic_count=3,
            potential=3.0,
        )
        assert score < 0.45

    def test_separation_between_cryptic_and_surface(self):
        scorer = PocketScorer()
        cryptic = scorer.calculate_composite_score(600, 18, 2, 6, 7.0)
        surface = scorer.calculate_composite_score(400, 6, 65, 3, 3.0)
        assert cryptic - surface > 0.20
