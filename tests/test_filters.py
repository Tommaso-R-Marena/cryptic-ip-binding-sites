"""Tests for cryptic-site candidate filtering and ranking."""

import pandas as pd
import pytest

from cryptic_ip.analysis.filters import CandidateFilter


def _sample_results() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "pocket_id": 1,
                "composite_score": 0.82,
                "basic_residues": 5,
                "sasa": 4.0,
                "volume": 450.0,
                "plddt_confidence": 85.0,
            },
            {
                "pocket_id": 2,
                "composite_score": 0.78,
                "basic_residues": 4,
                "sasa": 8.0,
                "volume": 520.0,
                "plddt_confidence": 72.0,
            },
            {
                "pocket_id": 3,
                "composite_score": 0.70,
                "basic_residues": 3,
                "sasa": 6.0,
                "volume": 400.0,
                "plddt_confidence": 80.0,
            },
            {
                "pocket_id": 4,
                "composite_score": 0.90,
                "basic_residues": 6,
                "sasa": 25.0,
                "volume": 500.0,
                "plddt_confidence": 90.0,
            },
        ]
    )


class TestCandidateFilter:
    def test_filter_by_score_uses_threshold(self):
        filt = CandidateFilter(min_score=0.75)
        kept = filt.filter_by_score(_sample_results())
        assert set(kept["pocket_id"]) == {1, 2, 4}

    def test_filter_by_criteria_enforces_burial_and_volume(self):
        filt = CandidateFilter()
        kept = filt.filter_by_criteria(_sample_results())
        assert set(kept["pocket_id"]) == {1, 2}

    def test_filter_cryptic_candidates_applies_all_gates(self):
        filt = CandidateFilter(min_score=0.75, min_plddt=70.0)
        hits = filt.filter_cryptic_candidates(_sample_results())
        assert list(hits["pocket_id"]) == [1, 2]
        assert "rank" in hits.columns
        assert "classification" in hits.columns
        assert hits.iloc[0]["rank"] == 1
        assert hits.iloc[0]["classification"] == "High confidence cryptic IP site"

    def test_rank_candidates_empty_frame(self):
        filt = CandidateFilter()
        ranked = filt.rank_candidates(pd.DataFrame())
        assert ranked.empty

    def test_get_top_candidates_limits_rows(self):
        filt = CandidateFilter(min_score=0.60)
        top = filt.get_top_candidates(_sample_results(), n=2)
        assert len(top) == 2
        assert top.iloc[0]["composite_score"] >= top.iloc[1]["composite_score"]

    def test_filter_by_confidence_uses_plddt_column(self):
        filt = CandidateFilter(min_plddt=75.0)
        kept = filt.filter_by_confidence(_sample_results())
        assert set(kept["pocket_id"]) == {1, 3, 4}
