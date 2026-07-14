"""Tests for high-level AnalysisPipeline and ScreeningPipeline wrappers."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from cryptic_ip.pipeline import AnalysisPipeline, ScreeningPipeline


def _scored_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "pocket_id": 1,
                "composite_score": 0.81,
                "basic_residues": 5,
                "sasa": 3.0,
                "volume": 420.0,
                "plddt_confidence": 88.0,
            },
            {
                "pocket_id": 2,
                "composite_score": 0.55,
                "basic_residues": 2,
                "sasa": 40.0,
                "volume": 900.0,
                "plddt_confidence": 60.0,
            },
        ]
    )


@patch("cryptic_ip.pipeline.ProteinAnalyzer")
def test_analysis_pipeline_returns_ranked_candidates(mock_analyzer_cls, tmp_path: Path):
    structure = tmp_path / "test.pdb"
    structure.write_text("HEADER\nEND\n", encoding="utf-8")

    mock_analyzer = MagicMock()
    mock_analyzer.run_pipeline.return_value = _scored_frame()
    mock_analyzer_cls.return_value = mock_analyzer

    pipeline = AnalysisPipeline(work_dir=tmp_path / "work", score_threshold=0.75)
    result = pipeline.analyze(structure)

    assert result["pockets_detected"] == 2
    assert len(result["candidates"]) == 1
    assert result["top_candidate"]["pocket_id"] == 1
    assert result["top_candidate"]["rank"] == 1


@patch("cryptic_ip.pipeline.AnalysisPipeline.analyze")
def test_screening_pipeline_concatenates_hits(mock_analyze, tmp_path: Path):
    structure_a = tmp_path / "a.pdb"
    structure_b = tmp_path / "b.pdb"
    structure_a.write_text("HEADER\nEND\n", encoding="utf-8")
    structure_b.write_text("HEADER\nEND\n", encoding="utf-8")

    candidate = _scored_frame().iloc[[0]].copy()
    mock_analyze.side_effect = [
        {"structure": str(structure_a), "pockets_detected": 2, "candidates": candidate, "top_candidate": None},
        {"structure": str(structure_b), "pockets_detected": 1, "candidates": candidate, "top_candidate": None},
    ]

    screen = ScreeningPipeline(output_dir=tmp_path / "screen", score_threshold=0.75)
    hits = screen.screen_structures([structure_a, structure_b])

    assert len(hits) == 2
    assert set(hits["structure_path"]) == {str(structure_a), str(structure_b)}
    assert (tmp_path / "screen" / "screen_hits.csv").exists()
