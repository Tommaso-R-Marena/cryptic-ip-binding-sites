"""End-to-end integration test for the ADAR2 -> ML -> MD -> figures workflow."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from cryptic_ip.analysis.comparative_analysis import ComparativeIPAnalysis
from cryptic_ip.analysis.ml_classifier import CrypticSiteMLClassifier, MLPocketScorer
from cryptic_ip.database.batch_processing import append_results_to_file
from cryptic_ip.validation import adar2 as adar2_module
from cryptic_ip.validation.md_validation import OpenMMMDValidationPipeline

try:
    import sklearn  # noqa: F401

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


@pytest.mark.integration
@pytest.mark.skipif(not SKLEARN_AVAILABLE, reason="scikit-learn not installed")
def test_end_to_end_pipeline(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, validation_feature_set):
    """Download ADAR2 -> analyze -> train ML -> score batch -> run mocked MD -> make figures."""

    fake_downloaded = tmp_path / "AF-P78563-F1-model_v4.pdb"
    fake_downloaded.write_text("ATOM\nEND\n", encoding="utf-8")

    def _fake_urlretrieve(url: str, out_path: Path):
        Path(out_path).write_text(fake_downloaded.read_text(encoding="utf-8"), encoding="utf-8")
        return str(out_path), None

    monkeypatch.setattr(adar2_module.urllib.request, "urlretrieve", _fake_urlretrieve)

    structures = adar2_module.download_adar2_structures(data_dir=tmp_path / "adar2_data")
    assert structures["alphafold"].exists()

    X, y = validation_feature_set
    classifier = CrypticSiteMLClassifier(model_type="random_forest", random_state=8)
    classifier.fit(X, y)
    scorer = MLPocketScorer(classifier)

    batch_candidates = [
        {
            "candidate_id": "cand_high",
            "organism": "human",
            "protein_id": "P11111",
            "structure_path": str(structures["alphafold"]),
            "features": {
                "volume": 640.0,
                "depth": 18.0,
                "sasa": 2.8,
                "basic_count": 6,
                "potential": 7.8,
                "plddt_confidence": 89.0,
            },
        },
        {
            "candidate_id": "cand_mid",
            "organism": "arabidopsis",
            "protein_id": "P33333",
            "structure_path": str(structures["alphafold"]),
            "features": {
                "volume": 520.0,
                "depth": 12.0,
                "sasa": 18.0,
                "basic_count": 3,
                "potential": 4.5,
                "plddt_confidence": 84.0,
            },
        },
        {
            "candidate_id": "cand_low",
            "organism": "yeast",
            "protein_id": "P22222",
            "structure_path": str(structures["alphafold"]),
            "features": {
                "volume": 360.0,
                "depth": 5.0,
                "sasa": 63.0,
                "basic_count": 1,
                "potential": 1.7,
                "plddt_confidence": 78.0,
            },
        },
    ]

    scored_rows = []
    for candidate in batch_candidates:
        feat = candidate["features"]
        score = scorer.calculate_composite_score(**feat)
        scored_rows.append(
            {
                "candidate_id": candidate["candidate_id"],
                "organism": candidate["organism"],
                "protein_id": candidate["protein_id"],
                "structure_path": candidate["structure_path"],
                "composite_score": score,
                "is_hit": score >= 0.6,
                "pocket_residues": "376,519,522",
                "pocket_center": "10.0,11.0,12.0",
            }
        )

    batch_csv = append_results_to_file(scored_rows, tmp_path / "batch" / "screening_results.csv")
    assert batch_csv.exists()

    md_pipeline = OpenMMMDValidationPipeline(output_dir=tmp_path / "md")

    def _fake_run_simulation(structure_path: Path, output_dir: Path):
        topology = output_dir / "prepared_system.pdb"
        trajectory = output_dir / "production.dcd"
        topology.write_text("ATOM\nEND\n", encoding="utf-8")
        trajectory.write_text("DUMMY", encoding="utf-8")
        return {"topology": str(topology), "trajectory": str(trajectory), "log": str(output_dir / "md.log")}

    def _fake_analyze_trajectory(*args, **kwargs):
        return {
            "avg_pocket_sasa_nm2": 0.8,
            "avg_pocket_rmsf_nm": 0.15,
            "avg_pocket_volume_nm3": 1.5,
            "std_pocket_volume_nm3": 0.2,
            "avg_waters_in_pocket": 0.5,
        }

    monkeypatch.setattr(md_pipeline, "run_simulation", _fake_run_simulation)
    monkeypatch.setattr(md_pipeline, "analyze_trajectory", _fake_analyze_trajectory)
    monkeypatch.setattr(md_pipeline, "generate_visualization_scripts", lambda *args, **kwargs: None)

    md_report = md_pipeline.validate_top_candidates(batch_csv, top_n=1)
    assert not md_report.empty
    assert set(md_report["classification"]) == {"stably buried"}
    assert (tmp_path / "md" / "stable_candidates.csv").exists()

    comp = ComparativeIPAnalysis()
    hits_df = pd.DataFrame(scored_rows)[["organism", "protein_id", "is_hit"]]
    hit_rates = comp.compute_hit_rates(hits_df, ip6_map={"human": 45.0, "arabidopsis": 12.0, "yeast": 7.5})
    ortholog_df = pd.DataFrame(
        {
            "orthogroup": ["OG1", "OG1"],
            "organism": ["human", "yeast"],
            "protein_id": ["P11111", "P22222"],
        }
    )
    ortholog_summary = comp.ortholog_conservation(hits_df, ortholog_df)
    go_stub = pd.DataFrame(
        {
            "organism": ["human", "yeast", "arabidopsis"],
            "go_term": ["phosphate binding", "RNA processing", "signal transduction"],
            "fdr_bh": [0.01, 0.02, 0.03],
            "enrichment": ["e", "e", "p"],
            "study_count": [20, 10, 8],
        }
    )

    figures_dir = tmp_path / "figures"
    venn_df = hits_df.merge(ortholog_df, on=["organism", "protein_id"], how="left")
    comp.generate_figures(hit_rates, venn_df, go_stub, output_dir=str(figures_dir))

    assert (figures_dir / "hit_rate_barplot.png").exists()
    assert (figures_dir / "ip6_vs_hit_rate_scatter.png").exists()
    assert (figures_dir / "ortholog_hits_venn.png").exists()
    assert (figures_dir / "functional_category_heatmap.png").exists()
