"""Integration tests for comparative analysis -> statistical validation -> figure generation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from cryptic_ip.analysis.comparative_analysis import ComparativeIPAnalysis
from cryptic_ip.analysis.statistical_validation import StatisticalValidation
import importlib.util



@pytest.mark.integration
def test_comparative_outputs_feed_statistical_validation_and_figure_generation(tmp_path: Path):
    """Comparative outputs should be consumable by statistical and publication figure modules."""
    hits_df = pd.DataFrame(
        {
            "organism": [
                "human", "human", "human", "human",
                "yeast", "yeast", "yeast", "yeast",
                "arabidopsis", "arabidopsis", "arabidopsis", "arabidopsis",
            ],
            "protein_id": ["H1", "H2", "H3", "H4", "Y1", "Y2", "Y3", "Y4", "A1", "A2", "A3", "A4"],
            "is_hit": [True, True, True, False, True, True, False, False, True, False, False, False],
        }
    )

    ip6_map = {"human": 45.0, "yeast": 9.0, "arabidopsis": 16.0}

    comp = ComparativeIPAnalysis()
    hit_rate_df = comp.compute_hit_rates(hits_df, ip6_map=ip6_map)
    correlation_df = comp.spearman_ip6_correlation(hit_rate_df)
    pairwise_df = comp.pairwise_organism_tests(hit_rate_df)

    power_df = StatisticalValidation.pairwise_hit_rate_power(
        {row.organism: row.hit_rate for row in hit_rate_df.itertuples(index=False)}
    )
    assert not pairwise_df.empty
    assert not power_df.empty

    fig2_inputs = tmp_path / "fig2_inputs"
    fig2_inputs.mkdir(parents=True, exist_ok=True)
    hit_rate_fig = hit_rate_df.rename(columns={"hit_rate": "hit_rate"})
    corr_fig = hit_rate_df[["organism", "ip6_uM", "hit_rate"]].rename(
        columns={"ip6_uM": "ip6_concentration"}
    )
    go_heatmap = pd.DataFrame(
        {
            "term": ["phosphate binding", "RNA binding"],
            "human": [2.1, 1.1],
            "yeast": [0.7, 1.9],
            "arabidopsis": [1.2, 0.8],
        }
    )
    phylo = pd.DataFrame(
        {
            "species": ["human", "yeast", "arabidopsis"],
            "candidate": ["cand1", "cand1", "cand1"],
            "conservation_score": [0.91, 0.42, 0.66],
        }
    )

    hit_path = fig2_inputs / "hit_rates.csv"
    corr_path = fig2_inputs / "correlation.csv"
    go_path = fig2_inputs / "go.csv"
    phylo_path = fig2_inputs / "phylo.csv"
    hit_rate_fig.to_csv(hit_path, index=False)
    corr_fig.to_csv(corr_path, index=False)
    go_heatmap.to_csv(go_path, index=False)
    phylo.to_csv(phylo_path, index=False)

    spec = importlib.util.spec_from_file_location(
        "generate_publication_figures", Path("scripts/generate_publication_figures.py")
    )
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)

    config = {
        "figure2": {
            "hit_rates_csv": str(hit_path),
            "correlation_csv": str(corr_path),
            "go_heatmap_csv": str(go_path),
            "phylo_conservation_csv": str(phylo_path),
        }
    }
    module.Image = None
    module.make_figure2(config, output_dir=tmp_path)

    assert (tmp_path / "Figure2_Comparative_Proteomics.png").exists()
    assert (tmp_path / "Figure2_Comparative_Proteomics.pdf").exists()
    assert not correlation_df.empty
