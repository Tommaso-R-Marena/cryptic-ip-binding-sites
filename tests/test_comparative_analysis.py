"""Tests for comparative multi-organism analysis workflows."""

import pandas as pd

from cryptic_ip.analysis.comparative_analysis import ComparativeIPAnalysis


def _mock_hits() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"organism": "S_cerevisiae", "protein_id": "Y1", "is_hit": True},
            {"organism": "S_cerevisiae", "protein_id": "Y2", "is_hit": False},
            {"organism": "S_cerevisiae", "protein_id": "Y3", "is_hit": False},
            {"organism": "H_sapiens", "protein_id": "H1", "is_hit": True},
            {"organism": "H_sapiens", "protein_id": "H2", "is_hit": True},
            {"organism": "H_sapiens", "protein_id": "H3", "is_hit": False},
            {"organism": "D_discoideum", "protein_id": "D1", "is_hit": True},
            {"organism": "D_discoideum", "protein_id": "D2", "is_hit": True},
            {"organism": "D_discoideum", "protein_id": "D3", "is_hit": True},
        ]
    )


def _mock_orthologs() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"orthogroup": "OG1", "organism": "S_cerevisiae", "protein_id": "Y1"},
            {"orthogroup": "OG1", "organism": "H_sapiens", "protein_id": "H1"},
            {"orthogroup": "OG1", "organism": "D_discoideum", "protein_id": "D1"},
            {"orthogroup": "OG2", "organism": "S_cerevisiae", "protein_id": "Y2"},
            {"orthogroup": "OG2", "organism": "H_sapiens", "protein_id": "H2"},
            {"orthogroup": "OG2", "organism": "D_discoideum", "protein_id": "D2"},
            {"orthogroup": "OG3", "organism": "S_cerevisiae", "protein_id": "Y3"},
            {"orthogroup": "OG3", "organism": "H_sapiens", "protein_id": "H3"},
            {"orthogroup": "OG3", "organism": "D_discoideum", "protein_id": "D3"},
        ]
    )


def test_hit_rate_correlation_and_pairwise_tests():
    analysis = ComparativeIPAnalysis()
    hit_rates = analysis.compute_hit_rates(
        _mock_hits(),
        {"S_cerevisiae": 20, "H_sapiens": 25, "D_discoideum": 520},
    )

    assert set(["hit_rate", "ci_low", "ci_high", "ip6_uM"]).issubset(hit_rates.columns)
    assert (hit_rates["ci_low"] <= hit_rates["hit_rate"]).all()
    assert (hit_rates["ci_high"] >= hit_rates["hit_rate"]).all()

    corr = analysis.spearman_ip6_correlation(hit_rates)
    assert corr.loc[0, "n_organisms"] == 3

    pairwise = analysis.pairwise_organism_tests(hit_rates)
    assert len(pairwise) == 3
    assert set(["chi2_p_value", "fisher_p_value", "fisher_fdr_q"]).issubset(pairwise.columns)


def test_ortholog_conservation_summarizes_groups():
    analysis = ComparativeIPAnalysis()
    summary = analysis.ortholog_conservation(_mock_hits(), _mock_orthologs())

    assert "orthogroup" in summary.columns
    assert "organisms_with_hits" in summary.columns
    assert summary["orthogroup"].nunique() == 3
    assert (summary["organisms_present"] == 3).all()
