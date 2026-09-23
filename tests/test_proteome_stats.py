"""Tests for proteome-screen hit calling and the organism comparison."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.analysis.proteome_stats import (
    HitCriteria,
    adjusted_comparison,
    clopper_pearson,
    hit_rates,
    keyword_enrichment,
    known_ip_annotation,
    pairwise_fisher,
    protein_table,
    rank_percentile,
    threshold_sweep,
)


def _pocket(uniprot_id, *, score=0.8, sasa=5.0, basic=6, volume=1000.0, plddt=90.0, organism="yeast"):
    return {
        "organism_key": organism,
        "uniprot_id": uniprot_id,
        "composite_score": score,
        "sasa": sasa,
        "basic_residues": basic,
        "volume": volume,
        "plddt_mean": plddt,
    }


def _screened(ids, organism="yeast", length=300):
    return pd.DataFrame(
        {
            "uniprot_id": ids,
            "organism_key": organism,
            "length": length,
            "mean_plddt": 85.0,
            "fraction_plddt_70": 0.9,
        }
    )


class TestHitCriteria:
    def test_volume_window_admits_the_adar2_cavity(self):
        """ADAR2's real InsP6 pocket measures ~1525 A^3; the old 800 cap rejected it."""
        pockets = pd.DataFrame([_pocket("P78563", volume=1525.0)])
        assert HitCriteria().passes(pockets).all()
        assert not HitCriteria(max_volume=800.0).passes(pockets).any()

    @pytest.mark.parametrize(
        "change, gate",
        [
            ({"score": 0.70}, "score"),
            ({"sasa": 12.0}, "sasa"),
            ({"basic": 3}, "basic"),
            ({"volume": 200.0}, "volume"),
            ({"plddt": 60.0}, "plddt"),
        ],
    )
    def test_each_gate(self, change, gate):
        gates = HitCriteria().gates(pd.DataFrame([_pocket("X", **change)]))
        assert not gates[gate].iloc[0]
        assert gates.drop(columns=gate).all(axis=1).iloc[0]


class TestProteinTable:
    def test_hit_eligibility_and_blocking_gate(self):
        pockets = pd.DataFrame(
            [
                _pocket("HIT"),
                _pocket("HIT", score=0.3),
                _pocket("LOWCONF", plddt=50.0),
                _pocket("NEAR", score=0.70),
                _pocket("NEAR", basic=1, sasa=20.0, score=0.9),
            ]
        )
        table = protein_table(pockets, _screened(["HIT", "LOWCONF", "NEAR", "NOPOCKET"]))
        table = table.set_index("uniprot_id")
        assert table.loc["HIT", "is_hit"] and table.loc["HIT", "blocking_gate"] == "passes"
        assert not table.loc["LOWCONF", "eligible"]
        assert table.loc["LOWCONF", "blocking_gate"] == "plddt"
        # NEAR's closest pocket fails only the score gate.
        assert table.loc["NEAR", "blocking_gate"] == "score"
        assert table.loc["NEAR", "best_confident_score"] == pytest.approx(0.9)
        assert table.loc["NOPOCKET", "blocking_gate"] == "no pocket"
        assert not table.loc["NOPOCKET", "is_hit"]

    def test_proteins_without_pockets_count_in_the_denominator(self):
        table = protein_table(pd.DataFrame([_pocket("A")]), _screened(["A", "B", "C", "D"]))
        rates = hit_rates(table)
        assert rates.loc[0, "screened"] == 4
        assert rates.loc[0, "hits"] == 1
        assert rates.loc[0, "eligible"] == 1
        assert rates.loc[0, "hit_rate_eligible"] == pytest.approx(1.0)


def test_clopper_pearson_matches_the_pilot_intervals():
    """The manuscript quotes these for the 499-protein pilot."""
    low, high = clopper_pearson(1, 499)
    assert low == pytest.approx(0.0000507, rel=0.01)
    assert high == pytest.approx(0.0111, rel=0.01)
    assert clopper_pearson(0, 499) == (0.0, pytest.approx(0.00737, rel=0.01))


def test_threshold_sweep_is_monotone():
    rng = np.random.default_rng(0)
    ids = [f"P{i}" for i in range(200)]
    pockets = pd.DataFrame([_pocket(i, score=float(rng.uniform(0.4, 1.0))) for i in ids])
    sweep = threshold_sweep(pockets, _screened(ids), [0.5, 0.6, 0.7, 0.8, 0.9])
    assert sweep["hits"].is_monotonic_decreasing


def _two_organisms(rate_a, rate_b, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    rows, screened = [], []
    for organism, rate in (("dictyostelium", rate_a), ("yeast", rate_b)):
        for i in range(n):
            uid = f"{organism[:2]}{i}"
            hit = rng.random() < rate
            rows.append(_pocket(uid, organism=organism, score=0.9 if hit else 0.5))
            screened.append(
                {
                    "uniprot_id": uid,
                    "organism_key": organism,
                    "length": int(rng.integers(100, 1500)),
                    "mean_plddt": 80.0,
                    "fraction_plddt_70": float(rng.uniform(0.3, 1.0)),
                }
            )
    return protein_table(pd.DataFrame(rows), pd.DataFrame(screened))


def test_fisher_detects_a_real_difference():
    proteins = _two_organisms(0.05, 0.005)
    fisher = pairwise_fisher(proteins)
    row = fisher[fisher["denominator"] == "screened"].iloc[0]
    assert row["odds_ratio"] > 3
    assert row["p_value"] < 1e-6


def test_adjusted_comparison_recovers_the_odds_ratio():
    pytest.importorskip("statsmodels")
    proteins = _two_organisms(0.05, 0.005)
    result = adjusted_comparison(proteins)
    term = result["hit_logit"]["terms"]["organism=yeast"]
    # yeast relative to Dictyostelium: odds ratio well below one.
    assert term["estimate"] < 0.3
    assert term["ci_high"] < 1.0
    assert "best_score_ols" in result


def test_adjusted_comparison_says_when_it_is_not_estimable():
    pytest.importorskip("statsmodels")
    proteins = _two_organisms(0.0, 0.0)
    result = adjusted_comparison(proteins)
    assert result["hit_logit"]["estimable"] is False


def test_keyword_enrichment():
    proteins = pd.DataFrame({"uniprot_id": [f"P{i}" for i in range(100)]})
    selected = pd.Series([i < 10 for i in range(100)])
    keywords = ["Nucleus; RNA-binding" if i < 9 else ("Nucleus" if i % 5 == 0 else "Membrane") for i in range(100)]
    annotations = pd.DataFrame({"uniprot_id": proteins["uniprot_id"], "keywords": keywords})
    table = keyword_enrichment(proteins, annotations, selected, ["RNA-binding", "Membrane"]).set_index("keyword")
    assert table.loc["RNA-binding", "in_selected"] == 9
    assert table.loc["RNA-binding", "p_value"] < 1e-6
    assert table.loc["Membrane", "odds_ratio"] < 1
    assert (table["p_holm"] >= table["p_value"]).all()


def test_known_ip_annotation():
    annotations = pd.DataFrame(
        {
            "uniprot_id": ["A", "B", "C"],
            "binding_site": [
                'BINDING 376; /ligand="1D-myo-inositol hexakisphosphate"',
                'BINDING 12; /ligand="ATP"',
                "",
            ],
            "function": ["", "", "Requires InsP6 for activity."],
        }
    )
    flags = known_ip_annotation(annotations)
    assert flags.to_dict() == {"A": True, "B": False, "C": True}


def test_rank_percentile():
    proteins = pd.DataFrame(
        {
            "uniprot_id": ["A", "B", "C", "D"],
            "organism_key": ["human"] * 4,
            "best_confident_score": [0.9, 0.5, 0.6, np.nan],
        }
    )
    assert rank_percentile(proteins, "A") == pytest.approx(2 / 3)
    assert rank_percentile(proteins, "B") == pytest.approx(0.0)
    assert rank_percentile(proteins, "D") is None
