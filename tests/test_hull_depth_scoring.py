"""Tests for hull depth as the rule-based scorer's depth measure.

The numbers pinned here are the Phase 1 report's measurements of each control's
site pocket on its AlphaFold model, scored as the proteome screen scores it.
They are the calibration panel, so these tests pin behaviour; whether hull depth
improves the score on data it was not chosen on is measured on the deposited
benchmark by the training workflow, which reports both depth measures.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.analysis.features import FEATURE_NAMES, PocketFeatureExtractor
from cryptic_ip.analysis.scorer import PocketScorer, ScoringParameters
from cryptic_ip.analysis.structure_arrays import load_structure_arrays
from cryptic_ip.testing.synthetic import default_benchmark_specs, write_synthetic_structure

#: (volume, burial depth, hull depth, enclosure, lining SASA, basic, Coulomb kT/e)
AF_SITES = {
    "ADAR2": (1493, 4.53, 17.30, 0.96, 30.7, 7, 10.02),
    "PLCd1_PH": (545, 2.71, 8.16, 0.82, 82.3, 6, 10.03),
    "Btk_PH": (306, 3.76, 7.78, 0.80, 37.1, 6, 11.28),
    "DAPP1_PH": (312, 3.15, 5.74, 0.77, 61.4, 5, 14.67),
    "Grp1_PH": (646, 3.00, 7.90, 0.82, 55.5, 6, 9.69),
}


def _score(scorer, site):
    volume, depth, hull, enclosure, sasa, basic, potential = site
    return scorer.calculate_composite_score(
        volume=volume, depth=depth, hull_depth=hull, enclosure=enclosure,
        sasa=sasa, basic_count=basic, potential=potential,
    )


class TestDepthMeasure:
    def test_hull_depth_is_used_when_available(self):
        scorer = PocketScorer()
        assert scorer.score_depth(depth=2.0, hull_depth=11.5) == pytest.approx(0.5)
        assert scorer.score_depth(depth=2.0, hull_depth=15.0) > 0.88
        assert scorer.score_depth(depth=2.0, hull_depth=8.0) < 0.12

    def test_falls_back_to_burial_depth(self):
        scorer = PocketScorer()
        params = ScoringParameters()
        expected = 1.0 / (1.0 + np.exp(-params.depth_slope * (6.0 - params.depth_midpoint)))
        assert scorer.score_depth(depth=6.0, hull_depth=None) == pytest.approx(expected)
        assert scorer.score_depth(depth=6.0, hull_depth=float("nan")) == pytest.approx(expected)

    def test_burial_measure_ignores_hull_depth(self):
        scorer = PocketScorer(parameters=ScoringParameters(depth_measure="burial"))
        assert scorer.score_depth(depth=6.0, hull_depth=30.0) == pytest.approx(
            scorer.score_depth(depth=6.0)
        )

    def test_score_frame_reads_hull_depth(self):
        frame = pd.DataFrame(
            [{"pocket_volume": 800, "burial_depth": 3.0, "hull_depth": 18.0, "enclosure": 0.95,
              "sasa_mean": 20.0, "n_basic_residues": 6, "coulomb_potential_kt": 8.0}]
        )
        hull = PocketScorer().score_frame(frame)[0]
        burial = PocketScorer(parameters=ScoringParameters(depth_measure="burial")).score_frame(frame)[0]
        assert hull > burial


def test_hull_depth_separates_adar2_from_the_ph_domains_more_widely():
    """On the AlphaFold models: ADAR2 above every PH domain, by a wider margin."""
    hull = PocketScorer()
    burial = PocketScorer(parameters=ScoringParameters(depth_measure="burial"))

    def margin(scorer):
        adar2 = _score(scorer, AF_SITES["ADAR2"])
        return adar2 - max(_score(scorer, s) for name, s in AF_SITES.items() if name != "ADAR2")

    assert margin(burial) > 0  # the narrow margin the Phase 1 report measured
    assert margin(hull) > margin(burial) + 0.1


class TestHullDepthDescriptor:
    def test_is_a_feature(self):
        assert "hull_depth" in FEATURE_NAMES
        assert FEATURE_NAMES[-1] == "hull_depth"  # appended: column order is a model contract

    def test_buried_site_is_deeper_than_surface_site(self, tmp_path):
        specs = default_benchmark_specs(n_buried=1, n_surface=1, n_decoy=0, seed=0)
        depths = {}
        for spec in specs:
            arrays = load_structure_arrays(write_synthetic_structure(spec, tmp_path))
            extractor = PocketFeatureExtractor(arrays, n_points=64)
            ligand = arrays.coords[arrays.is_hetero & ~arrays.is_solvent]
            depths[spec.expected_burial_class] = extractor.extract(1, ligand.mean(axis=0)).features["hull_depth"]
        assert set(depths) == {"cryptic", "surface"}
        assert all(np.isfinite(v) for v in depths.values())
        assert depths["cryptic"] > depths["surface"] + 5.0
