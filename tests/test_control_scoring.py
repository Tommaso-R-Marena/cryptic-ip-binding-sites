"""Unit tests for burial-aware Phase 1 control scoring."""

from cryptic_ip.validation.control_scoring import (
    burial_component,
    negative_passed,
    positive_passed,
    separation_quality,
    validation_score,
)


class TestBurialComponent:
    def test_cryptic_ligand(self):
        assert burial_component(0.0) == 1.0
        assert burial_component(5.0) == 1.0

    def test_surface_ligand(self):
        assert burial_component(50.0) == 0.0
        assert burial_component(100.0) == 0.0

    def test_intermediate_ligand(self):
        assert burial_component(27.5) == 0.5


class TestValidationScore:
    def test_positive_cryptic_site_scores_high(self):
        score = validation_score("positive", pocket_composite=0.8, ligand_sasa=0.0, basic_residues=5)
        assert score > 0.7

    def test_negative_surface_site_scores_low(self):
        score = validation_score("negative", pocket_composite=0.7, ligand_sasa=99.0, basic_residues=4)
        assert score < 0.35

    def test_positive_surface_site_scores_lower(self):
        score = validation_score("positive", pocket_composite=0.8, ligand_sasa=99.0, basic_residues=5)
        assert score < 0.5

    def test_tier1_separation_direction(self):
        adar2 = validation_score("positive", 0.637, 0.0, 8)
        plcd1 = validation_score("negative", 0.677, 99.1, 4)
        assert adar2 - plcd1 > 0.50


class TestPassCriteria:
    def test_adar2_like_positive_passes(self):
        val = validation_score("positive", 0.75, 0.0, 5)
        assert positive_passed(val, 0.0, 5, 0.75)

    def test_plc_like_negative_passes(self):
        val = validation_score("negative", 0.6, 99.0, 4)
        assert negative_passed(val, 99.0, 0.6)


class TestSeparationQuality:
    def test_clear_tier1_separation(self):
        sep = separation_quality([0.85], [0.25])
        assert sep["tier1_separation"] == 0.60
        assert sep["tier1_gate_passed"] is True
        assert sep["phase1_ready"] is True

    def test_poor_separation(self):
        sep = separation_quality([0.45, 0.50], [0.55, 0.48])
        assert sep["clear_separation"] is False
        assert sep["phase1_ready"] is False
