"""Tests for validation module."""

import pytest
from pathlib import Path
from cryptic_ip.validation.validation import ValidationSuite


class TestValidationSuite:
    """Test suite for ValidationSuite class."""

    def test_initialization(self):
        """Test validation suite initialization."""
        validator = ValidationSuite()
        assert validator is not None

    def test_adar2_validation(self, adar2_pdb, adar2_crystal, temp_output_dir):
        """Test ADAR2 structure validation."""
        if not adar2_pdb.exists() or not adar2_crystal.exists():
            pytest.skip("ADAR2 structures not available")

        validator = ValidationSuite()
        results = validator.validate_adar2(
            str(adar2_pdb),
            str(adar2_crystal),
            output_dir=str(temp_output_dir)
        )

        assert results is not None
        assert 'pocket_identified' in results
        assert 'rmsd' in results
        assert 'score' in results

    def test_structure_alignment(self, adar2_pdb, adar2_crystal):
        """Test structure alignment functionality."""
        if not adar2_pdb.exists() or not adar2_crystal.exists():
            pytest.skip("ADAR2 structures not available")

        validator = ValidationSuite()
        rmsd = validator.calculate_rmsd(
            str(adar2_pdb),
            str(adar2_crystal)
        )

        assert rmsd is not None
        assert rmsd >= 0
        # AlphaFold should be reasonably close to crystal
        assert rmsd < 5.0, f"RMSD too high: {rmsd}"

    def test_plddt_filtering(self, adar2_pdb):
        """Test pLDDT confidence score filtering."""
        if not adar2_pdb.exists():
            pytest.skip("ADAR2 structure not available")

        validator = ValidationSuite()
        plddt_scores = validator.extract_plddt_scores(str(adar2_pdb))

        assert plddt_scores is not None
        assert len(plddt_scores) > 0
        # Check scores are in valid range
        for score in plddt_scores:
            assert 0 <= score <= 100

    def test_negative_control_scoring(self, ph_domain_pdb, temp_output_dir):
        """Test that negative controls score appropriately low."""
        if not ph_domain_pdb.exists():
            pytest.skip("PH domain structure not available")

        validator = ValidationSuite()
        results = validator.validate_negative_control(
            str(ph_domain_pdb),
            output_dir=str(temp_output_dir)
        )

        assert results is not None
        if 'score' in results:
            assert results['score'] < 0.5, "Negative control scored too high"
