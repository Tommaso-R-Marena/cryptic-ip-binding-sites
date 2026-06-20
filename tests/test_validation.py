"""Tests for validation suite."""

import pytest

from cryptic_ip.validation.validation_suite import ValidationSuite


class TestValidationSuite:
    """Smoke tests for the current ValidationSuite API."""

    def test_initialization(self):
        suite = ValidationSuite(data_dir="data/validation", use_electrostatics=False)
        assert suite.data_dir.exists()

    @pytest.mark.slow
    def test_tier1_controls_offline(self):
        """Run tier-1 controls when structures are cached."""
        suite = ValidationSuite(data_dir="data/validation", use_electrostatics=False)
        pos = suite.validate_positive_control("ADAR2", suite.POSITIVE_CONTROLS["ADAR2"])
        neg = suite.validate_negative_control("PLCd1_PH", suite.NEGATIVE_CONTROLS["PLCd1_PH"])
        assert pos["passed"]
        assert neg["passed"]
        assert pos["score"] > neg["score"]
