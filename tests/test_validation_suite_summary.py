"""Regression tests for validation suite reporting."""

from cryptic_ip.validation.validation_suite import ValidationSuite


def test_print_summary_handles_missing_tier1_separation(capsys):
    suite = ValidationSuite(use_electrostatics=False)
    summary = {
        "separation_quality": {
            "positive_mean": 0.0,
            "negative_mean": 0.0,
            "separation": 0.0,
            "tier1_separation": None,
            "clear_separation": False,
            "phase1_ready": False,
            "all_positive_passed": False,
            "all_negative_passed": False,
        }
    }

    suite._print_summary(summary)
    captured = capsys.readouterr().out

    assert "Tier-1 separation (ADAR2 vs PLCδ1): n/a" in captured
