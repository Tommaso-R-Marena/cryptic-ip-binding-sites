"""Tests for Colab orchestrator script."""

from pathlib import Path
import importlib.util


def test_colab_run_all_presets():
    spec = importlib.util.spec_from_file_location(
        "colab_run_all", Path("scripts/colab_run_all.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert "quick" in module.PRESETS
    assert "pilot" in module.PRESETS
    assert module.PRESETS["quick"]["n_proteins"] == 25
    assert module.PRESETS["pilot"]["n_proteins"] == 500
    assert module.PRESETS["pilot"]["run_figures"] is True
