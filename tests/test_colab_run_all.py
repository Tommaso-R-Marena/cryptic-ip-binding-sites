"""Tests for Colab environment bootstrap."""

from pathlib import Path
import importlib.util


def _load_colab_env():
    spec = importlib.util.spec_from_file_location(
        "colab_env", Path("scripts/colab_env.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_colab_env_bootstrap_returns_python():
    colab_env = _load_colab_env()
    python = colab_env.bootstrap_colab_runtime()
    assert python
    assert colab_env.is_python_executable(python)


def test_colab_run_all_presets():
    spec = importlib.util.spec_from_file_location(
        "colab_run_all", Path("scripts/colab_run_all.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert "quick" in module.PRESETS
    assert module.PRESETS["quick"]["n_proteins"] == 25
    assert module.PRESETS["pilot"]["run_figures"] is True
