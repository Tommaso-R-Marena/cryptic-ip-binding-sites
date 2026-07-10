"""Tests for Colab runtime bootstrap helpers."""

import importlib.util
import os
import sys
from pathlib import Path


def _load_colab_env():
    spec = importlib.util.spec_from_file_location("colab_env", Path("scripts/colab_env.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bootstrap_returns_current_python_when_no_conda(tmp_path, monkeypatch):
    colab_env = _load_colab_env()
    monkeypatch.setattr(colab_env, "CONDA_PREFIXES", (tmp_path / "missing",))

    python = colab_env.bootstrap_colab_runtime()

    assert python == sys.executable


def test_write_and_read_runtime_marker(tmp_path, monkeypatch):
    colab_env = _load_colab_env()
    monkeypatch.setattr(colab_env, "CONDA_PREFIXES", (tmp_path / "missing",))

    marker = colab_env.write_runtime_marker(tmp_path)
    assert marker.exists()
    assert colab_env.read_runtime_python(tmp_path) == sys.executable


def test_bootstrap_prepends_conda_bin_to_path(tmp_path, monkeypatch):
    colab_env = _load_colab_env()
    prefix = tmp_path / "miniforge3"
    bin_dir = prefix / "bin"
    bin_dir.mkdir(parents=True)
    (bin_dir / "fpocket").write_text("", encoding="utf-8")
    (bin_dir / "python").write_text("#!/bin/sh\necho py\n", encoding="utf-8")
    monkeypatch.setattr(colab_env, "CONDA_PREFIXES", (prefix,))
    monkeypatch.setenv("PATH", "/usr/bin")

    python = colab_env.bootstrap_colab_runtime()

    assert python == str(bin_dir / "python")
    assert str(bin_dir) in os.environ["PATH"].split(os.pathsep)
