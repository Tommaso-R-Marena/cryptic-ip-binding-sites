"""Fast smoke tests for publication package scripts."""

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_run_publication_package_importable():
    spec = importlib.util.spec_from_file_location(
        "run_publication_package",
        ROOT / "scripts" / "run_publication_package.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert hasattr(module, "parse_args")
    assert hasattr(module, "main") or "if __name__" in (ROOT / "scripts" / "run_publication_package.py").read_text()


def test_publication_package_skip_flags(tmp_path):
    if shutil.which("fpocket") is None:
        pytest.skip("fpocket not installed")
    out = tmp_path / "publication_ci"
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "run_publication_package.py"),
        "--output-dir",
        str(out),
        "--skip-dataset-build",
        "--skip-ml-training",
        "--skip-figures",
    ]
    subprocess.run(cmd, check=True, cwd=ROOT)
    assert (out / "RESULTS_SUMMARY.md").exists()
    assert (out / "validation" / "control_benchmark.csv").exists()
