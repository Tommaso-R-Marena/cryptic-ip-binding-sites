"""Shared fixtures for integration tests."""

from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.analysis.ml_classifier import FEATURE_COLUMNS


@pytest.fixture(scope="session")
def cached_adar2_structure(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Download ADAR2 once when network is available, else provide a fallback stub."""
    cache_dir = tmp_path_factory.mktemp("integration_cache")
    out_path = cache_dir / "AF-P78563-F1-model_v4.pdb"

    if out_path.exists():
        return out_path

    url = "https://alphafold.ebi.ac.uk/files/AF-P78563-F1-model_v4.pdb"
    try:
        urllib.request.urlretrieve(url, out_path)
        return out_path
    except Exception:
        stub = (
            "ATOM      1  N   GLY A   1      10.154  12.345  15.678  1.00 20.00           N\n"
            "ATOM      2  CA  GLY A   1      11.154  12.845  16.078  1.00 20.00           C\n"
            "ATOM      3  C   GLY A   1      11.754  11.745  16.978  1.00 20.00           C\n"
            "TER\nEND\n"
        )
        out_path.write_text(stub, encoding="utf-8")
        return out_path


@pytest.fixture()
def validation_feature_set() -> tuple[pd.DataFrame, np.ndarray]:
    """Return separable features used across integration tests."""
    rows = []
    labels = []

    rng = np.random.default_rng(17)
    for _ in range(12):
        rows.append(
            {
                "pocket_depth": float(rng.normal(17.0, 0.8)),
                "sasa": float(rng.normal(3.0, 0.8)),
                "electrostatic_potential": float(rng.normal(7.0, 0.9)),
                "n_basic_residues": int(rng.integers(4, 8)),
                "pocket_volume": float(rng.normal(610.0, 50.0)),
                "plddt_confidence": float(rng.normal(86.0, 3.0)),
            }
        )
        labels.append(1)

    for _ in range(12):
        rows.append(
            {
                "pocket_depth": float(rng.normal(6.0, 1.0)),
                "sasa": float(rng.normal(58.0, 7.0)),
                "electrostatic_potential": float(rng.normal(2.4, 0.8)),
                "n_basic_residues": int(rng.integers(1, 3)),
                "pocket_volume": float(rng.normal(380.0, 40.0)),
                "plddt_confidence": float(rng.normal(82.0, 4.0)),
            }
        )
        labels.append(0)

    frame = pd.DataFrame(rows)
    return frame.loc[:, FEATURE_COLUMNS], np.asarray(labels)


@pytest.fixture()
def local_adar2_copy(tmp_path: Path, cached_adar2_structure: Path) -> Path:
    """Provide a writable ADAR2 path in a test-local directory."""
    target = tmp_path / cached_adar2_structure.name
    shutil.copyfile(cached_adar2_structure, target)
    return target
