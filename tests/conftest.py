"""Pytest configuration and fixtures."""

import pytest
import sys
from pathlib import Path

# Add package to path
package_dir = Path(__file__).parent.parent
sys.path.insert(0, str(package_dir))


def resolve_structure_path(pdb_id: str) -> Path:
    """Return the first available local path for a validation PDB."""
    for candidate in (
        Path("data/validation") / f"{pdb_id}.pdb",
        Path("tests/data/structures") / f"{pdb_id}.pdb",
    ):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(pdb_id)


@pytest.fixture(scope="session")
def adar2_structure_path() -> Path:
    try:
        return resolve_structure_path("1ZY7")
    except FileNotFoundError:
        pytest.skip("ADAR2 structure not available")


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory):
    """Create temporary data directory for test downloads."""
    return tmp_path_factory.mktemp("test_data")


@pytest.fixture(scope="session")
def alphafold_cache(tmp_path_factory):
    """AlphaFold cache directory."""
    return tmp_path_factory.mktemp("alphafold_cache")


@pytest.fixture(scope="session")
def pdb_cache(tmp_path_factory):
    """PDB cache directory."""
    return tmp_path_factory.mktemp("pdb_cache")


def pytest_configure(config):
    """Configure pytest."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "requires_network: marks tests that require internet"
    )
