"""Pytest configuration and fixtures."""

import functools
import os
import socket
import sys
import urllib.request
from pathlib import Path

import pytest

# Add package to path
package_dir = Path(__file__).parent.parent
sys.path.insert(0, str(package_dir))

#: Hosts the pipeline needs in order to fetch primary data.
DATA_SOURCE_HOSTS = (
    "https://files.rcsb.org",
    "https://alphafold.ebi.ac.uk",
)


@functools.lru_cache(maxsize=1)
def network_available(timeout: float = 5.0) -> bool:
    """Return whether the primary structural databases are reachable.

    Many environments - CI sandboxes, air-gapped clusters, and any network
    policy that allowlists only package registries - cannot reach the RCSB or
    AlphaFold. Tests that need those services should skip there rather than fail,
    so a red suite always means a real defect. Set ``CRYPTIC_IP_REQUIRE_NETWORK=1``
    to turn unreachable data sources back into failures, which is what a release
    check should do.

    Args:
        timeout: Per-host connection timeout in seconds.

    Returns:
        ``True`` when at least one data source responds.
    """
    if os.environ.get("CRYPTIC_IP_REQUIRE_NETWORK") == "1":
        return True
    for url in DATA_SOURCE_HOSTS:
        try:
            urllib.request.urlopen(url, timeout=timeout)  # noqa: S310 - fixed URLs
            return True
        except (urllib.error.URLError, socket.timeout, OSError):
            continue
    return False


#: Reusable marker for tests that need live database access.
requires_network = pytest.mark.skipif(
    not network_available(),
    reason=(
        "structural databases unreachable; set CRYPTIC_IP_REQUIRE_NETWORK=1 to fail instead"
    ),
)


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


def pytest_runtest_setup(item):
    """Skip network-marked tests when the structural databases are unreachable.

    The ``requires_network`` marker was previously registered but had no effect,
    so sandboxed runs reported connection errors as test failures. Honouring the
    marker means a red suite always indicates a real defect.
    """
    if item.get_closest_marker("requires_network") and not network_available():
        pytest.skip(
            "structural databases unreachable; "
            "set CRYPTIC_IP_REQUIRE_NETWORK=1 to fail instead of skipping"
        )


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
