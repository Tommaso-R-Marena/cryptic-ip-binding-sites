"""Pytest configuration and fixtures for test suite."""

import pytest
import os
from pathlib import Path
import tempfile
import shutil


@pytest.fixture(scope="session")
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "test_data"


@pytest.fixture(scope="function")
def temp_output_dir():
    """Create temporary output directory for tests."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)


@pytest.fixture(scope="session")
def adar2_pdb():
    """Path to ADAR2 test structure."""
    return Path(__file__).parent / "test_data" / "AF-P78563-F1-model_v4.pdb"


@pytest.fixture(scope="session")
def adar2_crystal():
    """Path to ADAR2 crystal structure."""
    return Path(__file__).parent / "test_data" / "1ZY7.pdb"


@pytest.fixture(scope="session")
def ph_domain_pdb():
    """Path to PH domain negative control."""
    return Path(__file__).parent / "test_data" / "1MAI.pdb"


@pytest.fixture
def mock_proteome_structures(temp_output_dir):
    """Create mock proteome structure files for testing."""
    structures = []
    for i in range(5):
        pdb_file = temp_output_dir / f"test_protein_{i}.pdb"
        # Create minimal PDB file
        with open(pdb_file, 'w') as f:
            f.write("HEADER    TEST PROTEIN\n")
            f.write(f"ATOM      1  CA  ALA A   1      {i}.000   {i}.000   {i}.000  1.00 90.00           C\n")
            f.write("END\n")
        structures.append(pdb_file)
    return structures
