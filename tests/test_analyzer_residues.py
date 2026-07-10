"""Regression tests for pocket residue detection and basic-residue counting."""

from pathlib import Path

import pandas as pd
import pytest

from cryptic_ip.analysis.analyzer import ProteinAnalyzer

MINI_PDB = """HEADER    TEST POCKET RESIDUE DETECTION
ATOM      1  N   ARG A  58       0.000   0.000   0.000  1.00 50.00           N
ATOM      2  CA  ARG A  58       1.000   0.000   0.000  1.00 50.00           C
ATOM      3  N   LYS A  59       2.000   0.000   0.000  1.00 50.00           N
ATOM      4  CA  LYS A  59       3.000   0.000   0.000  1.00 50.00           C
ATOM      5  N   HIS A  60       4.000   0.000   0.000  1.00 50.00           N
ATOM      6  CA  HIS A  60       5.000   0.000   0.000  1.00 50.00           C
ATOM      7  N   ALA A 128      50.000  50.000  50.000  1.00 50.00           N
ATOM      8  CA  ALA A 128      51.000  50.000  50.000  1.00 50.00           C
TER
END
"""


@pytest.fixture
def mini_structure(tmp_path: Path) -> Path:
    pdb = tmp_path / "mini.pdb"
    pdb.write_text(MINI_PDB, encoding="utf-8")
    return pdb


def _attach_pocket(
    analyzer: ProteinAnalyzer,
    *,
    fpocket_ids: str | None,
    center: tuple[float, float, float] = (50.0, 50.0, 50.0),
) -> None:
    row = {
        "pocket_id": 1,
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "volume": 400.0,
        "score": 0.5,
        "depth": 12.0,
    }
    if fpocket_ids is not None:
        row["fpocket_residue_ids"] = fpocket_ids
    analyzer.pockets = pd.DataFrame([row])


def test_get_pocket_residues_prefers_fpocket_ids(mini_structure: Path):
    analyzer = ProteinAnalyzer(str(mini_structure), skip_electrostatics=True)
    _attach_pocket(analyzer, fpocket_ids="58,59,60")

    residues = analyzer.get_pocket_residues(1)

    assert residues == [58, 59, 60]


def test_count_basic_residues_positive_with_fpocket_ids(mini_structure: Path):
    analyzer = ProteinAnalyzer(str(mini_structure), skip_electrostatics=True)
    _attach_pocket(analyzer, fpocket_ids="58,59,60")

    basic_count = analyzer.count_basic_residues(1)

    assert basic_count == 3


def test_centroid_fallback_undercounts_basic_residues(mini_structure: Path):
    """Centroid-only detection can miss lining Arg/Lys when the sphere center is offset."""
    analyzer = ProteinAnalyzer(str(mini_structure), skip_electrostatics=True)
    _attach_pocket(analyzer, fpocket_ids=None, center=(50.0, 50.0, 50.0))

    residues = analyzer.get_pocket_residues(1, distance_cutoff=5.0)
    basic_count = analyzer.count_basic_residues(1)

    assert residues == [128]
    assert basic_count == 0


def test_fpocket_ids_fix_basic_residue_undercount(mini_structure: Path):
    analyzer = ProteinAnalyzer(str(mini_structure), skip_electrostatics=True)
    _attach_pocket(analyzer, fpocket_ids="58,59,60", center=(50.0, 50.0, 50.0))

    assert analyzer.count_basic_residues(1) == 3
