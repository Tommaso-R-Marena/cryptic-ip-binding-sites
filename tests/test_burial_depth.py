"""Tests for geometric pocket burial-depth measurement and filtering."""

from pathlib import Path

import numpy as np
import pandas as pd

from cryptic_ip.analysis.analyzer import ProteinAnalyzer
from cryptic_ip.analysis.filters import CandidateFilter

# Three isolated (fully exposed) protein atoms plus one buried heteroatom ligand.
MINI_PDB = """HEADER    TEST BURIAL DEPTH
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 50.00           C
ATOM      2  CA  ALA A   2      10.000   0.000   0.000  1.00 50.00           C
ATOM      3  CA  ALA A   3       0.000  10.000   0.000  1.00 50.00           C
HETATM    4  P   IHP A 101       5.000   5.000   5.000  1.00 50.00           P
TER
END
"""


def _analyzer(tmp_path: Path) -> ProteinAnalyzer:
    pdb = tmp_path / "mini.pdb"
    pdb.write_text(MINI_PDB, encoding="utf-8")
    return ProteinAnalyzer(str(pdb), skip_electrostatics=True)


def test_surface_atoms_exclude_ligand_heteroatoms(tmp_path: Path):
    analyzer = _analyzer(tmp_path)
    coords = analyzer._compute_surface_atom_coords()

    # Only the three protein atoms define the surface; the ligand P is excluded.
    assert coords.shape == (3, 3)
    assert not np.any(np.all(np.isclose(coords, [5.0, 5.0, 5.0]), axis=1))


def test_burial_depth_is_distance_to_nearest_surface_atom(tmp_path: Path):
    analyzer = _analyzer(tmp_path)
    analyzer._surface_atom_coords = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])

    # Center sits 5 Å from both surface atoms.
    assert analyzer.calculate_pocket_burial_depth((5.0, 0.0, 0.0)) == 5.0
    # Center coincident with a surface atom has zero depth.
    assert analyzer.calculate_pocket_burial_depth((0.0, 0.0, 0.0)) == 0.0


def test_buried_center_is_deeper_than_surface_center(tmp_path: Path):
    analyzer = _analyzer(tmp_path)
    analyzer._surface_atom_coords = np.array(
        [[20.0, 0.0, 0.0], [-20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, -20.0, 0.0]]
    )

    buried = analyzer.calculate_pocket_burial_depth((0.0, 0.0, 0.0))
    surface = analyzer.calculate_pocket_burial_depth((19.0, 0.0, 0.0))

    assert buried > surface
    assert buried == 20.0


def test_burial_depth_zero_without_surface(tmp_path: Path):
    analyzer = _analyzer(tmp_path)
    analyzer._surface_atom_coords = np.empty((0, 3), dtype=float)
    assert analyzer.calculate_pocket_burial_depth((1.0, 2.0, 3.0)) == 0.0


def test_filter_by_burial_depth_keeps_deep_pockets():
    results = pd.DataFrame(
        {
            "pocket_id": [1, 2, 3],
            "composite_score": [0.8, 0.7, 0.9],
            "burial_depth": [20.0, 5.0, 16.0],
        }
    )
    filtered = CandidateFilter().filter_by_burial_depth(results, min_depth=15.0)
    assert sorted(filtered["pocket_id"].tolist()) == [1, 3]


def test_filter_by_burial_depth_missing_column_is_noop():
    results = pd.DataFrame({"pocket_id": [1, 2], "composite_score": [0.8, 0.7]})
    filtered = CandidateFilter().filter_by_burial_depth(results, min_depth=15.0)
    assert len(filtered) == 2


def test_analyze_pocket_reports_burial_depth(tmp_path: Path):
    analyzer = _analyzer(tmp_path)
    analyzer.pockets = pd.DataFrame(
        [
            {
                "pocket_id": 1,
                "center_x": 5.0,
                "center_y": 5.0,
                "center_z": 5.0,
                "volume": 400.0,
                "mean_local_hydrophobic_density": 12.0,
                "fpocket_residue_ids": "1,2,3",
            }
        ]
    )
    analyzer.calculate_sasa()

    metrics = analyzer.analyze_pocket(1)

    assert "burial_depth" in metrics
    assert metrics["burial_depth"] >= 0.0
