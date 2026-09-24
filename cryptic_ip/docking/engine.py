"""AutoDock Vina runs for the redocking benchmark (docs/REDOCKING_PLAN.md).

Parameters fixed by the plan: exhaustiveness 32, 20 poses, energy range 5
kcal/mol, three seeds for the primary arm; the box is centred on the site and
sized to the ligand's extent plus 8 Å on each side, at least 22 Å per side.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

EXHAUSTIVENESS = 32
N_POSES = 20
ENERGY_RANGE = 5.0
PADDING = 8.0
MIN_SIDE = 22.0
SEEDS = (1, 2, 3)
AD4_SPACING = 0.375


def box_size(ligand_xyz: np.ndarray, padding: float = PADDING, min_side: float = MIN_SIDE) -> List[float]:
    """Ligand extent plus ``padding`` per side, at least ``min_side``, per axis."""
    extent = np.ptp(np.asarray(ligand_xyz, dtype=float), axis=0)
    return [float(max(min_side, e + 2 * padding)) for e in extent]


@dataclass
class DockResult:
    scoring: str
    seed: int
    energies: np.ndarray  # (n_poses, k): column 0 is the total score
    poses_pdbqt: str
    centre: List[float]
    size: List[float]
    seconds: float = 0.0

    @property
    def scores(self) -> np.ndarray:
        return self.energies[:, 0] if len(self.energies) else np.array([])


@dataclass
class Receptor:
    pdbqt: Path
    ad4_maps: Optional[str] = None  # map prefix, when AD4 maps exist


def _vina(scoring: str, seed: int, cpu: int):
    from vina import Vina

    return Vina(sf_name=scoring, cpu=cpu, seed=seed, verbosity=0)


def dock(receptor: Receptor, ligand_pdbqt: str, centre: Sequence[float], size: Sequence[float], *,
         scoring: str = "vina", seed: int = 1, cpu: int = 0, exhaustiveness: int = EXHAUSTIVENESS,
         n_poses: int = N_POSES, energy_range: float = ENERGY_RANGE) -> DockResult:
    """One docking run; returns every pose within ``energy_range`` of the best."""
    import time

    start = time.time()
    v = _vina(scoring, seed, cpu)
    if scoring == "ad4":
        if not receptor.ad4_maps:
            raise RuntimeError("AD4 scoring needs autogrid4 maps")
        v.load_maps(receptor.ad4_maps)
    else:
        v.set_receptor(rigid_pdbqt_filename=str(receptor.pdbqt))
    v.set_ligand_from_string(ligand_pdbqt)
    if scoring != "ad4":
        v.compute_vina_maps(center=[float(c) for c in centre], box_size=[float(s) for s in size])
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    text = v.poses(n_poses=n_poses, energy_range=energy_range)
    energies = np.asarray(v.energies(n_poses=n_poses, energy_range=energy_range), dtype=float)
    return DockResult(scoring, seed, energies, text, [float(c) for c in centre], [float(s) for s in size],
                      time.time() - start)


def score_and_minimise(receptor: Receptor, ligand_pdbqt: str, centre: Sequence[float], size: Sequence[float],
                       *, scoring: str = "vina", cpu: int = 0) -> Tuple[float, float, str]:
    """Score a given pose, then locally minimise it; returns ``(score, minimised_score, minimised_pdbqt)``."""
    v = _vina(scoring, 1, cpu)
    if scoring == "ad4":
        v.load_maps(receptor.ad4_maps)
    else:
        v.set_receptor(rigid_pdbqt_filename=str(receptor.pdbqt))
    v.set_ligand_from_string(ligand_pdbqt)
    if scoring != "ad4":
        v.compute_vina_maps(center=[float(c) for c in centre], box_size=[float(s) for s in size])
    before = float(v.score()[0])
    minimised = v.optimize()
    after = float(minimised[0])
    return before, after, v.poses(n_poses=1) if False else _current_pose(v)


def _current_pose(v) -> str:
    import os
    import tempfile

    handle, path = tempfile.mkstemp(suffix=".pdbqt")
    os.close(handle)
    try:
        v.write_pose(path, overwrite=True)
        return Path(path).read_text()
    finally:
        os.unlink(path)


# ------------------------------------------------------------------ AD4
def ad4_available() -> bool:
    return shutil.which("autogrid4") is not None


def write_gpf(path: Path, receptor_pdbqt: Path, receptor_types: Sequence[str], ligand_types: Sequence[str],
              centre: Sequence[float], size: Sequence[float], spacing: float = AD4_SPACING) -> str:
    """An autogrid4 parameter file; returns the map prefix."""
    npts = [int(np.ceil(s / spacing / 2.0) * 2) for s in size]
    stem = receptor_pdbqt.stem
    lines = [
        f"npts {npts[0]} {npts[1]} {npts[2]}",
        f"gridfld {stem}.maps.fld",
        f"spacing {spacing}",
        "receptor_types " + " ".join(receptor_types),
        "ligand_types " + " ".join(ligand_types),
        f"receptor {receptor_pdbqt.name}",
        f"gridcenter {centre[0]:.3f} {centre[1]:.3f} {centre[2]:.3f}",
        "smooth 0.5",
    ]
    lines += [f"map {stem}.{t}.map" for t in ligand_types]
    lines += [f"elecmap {stem}.e.map", f"dsolvmap {stem}.d.map", "dielectric -0.1465"]
    path.write_text("\n".join(lines) + "\n")
    return str(receptor_pdbqt.parent / stem)


def run_autogrid(gpf: Path, timeout: int = 1800) -> None:
    proc = subprocess.run(["autogrid4", "-p", gpf.name, "-l", gpf.with_suffix(".glg").name],
                          cwd=gpf.parent, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(f"autogrid4 failed: {proc.stderr.strip()[-300:]}")


def pdbqt_types(text: str) -> List[str]:
    """Distinct AutoDock types in a PDBQT string, in first-seen order."""
    seen: Dict[str, None] = {}
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            seen[line[77:79].strip()] = None
    return list(seen)


# ----------------------------------------------------------- analysis
@dataclass
class PoseTable:
    """Per-pose RMSDs and scores for one docking run."""

    scores: List[float]
    rmsd: List[float]
    rmsd_p: List[float]
    centroid_distance: List[float] = field(default_factory=list)


def pose_table(result: DockResult, crystal, site_centroid: Optional[np.ndarray] = None) -> PoseTable:
    from .ligand import poses_from_pdbqt
    from .rmsd import symmetric_rmsd

    poses = poses_from_pdbqt(result.poses_pdbqt)
    scores = [float(s) for s in result.scores[: len(poses)]]
    rmsd, rmsd_p, dist = [], [], []
    for pose in poses:
        rmsd.append(symmetric_rmsd(crystal, pose) if crystal is not None else float("nan"))
        rmsd_p.append(symmetric_rmsd(crystal, pose, phosphorus_only=True) if crystal is not None else float("nan"))
        if site_centroid is not None:
            from rdkit import Chem

            xyz = Chem.RemoveHs(pose).GetConformer().GetPositions()
            dist.append(float(np.linalg.norm(xyz.mean(axis=0) - site_centroid)))
    return PoseTable(scores, rmsd, rmsd_p, dist)


def spearman(scores: Sequence[float], rmsd: Sequence[float]) -> float:
    """Rank correlation between pose score and RMSD within one run (nan with < 3 poses)."""
    from scipy.stats import spearmanr

    if len(scores) < 3 or np.ptp(scores) == 0 or np.ptp(rmsd) == 0:
        return float("nan")
    return float(spearmanr(scores, rmsd).correlation)


def greedy_clusters(rmsd_matrix: np.ndarray, cutoff: float = 2.0) -> List[List[int]]:
    """Leader clustering in the given (score) order: each pose joins the first leader within ``cutoff``."""
    leaders: List[int] = []
    clusters: List[List[int]] = []
    for i in range(len(rmsd_matrix)):
        for k, leader in enumerate(leaders):
            if rmsd_matrix[i, leader] <= cutoff:
                clusters[k].append(i)
                break
        else:
            leaders.append(i)
            clusters.append([i])
    return clusters
