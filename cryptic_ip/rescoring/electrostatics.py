"""A screened Coulomb energy between a docked ligand pose and its receptor (docs/RERANK_PLAN.md).

E = 332.0637 * sum_ij q_i q_j exp(-kappa r) / (eps(r) r), with eps(r) = 4 r,
kappa = 0.127 / Å (150 mM monovalent salt, 298 K), pairs within 12 Å, r floored
at 1.5 Å. Charges are the partial charges of the PDBQT files: PDB2PQR/AMBER for
the receptor, Meeko's Gasteiger charges for the ligand. Nothing here is fitted.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np

COULOMB = 332.0637
KAPPA = 0.127
CUTOFF = 12.0
R_FLOOR = 1.5
DIELECTRIC_SLOPE = 4.0


@dataclass
class ChargedAtoms:
    xyz: np.ndarray  # (n, 3)
    charge: np.ndarray  # (n,)

    @property
    def net(self) -> float:
        return float(self.charge.sum())


def _atom(line: str):
    return ([float(line[30:38]), float(line[38:46]), float(line[46:54])], float(line[70:76]))


def read_pdbqt(text: str) -> ChargedAtoms:
    """Coordinates and partial charges of every ATOM/HETATM record (the first MODEL only)."""
    xyz, q = [], []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            a, c = _atom(line)
            xyz.append(a)
            q.append(c)
    return ChargedAtoms(np.asarray(xyz, dtype=float).reshape(-1, 3), np.asarray(q, dtype=float))


def read_pdbqt_models(text: str) -> List[ChargedAtoms]:
    """Every MODEL of a Vina output PDBQT, in output (score) order."""
    models, xyz, q = [], [], []
    for line in text.splitlines():
        if line.startswith("MODEL"):
            xyz, q = [], []
        elif line.startswith("ENDMDL"):
            models.append(ChargedAtoms(np.asarray(xyz, dtype=float).reshape(-1, 3), np.asarray(q, dtype=float)))
            xyz, q = [], []
        elif line.startswith(("ATOM", "HETATM")):
            a, c = _atom(line)
            xyz.append(a)
            q.append(c)
    if not models and xyz:  # a single pose without MODEL records
        models.append(ChargedAtoms(np.asarray(xyz, dtype=float).reshape(-1, 3), np.asarray(q, dtype=float)))
    return models


class ReceptorField:
    """The receptor's charges, indexed once so that many poses are scored cheaply."""

    def __init__(self, receptor: ChargedAtoms) -> None:
        from scipy.spatial import cKDTree

        keep = receptor.charge != 0.0
        self.xyz = receptor.xyz[keep]
        self.charge = receptor.charge[keep]
        self.tree = cKDTree(self.xyz) if len(self.xyz) else None

    def energy(self, ligand: ChargedAtoms) -> float:
        """E_el of one pose in kcal/mol."""
        if self.tree is None or len(ligand.xyz) == 0:
            return 0.0
        total = 0.0
        for xyz, q in zip(ligand.xyz, ligand.charge):
            if q == 0.0:
                continue
            idx = self.tree.query_ball_point(xyz, CUTOFF)
            if not idx:
                continue
            r = np.linalg.norm(self.xyz[idx] - xyz, axis=1)
            r = np.maximum(r, R_FLOOR)
            total += q * float(np.sum(self.charge[idx] * np.exp(-KAPPA * r) / (DIELECTRIC_SLOPE * r * r)))
        return COULOMB * total
