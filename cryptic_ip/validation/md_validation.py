"""Molecular dynamics validation pipeline for cryptic pocket stability."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import json
import numpy as np
import pandas as pd


@dataclass
class MDSimulationConfig:
    """Configuration for OpenMM molecular dynamics runs."""

    temperature_k: float = 300.0
    pressure_bar: float = 1.0
    equilibration_ns: float = 1.0
    production_ns: float = 100.0
    timestep_fs: float = 2.0
    trajectory_interval_ps: float = 10.0
    friction_per_ps: float = 1.0
    water_padding_nm: float = 1.0
    ionic_strength_molar: float = 0.15


@dataclass
class PocketStabilityThresholds:
    """Thresholds for classifying pocket stability from MD metrics."""

    sasa_stably_buried: float = 1.0
    sasa_exposed: float = 5.0
    rmsf_stably_buried: float = 0.20
    rmsf_exposed: float = 0.45
    waters_stably_buried: float = 1.0
    waters_exposed: float = 5.0


class OpenMMMDValidationPipeline:
    """End-to-end MD validation pipeline for top cryptic-pocket candidates."""

    def __init__(
        self,
        output_dir: str | Path = "results/md_validation",
        config: Optional[MDSimulationConfig] = None,
        thresholds: Optional[PocketStabilityThresholds] = None,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config or MDSimulationConfig()
        self.thresholds = thresholds or PocketStabilityThresholds()

    def validate_top_candidates(
        self,
        candidates_csv: str | Path,
        top_n: int = 20,
    ) -> pd.DataFrame:
        """Run full MD validation for top candidates and write a report."""
        candidates = pd.read_csv(candidates_csv)
        if "composite_score" in candidates.columns:
            candidates = candidates.sort_values("composite_score", ascending=False)

        top_candidates = candidates.head(top_n).copy()
        results: List[Dict] = []

        for idx, row in top_candidates.iterrows():
            candidate_id = str(row.get("candidate_id", row.get("pocket_id", f"candidate_{idx}")))
            structure_path = Path(str(row["structure_path"]))
            pocket_residues = self._parse_residues(row.get("pocket_residues", ""))
            pocket_center = self._parse_center(row.get("pocket_center", ""))

            candidate_dir = self.output_dir / candidate_id
            candidate_dir.mkdir(parents=True, exist_ok=True)

            try:
                simulation_paths = self.run_simulation(
                    structure_path=structure_path,
                    output_dir=candidate_dir,
                )
                analysis = self.analyze_trajectory(
                    topology_path=simulation_paths["topology"],
                    trajectory_path=simulation_paths["trajectory"],
                    pocket_residues=pocket_residues,
                    pocket_center=pocket_center,
                )
                classification = self.classify_pocket_stability(analysis)
                self.generate_visualization_scripts(
                    topology_path=simulation_paths["topology"],
                    trajectory_path=simulation_paths["trajectory"],
                    pocket_residues=pocket_residues,
                    output_dir=candidate_dir,
                )

                result = {
                    "candidate_id": candidate_id,
                    "structure_path": str(structure_path),
                    "classification": classification,
                    **analysis,
                }
            except Exception as exc:  # pragma: no cover - runtime dependent on OpenMM/MDTraj.
                result = {
                    "candidate_id": candidate_id,
                    "structure_path": str(structure_path),
                    "classification": "failed",
                    "error": str(exc),
                }

            results.append(result)

        report_df = pd.DataFrame(results)
        report_path = self.output_dir / "md_validation_report.csv"
        report_df.to_csv(report_path, index=False)

        stable_df = report_df[report_df["classification"] == "stably buried"]
        stable_df.to_csv(self.output_dir / "stable_candidates.csv", index=False)

        with open(self.output_dir / "summary.json", "w", encoding="utf-8") as handle:
            json.dump(self._build_summary(report_df), handle, indent=2)

        return report_df

    def run_simulation(self, structure_path: Path, output_dir: Path) -> Dict[str, str]:
        """Prepare structure and run OpenMM equilibration + production MD."""
        try:
            from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, unit
            from openmm.app import (
                DCDReporter,
                ForceField,
                Modeller,
                PDBFile,
                PME,
                Simulation,
                StateDataReporter,
            )
        except ImportError as exc:  # pragma: no cover - dependency-gated.
            raise ImportError("OpenMM is required for MD simulation. Install `openmm`.") from exc

        pdb = PDBFile(str(structure_path))
        forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")

        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield)
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=self.config.water_padding_nm * unit.nanometer,
            ionicStrength=self.config.ionic_strength_molar * unit.molar,
            neutralize=True,
        )

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=None,
        )
        system.addForce(
            MonteCarloBarostat(
                self.config.pressure_bar * unit.bar,
                self.config.temperature_k * unit.kelvin,
            )
        )

        integrator = LangevinMiddleIntegrator(
            self.config.temperature_k * unit.kelvin,
            self.config.friction_per_ps / unit.picosecond,
            self.config.timestep_fs * unit.femtoseconds,
        )

        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        simulation.minimizeEnergy()

        topology_out = output_dir / "prepared_system.pdb"
        with open(topology_out, "w", encoding="utf-8") as handle:
            PDBFile.writeFile(modeller.topology, modeller.positions, handle)

        interval_steps = int(self.config.trajectory_interval_ps * 1000 / self.config.timestep_fs)
        equil_steps = int(self.config.equilibration_ns * 1_000_000 / self.config.timestep_fs)
        prod_steps = int(self.config.production_ns * 1_000_000 / self.config.timestep_fs)

        trajectory_out = output_dir / "production.dcd"
        simulation.reporters.append(DCDReporter(str(trajectory_out), interval_steps))
        simulation.reporters.append(
            StateDataReporter(
                str(output_dir / "md.log"),
                interval_steps,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                density=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=equil_steps + prod_steps,
                separator=",",
            )
        )

        simulation.step(equil_steps)
        simulation.step(prod_steps)

        return {
            "topology": str(topology_out),
            "trajectory": str(trajectory_out),
            "log": str(output_dir / "md.log"),
        }

    def analyze_trajectory(
        self,
        topology_path: str | Path,
        trajectory_path: str | Path,
        pocket_residues: Sequence[int],
        pocket_center: Optional[Tuple[float, float, float]] = None,
    ) -> Dict[str, float]:
        """Compute SASA, RMSF, pocket volume, and water penetration metrics."""
        try:
            import mdtraj as md
        except ImportError as exc:  # pragma: no cover - dependency-gated.
            raise ImportError("MDTraj is required for trajectory analysis. Install `mdtraj`.") from exc

        traj = md.load(str(trajectory_path), top=str(topology_path))

        residue_indices = self._resolve_residue_atom_indices(traj, pocket_residues)
        if not residue_indices:
            raise ValueError("No atoms selected for pocket residues; check residue numbering.")

        atom_sasa = md.shrake_rupley(traj, mode="atom")
        pocket_sasa = atom_sasa[:, residue_indices].sum(axis=1)

        pocket_traj = traj.atom_slice(residue_indices)
        ref_xyz = pocket_traj.xyz.mean(axis=0, keepdims=True)
        rmsf = md.rmsf(pocket_traj, reference=md.Trajectory(ref_xyz, pocket_traj.topology), frame=0)

        pocket_volume = self._estimate_pocket_volume(traj, residue_indices)
        water_count = self._count_waters_in_pocket(traj, residue_indices, pocket_center=pocket_center)

        return {
            "avg_pocket_sasa_nm2": float(np.mean(pocket_sasa)),
            "avg_pocket_rmsf_nm": float(np.mean(rmsf)),
            "avg_pocket_volume_nm3": float(np.mean(pocket_volume)),
            "std_pocket_volume_nm3": float(np.std(pocket_volume)),
            "avg_waters_in_pocket": float(np.mean(water_count)),
        }

    def classify_pocket_stability(self, analysis: Dict[str, float]) -> str:
        """Classify pocket dynamics into buried/transient/exposed classes."""
        buried_score = 0
        exposed_score = 0

        if analysis["avg_pocket_sasa_nm2"] <= self.thresholds.sasa_stably_buried:
            buried_score += 1
        elif analysis["avg_pocket_sasa_nm2"] >= self.thresholds.sasa_exposed:
            exposed_score += 1

        if analysis["avg_pocket_rmsf_nm"] <= self.thresholds.rmsf_stably_buried:
            buried_score += 1
        elif analysis["avg_pocket_rmsf_nm"] >= self.thresholds.rmsf_exposed:
            exposed_score += 1

        if analysis["avg_waters_in_pocket"] <= self.thresholds.waters_stably_buried:
            buried_score += 1
        elif analysis["avg_waters_in_pocket"] >= self.thresholds.waters_exposed:
            exposed_score += 1

        if buried_score >= 2:
            return "stably buried"
        if exposed_score >= 2:
            return "exposed"
        return "transiently accessible"

    def generate_visualization_scripts(
        self,
        topology_path: str | Path,
        trajectory_path: str | Path,
        pocket_residues: Sequence[int],
        output_dir: Path,
    ) -> None:
        """Generate PyMOL and VMD scripts for pocket-dynamics inspection."""
        selection = "+".join(str(r) for r in pocket_residues)

        pymol_script = output_dir / "view_pocket_dynamics.pml"
        pymol_script.write_text(
            "\n".join(
                [
                    f"load {Path(topology_path).name}, protein",
                    f"load_traj {Path(trajectory_path).name}, protein",
                    "hide everything, all",
                    "show cartoon, protein",
                    "color gray80, protein",
                    f"select pocket, resi {selection}",
                    "show sticks, pocket",
                    "color yellow, pocket",
                    "set cartoon_transparency, 0.2, protein",
                    "set sphere_scale, 0.2",
                    "show spheres, solvent within 4 of pocket",
                    "color cyan, solvent within 4 of pocket",
                ]
            ),
            encoding="utf-8",
        )

        vmd_script = output_dir / "view_pocket_dynamics.vmd.tcl"
        vmd_script.write_text(
            "\n".join(
                [
                    f"mol new {Path(topology_path).name}",
                    f"mol addfile {Path(trajectory_path).name} waitfor all",
                    "mol delrep 0 top",
                    "mol representation NewCartoon",
                    "mol color ColorID 8",
                    "mol selection protein",
                    "mol addrep top",
                    "mol representation Licorice",
                    f"mol selection resid {' '.join(str(r) for r in pocket_residues)}",
                    "mol color ColorID 4",
                    "mol addrep top",
                    "mol representation VDW 0.4",
                    f"mol selection water and within 4 of (resid {' '.join(str(r) for r in pocket_residues)})",
                    "mol color ColorID 0",
                    "mol addrep top",
                ]
            ),
            encoding="utf-8",
        )

    def _build_summary(self, report_df: pd.DataFrame) -> Dict[str, object]:
        counts = report_df["classification"].value_counts(dropna=False).to_dict()
        return {
            "total_candidates": int(len(report_df)),
            "classification_counts": counts,
            "stable_candidates": report_df[
                report_df["classification"] == "stably buried"
            ]["candidate_id"].tolist(),
            "filtered_out_candidates": report_df[
                report_df["classification"].isin(["exposed", "failed"])
            ]["candidate_id"].tolist(),
        }

    @staticmethod
    def _parse_residues(residue_text: str) -> List[int]:
        if pd.isna(residue_text) or str(residue_text).strip() == "":
            return []
        return [int(token.strip()) for token in str(residue_text).split(",") if token.strip()]

    @staticmethod
    def _parse_center(center_text: str) -> Optional[Tuple[float, float, float]]:
        if pd.isna(center_text) or str(center_text).strip() == "":
            return None
        parts = [float(token.strip()) for token in str(center_text).split(",") if token.strip()]
        if len(parts) != 3:
            raise ValueError("Pocket center must contain three comma-separated coordinates.")
        return parts[0], parts[1], parts[2]

    @staticmethod
    def _resolve_residue_atom_indices(traj, pocket_residues: Sequence[int]) -> List[int]:
        if not pocket_residues:
            return []
        residue_set = set(pocket_residues)
        atom_indices = [
            atom.index
            for atom in traj.topology.atoms
            if atom.residue.resSeq in residue_set and atom.element.symbol != "H"
        ]
        return atom_indices

    @staticmethod
    def _estimate_pocket_volume(traj, atom_indices: Sequence[int]) -> np.ndarray:
        pocket_xyz = traj.xyz[:, atom_indices, :]
        mins = pocket_xyz.min(axis=1)
        maxs = pocket_xyz.max(axis=1)
        lengths = np.maximum(maxs - mins, 1e-6)
        return np.prod(lengths, axis=1)

    @staticmethod
    def _count_waters_in_pocket(
        traj,
        atom_indices: Sequence[int],
        pocket_center: Optional[Tuple[float, float, float]] = None,
        radius_nm: float = 0.4,
    ) -> np.ndarray:
        pocket_coords = traj.xyz[:, atom_indices, :]
        if pocket_center is None:
            centers = pocket_coords.mean(axis=1)
        else:
            centers = np.tile(np.array(pocket_center), (traj.n_frames, 1))

        water_oxygen_indices = [
            atom.index
            for atom in traj.topology.atoms
            if atom.residue.is_water and atom.name.upper().startswith("O")
        ]

        if not water_oxygen_indices:
            return np.zeros(traj.n_frames)

        water_xyz = traj.xyz[:, water_oxygen_indices, :]
        diffs = water_xyz - centers[:, np.newaxis, :]
        dists = np.linalg.norm(diffs, axis=2)
        return (dists <= radius_nm).sum(axis=1)
