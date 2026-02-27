"""Electrostatics and pH-dependent APBS/PROPKA workflows."""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_PH_GRID: Tuple[float, ...] = (5.0, 6.0, 7.0, 7.4, 8.0)
BASIC_RESIDUES = {"ARG", "LYS", "HIS"}
PHYSIOLOGICAL_PH = 7.4


@dataclass
class PHAnalysisResult:
    """Container for pH-dependent electrostatics outputs."""

    propka_table: pd.DataFrame
    site_potentials: pd.DataFrame
    ph_sensitive_residues: pd.DataFrame
    optimal_binding_ph: Dict[str, float]
    profile_comparison: Optional[pd.DataFrame]
    plot_path: Optional[Path]


class ElectrostaticsCalculator:
    """Run APBS/PROPKA pipelines and derive pH-dependent IP-binding predictions."""

    def __init__(
        self,
        apbs_path: str = "apbs",
        pdb2pqr_path: str = "pdb2pqr",
        propka_path: str = "propka3",
    ):
        self.apbs_path = apbs_path
        self.pdb2pqr_path = pdb2pqr_path
        self.propka_path = propka_path

    def run_propka(self, pdb_path: Union[str, Path], output_dir: Union[str, Path]) -> pd.DataFrame:
        """Run PROPKA and parse pKa values for basic residues."""
        pdb_path = Path(pdb_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [self.propka_path, str(pdb_path)]
        try:
            subprocess.run(cmd, cwd=output_dir, capture_output=True, text=True, check=True)
        except FileNotFoundError as exc:
            raise RuntimeError(f"PROPKA executable not found: {self.propka_path}") from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"PROPKA failed: {exc.stderr.strip()}") from exc

        pka_file = output_dir / f"{pdb_path.stem}.pka"
        if not pka_file.exists():
            raise RuntimeError(f"PROPKA output file not found: {pka_file}")

        return self.parse_propka_output(pka_file)

    @staticmethod
    def parse_propka_output(pka_file: Union[str, Path]) -> pd.DataFrame:
        """Parse basic-residue pKa values from PROPKA output."""
        pattern = re.compile(
            r"^\s*(ARG|LYS|HIS)\s+(\d+)\s+([A-Za-z0-9])\s+([-+]?\d+(?:\.\d+)?)"
        )

        records: List[Dict[str, Union[str, int, float]]] = []
        with Path(pka_file).open("r", encoding="utf-8") as handle:
            for line in handle:
                match = pattern.match(line)
                if not match:
                    continue
                residue_name, residue_number, chain_id, pka = match.groups()
                records.append(
                    {
                        "residue_name": residue_name,
                        "residue_number": int(residue_number),
                        "chain_id": chain_id,
                        "pka": float(pka),
                    }
                )

        if not records:
            return pd.DataFrame(columns=["residue_name", "residue_number", "chain_id", "pka"])
        return pd.DataFrame(records)

    def generate_pqr(
        self,
        pdb_path: Union[str, Path],
        ph: float,
        output_dir: Union[str, Path],
    ) -> Path:
        """Generate a PQR file at a target pH using pdb2pqr."""
        pdb_path = Path(pdb_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        pqr_path = output_dir / f"{pdb_path.stem}_ph{ph:.1f}.pqr"

        cmd = [
            self.pdb2pqr_path,
            "--ff=PARSE",
            f"--with-ph={ph}",
            str(pdb_path),
            str(pqr_path),
        ]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except FileNotFoundError as exc:
            raise RuntimeError(f"pdb2pqr executable not found: {self.pdb2pqr_path}") from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"pdb2pqr failed at pH {ph}: {exc.stderr.strip()}") from exc

        if not pqr_path.exists():
            raise RuntimeError(f"pdb2pqr did not create expected output: {pqr_path}")
        return pqr_path

    def run_apbs(self, pqr_path: Union[str, Path], output_dir: Union[str, Path]) -> float:
        """Run APBS and return a scalar electrostatic proxy value from stdout."""
        pqr_path = Path(pqr_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        apbs_input = output_dir / f"{pqr_path.stem}.in"
        map_path = output_dir / f"{pqr_path.stem}.dx"
        apbs_input.write_text(
            "\n".join(
                [
                    "read",
                    f"    mol pqr {pqr_path}",
                    "end",
                    "elec",
                    "    mg-auto",
                    "    dime 65 65 65",
                    "    cglen 80 80 80",
                    "    fglen 60 60 60",
                    "    cgcent mol 1",
                    "    fgcent mol 1",
                    "    mol 1",
                    "    lpbe",
                    "    bcfl sdh",
                    "    pdie 2.0",
                    "    sdie 78.5",
                    "    chgm spl2",
                    "    sdens 10.0",
                    "    srfm smol",
                    "    srad 1.4",
                    "    swin 0.3",
                    "    temp 298.15",
                    "    calcenergy total",
                    "    calcforce no",
                    f"    write pot dx {map_path.with_suffix('')}",
                    "end",
                    "quit",
                ]
            ),
            encoding="utf-8",
        )

        try:
            result = subprocess.run(
                [self.apbs_path, str(apbs_input)],
                capture_output=True,
                text=True,
                check=True,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(f"APBS executable not found: {self.apbs_path}") from exc
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(f"APBS failed for {pqr_path.name}: {exc.stderr.strip()}") from exc

        energy_match = re.search(
            r"Global net ELEC energy\s*=\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)",
            result.stdout,
        )
        if not energy_match:
            raise RuntimeError("APBS output did not include 'Global net ELEC energy'.")

        return float(energy_match.group(1))

    def identify_ph_sensitive_residues(
        self,
        propka_table: pd.DataFrame,
        residue_numbers: Sequence[int],
        physiological_ph: float = PHYSIOLOGICAL_PH,
        window: float = 1.0,
    ) -> pd.DataFrame:
        """Return residues with pKa near physiological pH for the selected site."""
        site_table = propka_table[propka_table["residue_number"].isin(residue_numbers)].copy()
        delta = (site_table["pka"] - physiological_ph).abs()
        return site_table[delta <= window].reset_index(drop=True)

    def analyze_ph_dependent_binding(
        self,
        pdb_path: Union[str, Path],
        candidate_sites: Dict[str, Sequence[int]],
        output_dir: Union[str, Path],
        ph_values: Iterable[float] = DEFAULT_PH_GRID,
        site_types: Optional[Dict[str, str]] = None,
    ) -> PHAnalysisResult:
        """Full pH-dependent APBS/PROPKA analysis for candidate IP-binding sites."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        propka_table = self.run_propka(pdb_path, output_dir=output_dir)

        apbs_energy_by_ph: Dict[float, float] = {}
        for ph in ph_values:
            pqr_path = self.generate_pqr(pdb_path=pdb_path, ph=ph, output_dir=output_dir)
            apbs_energy_by_ph[ph] = self.run_apbs(pqr_path=pqr_path, output_dir=output_dir)

        potential_rows = []
        sensitive_frames = []
        optimal_binding_ph: Dict[str, float] = {}

        for site_name, residue_ids in candidate_sites.items():
            site_residues = propka_table[propka_table["residue_number"].isin(residue_ids)].copy()
            if site_residues.empty:
                site_residues = pd.DataFrame(
                    {
                        "residue_number": list(residue_ids),
                        "pka": [10.5] * len(residue_ids),
                    }
                )

            for ph, apbs_energy in apbs_energy_by_ph.items():
                protonation = 1.0 / (1.0 + 10.0 ** (ph - site_residues["pka"]))
                site_charge_factor = float(protonation.sum())
                potential_rows.append(
                    {
                        "site": site_name,
                        "ph": ph,
                        "apbs_energy": apbs_energy,
                        "site_charge_factor": site_charge_factor,
                        "site_potential": apbs_energy * site_charge_factor,
                    }
                )

            sensitive = self.identify_ph_sensitive_residues(propka_table, residue_ids)
            if not sensitive.empty:
                sensitive = sensitive.copy()
                sensitive["site"] = site_name
                sensitive_frames.append(sensitive)

        site_potentials = pd.DataFrame(potential_rows).sort_values(["site", "ph"])
        for site_name, site_df in site_potentials.groupby("site"):
            best_row = site_df.loc[site_df["site_potential"].idxmax()]
            optimal_binding_ph[site_name] = float(best_row["ph"])

        ph_sensitive_residues = (
            pd.concat(sensitive_frames, ignore_index=True)
            if sensitive_frames
            else pd.DataFrame(columns=["residue_name", "residue_number", "chain_id", "pka", "site"])
        )

        profile_comparison = self.compare_site_profiles(site_potentials, site_types)
        plot_path = output_dir / "ph_dependent_electrostatics.png"
        self.plot_potential_vs_ph(site_potentials, plot_path, site_types)

        return PHAnalysisResult(
            propka_table=propka_table,
            site_potentials=site_potentials,
            ph_sensitive_residues=ph_sensitive_residues,
            optimal_binding_ph=optimal_binding_ph,
            profile_comparison=profile_comparison,
            plot_path=plot_path,
        )

    @staticmethod
    def compare_site_profiles(
        site_potentials: pd.DataFrame,
        site_types: Optional[Dict[str, str]],
    ) -> Optional[pd.DataFrame]:
        """Compare electrostatic profiles between cryptic and surface sites."""
        if not site_types:
            return None

        typed = site_potentials.copy()
        typed["site_type"] = typed["site"].map(site_types).fillna("unlabeled")
        summary = (
            typed.groupby(["site_type", "ph"]) ["site_potential"]
            .mean()
            .reset_index(name="mean_site_potential")
        )
        return summary

    @staticmethod
    def plot_potential_vs_ph(
        site_potentials: pd.DataFrame,
        output_path: Union[str, Path],
        site_types: Optional[Dict[str, str]] = None,
    ) -> Path:
        """Generate publication-quality electrostatic potential vs pH curves."""
        output_path = Path(output_path)

        plt.style.use("seaborn-v0_8-whitegrid")
        fig, ax = plt.subplots(figsize=(8, 5), dpi=300)

        for site, data in site_potentials.groupby("site"):
            label = site
            if site_types and site in site_types:
                label = f"{site} ({site_types[site]})"
            ax.plot(data["ph"], data["site_potential"], marker="o", linewidth=2, label=label)

        ax.set_xlabel("pH")
        ax.set_ylabel("Site electrostatic potential (a.u.)")
        ax.set_title("pH-dependent electrostatic profiles of candidate IP-binding sites")
        ax.legend(frameon=True, fontsize=8)
        fig.tight_layout()
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)

        return output_path


def run_propka_wrapper(
    pdb_path: Union[str, Path], output_dir: Union[str, Path], propka_path: str = "propka3"
) -> pd.DataFrame:
    """Convenience wrapper for PROPKA with robust error handling."""
    return ElectrostaticsCalculator(propka_path=propka_path).run_propka(pdb_path, output_dir)


def run_apbs_wrapper(
    pqr_path: Union[str, Path], output_dir: Union[str, Path], apbs_path: str = "apbs"
) -> float:
    """Convenience wrapper for APBS with robust error handling."""
    return ElectrostaticsCalculator(apbs_path=apbs_path).run_apbs(pqr_path, output_dir)
