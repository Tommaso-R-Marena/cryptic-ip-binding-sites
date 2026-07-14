"""
Comprehensive validation suite for pipeline testing.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .adar2 import validate_adar2
from .burial_metrics import compute_burial_metrics
from .control_scoring import (
    negative_passed,
    positive_passed,
    separation_quality,
    validation_score,
)
from .structure_context import ligand_context
from ..analysis import ProteinAnalyzer


class ValidationSuite:
    """
    Complete validation suite for testing pipeline on known examples.

    Tier 1 controls gate Phase 1 readiness (ADAR2 + PLCδ1 PH).
    Tier 2 controls extend sensitivity/specificity testing.
    """

    POSITIVE_CONTROLS = {
        "ADAR2": {
            "pdb": "1ZY7",
            "ip_type": "IP6",
            "description": "Gold standard - completely buried IP6",
            "expected_score": 0.75,
            "use_adar2_validator": True,
            "tier": 1,
            "burial_class": "cryptic",
        },
        "Pds5B": {
            "pdb": "5HDT",
            "ip_type": "IP6",
            "description": "Cohesin regulator - crystal shows surface-exposed IP6 (artifact)",
            "expected_score": 0.65,
            "tier": 2,
            "burial_class": "crystal_artifact",
        },
        "HDAC1": {
            "pdb": "5ICN",
            "ip_type": "IP4",
            "description": "Histone deacetylase with semi-cryptic IP at interface",
            "expected_score": 0.60,
            "tier": 2,
            "burial_class": "semi_cryptic",
        },
    }

    NEGATIVE_CONTROLS = {
        "PLCd1_PH": {
            "pdb": "1MAI",
            "ip_type": "IP3",
            "description": "Classic surface-exposed PH domain",
            "expected_score": 0.30,
            "tier": 1,
        },
        "Btk_PH": {
            "pdb": "1BWN",
            "ip_type": "IP4",
            "description": "Kinase PH domain with surface-accessible IP4",
            "expected_score": 0.35,
            "tier": 2,
            "decoy_mode": True,
        },
    }

    def __init__(self, data_dir: str = "data/validation", use_electrostatics: bool = True):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.use_electrostatics = use_electrostatics
        self.results: List[Dict] = []

    def download_pdb(self, pdb_id: str) -> Path:
        """Download a PDB structure from RCSB if not already cached."""
        pdb_path = self.data_dir / f"{pdb_id.upper()}.pdb"
        if pdb_path.exists():
            return pdb_path

        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        print(f"  Downloading {pdb_id.upper()} from RCSB...")
        urllib.request.urlretrieve(url, pdb_path)
        return pdb_path

    def _select_best_pocket(
        self,
        scored: pd.DataFrame,
        analyzer: ProteinAnalyzer,
        site_residues: Set[int],
        ligand_centroid: Optional[Tuple[float, float, float]],
        *,
        decoy_mode: bool = False,
    ) -> pd.Series:
        if decoy_mode or ligand_centroid is None:
            return scored.iloc[0]

        best_row = scored.iloc[0]
        best_key = (-float("inf"), -1, -1.0)

        for _, row in scored.iterrows():
            pocket_id = int(row["pocket_id"])
            overlap = len(
                set(analyzer.get_pocket_residues(pocket_id, distance_cutoff=8.0)).intersection(site_residues)
            )
            distance = float("inf")
            if analyzer.pockets is not None:
                pocket = analyzer.pockets[analyzer.pockets["pocket_id"] == pocket_id].iloc[0]
                center = np.array([pocket["center_x"], pocket["center_y"], pocket["center_z"]], dtype=float)
                distance = float(np.linalg.norm(center - np.array(ligand_centroid, dtype=float)))

            key = (-distance, overlap, float(row["composite_score"]))
            if key > best_key:
                best_key = key
                best_row = row

        return best_row

    def score_structure_with_ligand_context(
        self,
        pdb_path: Path,
        *,
        control_type: str,
        decoy_mode: bool = False,
    ) -> Dict:
        site_residues, ligand_sasa, ligand_centroid = ligand_context(pdb_path)
        burial = compute_burial_metrics(pdb_path)

        analyzer = ProteinAnalyzer(
            str(pdb_path),
            work_dir=str(self.data_dir / "work" / pdb_path.stem),
            skip_electrostatics=not self.use_electrostatics,
        )
        scored = analyzer.run_pipeline(include_electrostatics=self.use_electrostatics)
        if scored.empty:
            raise RuntimeError(f"No pockets scored for {pdb_path.name}")

        top_pocket = self._select_best_pocket(
            scored, analyzer, site_residues, ligand_centroid, decoy_mode=decoy_mode
        )

        pocket_composite = float(top_pocket["composite_score"])
        pocket_potential = top_pocket.get("electrostatic_potential")
        basic_residues = len(site_residues) if site_residues else int(top_pocket["basic_residues"])
        sasa_metric = float(burial.ligand_sasa if burial.ligand_sasa is not None else top_pocket["sasa"])
        val_score = validation_score(
            control_type,
            pocket_composite,
            burial.ligand_sasa if burial.ligand_sasa is not None else ligand_sasa,
            basic_residues,
            pocket_potential=pocket_potential,
            use_electrostatics=self.use_electrostatics,
        )

        pocket_residues = analyzer.get_pocket_residues(int(top_pocket["pocket_id"]))
        plddt = top_pocket.get("plddt_confidence", float("nan"))

        return {
            "total_pockets": int(len(scored)),
            "top_pocket_id": int(top_pocket["pocket_id"]),
            "top_pocket_score": pocket_composite,
            "validation_score": val_score,
            "volume": float(top_pocket["volume"]),
            "sasa": sasa_metric,
            "phosphate_sasa": burial.phosphate_sasa,
            "delta_sasa": burial.delta_sasa,
            "burial_depth": burial.burial_depth,
            "burial_class": burial.burial_class,
            "pocket_residue_sasa": float(top_pocket["sasa"]),
            "ligand_sasa": burial.ligand_sasa if burial.ligand_sasa is not None else ligand_sasa,
            "basic_residues": basic_residues,
            "depth": float(top_pocket["depth"]),
            "electrostatic_potential": pocket_potential,
            "plddt_confidence": plddt,
            "pocket_residues": ",".join(str(r) for r in pocket_residues),
            "structure_path": str(pdb_path),
            "decoy_mode": decoy_mode,
        }

    def validate_positive_control(self, name: str, info: Dict) -> Dict:
        """Validate a single positive control structure."""
        burial_class = info.get("burial_class", "cryptic")

        if info.get("use_adar2_validator"):
            result = validate_adar2(use_alphafold=False, use_electrostatics=self.use_electrostatics)
            pocket_composite = float(result["top_pocket_score"])
            ligand_sasa = result.get("ligand_sasa")
            pocket_potential = result.get("electrostatic_potential")
            val_score = validation_score(
                "positive",
                pocket_composite,
                ligand_sasa,
                int(result["basic_residues"]),
                pocket_potential=pocket_potential,
                use_electrostatics=self.use_electrostatics,
            )
            passed = positive_passed(
                val_score, ligand_sasa, int(result["basic_residues"]), pocket_composite, burial_class=burial_class
            )
            return {
                "protein": name,
                "control_type": "positive",
                "tier": info.get("tier", 2),
                "pdb_id": info["pdb"],
                "ip_type": info["ip_type"],
                "score": val_score,
                "pocket_score": pocket_composite,
                "expected": info["expected_score"],
                "passed": passed,
                "burial_class": burial_class,
                "total_pockets": result["total_pockets"],
                "sasa": result["sasa"],
                "basic_residues": result["basic_residues"],
                "electrostatic_potential": pocket_potential,
                "description": info["description"],
            }

        pdb_path = self.download_pdb(info["pdb"])
        metrics = self.score_structure_with_ligand_context(pdb_path, control_type="positive")
        if metrics["ligand_sasa"] is not None and metrics["ligand_sasa"] > 100:
            burial_class = "crystal_artifact"
            metrics["burial_class"] = burial_class

        passed = positive_passed(
            metrics["validation_score"],
            metrics.get("ligand_sasa"),
            metrics["basic_residues"],
            metrics["top_pocket_score"],
            burial_class=burial_class,
        )
        return self._result_row(name, "positive", info, metrics, passed, burial_class=burial_class)

    def validate_negative_control(self, name: str, info: Dict) -> Dict:
        """Validate a single negative control structure."""
        pdb_path = self.download_pdb(info["pdb"])
        decoy_mode = bool(info.get("decoy_mode", False))
        metrics = self.score_structure_with_ligand_context(
            pdb_path, control_type="negative", decoy_mode=decoy_mode
        )
        passed = negative_passed(
            metrics["validation_score"],
            metrics.get("ligand_sasa"),
            metrics["top_pocket_score"],
            decoy_mode=decoy_mode,
        )
        return self._result_row(name, "negative", info, metrics, passed)

    def _result_row(
        self,
        name: str,
        control_type: str,
        info: Dict,
        metrics: Dict,
        passed: bool,
        burial_class: Optional[str] = None,
    ) -> Dict:
        return {
            "protein": name,
            "control_type": control_type,
            "tier": info.get("tier", 2),
            "pdb_id": info["pdb"],
            "ip_type": info["ip_type"],
            "score": metrics["validation_score"],
            "pocket_score": metrics["top_pocket_score"],
            "expected": info["expected_score"],
            "passed": passed,
            "burial_class": burial_class or metrics.get("burial_class"),
            "total_pockets": metrics["total_pockets"],
            "sasa": metrics["sasa"],
            "phosphate_sasa": metrics.get("phosphate_sasa"),
            "burial_depth": metrics.get("burial_depth"),
            "basic_residues": metrics["basic_residues"],
            "electrostatic_potential": metrics.get("electrostatic_potential"),
            "plddt_confidence": metrics.get("plddt_confidence"),
            "decoy_mode": metrics.get("decoy_mode", False),
            "description": info["description"],
        }

    def run_positive_controls(self) -> pd.DataFrame:
        """Test pipeline on positive controls (known buried sites)."""
        print("\n" + "=" * 60)
        print("TESTING POSITIVE CONTROLS")
        print("=" * 60)

        results = []
        for name, info in self.POSITIVE_CONTROLS.items():
            print(f"\nTesting {name}...")
            try:
                results.append(self.validate_positive_control(name, info))
                status = "PASS" if results[-1]["passed"] else "FAIL"
                print(
                    f"  Validation score: {results[-1]['score']:.3f} "
                    f"(pocket={results[-1]['pocket_score']:.3f}, SASA={results[-1]['sasa']:.1f}, "
                    f"class={results[-1].get('burial_class')}) ({status})"
                )
            except Exception as exc:
                print(f"  Error: {exc}")
                results.append(
                    {
                        "protein": name,
                        "control_type": "positive",
                        "tier": info.get("tier", 2),
                        "pdb_id": info.get("pdb"),
                        "ip_type": info["ip_type"],
                        "score": float("nan"),
                        "pocket_score": float("nan"),
                        "expected": info["expected_score"],
                        "passed": False,
                        "burial_class": info.get("burial_class"),
                        "description": info["description"],
                        "error": str(exc),
                    }
                )

        return pd.DataFrame(results)

    def run_negative_controls(self) -> pd.DataFrame:
        """Test pipeline on negative controls (surface binding)."""
        print("\n" + "=" * 60)
        print("TESTING NEGATIVE CONTROLS")
        print("=" * 60)

        results = []
        for name, info in self.NEGATIVE_CONTROLS.items():
            print(f"\nTesting {name}...")
            try:
                results.append(self.validate_negative_control(name, info))
                status = "PASS" if results[-1]["passed"] else "FAIL"
                print(
                    f"  Validation score: {results[-1]['score']:.3f} "
                    f"(pocket={results[-1]['pocket_score']:.3f}, SASA={results[-1]['sasa']:.1f}) ({status})"
                )
            except Exception as exc:
                print(f"  Error: {exc}")
                results.append(
                    {
                        "protein": name,
                        "control_type": "negative",
                        "tier": info.get("tier", 2),
                        "pdb_id": info.get("pdb"),
                        "ip_type": info["ip_type"],
                        "score": float("nan"),
                        "pocket_score": float("nan"),
                        "expected": info["expected_score"],
                        "passed": False,
                        "description": info["description"],
                        "error": str(exc),
                    }
                )

        return pd.DataFrame(results)

    def run_full_validation(self, output_dir: Optional[Path] = None) -> Dict:
        """Run complete validation suite and optionally persist CSV summaries."""
        positive_results = self.run_positive_controls()
        negative_results = self.run_negative_controls()

        pos_classes = positive_results.get("burial_class", pd.Series(dtype=str)).tolist()
        sep = separation_quality(
            positive_results["score"],
            negative_results["score"],
            positive_classes=pos_classes,
        )
        tier1_pos = positive_results[positive_results["tier"] == 1]
        tier1_neg = negative_results[negative_results["tier"] == 1]
        if not tier1_pos.empty and not tier1_neg.empty:
            tier_sep = separation_quality(tier1_pos["score"], tier1_neg["score"])
            sep["tier1_separation"] = tier_sep.get("separation", sep.get("tier1_separation"))
            sep["tier1_gate_passed"] = bool(
                tier1_pos["passed"].all() and tier1_neg["passed"].all() and tier_sep.get("separation", 0) > 0.50
            )
            sep["phase1_ready"] = sep["tier1_gate_passed"]

        scored_pos = positive_results[positive_results.get("burial_class", "") != "crystal_artifact"]
        sep["all_positive_passed"] = bool(scored_pos["passed"].fillna(False).all())
        sep["all_negative_passed"] = bool(negative_results["passed"].fillna(False).all())

        summary = {
            "positive_controls": positive_results,
            "negative_controls": negative_results,
            "separation_quality": sep,
        }

        if output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            positive_results.to_csv(output_dir / "positive_controls.csv", index=False)
            negative_results.to_csv(output_dir / "negative_controls.csv", index=False)
            combined = pd.concat([positive_results, negative_results], ignore_index=True)
            combined.to_csv(output_dir / "control_benchmark.csv", index=False)

        self._print_summary(summary)
        return summary

    def _print_summary(self, summary: Dict) -> None:
        """Print validation summary."""
        print("\n" + "=" * 60)
        print("VALIDATION SUMMARY")
        print("=" * 60)

        sep = summary.get("separation_quality") or {}
        if sep:
            print("\nScore Separation (burial-aware validation scores):")
            print(f"  Positive controls: {sep.get('positive_mean', 0):.3f}")
            print(f"  Negative controls: {sep.get('negative_mean', 0):.3f}")
            print(f"  Separation: {sep.get('separation', 0):.3f}")
            tier1 = sep.get("tier1_separation")
            tier1_display = f"{float(tier1):.3f}" if tier1 is not None else "n/a"
            print(f"  Tier-1 separation (ADAR2 vs PLCδ1): {tier1_display}")
            print(f"  Clear separation (>0.30): {sep.get('clear_separation', False)}")
            print(f"  Phase 1 gate passed: {sep.get('phase1_ready', False)}")
            print(f"  All positives passed: {sep.get('all_positive_passed', False)}")
            print(f"  All negatives passed: {sep.get('all_negative_passed', False)}")

        print("\n" + "=" * 60 + "\n")
