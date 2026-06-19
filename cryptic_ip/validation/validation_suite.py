"""
Comprehensive validation suite for pipeline testing.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .adar2 import validate_adar2
from ..analysis import ProteinAnalyzer


class ValidationSuite:
    """
    Complete validation suite for testing pipeline on known examples.

    Positive controls: ADAR2, Pds5B, HDAC1 (should detect)
    Negative controls: PH domains (should reject)
    """

    POSITIVE_CONTROLS = {
        "ADAR2": {
            "pdb": "1ZY7",
            "ip_type": "IP6",
            "description": "Gold standard - completely buried IP6",
            "expected_score": 0.75,
            "use_adar2_validator": True,
        },
        "Pds5B": {
            "pdb": "5HDT",
            "ip_type": "IP6",
            "description": "Cohesin regulator with buried IP6",
            "expected_score": 0.65,
        },
        "HDAC1": {
            "pdb": "5ICN",
            "ip_type": "IP4",
            "description": "Histone deacetylase with IP4 at interface",
            "expected_score": 0.60,
        },
    }

    NEGATIVE_CONTROLS = {
        "PLCd1_PH": {
            "pdb": "1MAI",
            "ip_type": "IP3",
            "description": "Classic surface-exposed PH domain",
            "expected_score": 0.30,
        },
        "Btk_PH": {
            "pdb": "1BTK",
            "ip_type": "IP4",
            "description": "Kinase PH domain - membrane targeting",
            "expected_score": 0.35,
        },
    }

    POSITIVE_SCORE_THRESHOLD = 0.55
    NEGATIVE_SCORE_THRESHOLD = 0.50

    def __init__(self, data_dir: str = "data/validation"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
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

    def score_structure(
        self,
        pdb_path: Path,
        *,
        skip_electrostatics: bool = True,
    ) -> Dict:
        """Run pocket detection and scoring on a single structure."""
        analyzer = ProteinAnalyzer(
            str(pdb_path),
            work_dir=str(self.data_dir / "work" / pdb_path.stem),
            skip_electrostatics=skip_electrostatics,
        )
        scored = analyzer.run_pipeline(include_electrostatics=not skip_electrostatics)
        if scored.empty:
            raise RuntimeError(f"No pockets scored for {pdb_path.name}")

        top_pocket = scored.iloc[0]
        return {
            "total_pockets": int(len(scored)),
            "top_pocket_id": int(top_pocket["pocket_id"]),
            "top_pocket_score": float(top_pocket["composite_score"]),
            "volume": float(top_pocket["volume"]),
            "sasa": float(top_pocket["sasa"]),
            "basic_residues": int(top_pocket["basic_residues"]),
            "depth": float(top_pocket["depth"]),
            "structure_path": str(pdb_path),
        }

    def _positive_passed(self, score: float, metrics: Dict) -> bool:
        return score >= self.POSITIVE_SCORE_THRESHOLD and metrics["basic_residues"] >= 3

    def _negative_passed(self, score: float) -> bool:
        return score < self.NEGATIVE_SCORE_THRESHOLD

    def validate_positive_control(self, name: str, info: Dict) -> Dict:
        """Validate a single positive control structure."""
        if info.get("use_adar2_validator"):
            result = validate_adar2(use_alphafold=False)
            return {
                "protein": name,
                "control_type": "positive",
                "pdb_id": info["pdb"],
                "ip_type": info["ip_type"],
                "score": result["top_pocket_score"],
                "expected": info["expected_score"],
                "passed": result["validation_passed"],
                "total_pockets": result["total_pockets"],
                "sasa": result["sasa"],
                "basic_residues": result["basic_residues"],
                "description": info["description"],
            }

        pdb_path = self.download_pdb(info["pdb"])
        metrics = self.score_structure(pdb_path)
        passed = self._positive_passed(metrics["top_pocket_score"], metrics)
        return {
            "protein": name,
            "control_type": "positive",
            "pdb_id": info["pdb"],
            "ip_type": info["ip_type"],
            "score": metrics["top_pocket_score"],
            "expected": info["expected_score"],
            "passed": passed,
            "total_pockets": metrics["total_pockets"],
            "sasa": metrics["sasa"],
            "basic_residues": metrics["basic_residues"],
            "description": info["description"],
        }

    def validate_negative_control(self, name: str, info: Dict) -> Dict:
        """Validate a single negative control structure."""
        pdb_path = self.download_pdb(info["pdb"])
        metrics = self.score_structure(pdb_path)
        passed = self._negative_passed(metrics["top_pocket_score"])
        return {
            "protein": name,
            "control_type": "negative",
            "pdb_id": info["pdb"],
            "ip_type": info["ip_type"],
            "score": metrics["top_pocket_score"],
            "expected": info["expected_score"],
            "passed": passed,
            "total_pockets": metrics["total_pockets"],
            "sasa": metrics["sasa"],
            "basic_residues": metrics["basic_residues"],
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
                print(f"  Score: {results[-1]['score']:.3f} ({status})")
            except Exception as exc:
                print(f"  Error: {exc}")
                results.append(
                    {
                        "protein": name,
                        "control_type": "positive",
                        "pdb_id": info.get("pdb"),
                        "ip_type": info["ip_type"],
                        "score": float("nan"),
                        "expected": info["expected_score"],
                        "passed": False,
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
                print(f"  Score: {results[-1]['score']:.3f} ({status})")
            except Exception as exc:
                print(f"  Error: {exc}")
                results.append(
                    {
                        "protein": name,
                        "control_type": "negative",
                        "pdb_id": info.get("pdb"),
                        "ip_type": info["ip_type"],
                        "score": float("nan"),
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

        summary = {
            "positive_controls": positive_results,
            "negative_controls": negative_results,
            "separation_quality": self._calculate_separation(positive_results, negative_results),
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

    def _calculate_separation(self, positive: pd.DataFrame, negative: pd.DataFrame) -> Dict:
        """Calculate score separation between positive and negative controls."""
        pos_scores = positive["score"].dropna()
        neg_scores = negative["score"].dropna()
        if len(pos_scores) == 0 or len(neg_scores) == 0:
            return {}

        pos_mean = float(pos_scores.mean())
        neg_mean = float(neg_scores.mean())
        separation = pos_mean - neg_mean

        return {
            "positive_mean": pos_mean,
            "negative_mean": neg_mean,
            "separation": separation,
            "clear_separation": separation > 0.3,
            "all_positive_passed": bool(positive["passed"].fillna(False).all()),
            "all_negative_passed": bool(negative["passed"].fillna(False).all()),
        }

    def _print_summary(self, summary: Dict) -> None:
        """Print validation summary."""
        print("\n" + "=" * 60)
        print("VALIDATION SUMMARY")
        print("=" * 60)

        sep = summary.get("separation_quality") or {}
        if sep:
            print("\nScore Separation:")
            print(f"  Positive controls: {sep.get('positive_mean', 0):.3f}")
            print(f"  Negative controls: {sep.get('negative_mean', 0):.3f}")
            print(f"  Separation: {sep.get('separation', 0):.3f}")
            print(f"  Clear separation: {sep.get('clear_separation', False)}")
            print(f"  All positives passed: {sep.get('all_positive_passed', False)}")
            print(f"  All negatives passed: {sep.get('all_negative_passed', False)}")

        print("\n" + "=" * 60 + "\n")
