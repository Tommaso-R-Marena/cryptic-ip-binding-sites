"""
Protein structure analysis and pocket detection.
"""

import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from prody import parsePDB, calcSASA

from .electrostatics import ElectrostaticsCalculator
from .features import FEATURE_NAMES, PocketFeatureExtractor, legacy_feature_row
from .fpocket_parser import FpocketParser
from .ml_classifier import CrypticSiteMLClassifier, MLPocketScorer
from .scorer import PocketScorer
from .structure_arrays import load_structure_arrays
from ..utils.profiling import timed
from ..utils.resources import cleanup_files, set_memory_limit


class ProteinAnalyzer:
    """
    Main analyzer class for detecting and characterizing pockets in protein structures.

    Integrates fpocket for pocket detection, FreeSASA for accessibility,
    and APBS for electrostatics.
    """

    def __init__(
        self,
        pdb_path: str,
        work_dir: Optional[str] = None,
        scorer: Optional[Any] = None,
        use_ml_model: bool = False,
        model_path: Optional[str] = None,
        skip_electrostatics: bool = False,
        memory_limit_gb: Optional[float] = None,
    ):
        """
        Initialize analyzer with a protein structure.

        Args:
            pdb_path: Path to PDB file
            work_dir: Working directory for outputs (temp if None)
            scorer: Optional scorer object implementing
                ``calculate_composite_score`` and ``classify_site``.
                Uses threshold-based :class:`PocketScorer` when omitted.
            use_ml_model: When True, attempt to load a pre-trained ML scorer.
            model_path: Optional path to serialized ML model.
        """
        self.pdb_path = Path(pdb_path).resolve()
        self.work_dir = Path(work_dir) if work_dir else Path(tempfile.mkdtemp())
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # Parse structure
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure("protein", str(self.pdb_path))

        # Initialize components
        self.fpocket_parser = FpocketParser()
        self.model_path = self._resolve_model_path(model_path)
        self.scorer = scorer or self._build_default_scorer(use_ml_model)
        self.skip_electrostatics = skip_electrostatics
        if memory_limit_gb:
            set_memory_limit(memory_limit_gb)

        # Storage for results
        self.pockets = None
        self.sasa_data = None
        self.electrostatic_data = None
        self.electrostatic_map_path: Optional[Path] = None
        self._surface_atom_coords: Optional[np.ndarray] = None
        self._arrays = None
        self._feature_extractor: Optional[PocketFeatureExtractor] = None
        self._pocket_residue_cache: Dict[Tuple[int, float], List[int]] = {}

    @property
    def feature_extractor(self) -> PocketFeatureExtractor:
        """Chain-aware descriptor extractor for this structure.

        Built lazily and cached: it computes whole-structure SASA once, which
        dominates the cost and is shared by every pocket.
        """
        if self._feature_extractor is None:
            if self._arrays is None:
                self._arrays = load_structure_arrays(self.pdb_path)
            self._feature_extractor = PocketFeatureExtractor(self._arrays)
        return self._feature_extractor

    @timed()
    def run_pipeline(self, include_electrostatics: Optional[bool] = None) -> pd.DataFrame:
        """Run the full single-protein pipeline with intra-protein parallelization."""
        self.detect_pockets()
        use_electrostatics = (not self.skip_electrostatics) if include_electrostatics is None else include_electrostatics

        jobs = [self.calculate_sasa]
        if use_electrostatics:
            jobs.append(self.calculate_electrostatics)

        with ThreadPoolExecutor(max_workers=len(jobs)) as executor:
            futures = [executor.submit(job) for job in jobs]
            for future in futures:
                future.result()

        return self.score_all_pockets()

    def _resolve_model_path(self, model_path: Optional[str]) -> Path:
        if model_path:
            return Path(model_path)
        return Path(__file__).resolve().parents[2] / "models" / "cryptic_ip_classifier_v1.pkl"

    #: A model whose recorded out-of-fold AUROC is at or below this value has no
    #: demonstrated skill, and using it is worse than the transparent rule-based
    #: score, which at least means something.
    MIN_USABLE_ROC_AUC = 0.55

    def _recorded_model_roc_auc(self) -> Optional[float]:
        """Read the out-of-fold AUROC recorded alongside the model, if any.

        Returns:
            The recorded AUROC, or ``None`` when no metadata is available.
        """
        metadata_path = self.model_path.with_suffix("").with_suffix(".metadata.json")
        if not metadata_path.exists():
            metadata_path = self.model_path.parent / f"{self.model_path.stem}.metadata.json"
        if not metadata_path.exists():
            return None
        try:
            import json

            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return None

        selected = metadata.get("selected_model_type")
        candidates = metadata.get("model_candidates") or {}
        entry = candidates.get(selected, {})
        for source in (entry.get("metrics", {}), entry, metadata):
            value = source.get("roc_auc") if isinstance(source, dict) else None
            if isinstance(value, (int, float)):
                return float(value)
        return None

    def _build_default_scorer(self, use_ml_model: bool) -> Any:
        """Choose the scorer, refusing to deploy a model with no demonstrated skill.

        Silently loading a model that performed at chance during validation is the
        worst outcome available: it produces confident-looking numbers with no
        information behind them. When the recorded metrics show no skill, this
        falls back to the interpretable rule-based score and says why.
        """
        if not use_ml_model:
            return PocketScorer()

        try:
            classifier = CrypticSiteMLClassifier.load(str(self.model_path))
        except Exception as exc:  # noqa: BLE001 - fallback by design
            print(
                f"Warning: Unable to load ML model from {self.model_path} ({exc}). "
                "Falling back to rule-based scoring."
            )
            return PocketScorer()

        recorded_auc = self._recorded_model_roc_auc()
        if recorded_auc is not None and recorded_auc <= self.MIN_USABLE_ROC_AUC:
            print(
                f"Warning: the model at {self.model_path} recorded a validation AUROC of "
                f"{recorded_auc:.3f}, which is not distinguishable from chance. Falling back "
                "to rule-based scoring. Retrain with scripts/train_ml_classifier.py, or pass "
                "an explicit --model-path to override this check."
            )
            return PocketScorer()

        return MLPocketScorer(classifier)

    @timed()
    def detect_pockets(self, min_alpha_sphere: int = 3) -> pd.DataFrame:
        """
        Detect pockets using fpocket.

        Args:
            min_alpha_sphere: Minimum number of alpha spheres for a pocket

        Returns:
            DataFrame with pocket properties
        """
        # Check if fpocket is available
        try:
            subprocess.run(["fpocket", "-h"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError(
                "fpocket not found. Please install from https://github.com/Discngine/fpocket"
            )

        # fpocket writes output adjacent to the input structure file
        output_dir = self.pdb_path.parent / f"{self.pdb_path.stem}_out"
        cmd = ["fpocket", "-f", str(self.pdb_path), "-m", str(min_alpha_sphere)]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(f"fpocket failed: {result.stderr}")

        # Parse fpocket output
        info_file = output_dir / f"{self.pdb_path.stem}_info.txt"
        if not info_file.exists():
            raise RuntimeError(f"fpocket output not found: {info_file}")

        self.pockets = self.fpocket_parser.parse_info_file(info_file)
        self.pockets["pdb_path"] = str(self.pdb_path)

        return self.pockets

    @timed()
    def calculate_sasa(self) -> Dict[int, float]:
        """
        Calculate solvent accessible surface area for all residues.

        Returns:
            Dictionary mapping residue number to SASA value
        """
        try:
            from Bio.PDB.SASA import ShrakeRupley

            ShrakeRupley().compute(self.structure, level="R")
            residue_sasa: Dict[int, float] = {}
            for model in self.structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] != " ":
                            continue
                        residue_sasa[residue.id[1]] = float(getattr(residue, "sasa", 0.0))

            self.sasa_data = residue_sasa
            return residue_sasa

        except Exception as e:
            raise RuntimeError(f"SASA calculation failed: {e}")

    def _compute_surface_atom_coords(self, sasa_threshold: float = 1.0) -> np.ndarray:
        """
        Return coordinates of solvent-exposed protein atoms.

        An atom is considered part of the molecular surface when its
        Shrake-Rupley solvent accessible surface area exceeds ``sasa_threshold``
        Å². Heteroatoms/ligands (e.g. bound inositol phosphates) are excluded so
        that burial is measured relative to the protein surface only.

        Args:
            sasa_threshold: Minimum per-atom SASA (Å²) to count as exposed.

        Returns:
            Array of shape ``(n_surface_atoms, 3)`` with atom coordinates. Falls
            back to all protein atom coordinates if no atom clears the threshold.
        """
        from Bio.PDB.SASA import ShrakeRupley

        ShrakeRupley().compute(self.structure, level="A")

        surface_coords: List[np.ndarray] = []
        all_protein_coords: List[np.ndarray] = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    # Skip heteroatoms/ligands; only protein residues define the surface.
                    if residue.id[0] != " ":
                        continue
                    for atom in residue:
                        all_protein_coords.append(atom.get_coord())
                        if float(getattr(atom, "sasa", 0.0)) > sasa_threshold:
                            surface_coords.append(atom.get_coord())

        chosen = surface_coords or all_protein_coords
        if not chosen:
            return np.empty((0, 3), dtype=float)
        return np.asarray(chosen, dtype=float)

    @timed()
    def calculate_pocket_burial_depth(
        self, pocket_center: Tuple[float, float, float]
    ) -> float:
        """
        Geometric burial depth of a pocket center below the protein surface.

        Depth is defined as the Euclidean distance (Å) from the pocket center to
        the nearest solvent-exposed protein atom. Deeply buried structural pockets
        (e.g. the ADAR2 IP6 cavity) sit far from any exposed atom and yield large
        depths (>15 Å), whereas surface signaling pockets sit adjacent to exposed
        atoms and yield small depths. This implements the pipeline's
        "pocket depth from surface" criterion as a true geometric measurement,
        rather than fpocket's local hydrophobic density proxy.

        Args:
            pocket_center: (x, y, z) coordinate of the pocket center.

        Returns:
            Burial depth in Angstroms (0.0 when no surface atoms are available).
        """
        if self._surface_atom_coords is None:
            self._surface_atom_coords = self._compute_surface_atom_coords()

        if self._surface_atom_coords.size == 0:
            return 0.0

        center = np.asarray(pocket_center, dtype=float)
        distances = np.linalg.norm(self._surface_atom_coords - center, axis=1)
        return float(np.min(distances))

    @timed()
    def calculate_electrostatics(self, ph: float = 7.4) -> Optional[float]:
        """
        Calculate electrostatic potential using pdb2pqr + APBS.

        Args:
            ph: Solution pH used during protonation assignment

        Returns:
            Scalar APBS electrostatic energy proxy, or ``None`` on failure.
        """
        calculator = ElectrostaticsCalculator()
        electro_dir = self.work_dir / "electrostatics"

        try:
            pqr_path = calculator.generate_pqr(self.pdb_path, ph=ph, output_dir=electro_dir)
            potential, dx_path = calculator.run_apbs_with_map(pqr_path=pqr_path, output_dir=electro_dir)
            self.electrostatic_map_path = dx_path
        except RuntimeError as exc:
            print(f"Warning: Electrostatics calculation failed: {exc}")
            return None

        self.electrostatic_data = potential
        return potential

    def pocket_electrostatic_potential(self, pocket_center: Tuple[float, float, float]) -> Optional[float]:
        """Sample APBS potential at a pocket center when a map is available."""
        if self.electrostatic_map_path is None or not self.electrostatic_map_path.exists():
            return self.electrostatic_data
        calculator = ElectrostaticsCalculator()
        try:
            return calculator.sample_potential_at_point(self.electrostatic_map_path, pocket_center)
        except Exception:
            return self.electrostatic_data

    def pocket_plddt_confidence(self, pocket_residues: List[int]) -> float:
        """Mean pLDDT (B-factor) for pocket-lining residues in AlphaFold models."""
        from ..validation.plddt import pocket_plddt_confidence

        summary = pocket_plddt_confidence(self.pdb_path, pocket_residues)
        return float(summary["plddt_mean"])

    @timed()
    def get_pocket_residues(self, pocket_id: int, distance_cutoff: float = 8.0) -> List[int]:
        """
        Get residue numbers lining a pocket.

        Uses fpocket pocket atom residue IDs when available, otherwise falls back
        to CA atoms within ``distance_cutoff`` of the pocket alpha-sphere centroid.

        Results are memoised per analyzer instance. An ``lru_cache`` on the method
        was previously used, which had two defects: the cache was keyed on ``self``
        and so kept every analyzer (and its parsed structure) alive for the life of
        the process, and stacking it under ``@timed`` hid ``cache_clear`` from
        :meth:`cleanup`, which therefore raised ``AttributeError`` on every call.

        Note the residue numbers here are chain-blind, and are retained for
        callers that expect that. Descriptor extraction uses chain-aware residue
        keys via :mod:`cryptic_ip.analysis.features`.
        """
        if self.pockets is None:
            raise ValueError("Run detect_pockets() first")

        cache_key = (int(pocket_id), float(distance_cutoff))
        cached = self._pocket_residue_cache.get(cache_key)
        if cached is not None:
            return list(cached)

        pocket = self.pockets[self.pockets["pocket_id"] == pocket_id].iloc[0]
        if "fpocket_residue_ids" in pocket and pd.notna(pocket["fpocket_residue_ids"]):
            fpocket_ids = [
                int(token)
                for token in str(pocket["fpocket_residue_ids"]).split(",")
                if token.strip()
            ]
            if fpocket_ids:
                self._pocket_residue_cache[cache_key] = fpocket_ids
                return list(fpocket_ids)

        center = np.array([pocket["center_x"], pocket["center_y"], pocket["center_z"]])
        residue_numbers = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if "CA" not in residue:
                        continue
                    ca_coord = residue["CA"].get_coord()
                    distance = np.linalg.norm(ca_coord - center)
                    if distance <= distance_cutoff:
                        residue_numbers.append(residue.id[1])

        resolved = sorted(set(residue_numbers))
        self._pocket_residue_cache[cache_key] = resolved
        return list(resolved)

    def count_basic_residues(self, pocket_id: int, distance_cutoff: float = 5.0) -> int:
        """
        Count basic residues (Arg, Lys, His) near pocket.

        Args:
            pocket_id: Pocket identifier
            distance_cutoff: Distance threshold in Angstroms

        Returns:
            Number of basic residues
        """
        pocket_residues = self.get_pocket_residues(pocket_id, distance_cutoff)
        basic_count = 0

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.id[1] in pocket_residues:
                        if residue.get_resname() in ["ARG", "LYS", "HIS"]:
                            basic_count += 1

        return basic_count

    def _pocket_alpha_spheres(self, pocket_id: int) -> Optional[np.ndarray]:
        """Return alpha-sphere centres for a pocket, if fpocket wrote them.

        Alpha-sphere centres describe the cavity itself, whereas the pocket atom
        file describes the residues lining it. Shape and hull-volume descriptors
        are computed from the cavity, so the vertex file is preferred.

        Args:
            pocket_id: Pocket identifier.

        Returns:
            Coordinates, or ``None`` when no pocket file is available.
        """
        pockets_dir = self.pdb_path.parent / f"{self.pdb_path.stem}_out" / "pockets"
        for filename in (f"pocket{int(pocket_id)}_vert.pqr", f"pocket{int(pocket_id)}_atm.pdb"):
            path = pockets_dir / filename
            if path.exists():
                atoms = self.fpocket_parser.parse_pocket_atoms(path)
                if atoms:
                    return np.asarray(
                        [[atom["x"], atom["y"], atom["z"]] for atom in atoms], dtype=float
                    )
        return None

    @timed()
    def analyze_pocket(self, pocket_id: int) -> Dict:
        """
        Complete analysis of a single pocket.

        Returns the full descriptor suite from
        :mod:`cryptic_ip.analysis.features` alongside the legacy keys.

        Note on ``depth``: it now carries the **geometric burial depth** -
        distance from the pocket centre to the nearest solvent-exposed atom.
        Earlier versions populated it from fpocket's mean local hydrophobic
        density, a composition statistic unrelated to depth, while computing the
        real depth and discarding it. fpocket's value is still reported, under
        the accurate name ``mean_local_hydrophobic_density``.

        Args:
            pocket_id: Pocket identifier

        Returns:
            Dictionary with all pocket metrics
        """
        if self.pockets is None:
            self.detect_pockets()

        if self.sasa_data is None:
            self.calculate_sasa()

        pocket = self.pockets[self.pockets["pocket_id"] == pocket_id].iloc[0]
        pocket_residues = self.get_pocket_residues(pocket_id)
        pocket_center = (
            float(pocket["center_x"]),
            float(pocket["center_y"]),
            float(pocket["center_z"]),
        )

        pocket_potential = self.pocket_electrostatic_potential(pocket_center)
        fpocket_volume = pocket.get("volume", None)
        if fpocket_volume is not None and not pd.isna(fpocket_volume):
            fpocket_volume = float(fpocket_volume)
        else:
            fpocket_volume = None

        descriptors = self.feature_extractor.extract(
            int(pocket_id),
            pocket_center,
            alpha_sphere_coords=self._pocket_alpha_spheres(pocket_id),
            fpocket_volume=fpocket_volume,
            apbs_potential=pocket_potential,
        )
        features = dict(descriptors.features)

        result: Dict[str, Any] = {"pocket_id": int(pocket_id)}
        result.update(features)
        # Legacy keys retained for existing callers and stored result schemas.
        result.update(
            {
                "volume": fpocket_volume if fpocket_volume is not None else np.nan,
                "depth": features.get("burial_depth", np.nan),
                "burial_depth": features.get("burial_depth", np.nan),
                "sasa": features.get("sasa_mean", np.nan),
                "basic_residues": features.get("n_basic_residues", np.nan),
                "residue_count": len(descriptors.residue_keys),
                "electrostatic_potential": pocket_potential,
                "plddt_confidence": features.get("plddt_mean", np.nan),
                "mean_local_hydrophobic_density": pocket.get(
                    "mean_local_hydrophobic_density", np.nan
                ),
                "center": pocket_center,
                "pocket_residue_numbers": pocket_residues,
            }
        )
        return result

    @timed()
    def score_all_pockets(self) -> pd.DataFrame:
        """
        Score all detected pockets for cryptic IP binding.

        Returns:
            DataFrame with scores for all pockets
        """
        if self.pockets is None:
            self.detect_pockets()

        if self.sasa_data is None:
            self.calculate_sasa()

        results = []
        for pocket_id in self.pockets["pocket_id"]:
            try:
                results.append(self.analyze_pocket(pocket_id))
            except Exception as e:
                print(f"Warning: Failed to analyze pocket {pocket_id}: {e}")
                continue

        if not results:
            return pd.DataFrame(results)

        if hasattr(self.scorer, "calculate_composite_scores"):
            # Supply both the full descriptor schema and the legacy column names,
            # so a model trained on either schema can be scored without the
            # caller needing to know which one it was.
            samples = pd.DataFrame(
                [
                    {
                        **{name: row.get(name, np.nan) for name in FEATURE_NAMES},
                        **legacy_feature_row(row),
                    }
                    for row in results
                ]
            )
            scores = self.scorer.calculate_composite_scores(samples)
            for row, score in zip(results, scores):
                row["composite_score"] = float(score)
        else:
            for row in results:
                row["composite_score"] = self.scorer.calculate_composite_score(
                    volume=row.get("pocket_volume", row.get("volume")),
                    depth=row.get("burial_depth"),
                    sasa=row.get("sasa_mean", row.get("sasa")),
                    basic_count=row.get("n_basic_residues", row.get("basic_residues")),
                    potential=(
                        row.get("electrostatic_potential")
                        if row.get("electrostatic_potential") is not None
                        and np.isfinite(
                            np.asarray(row.get("electrostatic_potential", np.nan), dtype=float)
                        )
                        else row.get("coulomb_potential_kt")
                    ),
                    enclosure=row.get("enclosure"),
                )

        return pd.DataFrame(results).sort_values("composite_score", ascending=False)

    @timed()
    def cleanup(self) -> None:
        """Clear analyzer intermediates from disk and reset in-memory caches."""
        cleanup_files([self.pdb_path.parent / f"{self.pdb_path.stem}_out", self.work_dir / "electrostatics"])
        self._pocket_residue_cache.clear()
