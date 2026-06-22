"""High-level pipeline wrappers documented in README and tutorials."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .analysis import ProteinAnalyzer
from .analysis.filters import CandidateFilter
from .database.batch_processing import AlphaFoldBatchDownloader


class AnalysisPipeline:
    """Single-structure analysis pipeline (README-compatible API)."""

    def __init__(
        self,
        work_dir: Optional[Union[str, Path]] = None,
        score_threshold: float = 0.75,
        use_ml_model: bool = False,
        model_path: Optional[str] = None,
        skip_electrostatics: bool = True,
    ) -> None:
        self.work_dir = Path(work_dir) if work_dir else Path("results/pipeline")
        self.score_threshold = score_threshold
        self.use_ml_model = use_ml_model
        self.model_path = model_path
        self.skip_electrostatics = skip_electrostatics
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def analyze(
        self,
        structure_path: Union[str, Path],
        *,
        include_electrostatics: bool = False,
        min_plddt: float = 70.0,
    ) -> Dict[str, Any]:
        """Run pocket detection, scoring, and cryptic-site filtering on one structure."""
        structure_path = Path(structure_path)
        analyzer = ProteinAnalyzer(
            str(structure_path),
            work_dir=str(self.work_dir / structure_path.stem),
            use_ml_model=self.use_ml_model,
            model_path=self.model_path,
            skip_electrostatics=self.skip_electrostatics and not include_electrostatics,
        )
        scored = analyzer.run_pipeline(include_electrostatics=include_electrostatics)
        filt = CandidateFilter(min_score=self.score_threshold, min_plddt=min_plddt)
        candidates = filt.filter_cryptic_candidates(scored, structure_path=str(structure_path))
        return {
            "structure": str(structure_path),
            "pockets_detected": int(len(scored)),
            "candidates": candidates,
            "top_candidate": candidates.iloc[0].to_dict() if not candidates.empty else None,
        }


class ScreeningPipeline:
    """Proteome-scale screening wrapper around structure download and batch analysis."""

    def __init__(
        self,
        output_dir: Union[str, Path] = "results/screen",
        score_threshold: float = 0.75,
        min_plddt: float = 70.0,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.score_threshold = score_threshold
        self.min_plddt = min_plddt
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._analysis = AnalysisPipeline(
            work_dir=self.output_dir / "work",
            score_threshold=score_threshold,
            skip_electrostatics=True,
        )

    def download_proteome(self, proteome_id: str, limit: Optional[int] = None) -> List[Path]:
        """Download AlphaFold structures for a UniProt proteome ID."""
        structures_dir = self.output_dir / "structures"
        downloader = AlphaFoldBatchDownloader(output_dir=structures_dir)
        uniprot_ids = downloader.fetch_proteome_uniprot_ids(proteome_id)
        if limit is not None:
            uniprot_ids = uniprot_ids[:limit]
        paths: List[Path] = []
        for uniprot_id in uniprot_ids:
            try:
                paths.append(Path(downloader.af_client.fetch_structure(uniprot_id)))
            except Exception:
                continue
        return paths

    def screen_structures(self, structure_paths: List[Union[str, Path]]) -> pd.DataFrame:
        """Screen a list of structures and return concatenated hit table."""
        rows: List[Dict[str, Any]] = []
        for path in structure_paths:
            result = self._analysis.analyze(path, min_plddt=self.min_plddt)
            for hit in result["candidates"].to_dict(orient="records"):
                hit["structure_path"] = str(path)
                rows.append(hit)
        hits = pd.DataFrame(rows)
        hits.to_csv(self.output_dir / "screen_hits.csv", index=False)
        return hits
