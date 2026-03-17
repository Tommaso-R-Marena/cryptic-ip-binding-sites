"""Database module for proteome data management."""

from __future__ import annotations

__all__ = [
    "ProteomeDownloader",
    "ProteomeManager",
    "AnalysisCache",
    "AlphaFoldBatchDownloader",
    "ParallelProcessor",
    "append_results_to_file",
    "DatabaseIntegrityChecker",
]


def __getattr__(name: str):
    if name == "ProteomeDownloader":
        from .downloader import ProteomeDownloader

        return ProteomeDownloader
    if name == "ProteomeManager":
        from .manager import ProteomeManager

        return ProteomeManager
    if name in {"AnalysisCache", "AlphaFoldBatchDownloader", "ParallelProcessor", "append_results_to_file"}:
        from .batch_processing import (
            AlphaFoldBatchDownloader,
            AnalysisCache,
            ParallelProcessor,
            append_results_to_file,
        )

        return {
            "AnalysisCache": AnalysisCache,
            "AlphaFoldBatchDownloader": AlphaFoldBatchDownloader,
            "ParallelProcessor": ParallelProcessor,
            "append_results_to_file": append_results_to_file,
        }[name]
    if name == "DatabaseIntegrityChecker":
        from .integrity_checker import DatabaseIntegrityChecker

        return DatabaseIntegrityChecker
    raise AttributeError(name)
