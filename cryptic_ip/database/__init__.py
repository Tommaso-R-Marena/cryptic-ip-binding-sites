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
    "RcsbClient",
    "RcsbUnavailableError",
    "IPLigand",
    "discover_ip_ligands",
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
    if name in {"RcsbClient", "RcsbUnavailableError"}:
        from .rcsb_client import RcsbClient, RcsbUnavailableError

        return {"RcsbClient": RcsbClient, "RcsbUnavailableError": RcsbUnavailableError}[name]
    if name in {"IPLigand", "discover_ip_ligands"}:
        from .ip_ligands import IPLigand, discover_ip_ligands

        return {"IPLigand": IPLigand, "discover_ip_ligands": discover_ip_ligands}[name]
    raise AttributeError(name)
