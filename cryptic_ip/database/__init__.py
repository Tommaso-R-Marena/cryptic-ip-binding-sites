"""
Database module for proteome data management.
"""

from .downloader import ProteomeDownloader
from .manager import ProteomeManager
from .batch_processing import (
    AnalysisCache,
    AlphaFoldBatchDownloader,
    ParallelProcessor,
    append_results_to_file,
)
from .integrity_checker import DatabaseIntegrityChecker

__all__ = [
    "ProteomeDownloader",
    "ProteomeManager",
    "AnalysisCache",
    "AlphaFoldBatchDownloader",
    "ParallelProcessor",
    "append_results_to_file",
    "DatabaseIntegrityChecker",
]
