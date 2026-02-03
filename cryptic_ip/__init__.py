"""
Cryptic IP Binding Sites
========================

Computational pipeline for identifying buried inositol phosphate binding sites in protein structures.

Modules:
    analysis: Pocket detection, scoring, and filtering
    database: Proteome data management and downloads
    validation: Known site validation and benchmarking
    visualization: Structure visualization and plotting
"""

__version__ = "0.1.0"
__author__ = "Tommaso R. Marena"

from .analysis import ProteinAnalyzer, PocketScorer
from .validation import validate_adar2, ValidationSuite

__all__ = [
    "ProteinAnalyzer",
    "PocketScorer",
    "validate_adar2",
    "ValidationSuite",
]
