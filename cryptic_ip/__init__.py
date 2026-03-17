"""Top-level package for cryptic IP binding site pipeline."""

from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Tommaso R. Marena"

__all__ = [
    "ProteinAnalyzer",
    "PocketScorer",
    "CrypticSiteMLClassifier",
    "MLPocketScorer",
    "FEATURE_COLUMNS",
    "validate_adar2",
    "ValidationSuite",
]


def __getattr__(name: str):
    if name in {"ProteinAnalyzer", "PocketScorer", "CrypticSiteMLClassifier", "MLPocketScorer", "FEATURE_COLUMNS"}:
        from .analysis import (
            ProteinAnalyzer,
            PocketScorer,
            CrypticSiteMLClassifier,
            MLPocketScorer,
            FEATURE_COLUMNS,
        )

        return {
            "ProteinAnalyzer": ProteinAnalyzer,
            "PocketScorer": PocketScorer,
            "CrypticSiteMLClassifier": CrypticSiteMLClassifier,
            "MLPocketScorer": MLPocketScorer,
            "FEATURE_COLUMNS": FEATURE_COLUMNS,
        }[name]
    if name in {"validate_adar2", "ValidationSuite"}:
        from .validation import ValidationSuite, validate_adar2

        return {"validate_adar2": validate_adar2, "ValidationSuite": ValidationSuite}[name]
    raise AttributeError(name)
