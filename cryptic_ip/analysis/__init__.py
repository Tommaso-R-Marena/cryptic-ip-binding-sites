"""
Analysis module for pocket detection and scoring.
"""

from .analyzer import ProteinAnalyzer
from .scorer import CalibratedRuleScorer, PocketScorer, ScoringParameters
from .filters import CandidateFilter
from .features import (
    FEATURE_NAMES,
    LEGACY_FEATURE_NAMES,
    PocketFeatureExtractor,
    PocketFeatures,
    feature_documentation,
)
from .labeling import (
    LigandSite,
    PocketAssignment,
    PocketLabel,
    assign_pocket_labels,
    summarise_labels,
)
from .ml_classifier import (
    DEFAULT_FEATURE_COLUMNS,
    FEATURE_COLUMNS,
    CrypticSiteMLClassifier,
    MLPocketScorer,
    classification_metrics,
    compare_models,
    delong_roc_test,
    model_comparison_table,
)
from .statistical_validation import StatisticalValidation, BootstrapCurveResult
from .comparative_analysis import ComparativeIPAnalysis, ComparativeResult
from .electrostatics import ElectrostaticsCalculator, PHAnalysisResult, run_apbs_wrapper, run_propka_wrapper

__all__ = [
    "ProteinAnalyzer",
    "PocketScorer",
    "CalibratedRuleScorer",
    "ScoringParameters",
    "CandidateFilter",
    "PocketFeatureExtractor",
    "PocketFeatures",
    "FEATURE_NAMES",
    "LEGACY_FEATURE_NAMES",
    "feature_documentation",
    "LigandSite",
    "PocketAssignment",
    "PocketLabel",
    "assign_pocket_labels",
    "summarise_labels",
    "CrypticSiteMLClassifier",
    "MLPocketScorer",
    "FEATURE_COLUMNS",
    "DEFAULT_FEATURE_COLUMNS",
    "classification_metrics",
    "compare_models",
    "delong_roc_test",
    "model_comparison_table",
    "StatisticalValidation",
    "BootstrapCurveResult",
    "ComparativeIPAnalysis",
    "ComparativeResult",
    "ElectrostaticsCalculator",
    "PHAnalysisResult",
    "run_propka_wrapper",
    "run_apbs_wrapper",
]
