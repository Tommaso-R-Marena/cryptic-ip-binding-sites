"""Validation utilities for pocket scoring outputs."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from scipy.stats import normaltest

from .structure_validator import ValidationIssue, ValidationReport


class ResultsValidator:
    """Validate analysis output tables and serialized result files."""

    REQUIRED_FEATURES = {
        "pocket_id",
        "pocket_depth",
        "sasa",
        "volume",
        "composite_score",
    }

    NORMALIZED_SCORE_COLUMNS = {"composite_score", "normalized_score", "confidence_score"}

    def validate(self, results_path: str) -> ValidationReport:
        path = Path(results_path)
        issues: List[ValidationIssue] = []
        metrics: Dict[str, float] = {}

        frame = self._load_table(path, issues)
        if frame is None:
            return ValidationReport(valid=False, issues=issues, metrics=metrics)

        metrics["rows"] = float(len(frame))
        metrics["columns"] = float(len(frame.columns))

        self._check_required_features(frame, issues)
        self._check_score_ranges(frame, issues)
        self._check_suspicious_geometry(frame, issues)
        self._check_statistical_assumptions(frame, issues, metrics)

        is_valid = not any(issue.level == "error" for issue in issues)
        return ValidationReport(valid=is_valid, issues=issues, metrics=metrics)

    def validate_schema_consistency(self, paths: Sequence[str]) -> ValidationReport:
        issues: List[ValidationIssue] = []
        metrics: Dict[str, float] = {}
        schemas = {}

        for path in paths:
            loaded = self._load_table(Path(path), issues)
            if loaded is not None:
                schemas[path] = set(loaded.columns)

        if not schemas:
            return ValidationReport(valid=False, issues=issues, metrics=metrics)

        canonical_path = next(iter(schemas))
        canonical_schema = schemas[canonical_path]
        for path, schema in schemas.items():
            if schema != canonical_schema:
                missing = sorted(canonical_schema - schema)
                extra = sorted(schema - canonical_schema)
                issues.append(
                    ValidationIssue(
                        level="error",
                        check="schema_consistency",
                        message=f"Schema mismatch in {Path(path).name}: missing={missing}, extra={extra}",
                        suggested_fix="Regenerate outputs from a single pipeline version to keep schemas aligned.",
                    )
                )

        metrics["files_checked"] = float(len(schemas))
        is_valid = not any(issue.level == "error" for issue in issues)
        return ValidationReport(valid=is_valid, issues=issues, metrics=metrics)

    def _load_table(self, path: Path, issues: List[ValidationIssue]):
        try:
            if path.suffix.lower() == ".csv":
                return pd.read_csv(path)
            if path.suffix.lower() == ".json":
                with path.open("r", encoding="utf-8") as handle:
                    payload = json.load(handle)
                if isinstance(payload, list):
                    return pd.DataFrame(payload)
                if isinstance(payload, dict) and "results" in payload:
                    return pd.DataFrame(payload["results"])
                return pd.DataFrame(payload)
        except Exception as exc:
            issues.append(
                ValidationIssue(
                    level="error",
                    check="file_read",
                    message=f"Failed to parse {path.name}: {exc}",
                    suggested_fix="Ensure the output is valid CSV/JSON and not truncated.",
                )
            )
            return None

        issues.append(
            ValidationIssue(
                level="error",
                check="file_format",
                message=f"Unsupported results format: {path.suffix}",
                suggested_fix="Use .csv or .json output files for validation.",
            )
        )
        return None

    def _check_required_features(self, frame: pd.DataFrame, issues: List[ValidationIssue]) -> None:
        missing = sorted(self.REQUIRED_FEATURES - set(frame.columns))
        if missing:
            issues.append(
                ValidationIssue(
                    level="error",
                    check="required_features",
                    message=f"Missing required feature columns: {missing}",
                    suggested_fix="Re-run feature extraction to populate all pocket descriptors.",
                )
            )

    def _check_score_ranges(self, frame: pd.DataFrame, issues: List[ValidationIssue]) -> None:
        for col in self.NORMALIZED_SCORE_COLUMNS:
            if col not in frame.columns:
                continue
            bad = frame[(frame[col] < 0) | (frame[col] > 1)]
            if not bad.empty:
                issues.append(
                    ValidationIssue(
                        level="error",
                        check="score_range",
                        message=f"Column '{col}' has {len(bad)} values outside [0, 1].",
                        suggested_fix="Check score normalization logic and clamp/scale values to [0,1].",
                    )
                )

    def _check_suspicious_geometry(
        self, frame: pd.DataFrame, issues: List[ValidationIssue]
    ) -> None:
        for col in ("volume", "sasa", "pocket_depth"):
            if col in frame.columns and (frame[col] <= 0).any():
                issues.append(
                    ValidationIssue(
                        level="error",
                        check="impossible_geometry",
                        message=f"Found non-positive values in '{col}'.",
                        suggested_fix="Recompute geometric features; pocket volumes/depth/SASA must be positive.",
                    )
                )

        overlap_cols = {"overlap_fraction", "overlap_with_other_pockets"}
        for col in overlap_cols.intersection(frame.columns):
            high_overlap = frame[frame[col] > 0.7]
            if not high_overlap.empty:
                issues.append(
                    ValidationIssue(
                        level="warning",
                        check="pocket_overlap",
                        message=f"{len(high_overlap)} pockets have overlap > 0.7 in '{col}'.",
                        suggested_fix="Deduplicate or merge overlapping pockets before downstream ranking.",
                    )
                )

    def _check_statistical_assumptions(
        self,
        frame: pd.DataFrame,
        issues: List[ValidationIssue],
        metrics: Dict[str, float],
    ) -> None:
        if "composite_score" in frame.columns and len(frame) >= 8:
            stat, pvalue = normaltest(frame["composite_score"].dropna())
            metrics["normality_pvalue"] = float(pvalue)
            if pvalue < 0.05:
                issues.append(
                    ValidationIssue(
                        level="warning",
                        check="normality_assumption",
                        message="Composite scores significantly deviate from normality (p < 0.05).",
                        suggested_fix="Use non-parametric tests or transform scores before parametric testing.",
                    )
                )

        if "protein_file" in frame.columns and "pocket_id" in frame.columns:
            duplicates = frame.duplicated(subset=["protein_file", "pocket_id"]).sum()
            metrics["duplicate_pairs"] = float(duplicates)
            if duplicates > 0:
                issues.append(
                    ValidationIssue(
                        level="warning",
                        check="independence_assumption",
                        message=f"Found {duplicates} duplicated protein_file/pocket_id pairs.",
                        suggested_fix="Remove duplicate pocket observations before hypothesis testing.",
                    )
                )
