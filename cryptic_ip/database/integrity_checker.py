"""Database and cache integrity checks for Cryptic-IP pipeline."""

from __future__ import annotations

import hashlib
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

from ..validation.structure_validator import ValidationIssue, ValidationReport
from .uniprot_client import UniProtClient


class DatabaseIntegrityChecker:
    """Validate SQLite cache consistency and external data references."""

    def __init__(self, max_alphafold_age_days: int = 365) -> None:
        self.max_alphafold_age_days = max_alphafold_age_days
        self.uniprot = UniProtClient()

    def validate(
        self,
        db_path: str,
        expected_checksums: Optional[Dict[str, str]] = None,
    ) -> ValidationReport:
        issues: List[ValidationIssue] = []
        metrics: Dict[str, float] = {}

        conn = self._open_database(db_path, issues)
        if conn is None:
            return ValidationReport(valid=False, issues=issues, metrics=metrics)

        with conn:
            self._check_required_tables(conn, issues)
            self._check_uniprot_ids(conn, issues, metrics)
            self._check_checksums(conn, issues, expected_checksums or {})
            self._check_alphafold_dates(conn, issues, metrics)

        is_valid = not any(issue.level == "error" for issue in issues)
        return ValidationReport(valid=is_valid, issues=issues, metrics=metrics)

    def _open_database(self, db_path: str, issues: List[ValidationIssue]):
        try:
            conn = sqlite3.connect(db_path)
            conn.execute("PRAGMA integrity_check")
            integrity = conn.execute("PRAGMA integrity_check").fetchone()[0]
            if integrity != "ok":
                issues.append(
                    ValidationIssue(
                        level="error",
                        check="sqlite_integrity",
                        message=f"SQLite integrity check failed: {integrity}",
                        suggested_fix="Restore the cache from backup or rebuild the SQLite database.",
                    )
                )
            return conn
        except sqlite3.DatabaseError as exc:
            issues.append(
                ValidationIssue(
                    level="error",
                    check="sqlite_open",
                    message=f"Could not open SQLite database: {exc}",
                    suggested_fix="Ensure the DB path is correct and file is a valid SQLite database.",
                )
            )
            return None

    def _check_required_tables(self, conn, issues: List[ValidationIssue]) -> None:
        existing = {
            row[0]
            for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        }
        required = {"structures"}
        missing = sorted(required - existing)
        if missing:
            issues.append(
                ValidationIssue(
                    level="error",
                    check="required_tables",
                    message=f"Missing required tables: {missing}",
                    suggested_fix="Run database migration/initialization to create required cache tables.",
                )
            )

    def _check_uniprot_ids(
        self, conn, issues: List[ValidationIssue], metrics: Dict[str, float]
    ) -> None:
        if not self._table_exists(conn, "structures"):
            return

        rows = conn.execute(
            "SELECT DISTINCT uniprot_id FROM structures WHERE uniprot_id IS NOT NULL"
        ).fetchall()
        ids = [row[0] for row in rows if row[0]]
        invalid = []
        for uniprot_id in ids[:20]:
            try:
                self.uniprot.get_protein_info(uniprot_id)
            except Exception:
                invalid.append(uniprot_id)

        metrics["uniprot_checked"] = float(min(len(ids), 20))
        if invalid:
            issues.append(
                ValidationIssue(
                    level="warning",
                    check="uniprot_mapping",
                    message=f"Invalid or unavailable UniProt IDs detected: {invalid[:5]}",
                    suggested_fix="Correct stale accessions in the cache and re-sync metadata.",
                )
            )

    def _check_checksums(
        self,
        conn,
        issues: List[ValidationIssue],
        expected_checksums: Dict[str, str],
    ) -> None:
        if not expected_checksums:
            return

        for file_path, expected in expected_checksums.items():
            path = Path(file_path)
            if not path.exists():
                issues.append(
                    ValidationIssue(
                        level="error",
                        check="checksum_file_exists",
                        message=f"Missing structure for checksum validation: {file_path}",
                        suggested_fix="Re-download the missing structure and update cache entries.",
                    )
                )
                continue
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest != expected:
                issues.append(
                    ValidationIssue(
                        level="error",
                        check="checksum_match",
                        message=f"Checksum mismatch for {path.name}",
                        suggested_fix="Delete and re-download structure to eliminate partial/corrupt files.",
                    )
                )

    def _check_alphafold_dates(
        self, conn, issues: List[ValidationIssue], metrics: Dict[str, float]
    ) -> None:
        if not self._table_exists(conn, "structures"):
            return

        columns = {row[1] for row in conn.execute("PRAGMA table_info(structures)").fetchall()}
        date_columns = [
            col for col in ("model_date", "download_date", "updated_at") if col in columns
        ]
        if not date_columns:
            return

        col = date_columns[0]
        rows = conn.execute(
            f"SELECT uniprot_id, {col} FROM structures WHERE {col} IS NOT NULL"
        ).fetchall()
        stale = []
        now = datetime.now(timezone.utc)
        for uniprot_id, raw_date in rows:
            try:
                dt = datetime.fromisoformat(str(raw_date).replace("Z", "+00:00"))
                age_days = (now - dt.astimezone(timezone.utc)).days
                if age_days > self.max_alphafold_age_days:
                    stale.append((uniprot_id, age_days))
            except ValueError:
                continue

        metrics["alphafold_rows_checked"] = float(len(rows))
        if stale:
            sample = [f"{uid} ({age} days)" for uid, age in stale[:5]]
            issues.append(
                ValidationIssue(
                    level="warning",
                    check="alphafold_freshness",
                    message=f"Outdated AlphaFold entries detected: {sample}",
                    suggested_fix="Refresh stale structures from the latest AlphaFold release.",
                )
            )

    def _table_exists(self, conn, name: str) -> bool:
        row = conn.execute(
            "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?", (name,)
        ).fetchone()
        return bool(row and row[0])
