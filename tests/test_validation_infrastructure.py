"""Tests for validation infrastructure and CLI entry points."""

import sqlite3
from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from cryptic_ip.cli import main
from cryptic_ip.database.integrity_checker import DatabaseIntegrityChecker
from cryptic_ip.validation.results_validator import ResultsValidator
from cryptic_ip.validation.structure_validator import StructureValidator


def _write_minimal_pdb(path: Path) -> None:
    path.write_text(
        "HEADER    TEST STRUCTURE\n"
        "REMARK   2 RESOLUTION.    2.00 ANGSTROMS.\n"
        "ATOM      1  N   ALA A   1      11.104  13.207  14.503  1.00 82.00           N\n"
        "ATOM      2  CA  ALA A   1      11.551  11.823  14.240  1.00 82.00           C\n"
        "ATOM      3  C   ALA A   1      12.001  11.089  15.505  1.00 82.00           C\n"
        "ATOM      4  N   GLY A   2      13.104  10.207  16.503  1.00 80.00           N\n"
        "ATOM      5  CA  GLY A   2      13.551   8.823  16.240  1.00 80.00           C\n"
        "TER\nEND\n"
    )


def test_structure_validator_accepts_valid_pdb(tmp_path):
    pdb_path = tmp_path / "ok.pdb"
    _write_minimal_pdb(pdb_path)

    report = StructureValidator(min_chain_length=1).validate(str(pdb_path))

    assert report.valid
    assert not [i for i in report.issues if i.level == "error"]


def test_results_validator_flags_missing_required_columns(tmp_path):
    csv_path = tmp_path / "results.csv"
    pd.DataFrame({"pocket_id": [1], "volume": [200.0]}).to_csv(csv_path, index=False)

    report = ResultsValidator().validate(str(csv_path))

    assert not report.valid
    assert any(issue.check == "required_features" for issue in report.issues)


def test_database_integrity_checker_detects_missing_table(tmp_path):
    db_path = tmp_path / "cache.db"
    conn = sqlite3.connect(db_path)
    conn.execute("CREATE TABLE other_table (id INTEGER)")
    conn.commit()
    conn.close()

    report = DatabaseIntegrityChecker().validate(str(db_path))

    assert not report.valid
    assert any(issue.check == "required_tables" for issue in report.issues)


def test_cli_validate_structure_command(tmp_path):
    pdb_path = tmp_path / "ok.pdb"
    _write_minimal_pdb(pdb_path)

    runner = CliRunner()
    result = runner.invoke(main, ["validate", "--structure", str(pdb_path)])

    assert result.exit_code == 0
    assert "validation" in result.output.lower()


def test_cli_validate_database_command(tmp_path):
    db_path = tmp_path / "cache.db"
    conn = sqlite3.connect(db_path)
    conn.execute("CREATE TABLE structures (uniprot_id TEXT, model_date TEXT, file_path TEXT)")
    conn.execute(
        "INSERT INTO structures (uniprot_id, model_date, file_path) VALUES (?, ?, ?)",
        (None, "2000-01-01T00:00:00+00:00", "missing.pdb"),
    )
    conn.commit()
    conn.close()

    runner = CliRunner()
    result = runner.invoke(main, ["validate", "--database", str(db_path)])

    assert result.exit_code == 0
    assert "database validation" in result.output.lower()
