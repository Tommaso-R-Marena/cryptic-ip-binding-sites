"""Tests for pH-dependent electrostatics workflows."""

from pathlib import Path

import pandas as pd
import pytest

from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator


def test_parse_propka_output_extracts_basic_residues(tmp_path):
    pka_text = """
 ARG   42 A    12.20
 LYS   51 A    10.45
 ASP   99 A     4.00
 HIS  120 B     6.75
"""
    pka_path = tmp_path / "sample.pka"
    pka_path.write_text(pka_text, encoding="utf-8")

    parsed = ElectrostaticsCalculator.parse_propka_output(pka_path)

    assert list(parsed["residue_name"]) == ["ARG", "LYS", "HIS"]
    assert list(parsed["residue_number"]) == [42, 51, 120]


def test_generate_pqr_raises_helpful_error(monkeypatch, tmp_path):
    calculator = ElectrostaticsCalculator(pdb2pqr_path="missing-pdb2pqr")
    pdb_path = tmp_path / "protein.pdb"
    pdb_path.write_text("ATOM", encoding="utf-8")

    def _fake_run(*args, **kwargs):
        raise FileNotFoundError("missing")

    monkeypatch.setattr("subprocess.run", _fake_run)

    with pytest.raises(RuntimeError, match="pdb2pqr executable not found"):
        calculator.generate_pqr(pdb_path=pdb_path, ph=7.4, output_dir=tmp_path)


def test_analyze_ph_dependent_binding(monkeypatch, tmp_path):
    calculator = ElectrostaticsCalculator()
    pdb_path = tmp_path / "protein.pdb"
    pdb_path.write_text("ATOM", encoding="utf-8")

    propka_df = pd.DataFrame(
        [
            {"residue_name": "ARG", "residue_number": 10, "chain_id": "A", "pka": 7.8},
            {"residue_name": "LYS", "residue_number": 20, "chain_id": "A", "pka": 10.5},
            {"residue_name": "HIS", "residue_number": 30, "chain_id": "A", "pka": 6.9},
        ]
    )

    monkeypatch.setattr(calculator, "run_propka", lambda *args, **kwargs: propka_df)
    monkeypatch.setattr(
        calculator,
        "generate_pqr",
        lambda pdb_path, ph, output_dir: Path(output_dir) / f"fake_{ph:.1f}.pqr",
    )
    monkeypatch.setattr(
        calculator,
        "run_apbs",
        lambda pqr_path, output_dir: 2.0 + float(pqr_path.stem.split("_")[-1]),
    )

    result = calculator.analyze_ph_dependent_binding(
        pdb_path=pdb_path,
        candidate_sites={"cryptic_1": [10, 30], "surface_1": [20]},
        output_dir=tmp_path,
        site_types={"cryptic_1": "cryptic", "surface_1": "surface"},
    )

    assert set(result.optimal_binding_ph.keys()) == {"cryptic_1", "surface_1"}
    assert not result.ph_sensitive_residues.empty
    assert result.profile_comparison is not None
    assert result.plot_path.exists()
    assert set(result.profile_comparison["site_type"]) == {"cryptic", "surface"}
