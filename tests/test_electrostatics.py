"""Tests for pH-dependent electrostatics workflows."""

from pathlib import Path

import numpy as np
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


def _write_apbs_dx(path, origin, delta, counts, func):
    """Write an OpenDX map the way APBS does: z varies fastest, x slowest."""
    nx, ny, nz = counts
    values = [
        func(origin[0] + i * delta[0], origin[1] + j * delta[1], origin[2] + k * delta[2])
        for i in range(nx)
        for j in range(ny)
        for k in range(nz)
    ]
    lines = [
        f"object 1 class gridpositions counts {nx} {ny} {nz}",
        f"origin {origin[0]} {origin[1]} {origin[2]}",
        f"delta {delta[0]} 0 0",
        f"delta 0 {delta[1]} 0",
        f"delta 0 0 {delta[2]}",
        f"object 2 class gridconnections counts {nx} {ny} {nz}",
        f"object 3 class array type double rank 0 items {nx * ny * nz} data follows",
    ]
    for start in range(0, len(values), 3):
        lines.append(" ".join(f"{v:.6e}" for v in values[start : start + 3]))
    lines.append('attribute "dep" string "positions"')
    lines.append('object "regular positions regular connections" class field')
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize("counts", [(5, 6, 7), (6, 6, 6)])
def test_dx_sampling_reads_the_right_point(tmp_path, counts):
    """A linear field is reproduced exactly by trilinear interpolation.

    With x and z transposed the sample comes from the mirror point, which a
    field that weights the axes differently exposes - including on a cubic
    grid, where the transposed array still has the right shape.
    """
    from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator

    def field(x, y, z):
        return 1.0 * x + 10.0 * y + 100.0 * z

    dx = tmp_path / "pot.dx"
    _write_apbs_dx(dx, origin=(-2.0, 1.0, 3.0), delta=(0.5, 0.75, 1.0), counts=counts, func=field)
    point = (-0.9, 2.3, 5.4)
    sampled = ElectrostaticsCalculator.__new__(ElectrostaticsCalculator).sample_potential_at_point(dx, point)
    assert sampled == pytest.approx(field(*point), rel=1e-6)


def test_sampling_outside_the_map_is_refused(tmp_path):
    from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator

    dx = tmp_path / "pot.dx"
    _write_apbs_dx(dx, origin=(0.0, 0.0, 0.0), delta=(1.0, 1.0, 1.0), counts=(4, 4, 4), func=lambda x, y, z: x)
    calculator = ElectrostaticsCalculator.__new__(ElectrostaticsCalculator)
    with pytest.raises(ValueError, match="outside the potential map"):
        calculator.sample_potential_at_point(dx, (10.0, 1.0, 1.0))


def test_grid_contains_the_molecule_and_focuses_on_the_point(tmp_path):
    from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator

    pqr = tmp_path / "big.pqr"
    lines = []
    for i, x in enumerate(np.linspace(-60.0, 60.0, 50)):
        lines.append(f"ATOM  {i + 1:5d}  CA  ALA A{i + 1:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  0.0000 1.8500")
    pqr.write_text("\n".join(lines) + "\n")
    cglen, cgcent, fglen, fgcent = ElectrostaticsCalculator._grid_geometry(pqr, (50.0, 0.0, 0.0))
    coarse = [float(v) for v in cglen.split()]
    fine = [float(v) for v in fglen.split()]
    # The coarse box spans the 120 A molecule and holds a fine box centred at its end.
    assert coarse[0] >= 120.0 + fine[0]
    assert fgcent == "50.000 0.000 0.000"
    assert cgcent.startswith("0.000")


def test_apbs_without_energy_line_still_returns_the_map(tmp_path, monkeypatch):
    """The energy is a by-product; a written map is what matters."""
    import math
    import subprocess

    from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator

    pqr = tmp_path / "x.pqr"
    pqr.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  0.0000 1.8500\n")

    def fake_run(cmd, **kwargs):
        (tmp_path / "out" / "x.dx").write_text("map")
        return subprocess.CompletedProcess(cmd, 0, "no energy printed here", "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    energy, dx = ElectrostaticsCalculator().run_apbs_with_map(pqr, tmp_path / "out")
    assert math.isnan(energy)
    assert dx.name == "x.dx"


def test_apbs_with_neither_energy_nor_map_raises(tmp_path, monkeypatch):
    import subprocess

    from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator

    pqr = tmp_path / "x.pqr"
    pqr.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  0.0000 1.8500\n")
    monkeypatch.setattr(
        subprocess, "run", lambda cmd, **k: subprocess.CompletedProcess(cmd, 0, "error: grid", "")
    )
    with pytest.raises(RuntimeError, match="neither an energy nor a potential map"):
        ElectrostaticsCalculator().run_apbs(pqr, tmp_path / "out")
