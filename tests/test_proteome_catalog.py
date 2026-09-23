"""Tests for the Phase 2 proteome catalogue and its QC checklist."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.database.proteome_catalog import (
    build_catalog,
    ca_geometry_fraction,
    measure_model,
    qc_report,
    shard,
)


def _model_text(n_residues: int, *, spacing: float = 3.8, plddt: float = 90.0, end: bool = True) -> str:
    lines = []
    for i in range(n_residues):
        # A zig-zag keeps consecutive C-alpha atoms at ``spacing`` apart.
        x = i * spacing * np.cos(np.radians(20))
        y = (i % 2) * spacing * np.sin(np.radians(20))
        lines.append(
            f"ATOM  {i + 1:5d}  CA  ALA A{i + 1:4d}    {x:8.3f}{y:8.3f}{0.0:8.3f}"
            f"  1.00{plddt:6.2f}           C"
        )
    if end:
        lines.append("END")
    return "\n".join(lines) + "\n"


def _write(directory: Path, name: str, text: str, *, gz: bool = True) -> Path:
    path = directory / (name + (".gz" if gz else ""))
    if gz:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


class TestMeasureModel:
    def test_good_model(self, tmp_path):
        path = _write(tmp_path, "AF-P12345-F1-model_v4.pdb", _model_text(400, plddt=85.0))
        record = measure_model(path)
        assert record.uniprot_id == "P12345"
        assert record.fragment == 1 and record.version == 4
        assert record.length == 400
        assert record.mean_plddt == pytest.approx(85.0)
        assert record.fraction_plddt_70 == pytest.approx(1.0)
        assert record.ca_geometry_fraction == pytest.approx(1.0)
        assert record.qc_pass

    def test_truncated_model_fails(self, tmp_path):
        path = _write(tmp_path, "AF-P12345-F1-model_v4.pdb", _model_text(60, end=False))
        record = measure_model(path)
        assert not record.complete
        assert not record.qc_pass
        assert "truncated" in record.error

    def test_cut_gzip_stream_is_recorded_not_raised(self, tmp_path):
        good = _write(tmp_path, "AF-P1-F1-model_v4.pdb", _model_text(200))
        cut = tmp_path / "AF-P2-F1-model_v4.pdb.gz"
        cut.write_bytes(good.read_bytes()[: len(good.read_bytes()) // 2])
        record = measure_model(cut)
        assert not record.qc_pass
        assert record.error

    def test_zero_byte_file(self, tmp_path):
        path = tmp_path / "AF-P12345-F1-model_v4.pdb.gz"
        path.write_bytes(b"")
        record = measure_model(path)
        assert record.error == "zero-byte file"
        assert not record.qc_pass

    def test_nonconforming_name_fails(self, tmp_path):
        path = _write(tmp_path, "P12345_model.pdb", _model_text(60), gz=False)
        record = measure_model(path)
        assert not record.name_ok
        assert not record.qc_pass

    def test_wrong_geometry_fails(self, tmp_path):
        """Coordinates in nm rather than A put neighbours 0.38 apart."""
        path = _write(tmp_path, "AF-P12345-F1-model_v4.pdb", _model_text(60, spacing=0.38))
        record = measure_model(path)
        assert record.ca_geometry_fraction == pytest.approx(0.0)
        assert not record.qc_pass


def test_geometry_ignores_chain_breaks():
    coords = np.array([[0, 0, 0], [3.8, 0, 0], [50, 0, 0], [53.8, 0, 0]], dtype=float)
    resseqs = np.array([1, 2, 10, 11])
    chains = np.array(["A", "A", "A", "A"])
    assert ca_geometry_fraction(coords, resseqs, chains) == pytest.approx(1.0)


class TestCatalogAndQc:
    @pytest.fixture
    def catalog(self, tmp_path):
        _write(tmp_path, "AF-P00001-F1-model_v4.pdb", _model_text(300))
        _write(tmp_path, "AF-P00002-F1-model_v4.pdb", _model_text(250, plddt=50.0))
        _write(tmp_path, "AF-P00003-F1-model_v4.pdb", _model_text(2700))
        _write(tmp_path, "AF-P00003-F2-model_v4.pdb", _model_text(2700))
        _write(tmp_path, "AF-P00004-F1-model_v4.pdb", _model_text(80, end=False))
        return build_catalog(sorted(tmp_path.glob("*.gz")), "yeast")

    def test_master_csv_columns(self, catalog):
        for column in ("uniprot_id", "filename", "organism", "length", "mean_plddt", "qc_pass"):
            assert column in catalog.columns
        assert set(catalog["organism"]) == {"Saccharomyces cerevisiae"}

    def test_report(self, catalog):
        report = qc_report(catalog, "yeast")
        items = report["checklist"]
        assert items["file_count"]["observed_proteins"] == 4
        assert items["file_count"]["multi_fragment_proteins"] == 1
        assert not items["file_count"]["passes"]  # 4 is not ~6,049
        assert items["no_zero_byte_or_truncated"]["truncated"] == 1
        assert not items["no_zero_byte_or_truncated"]["passes"]
        assert items["consistent_naming"]["passes"]
        # P00001, P00002, P00003-F1 are screenable; the truncated one is not.
        assert report["screenable_f1_models"] == 3

    def test_shards_partition_the_screenable_f1_models(self, catalog):
        parts = [shard(catalog, i, 2) for i in range(2)]
        ids = pd.concat(parts)["uniprot_id"].tolist()
        assert sorted(ids) == ["P00001", "P00002", "P00003"]
        assert len(ids) == len(set(ids))
        # Longest first, dealt round-robin: the two longest land on different shards.
        assert parts[0]["uniprot_id"].iloc[0] == "P00003"
        assert parts[1]["uniprot_id"].iloc[0] == "P00001"

    def test_bad_shard_index(self, catalog):
        with pytest.raises(ValueError):
            shard(catalog, 2, 2)


def test_local_distortion_is_flagged_not_excluded(tmp_path):
    """Stretched steps in a low-confidence region flag a model; they do not drop it."""
    text = _model_text(400)
    lines = text.splitlines()
    # Stretch 20 % of the chain: shift every residue from 320 on by 1.5 A in z,
    # making the 319->320 step long, and give 60 steps extra length by spacing.
    stretched = []
    for i, line in enumerate(lines[:-1]):
        if i >= 320:
            z = float(line[46:54]) + 0.0
            x = float(line[30:38]) + (i - 319) * 1.2
            line = f"{line[:30]}{x:8.3f}{line[38:46]}{z:8.3f}{line[54:]}"
        stretched.append(line)
    stretched.append("END")
    path = _write(tmp_path, "AF-P77777-F1-model_v6.pdb", "\n".join(stretched) + "\n")
    record = measure_model(path)
    assert 0.5 <= record.ca_geometry_fraction < 0.95
    assert record.ca_pairs_long > 0
    assert record.qc_pass
    report = qc_report(build_catalog([path], "yeast"), "yeast")
    geometry = report["checklist"]["geometry_is_protein"]
    assert geometry["flagged_below_95pct"] == 1
    assert geometry["failing_models"] == 0
