"""scripts/redocking.py, scripts/redocking_report.py and scripts/extract_log_block.py on synthetic data."""

from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
MYO_IP6 = ("O[P](=O)(O)O[C@H]1[C@H](O[P](=O)(O)O)[C@@H](O[P](=O)(O)O)[C@H](O[P](=O)(O)O)"
           "[C@@H](O[P](=O)(O)O)[C@@H]1O[P](=O)(O)O")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# ------------------------------------------------------------ log blocks
def test_log_block_round_trip_through_a_timestamped_log(tmp_path):
    blocks = _load("extract_log_block")
    data = tmp_path / "r.json"
    data.write_text(json.dumps({"a": 1.5, "b": [1, 2]}))
    table = tmp_path / "t.csv"
    table.write_text("x,y\n" + "\n".join(f"{i},{i * i}" for i in range(500)) + "\n")
    log = "\n".join(
        f"2026-09-24T12:00:{i % 60:02d}.1234567Z {line}"
        for i, line in enumerate(["noise", blocks.emit("R_JSON", data), "more", blocks.emit("T_B64", table)])
    )
    wrapped = tmp_path / "saved.txt"
    wrapped.write_text(json.dumps({"job_id": 1, "logs_content": "﻿" + log}))
    text = blocks.log_text(wrapped.read_text())
    assert json.loads(blocks.extract(text, "R_JSON")) == {"a": 1.5, "b": [1, 2]}
    assert blocks.extract(text, "T_B64") == table.read_bytes()
    with pytest.raises(SystemExit):
        blocks.extract(text, "MISSING_JSON")


# ------------------------------------------------------------- selection
def test_selection_one_copy_per_entry_and_class_complete_first():
    redock = _load("redocking")
    rows = pd.DataFrame([
        {"pdb_id": "1AAA", "burial_class": "cryptic", "status": "eligible", "complete": False, "chain": "A",
         "resseq": 1, "icode": ""},
        {"pdb_id": "1AAA", "burial_class": "cryptic", "status": "eligible", "complete": True, "chain": "B",
         "resseq": 1, "icode": ""},
        {"pdb_id": "1AAA", "burial_class": "surface", "status": "eligible", "complete": True, "chain": "C",
         "resseq": 1, "icode": ""},
        {"pdb_id": "1AAA", "burial_class": "surface", "status": "eligible", "complete": True, "chain": "D",
         "resseq": 1, "icode": ""},
        {"pdb_id": "1AAA", "burial_class": "crystal_artifact", "status": "excluded: crystal artefact",
         "complete": True, "chain": "E", "resseq": 1, "icode": ""},
        {"pdb_id": "2BBB", "burial_class": "surface", "status": "excluded: configuration differs from the CCD",
         "complete": True, "chain": "A", "resseq": 1, "icode": ""},
    ])
    selected = redock.select_copies(rows)
    assert rows.loc[selected, "chain"].tolist() == ["B", "C"]


def test_shards_partition_entries():
    redock = _load("redocking")
    entries = [f"{i}ABC" for i in range(37)]
    shards = [redock.shard_entries(entries, i, 5) for i in range(5)]
    flat = [e for s in shards for e in s]
    assert sorted(flat) == sorted(entries) and len(flat) == len(set(flat))


def test_kabsch_recovers_a_rigid_motion():
    from scipy.spatial.transform import Rotation

    redock = _load("redocking")
    x = np.random.default_rng(1).normal(size=(10, 3))
    rot = Rotation.random(random_state=2).as_matrix()
    y = x @ rot.T + np.array([1.0, -2.0, 3.0])
    r, t = redock.kabsch(x, y)
    assert np.allclose(x @ r.T + t, y, atol=1e-8)


def test_sequence_alignment_identity_and_coverage():
    redock = _load("redocking")
    pairs, identity, coverage = redock.align_sequences("MKTAYIAKQR", "MKTAYIAKQRQISFVKSHFSRQ")
    assert identity == 1.0 and coverage == 1.0 and len(pairs) == 10


# ---------------------------------------------------------- end to end
def _synthetic_entry(out: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    pep = Chem.AddHs(Chem.MolFromSequence("KRKHGSKRWKHKRSGKRKYK"))
    params = AllChem.ETKDGv3()
    params.randomSeed, params.useRandomCoords = 7, True
    AllChem.EmbedMolecule(pep, params)
    AllChem.MMFFOptimizeMolecule(pep, maxIters=300)
    pep = Chem.RemoveHs(pep)
    lig = Chem.AddHs(Chem.MolFromSmiles(MYO_IP6))
    AllChem.EmbedMolecule(lig, randomSeed=11)
    lig = Chem.RemoveHs(lig)
    shift = pep.GetConformer().GetPositions().mean(axis=0) - lig.GetConformer().GetPositions().mean(axis=0)
    shift = shift + np.array([4.0, 0, 0])
    lines = [ln for ln in Chem.MolToPDBBlock(pep, flavor=4).splitlines() if ln.startswith("ATOM")]
    counts = {}
    for i, atom in enumerate(lig.GetAtoms()):
        el = atom.GetSymbol()
        counts[el] = counts.get(el, 0) + 1
        x, y, z = lig.GetConformer().GetPositions()[i] + shift
        lines.append(f"HETATM{9000 + i:5d} {el + str(counts[el]):<4s} IHP B 501    {x:8.3f}{y:8.3f}{z:8.3f}"
                     f"  1.00 20.00          {el:>2s}")
    (out / "SYN1.pdb").write_text("CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n"
                                  + "\n".join(lines) + "\nEND\n")
    (out / "ccd").mkdir()
    (out / "ccd" / "IHP.cif").write_text(
        "data_IHP\n_chem_comp.id IHP\nloop_\n_pdbx_chem_comp_descriptor.comp_id\n_pdbx_chem_comp_descriptor.type\n"
        "_pdbx_chem_comp_descriptor.program\n_pdbx_chem_comp_descriptor.program_version\n"
        f'_pdbx_chem_comp_descriptor.descriptor\nIHP SMILES_CANONICAL CACTVS 3.385 "{MYO_IP6}"\n')
    pd.DataFrame([{"pdb_id": "SYN1", "uniprot_ids": "P00000", "resolution": "2.00",
                   "experimental_method": "X-RAY DIFFRACTION", "release_date": "2020-01-01",
                   "homology_group": "G1", "homology_group_strict": "S1"}]).to_csv(out / "entries.csv", index=False)
    pd.DataFrame([{"structure_id": "SYN1", "pocket_id": 1}]).to_csv(out / "table.csv.gz", index=False)


@pytest.mark.skipif(shutil.which("pdb2pqr") is None and shutil.which("pdb2pqr30") is None,
                    reason="pdb2pqr not installed")
def test_census_dock_and_report_end_to_end(tmp_path, monkeypatch):
    pytest.importorskip("vina")
    pytest.importorskip("meeko")
    pytest.importorskip("gemmi")
    redock = _load("redocking")
    report = _load("redocking_report")
    _synthetic_entry(tmp_path)
    census = tmp_path / "census.csv"
    assert redock.main(["census", "--table", str(tmp_path / "table.csv.gz"), "--entries",
                        str(tmp_path / "entries.csv"), "--structures-dir", str(tmp_path), "--ccd-dir",
                        str(tmp_path / "ccd"), "--output", str(census)]) == 0
    frame = pd.read_csv(census)
    assert frame["status"].tolist() == ["eligible"] and bool(frame["selected"].iloc[0])
    assert bool(frame["configuration_match"].iloc[0]) and frame["species"].iloc[0] == "InsP6"
    assert redock.main(["dock", "--arm", "primary", "--census", str(census), "--structures-dir", str(tmp_path),
                        "--ccd-dir", str(tmp_path / "ccd"), "--out-dir", str(tmp_path / "arms"),
                        "--work-dir", str(tmp_path / "work"), "--exhaustiveness", "1"]) == 0
    record = json.loads((tmp_path / "arms" / "primary_0.jsonl").read_text().splitlines()[0])
    assert "error" not in record, record
    arms = {r["arm"] for r in record["runs"]}
    assert {"vina_s1", "vina_s2", "vina_s3", "crystal_control"} <= arms
    run = record["runs"][0]
    assert run["start_rmsd"] > 2.0 and run["n_poses"] >= 1 and run["box_side"] >= 22.0
    out = tmp_path / "rep"
    assert report.main(["--census", str(census), "--results-dir", str(tmp_path / "arms"), "--out-dir", str(out),
                        "--n-bootstrap", "20"]) == 0
    rep = json.loads((out / "redocking.json").read_text())
    assert rep["accounting"]["primary_set"] == 1
    assert rep["decisions"]["R1"]["decision"].startswith("not evaluable")  # one group only
    assert (out / "report.html").read_text().startswith("<!DOCTYPE html>")
