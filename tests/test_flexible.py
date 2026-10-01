"""scripts/flexible.py (docs/FLEXIBLE_PLAN.md + amendment 1): study J."""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def fx():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("flexible", ROOT / "scripts" / "flexible.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["flexible"] = module
    spec.loader.exec_module(module)
    return module


def _arrays(atoms):
    """atoms: list of (resname, chain, resseq, atom_name, element, xyz)."""
    keys, order = [], []
    for resname, chain, resseq, *_ in atoms:
        key = f"{chain}:{resseq}"
        keys.append(key)
        if key not in order:
            order.append(key)
    return types.SimpleNamespace(
        resnames=np.array([a[0] for a in atoms]),
        chain_ids=np.array([a[1] for a in atoms]),
        resseqs=np.array([a[2] for a in atoms]),
        icodes=np.array(["" for _ in atoms]),
        atom_names=np.array([a[3] for a in atoms]),
        elements=np.array([a[4] for a in atoms]),
        coords=np.array([a[5] for a in atoms], dtype=float),
        is_polymer=np.full(len(atoms), True),
        residue_index=np.array([order.index(k) for k in keys]),
    )


LIGAND = np.array([[0.0, 0.0, 0.0]])


# ------------------------------------------------------- residue selection
def test_only_residues_within_the_radius_are_flexible(fx):
    arrays = _arrays([
        ("LYS", "A", 1, "NZ", "N", [3.0, 0, 0]),      # inside 4 A
        ("ARG", "A", 2, "NH1", "N", [9.0, 0, 0]),     # outside
    ])
    picked = fx.flexible_residues(arrays, LIGAND)
    assert [r["resseq"] for r in picked] == [1]


def test_backbone_proximity_alone_does_not_make_a_residue_flexible(fx):
    """The plan says a *side-chain* heavy atom must be close."""
    arrays = _arrays([("LEU", "A", 1, "CA", "C", [2.0, 0, 0]),
                      ("LEU", "A", 1, "N", "N", [2.2, 0, 0]),
                      ("LEU", "A", 1, "CD1", "C", [9.0, 0, 0])])
    assert fx.flexible_residues(arrays, LIGAND) == []


def test_residues_without_a_rotatable_side_chain_are_excluded(fx):
    arrays = _arrays([("GLY", "A", 1, "CA", "C", [1.0, 0, 0]),
                      ("ALA", "A", 2, "CB", "C", [1.5, 0, 0]),
                      ("PRO", "A", 3, "CG", "C", [2.0, 0, 0]),
                      ("SER", "A", 4, "OG", "O", [2.5, 0, 0])])
    assert [r["resname"] for r in fx.flexible_residues(arrays, LIGAND)] == ["SER"]


def test_a_disulphide_cysteine_stays_fixed_but_a_free_one_does_not(fx):
    bonded = _arrays([("CYS", "A", 1, "SG", "S", [2.0, 0, 0]),
                      ("CYS", "A", 2, "SG", "S", [2.0, 2.0, 0])])   # 2.0 A apart
    assert fx.flexible_residues(bonded, LIGAND) == []
    free = _arrays([("CYS", "A", 1, "SG", "S", [2.0, 0, 0]),
                    ("CYS", "A", 2, "SG", "S", [2.0, 8.0, 0])])     # far apart
    assert [r["resseq"] for r in fx.flexible_residues(free, LIGAND)] == [1]


def test_selection_is_ordered_by_minimum_distance_to_the_ligand(fx):
    arrays = _arrays([("LYS", "A", 1, "NZ", "N", [3.5, 0, 0]),
                      ("ARG", "A", 2, "NH1", "N", [1.0, 0, 0]),
                      ("SER", "A", 3, "OG", "O", [2.0, 0, 0])])
    assert [r["resseq"] for r in fx.flexible_residues(arrays, LIGAND)] == [2, 3, 1]


def test_the_cap_is_honoured(fx):
    atoms = [("LYS", "A", i, "NZ", "N", [1.0 + 0.1 * i, 0, 0]) for i in range(1, 15)]
    picked = fx.flexible_residues(_arrays(atoms), LIGAND)
    assert len(picked) == fx.MAX_FLEX == 8


def test_a_residue_contributes_its_closest_side_chain_atom_only_once(fx):
    arrays = _arrays([("LYS", "A", 1, "CB", "C", [3.9, 0, 0]),
                      ("LYS", "A", 1, "NZ", "N", [1.1, 0, 0])])
    picked = fx.flexible_residues(arrays, LIGAND)
    assert len(picked) == 1
    assert picked[0]["min_distance"] == pytest.approx(1.1)


def test_meeko_residue_id_matches_the_monomer_key_format(fx):
    assert fx.meeko_residue_id({"chain": "A", "resseq": 42, "icode": ""}) == "A:42"
    assert fx.meeko_residue_id({"chain": "B", "resseq": 7, "icode": "A"}) == "B:7A"


# ------------------------------------------------------- the pose-parsing trap
def test_flex_residue_blocks_are_stripped_from_poses(fx):
    """Vina puts the moved side chains inside each pose; they are not ligand atoms."""
    pdbqt = "\n".join([
        "MODEL 1",
        "ATOM      1  P   IHP X 999       1.000   2.000   3.000  1.00  0.00    +0.000 P",
        "BEGIN_RES LYS A 2",
        "ROOT",
        "ATOM      2  CA  LYS A   2       9.000   9.000   9.000  1.00  0.00    +0.171 C",
        "ENDROOT",
        "END_RES LYS A 2",
        "ENDMDL",
    ])
    out = fx.strip_flex_residues(pdbqt)
    assert "LYS" not in out
    assert "BEGIN_RES" not in out and "END_RES" not in out and "ROOT" not in out
    assert "IHP" in out and "MODEL 1" in out and "ENDMDL" in out


def test_stripping_leaves_a_rigid_pose_file_untouched(fx):
    pdbqt = "MODEL 1\nATOM      1  P   IHP X 999       1.0   2.0   3.0  1.00  0.00    +0.0 P\nENDMDL"
    assert fx.strip_flex_residues(pdbqt).strip() == pdbqt.strip()


def test_stripping_handles_several_flexible_residues_in_one_model(fx):
    body = ["MODEL 1", "ATOM      1  P   IHP X 999       1.0   2.0   3.0  1.00  0.00    +0.0 P"]
    for res in ("LYS A 2", "ARG A 5"):
        body += [f"BEGIN_RES {res}", "ROOT",
                 "ATOM      9  CA  XXX A   9       9.0   9.0   9.0  1.00  0.00    +0.0 C",
                 "ENDROOT", f"END_RES {res}"]
    body.append("ENDMDL")
    out = fx.strip_flex_residues("\n".join(body))
    assert out.count("ATOM") == 1 and "XXX" not in out


# ------------------------------------------------------- outcomes
def test_success_uses_the_top_pose_and_the_ceiling_uses_the_best(fx):
    assert fx.success([1.2, 5.0]) is True
    assert fx.success([3.0, 0.5]) is False          # top pose misses, even though one hits
    assert fx.best_of_list([3.0, 0.5]) is True
    assert fx.success([]) is None and fx.best_of_list([]) is None
    assert fx.success([float("nan")]) is None


def _records(flex_rate, rigid_rate, n=12, ceiling=None):
    """n copies, each in its own strict group, with per-arm success rates realised exactly."""
    out = []
    for i in range(n):
        runs = {}
        for arm, rate in (("flex", flex_rate), ("rigid", rigid_rate)):
            hit = i < round(rate * n)
            top = 1.0 if hit else 5.0
            best = top if ceiling is None else (1.0 if i < round(ceiling * n) else 5.0)
            runs[arm] = [{"seed": s, "rmsd": [top, best], "vina": [-8.0, -7.0]} for s in (1, 2, 3)]
        out.append({"copy_key": f"X{i}", "homology_group_strict": f"G{i}",
                    "burial_class": "cryptic" if i % 2 else "surface", "n_flex": 4, "runs": runs})
    return out


def test_per_copy_collapses_seeds_into_one_row_per_arm(fx):
    frame = fx.per_copy(_records(0.5, 0.25, n=8))
    assert set(frame["arm"]) == {"flex", "rigid"}
    assert len(frame) == 16
    assert frame[frame.arm == "flex"]["success"].mean() == pytest.approx(0.5)


def test_j1_reports_better_when_flex_wins(fx):
    r = fx.build(_records(1.0, 0.0, n=12), n_bootstrap=300)
    assert r["J1"]["decision"] == "better"
    assert r["J1"]["estimate"]["per_group"]["point"] == pytest.approx(1.0)
    assert set(r["holm"]) == {"J1", "J2"}


def test_j1_reports_no_material_change_within_the_margin(fx):
    r = fx.build(_records(0.0, 0.0, n=12), n_bootstrap=300)
    assert r["J1"]["decision"] == "no material change"


def test_j1_is_not_evaluable_below_the_group_floor(fx):
    r = fx.build(_records(1.0, 0.0, n=3), n_bootstrap=300)
    assert r["J1"]["decision"] == "not evaluable"


def test_the_ceiling_is_reported_separately_from_the_top_pose(fx):
    """A receptor that lifts the ceiling but not the top pose is a scoring result."""
    records = _records(0.0, 0.0, n=12, ceiling=1.0)
    r = fx.build(records, n_bootstrap=300)
    assert r["J1"]["decision"] == "no material change"
    assert r["J2"]["flex_mean"] == pytest.approx(1.0)


def test_the_rigid_arm_is_audited_against_study_f(fx):
    """Amendment 1: the receptor pipeline change must be visible, not folded into J1."""
    r = fx.build(_records(1.0, 0.0, n=12), n_bootstrap=300)
    check = r["rigid_arm_check"]
    assert check["study_f_success"] == 0.110
    assert check["meeko_receptor_success"] == pytest.approx(0.0)
    assert "receptor pipeline, audited" in fx.markdown(r)


def test_j3_stratifies_by_burial_and_is_not_holm_corrected(fx):
    r = fx.build(_records(1.0, 0.0, n=12), n_bootstrap=300)
    assert set(r["J3"]) == {"cryptic", "surface"}
    assert set(r["holm"]) == {"J1", "J2"}          # J3 is exploratory


def test_errors_are_reported_rather_than_dropped_silently(fx):
    records = _records(1.0, 0.0, n=12) + [{"copy_key": "BAD", "error": "ReceptorError: no polymer atoms"}]
    r = fx.build(records, n_bootstrap=300)
    assert r["errors"] == ["ReceptorError: no polymer atoms"]
    assert r["n_copies_scored"] == 12
    assert "Copies not docked" in fx.markdown(r)


def test_markdown_states_the_flexible_residue_budget(fx):
    text = fx.markdown(fx.build(_records(1.0, 0.0, n=12), n_bootstrap=300))
    assert f"within {fx.FLEX_RADIUS}" in text and f"at most {fx.MAX_FLEX}" in text


# ------------------------------------------------------- real receptor preparation
def _external_tools_present() -> bool:
    import shutil

    if shutil.which("pdb2pqr") is None and shutil.which("pdb2pqr30") is None:
        return False
    try:
        import meeko  # noqa: F401
    except ImportError:
        return False
    return True


requires_tools = pytest.mark.skipif(not _external_tools_present(),
                                    reason="needs pdb2pqr and meeko")


@pytest.fixture(scope="module")
def peptide(tmp_path_factory):
    """A real protonated receptor pair, through the project's own PDB2PQR path."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from cryptic_ip.docking.receptor import protonate

    mol = Chem.AddHs(Chem.MolFromSequence("AKA"), addCoords=True)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(mol)
    block = Chem.MolToPDBBlock(mol)
    work = tmp_path_factory.mktemp("peptide")
    src = work / "pep.pdb"
    # RDKit labels the N-terminal fragment UNL; PDB2PQR needs standard residues only.
    src.write_text("".join(line + "\n" for line in block.splitlines()
                           if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() != "UNL") + "END\n")
    return protonate(src, work), work


@requires_tools
def test_no_flexible_residues_yields_no_flex_file(fx, peptide):
    """The plan's pinning: with nothing flexible, "rigid" must mean rigid."""
    protonated, work = peptide
    rigid, flex, applied = fx.prepare_pair(protonated, work / "rigid_only", ())
    assert flex is None
    assert applied == []
    assert rigid.exists() and rigid.read_text().strip()


@requires_tools
def test_a_flexible_lysine_moves_out_of_the_rigid_file_into_a_torsion_tree(fx, peptide):
    protonated, work = peptide
    _pdb, pqr = protonated
    lysines = sorted({(line[21:22].strip() or "A", int(line[22:26]))
                      for line in Path(pqr).read_text().splitlines()
                      if line.startswith("ATOM") and line[17:20].strip() == "LYS"})
    assert lysines, "the fixture peptide should contain a lysine"
    record = {"chain": lysines[0][0], "resseq": lysines[0][1], "icode": ""}

    rigid_only, _none, _ = fx.prepare_pair(protonated, work / "baseline", ())
    rigid, flex, applied = fx.prepare_pair(protonated, work / "flexed", [record])

    assert applied == [fx.meeko_residue_id(record)]
    assert flex is not None
    tree = flex.read_text()
    # A real AutoDock flexible residue: a named block with a rooted, branched torsion tree.
    for token in ("BEGIN_RES", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH"):
        assert token in tree, token
    # The side chain must leave the rigid file, or it would be present twice.
    assert len(rigid.read_text().splitlines()) < len(rigid_only.read_text().splitlines())


@requires_tools
def test_an_unknown_residue_id_is_skipped_rather_than_crashing_the_copy(fx, peptide):
    protonated, work = peptide
    _rigid, flex, applied = fx.prepare_pair(
        protonated, work / "bogus", [{"chain": "Z", "resseq": 9999, "icode": ""}])
    assert applied == [] and flex is None


# ------------------------------------------------------- the residues census
def test_nan_text_reads_as_empty_not_as_the_string_nan(fx):
    """`NaN or ""` is NaN, so formatting it gives "nan" and matches no copy.

    This cost a whole study J shard: every census row without an insertion code carries
    NaN there, so find_copy was asked for "A:281nan".
    """
    assert fx._text(float("nan")) == ""
    assert fx._text(None) == ""
    assert fx._text("  A ") == "A"


def test_the_residues_census_records_a_bad_copy_instead_of_aborting_the_shard(fx, tmp_path, monkeypatch):
    """find_copy raises rather than returning None, which killed the whole step."""
    import redocking

    census = tmp_path / "census.csv"
    pd_mod = pytest.importorskip("pandas")
    pd_mod.DataFrame([
        {"copy_key": "1AAA:A:1", "pdb_id": "1AAA", "chain": "A", "resseq": 1,
         "icode": float("nan"), "selected": True},
        {"copy_key": "1BBB:A:2", "pdb_id": "1BBB", "chain": "A", "resseq": 2,
         "icode": float("nan"), "selected": True},
    ]).to_csv(census, index=False)

    structures = tmp_path / "structures"
    structures.mkdir()
    (structures / "1AAA.pdb").write_text("END\n")
    (structures / "1BBB.pdb").write_text("END\n")

    seen = {}

    def fake_find_copy(arrays, chain, resseq, icode):
        seen[resseq] = icode
        if resseq == 1:
            raise LookupError(f"copy {chain}:{resseq}{icode} not found")
        return ("k", "IHP", np.array([0]))

    monkeypatch.setattr(redocking, "find_copy", fake_find_copy)
    monkeypatch.setattr(redocking, "structure_path", lambda d, p: Path(d) / f"{p}.pdb")
    monkeypatch.setattr(fx, "flexible_residues", lambda arrays, xyz: [{"resseq": 5}])
    monkeypatch.setattr(
        "cryptic_ip.analysis.structure_arrays.load_structure_arrays",
        lambda path: types.SimpleNamespace(elements=np.array(["P"]), coords=np.zeros((1, 3))))

    out = tmp_path / "residues.json"
    assert fx.main(["residues", "--census", str(census), "--structures", str(structures),
                    "--out", str(out)]) == 0
    import json
    records = json.loads(out.read_text())
    # The raising copy is recorded as its own failure; the good copy still gets counted.
    assert records[0]["error"].startswith("LookupError")
    assert records[1]["n_flex"] == 1
    # And the icode reached find_copy as empty, never as "nan".
    assert set(seen.values()) == {""}


def test_the_report_survives_a_run_whose_shards_all_failed(fx, tmp_path):
    """The report job runs on always() so a dead run still reports; it must not crash.

    pivot_table raises KeyError on an empty frame, so the emptiness guard has to come
    before the pivot rather than after it.
    """
    import json

    records = tmp_path / "records.json"
    records.write_text("[]")
    out, md = tmp_path / "flexible.json", tmp_path / "FLEXIBLE.md"
    assert fx.main(["report", "--records", str(records), "--out", str(out),
                    "--markdown", str(md)]) == 0
    result = json.loads(out.read_text())
    assert result["J1"]["decision"] == "not evaluable"
    assert result["J2"]["decision"] == "not evaluable"
    assert result["n_copies_scored"] == 0
    assert "holm" not in result
    assert "not evaluable" in md.read_text()


def test_the_report_survives_records_that_are_all_errors(fx, tmp_path):
    import json

    records = tmp_path / "records.json"
    records.write_text(json.dumps([{"copy_key": "X", "error": "ReceptorError: no polymer atoms"}]))
    out, md = tmp_path / "flexible.json", tmp_path / "FLEXIBLE.md"
    assert fx.main(["report", "--records", str(records), "--out", str(out),
                    "--markdown", str(md)]) == 0
    result = json.loads(out.read_text())
    assert result["J1"]["decision"] == "not evaluable"
    assert result["errors"] == ["ReceptorError: no polymer atoms"]


# ------------------------------------------------- the shard's own protections
#
# Run 36828409167 lost thirteen shards' docking work: one copy overran the per-copy
# timeout, subprocess.TimeoutExpired raised out of the loop, and the output file - written
# only after the loop - was never written at all. Each of these pins one of the three
# protections added in response.
def _shard_census(tmp_path, keys):
    import pandas as pd

    census = tmp_path / "census.csv"
    pd.DataFrame([{"copy_key": key, "pdb_id": key.split(":")[0], "chain": "A", "resseq": 1,
                   "icode": None, "selected": "True", "homology_group_strict": "g1",
                   "burial_class": "buried"} for key in keys]).to_csv(census, index=False)
    return census


def _dock_argv(tmp_path, census, out, extra=()):
    return ["dock", "--census", str(census), "--structures-dir", str(tmp_path),
            "--ccd-dir", str(tmp_path), "--work-dir", str(tmp_path / "work"),
            "--shard", "0", "--shards", "1", "--out", str(out), *extra]


def _fake_run(monkeypatch, behaviour):
    """behaviour(call_index, out_json) -> None to succeed, or raises."""
    import subprocess
    import types as _types

    calls = {"n": 0}

    def run(cmd, **kwargs):
        index = calls["n"]
        calls["n"] += 1
        out_json = Path(cmd[cmd.index("--out-json") + 1])
        behaviour(index, out_json, kwargs)
        return _types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", run)
    return calls


def test_a_copy_that_overruns_its_timeout_does_not_take_the_shard(fx, tmp_path, monkeypatch):
    import json
    import subprocess

    census = _shard_census(tmp_path, ["1AAA:A:1", "2BBB:A:1", "3CCC:A:1"])
    out = tmp_path / "flexible-0.json"

    def behaviour(index, out_json, kwargs):
        if index == 1:
            raise subprocess.TimeoutExpired(cmd="dock-one", timeout=5400)
        out_json.write_text(json.dumps({"copy_key": "k", "runs": {}}))

    calls = _fake_run(monkeypatch, behaviour)
    assert fx.main(_dock_argv(tmp_path, census, out)) == 0
    assert calls["n"] == 3, "the shard stopped at the overrunning copy"
    records = json.loads(out.read_text())
    assert len(records) == 3
    assert "timed out after 5400 s" == records[1]["error"]
    assert [bool(r.get("error")) for r in records] == [False, True, False]


def test_the_output_file_is_rewritten_after_every_copy(fx, tmp_path, monkeypatch):
    import json

    census = _shard_census(tmp_path, ["1AAA:A:1", "2BBB:A:1"])
    out = tmp_path / "flexible-0.json"
    seen = []

    def behaviour(index, out_json, kwargs):
        # What the file holds at the moment this copy starts, before it has written anything.
        seen.append(len(json.loads(out.read_text())))
        out_json.write_text(json.dumps({"copy_key": "k", "runs": {}}))

    _fake_run(monkeypatch, behaviour)
    assert fx.main(_dock_argv(tmp_path, census, out)) == 0
    assert seen == [0, 1], "a shard killed mid-copy would have uploaded nothing"
    assert len(json.loads(out.read_text())) == 2


def test_copies_the_budget_never_reached_say_so(fx, tmp_path, monkeypatch):
    import json

    census = _shard_census(tmp_path, ["1AAA:A:1", "2BBB:A:1", "3CCC:A:1"])
    out = tmp_path / "flexible-0.json"

    def behaviour(index, out_json, kwargs):
        out_json.write_text(json.dumps({"copy_key": "k", "runs": {}}))

    calls = _fake_run(monkeypatch, behaviour)
    # A budget under --min-copy-seconds: no copy may be started at all.
    assert fx.main(_dock_argv(tmp_path, census, out,
                              ("--shard-budget", "0", "--min-copy-seconds", "600"))) == 0
    assert calls["n"] == 0
    records = json.loads(out.read_text())
    assert len(records) == 3
    assert all(r["error"] == "not reached: the shard's budget ran out" for r in records)


def test_the_per_copy_timeout_never_exceeds_the_budget_left(fx, tmp_path, monkeypatch):
    import json

    census = _shard_census(tmp_path, ["1AAA:A:1"])
    out = tmp_path / "flexible-0.json"
    timeouts = []

    def behaviour(index, out_json, kwargs):
        timeouts.append(kwargs["timeout"])
        out_json.write_text(json.dumps({"copy_key": "k", "runs": {}}))

    _fake_run(monkeypatch, behaviour)
    assert fx.main(_dock_argv(tmp_path, census, out,
                              ("--shard-budget", "900", "--per-copy-timeout", "5400"))) == 0
    assert 0 < timeouts[0] <= 900


def test_error_classes_count_causes_and_keep_one_traceback_each(fx):
    records = [
        {"copy_key": "A", "error": "PolymerCreationError: template matching failed",
         "traceback": "Traceback\n  polymer.py line 1\nPolymerCreationError"},
        {"copy_key": "B", "error": "PolymerCreationError: H discrepancy"},
        {"copy_key": "C", "error": "ValueError: invalid literal for int()",
         "traceback": "Traceback\n  pdbutils.py line 9\nValueError"},
        {"copy_key": "D", "runs": {}},
    ]
    classes = fx.error_classes(records)
    assert classes["PolymerCreationError"]["n"] == 2
    assert classes["PolymerCreationError"]["copies"] == ["A", "B"]
    assert "polymer.py" in classes["PolymerCreationError"]["traceback"]
    assert classes["ValueError"]["n"] == 1
    assert "D" not in str(classes)


def test_the_report_tables_the_causes_a_run_failed_for(fx, tmp_path):
    import json

    records = tmp_path / "records.json"
    records.write_text(json.dumps([
        {"copy_key": "A", "error": "PolymerCreationError: template matching failed"},
        {"copy_key": "B", "error": "PolymerCreationError: H discrepancy"},
        {"copy_key": "C", "error": "ValueError: invalid literal for int()"}]))
    out, md = tmp_path / "flexible.json", tmp_path / "FLEXIBLE.md"
    assert fx.main(["report", "--records", str(records), "--out", str(out),
                    "--markdown", str(md)]) == 0
    text = md.read_text()
    assert "Copies not docked, by cause" in text
    assert "| PolymerCreationError | 2 |" in text
    assert "| ValueError | 1 |" in text
