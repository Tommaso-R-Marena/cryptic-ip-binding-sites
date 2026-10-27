"""The ground-truth arm's data step has to survive not finishing.

``within_protein.py sites`` runs fpocket over every deposited structure in the
census, which took four hours and fifteen minutes on a hosted runner. That is
long enough that a job will sometimes be killed at its budget, and the first
version of this command wrote nothing until the final copy was done - so a run
that did not fit discarded every hour it had spent and the next attempt began
again from zero.

Two properties fix that, and each has tests here:

* the tables are on disk after every copy, not after the last one;
* ``--resume`` skips what an earlier attempt's report already records, counting
  the copies it failed on as well as the ones it measured.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "within_protein_resume", ROOT / "scripts" / "within_protein.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["within_protein_resume"] = script
SPEC.loader.exec_module(script)


CENSUS = pd.DataFrame([
    {"comp_id": "IHP", "status": "eligible", "complete": True, "selected": True,
     "symmetry_contact": False, "pdb_id": "1aaa", "chain": "A", "resseq": 301,
     "homology_group_strict": "g1", "burial_class": "buried"},
    {"comp_id": "IHP", "status": "eligible", "complete": True, "selected": True,
     "symmetry_contact": False, "pdb_id": "2bbb", "chain": "B", "resseq": 302,
     "homology_group_strict": "g2", "burial_class": "cryptic"},
    {"comp_id": "IHP", "status": "eligible", "complete": True, "selected": True,
     "symmetry_contact": False, "pdb_id": "3ccc", "chain": "C", "resseq": 303,
     "homology_group_strict": "g3", "burial_class": "exposed"},
])


def _args(tmp_path: Path, *, resume: bool = False) -> argparse.Namespace:
    census = tmp_path / "census.csv"
    CENSUS.to_csv(census, index=False)
    structures = tmp_path / "structures"
    structures.mkdir(exist_ok=True)
    for pdb_id in CENSUS["pdb_id"]:
        (structures / f"{pdb_id}.pdb").write_text("ATOM\n")
    return argparse.Namespace(
        census=census, structures=structures, comp_id="IHP",
        out=tmp_path / "site_pockets.csv", sites=tmp_path / "deposited_sites.csv",
        report=tmp_path / "report.json", summary=tmp_path / "summary.json",
        fpocket_timeout=900.0, limit=0, resume=resume)


def _record(n_pockets: int = 2, error: str = "") -> dict:
    return {"site": [10, 11, 12], "error": error, "chain": "A",
            "n_chains_in_contact": 1,
            "pockets": [{"pocket_id": i, "residues": f"{i},{i + 1}"}
                        for i in range(1, n_pockets + 1)]}


def _patch(monkeypatch, fake):
    monkeypatch.setattr(script, "site_one", fake)


# ----------------------------------------------- written after every copy

def test_the_tables_are_on_disk_before_the_last_copy_is_measured(tmp_path, monkeypatch):
    """A job killed mid-run must leave what it already measured behind.

    The detector itself reads the table back on its second call, which is the
    only way to observe a write that happens between two copies.
    """
    args = _args(tmp_path)
    seen = {}

    def fake(path, chain, resseq, icode, work, comp_id="IHP", timeout_s=900.0):
        if resseq == 302:
            seen["rows"] = len(pd.read_csv(args.out))
            seen["sites"] = len(pd.read_csv(args.sites))
            seen["report"] = json.loads(args.report.read_text())
        return _record()

    _patch(monkeypatch, fake)
    script.cmd_sites(args)
    assert seen["rows"] == 2, "the first copy's pockets were not written before the second ran"
    assert seen["sites"] == 1
    assert [e["copy"] for e in seen["report"]["status"]] == ["1aaa_A_301"]


def test_a_failed_copy_is_recorded_on_disk_too(tmp_path, monkeypatch):
    """Otherwise a structure that cannot be read is retried on every attempt."""
    args = _args(tmp_path)

    def fake(path, chain, resseq, icode, work, comp_id="IHP", timeout_s=900.0):
        return _record(error="ligand copy not present in the structure") if resseq == 301 \
            else _record()

    _patch(monkeypatch, fake)
    script.cmd_sites(args)
    report = json.loads(args.report.read_text())
    assert {e["copy"] for e in report["status"]} == {"1aaa_A_301", "2bbb_B_302", "3ccc_C_303"}
    assert report["excluded"] == {"ligand copy not present in the structure": 1}
    assert report["n_measured"] == 2


def test_a_missing_structure_is_recorded_on_disk_too(tmp_path, monkeypatch):
    args = _args(tmp_path)
    (args.structures / "2bbb.pdb").unlink()
    _patch(monkeypatch, lambda *a, **k: _record())
    script.cmd_sites(args)
    report = json.loads(args.report.read_text())
    assert report["excluded"] == {"structure not downloaded": 1}
    assert {e["copy"] for e in report["status"]} == {"1aaa_A_301", "2bbb_B_302", "3ccc_C_303"}


# ------------------------------------------------------------- resuming

def test_resume_skips_the_copies_an_earlier_attempt_measured(tmp_path, monkeypatch):
    first = _args(tmp_path)

    def only_two(path, chain, resseq, icode, work, comp_id="IHP", timeout_s=900.0):
        if resseq == 303:
            raise KeyboardInterrupt("the job hit its budget")
        return _record()

    _patch(monkeypatch, only_two)
    with pytest.raises(KeyboardInterrupt):
        script.cmd_sites(first)

    measured = []

    def record_calls(path, chain, resseq, icode, work, comp_id="IHP", timeout_s=900.0):
        measured.append(resseq)
        return _record()

    second = _args(tmp_path, resume=True)
    _patch(monkeypatch, record_calls)
    script.cmd_sites(second)

    assert measured == [303], "resume re-ran copies the first attempt had already done"
    report = json.loads(second.report.read_text())
    assert report["n_measured"] == 3
    assert report["n_pockets"] == 6
    assert len(pd.read_csv(second.sites)) == 3


def test_resume_does_not_retry_a_copy_that_failed(tmp_path, monkeypatch):
    first = _args(tmp_path)

    def fail_then_stop(path, chain, resseq, icode, work, comp_id="IHP", timeout_s=900.0):
        if resseq == 303:
            raise KeyboardInterrupt("the job hit its budget")
        return _record(error="no contact with the ligand's own chain") if resseq == 301 \
            else _record()

    _patch(monkeypatch, fail_then_stop)
    with pytest.raises(KeyboardInterrupt):
        script.cmd_sites(first)

    measured = []
    second = _args(tmp_path, resume=True)
    _patch(monkeypatch, lambda p, c, r, i, w, comp_id="IHP", timeout_s=900.0:
           (measured.append(r), _record())[1])
    script.cmd_sites(second)

    assert measured == [303]
    report = json.loads(second.report.read_text())
    assert report["excluded"] == {"no contact with the ligand's own chain": 1}
    assert report["n_measured"] == 2


def test_without_resume_everything_is_measured_again(tmp_path, monkeypatch):
    """The flag is opt-in: the default stays a clean run from nothing."""
    _patch(monkeypatch, lambda *a, **k: _record())
    script.cmd_sites(_args(tmp_path))

    measured = []
    again = _args(tmp_path, resume=False)
    _patch(monkeypatch, lambda p, c, r, i, w, comp_id="IHP", timeout_s=900.0:
           (measured.append(r), _record())[1])
    script.cmd_sites(again)
    assert measured == [301, 302, 303]
    assert len(pd.read_csv(again.sites)) == 3


def test_resume_with_nothing_to_resume_from_is_an_ordinary_run(tmp_path, monkeypatch):
    measured = []
    args = _args(tmp_path, resume=True)
    _patch(monkeypatch, lambda p, c, r, i, w, comp_id="IHP", timeout_s=900.0:
           (measured.append(r), _record())[1])
    script.cmd_sites(args)
    assert measured == [301, 302, 303]


def test_a_resumed_table_keeps_the_rows_it_read_back(tmp_path, monkeypatch):
    """The pockets of the first attempt must appear in the final table unchanged."""
    first = _args(tmp_path)

    def one_then_stop(path, chain, resseq, icode, work, comp_id="IHP", timeout_s=900.0):
        if resseq != 301:
            raise KeyboardInterrupt("the job hit its budget")
        return _record(n_pockets=3)

    _patch(monkeypatch, one_then_stop)
    with pytest.raises(KeyboardInterrupt):
        script.cmd_sites(first)
    before = pd.read_csv(first.out)

    second = _args(tmp_path, resume=True)
    _patch(monkeypatch, lambda *a, **k: _record(n_pockets=1))
    script.cmd_sites(second)

    after = pd.read_csv(second.out)
    kept = after[after["uniprot_id"] == "1aaa_A_301"].reset_index(drop=True)
    pd.testing.assert_frame_equal(before.reset_index(drop=True), kept)
    assert len(after) == 3 + 1 + 1


def test_the_summary_is_written_at_the_end_of_a_resumed_run(tmp_path, monkeypatch):
    first = _args(tmp_path)
    _patch(monkeypatch, lambda p, c, r, i, w, comp_id="IHP", timeout_s=900.0:
           _record() if r == 301 else (_ for _ in ()).throw(KeyboardInterrupt("budget")))
    with pytest.raises(KeyboardInterrupt):
        script.cmd_sites(first)

    second = _args(tmp_path, resume=True)
    _patch(monkeypatch, lambda *a, **k: _record())
    assert script.cmd_sites(second) == 0
    summary = json.loads(second.summary.read_text())
    assert summary["n_measured"] == 3
    assert "status" not in summary, "the summary is the report without its 137 entries"
