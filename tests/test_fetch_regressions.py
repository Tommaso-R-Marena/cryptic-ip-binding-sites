"""Regression guards for data fetching and the analyzer's pipeline plumbing.

Each test pins a defect that was found and fixed, so it cannot come back
unnoticed:

* downloads that bypassed the verified fetch layer (no retry, no validation,
  non-atomic writes);
* the unused whole-structure Biopython SASA computed on every structure;
* fpocket without a time limit;
* the APBS *energy* returned in place of a pocket *potential*;
* the fetch CLI's exit status, which CI relies on to fail loudly.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]

#: Modules allowed to open network connections themselves. Everything else
#: fetches through cryptic_ip.database.async_fetch.
NETWORK_MODULES = {
    "cryptic_ip/database/async_fetch.py",
    "cryptic_ip/database/rcsb_client.py",  # its own retrying, validating client
    "cryptic_ip/database/batch_processing.py",  # UniProt listing via a retried session
    "cryptic_ip/database/alphafold_client.py",  # metadata endpoint
    "cryptic_ip/database/pdb_client.py",  # metadata endpoint
    "cryptic_ip/database/uniprot_client.py",  # annotation endpoint
}
#: A download call, not a mention of one in prose.
RAW_DOWNLOAD = re.compile(r"urllib\.request\.url(open|retrieve)\(|[\"']wget[\"']|^\s*wget\s", re.MULTILINE)


def _python_sources():
    for folder in ("cryptic_ip", "scripts"):
        for path in (ROOT / folder).rglob("*.py"):
            yield path.relative_to(ROOT).as_posix(), path.read_text(encoding="utf-8")


@pytest.mark.parametrize("name, text", list(_python_sources()), ids=lambda v: v if isinstance(v, str) and v.endswith(".py") else "")
def test_no_unverified_downloads_in_code(name, text):
    """Structure downloads go through the verified, retried fetch layer."""
    if name in NETWORK_MODULES:
        return
    assert not RAW_DOWNLOAD.search(text), f"{name} downloads without verification or retry"


@pytest.mark.parametrize("path", sorted((ROOT / ".github" / "workflows").glob("*.yml")), ids=lambda p: p.name)
def test_workflows_parse_and_do_not_download_raw(path):
    text = path.read_text(encoding="utf-8")
    workflow = yaml.safe_load(text)
    assert "jobs" in workflow
    if path.name == "data-validation.yml":  # endpoint health probes only
        return
    assert not RAW_DOWNLOAD.search(text), f"{path.name} downloads outside scripts/fetch_structures.py"


def test_workflow_actions_are_off_node20():
    """Majors that still declare node20 must not creep back in."""
    stale = {
        "actions/checkout@v4", "actions/cache@v4", "actions/upload-artifact@v4",
        "actions/upload-artifact@v5", "actions/download-artifact@v4",
        "actions/download-artifact@v5", "actions/download-artifact@v6",
        "actions/setup-python@v5", "conda-incubator/setup-miniconda@v3",
    }
    for path in (ROOT / ".github" / "workflows").glob("*.yml"):
        used = set(re.findall(r"uses:\s*([\w./-]+@v\d+)", path.read_text()))
        assert not used & stale, f"{path.name}: {sorted(used & stale)}"


# ------------------------------------------------------------------ analyzer
def _analyzer(tmp_path: Path):
    from cryptic_ip.analysis.analyzer import ProteinAnalyzer

    pdb = tmp_path / "tiny.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n"
        "ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00 90.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    return ProteinAnalyzer(str(pdb), work_dir=str(tmp_path / "work"), skip_electrostatics=True)


def test_run_pipeline_does_not_compute_the_unused_biopython_sasa(tmp_path, monkeypatch):
    import pandas as pd

    analyzer = _analyzer(tmp_path)
    calls = []
    monkeypatch.setattr(analyzer, "detect_pockets", lambda *a, **k: calls.append("detect"))
    monkeypatch.setattr(analyzer, "calculate_sasa", lambda: calls.append("biopython_sasa"))
    monkeypatch.setattr(analyzer, "score_all_pockets", lambda: pd.DataFrame())
    analyzer.run_pipeline(include_electrostatics=False)
    assert calls == ["detect"]


def test_fpocket_timeout_is_an_error_for_that_structure(tmp_path, monkeypatch):
    analyzer = _analyzer(tmp_path)
    analyzer.fpocket_timeout_s = 1.0

    def fake_run(cmd, **kwargs):
        if "-h" in cmd:
            return subprocess.CompletedProcess(cmd, 0, "", "")
        assert kwargs.get("timeout") == 1.0
        raise subprocess.TimeoutExpired(cmd, kwargs["timeout"])

    monkeypatch.setattr(subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="fpocket exceeded"):
        analyzer.detect_pockets()


def test_no_potential_map_means_no_potential(tmp_path):
    """The whole-structure APBS energy must never stand in for a pocket potential."""
    analyzer = _analyzer(tmp_path)
    analyzer.electrostatic_data = -12345.6  # an energy, in kJ/mol
    analyzer.electrostatic_map_path = None
    assert analyzer.pocket_electrostatic_potential((0.0, 0.0, 0.0)) is None
    analyzer.electrostatic_map_path = tmp_path / "absent.dx"
    assert analyzer.pocket_electrostatic_potential((0.0, 0.0, 0.0)) is None


def test_unreadable_potential_map_means_no_potential(tmp_path):
    analyzer = _analyzer(tmp_path)
    analyzer.electrostatic_data = -12345.6
    broken = tmp_path / "broken.dx"
    broken.write_text("not a dx map\n")
    analyzer.electrostatic_map_path = broken
    assert analyzer.pocket_electrostatic_potential((0.0, 0.0, 0.0)) is None
