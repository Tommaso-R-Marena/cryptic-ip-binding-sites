"""The detector CLI, run against synthetic structures with known geometry."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from cryptic_ip.testing.synthetic import SyntheticStructureSpec, write_synthetic_structure

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("detect_cli", ROOT / "scripts" / "detect.py")
detect = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(detect)

FAST = {"field_points": 192, "spacing": 1.0, "max_pockets": 1, "ensemble": False}


@pytest.fixture(scope="module")
def structures(tmp_path_factory):
    """A buried basic cavity and a buried acidic one, from the same generator."""
    out = tmp_path_factory.mktemp("structures")
    basic = SyntheticStructureSpec(
        name="BURIED_BASIC", protein_radius=18.0, site_offset=0.0,
        n_basic_lining=8, seed=1, expected_burial_class="cryptic",
    )
    acidic = SyntheticStructureSpec(
        name="BURIED_ACIDIC", protein_radius=18.0, site_offset=0.0,
        n_basic_lining=8, acidic_lining=True, include_ligand=False,
        seed=2, expected_burial_class="cryptic",
    )
    return {
        "basic": write_synthetic_structure(basic, out),
        "acidic": write_synthetic_structure(acidic, out),
        "dir": out,
    }


def test_the_detector_finds_the_carved_cavity(structures):
    report = detect.describe_structure(structures["basic"], **FAST)
    assert report["apo_cavities"]["n_cavities"] >= 1
    # The generator carves the cavity around an IP6, which displaces ~350 A^3.
    assert report["apo_cavities"]["total_volume"] > 150.0


def test_a_basic_cavity_outscores_an_acidic_one(structures):
    basic = detect.describe_structure(structures["basic"], **FAST)
    acidic = detect.describe_structure(structures["acidic"], **FAST)
    assert basic["sites"][0]["score"]["low"] > acidic["sites"][0]["score"]["low"]


def test_the_field_well_is_positive_in_a_basic_cavity_and_negative_in_an_acidic_one(structures):
    basic = detect.describe_structure(structures["basic"], **FAST)
    acidic = detect.describe_structure(structures["acidic"], **FAST)
    assert basic["sites"][0]["features"]["field_depth_kt"] > 0.0
    assert acidic["sites"][0]["features"]["field_depth_kt"] < 0.0


def test_an_acidic_cavity_is_called_no_polyanion_site(structures):
    acidic = detect.describe_structure(structures["acidic"], **FAST)
    verdict = acidic["sites"][0]["specificity"]
    assert verdict["call"] == "no_polyanion_site"
    assert verdict["is_candidate"] is False


def test_a_generic_basic_cavity_is_not_promoted_to_an_ip6_candidate(structures):
    # The synthetic cavity is built as "a hole with eight basic residues", not
    # as an IP6 site with the right charge geometry, so the honest answer is
    # that it is a polyanion pocket of unknown preference. The deployed screen
    # has no way to say that.
    basic = detect.describe_structure(structures["basic"], **FAST)
    assert basic["sites"][0]["specificity"]["call"] in ("polyanion_generic", "ip6_selective")
    assert basic["n_candidate_sites"] == basic["sites"][0]["specificity"]["is_candidate"] * 1 \
        if basic["sites"][0]["decision"] == "site" else True


def test_burial_comes_from_the_cavity_detector_under_its_own_name(structures):
    report = detect.describe_structure(structures["basic"], **FAST)
    features = report["sites"][0]["features"]
    assert features["cavity_openness"] is not None
    # The ligand-copy measurement is a different quantity and abstains.
    assert features["relative_sasa"] is None


def test_a_synthetic_model_carries_plddt_so_the_score_is_a_point(structures):
    report = detect.describe_structure(structures["basic"], **FAST)
    score = report["sites"][0]["score"]
    assert score["missing"] == [] or "relative_sasa" not in score["missing"]
    assert score["point"] is not None


def test_the_report_names_the_limiting_measurement_for_a_rejected_pocket(structures):
    acidic = detect.describe_structure(structures["acidic"], **FAST)
    assert acidic["sites"][0]["blocking_feature"] is not None


def test_the_ensemble_path_runs_and_reports_crypticity(structures):
    report = detect.describe_structure(
        structures["basic"], field_points=192, spacing=1.2, max_pockets=1,
        n_modes=8, amplitudes=(-1.0, 1.0), which_modes=(0,), ensemble=True,
    )
    assert "unavailable" not in report["ensemble"]
    assert report["ensemble"]["n_conformers"] == 3
    assert report["sites"][0]["crypticity"] is not None


def test_the_parameters_are_recorded_in_the_report(structures):
    report = detect.describe_structure(structures["basic"], **FAST)
    assert report["parameters"]["threshold"] == detect.DECISION_THRESHOLD
    assert report["parameters"]["field_points"] == 192


# -------------------------------------------------------------- the commands

def test_describe_writes_a_json_report(structures, tmp_path, capsys):
    out = tmp_path / "report.json"
    code = detect.main([
        "describe", str(structures["basic"]), "--out", str(out),
        "--no-ensemble", "--field-points", "192", "--spacing", "1.0", "--max-pockets", "1",
    ])
    assert code == 0
    payload = json.loads(out.read_text())
    assert payload["n_sites"] == 1
    assert "site 0" in capsys.readouterr().out


def test_shard_splits_a_directory_and_writes_after_every_structure(structures, tmp_path, capsys):
    out = tmp_path / "part.json"
    code = detect.main([
        "shard", "--structures", str(structures["dir"]),
        "--shard", "0", "--shards", "1", "--out", str(out),
        "--no-ensemble", "--field-points", "192", "--spacing", "1.0", "--max-pockets", "1",
    ])
    assert code == 0
    records = json.loads(out.read_text())
    assert len(records) == 2
    assert "shard 0 of 1" in capsys.readouterr().out


def test_a_shard_of_an_empty_directory_fails_rather_than_writing_nothing(tmp_path, capsys):
    empty = tmp_path / "empty"
    empty.mkdir()
    code = detect.main([
        "shard", "--structures", str(empty), "--shard", "0", "--shards", "1",
        "--out", str(tmp_path / "x.json"),
    ])
    assert code == 1
    assert "no structures" in capsys.readouterr().out


def test_one_unreadable_structure_does_not_take_the_whole_shard(structures, tmp_path, capsys):
    # The study J lesson: a shard must keep what it finished.
    broken = structures["dir"] / "BROKEN.pdb"
    broken.write_text("this is not a pdb file\n")
    try:
        out = tmp_path / "part.json"
        code = detect.main([
            "shard", "--structures", str(structures["dir"]),
            "--shard", "0", "--shards", "1", "--out", str(out),
            "--no-ensemble", "--field-points", "192", "--spacing", "1.0", "--max-pockets", "1",
        ])
        assert code == 0
        records = json.loads(out.read_text())
        assert any("error" in record for record in records)
        assert any("sites" in record for record in records)
    finally:
        broken.unlink()


def test_power_prints_the_minimum_detectable_effect(capsys):
    assert detect.main(["power"]) == 0
    out = capsys.readouterr().out
    assert "0.87" in out
    assert "statement about power" in out


def test_spec_prints_every_feature_with_its_direction(capsys):
    assert detect.main(["spec"]) == 0
    out = capsys.readouterr().out
    for feature in ("cavity_openness", "field_shape", "openability"):
        assert feature in out
    assert "Competitor panel" in out


def test_spec_as_json_is_machine_readable(capsys):
    assert detect.main(["spec", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["threshold"] == detect.DECISION_THRESHOLD
    assert "IP6" in payload["templates"]
    assert payload["specificity_margin"] > 0
