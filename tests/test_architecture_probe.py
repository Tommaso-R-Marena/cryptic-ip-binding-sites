"""Tests for the architecture probe driver.

The probes are the evidence behind every constant in `architecture.py`, so the
checks here are that each one still reports what the module's docstrings and
`docs/ARCHITECTURE_PLAN.md` claim it reports. A probe that silently stopped
measuring would leave the constants unjustified.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "architecture_probe", ROOT / "scripts" / "architecture_probe.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(script)

ARTEFACT = ROOT / "results" / "architecture" / "architecture.json"


def rng():
    return np.random.default_rng(7)


# --------------------------------------------------------------- the barrel

def test_the_barrel_has_a_face_and_it_is_a_minority_of_residues():
    inward, heights = script.barrel_inward()
    assert 30 < len(inward) < 60
    assert len(heights) == len(inward)
    assert set(heights) == set(inward)


def test_linkers_renumber_without_changing_the_count():
    inward, _ = script.barrel_inward()
    renumbered = script.with_linkers(inward)
    assert len(renumbered) == len(inward)
    assert renumbered == sorted(renumbered)
    assert max(renumbered) > max(inward)


# ------------------------------------------------------------- the probes

def test_the_barcode_identity_probe_finds_no_mismatch():
    out = script.probe_identities(rng())
    assert out["barcode_mismatches"] == 0
    assert out["barcode_identity_holds"] is True


def test_the_axial_probe_still_finds_the_degeneracy():
    """If this ever passes, the statistic became usable and the docs are wrong."""
    out = script.probe_axial(rng())
    assert out["max_abs_deviation_at_period_2"] < 1e-12
    assert out["axial_statistic_usable"] is False


def test_the_scan_inflation_probe_shows_the_full_scan_is_worst():
    out = script.probe_scan_inflation(rng())
    by_name = {r["statistic"]: r for r in out["rows"]}
    single = by_name["single period 3.6"]["true_alpha_at_analytic_threshold"]
    full = by_name["full scan 2-12"]["true_alpha_at_analytic_threshold"]
    assert full > single
    assert full > 0.5, "a full scan at the analytic cutoff should admit most linings"


def test_the_peak_width_probe_recovers_a_flat_ratio():
    out = script.probe_peak_width()
    assert out["ratio_spread"] < 0.1, "the law should hold across spans"
    assert 0.4 < out["ratio_mean"] < 0.8
    assert out["worst_loss"] > 0.2, "a fixed scan should visibly lose a face"


def test_the_orderings_probe_recovers_the_one_over_k_floor():
    out = script.probe_orderings_floor()
    for row in out["rows"]:
        assert row["min_attainable_p"] >= 1.0 / row["orderings"] - 1e-9
    small = [r for r in out["rows"] if r["orderings"] < 20]
    assert small and all(r["min_attainable_p"] > 0.05 for r in small)


def test_the_calibration_probe_shows_the_family_needs_correcting():
    out = script.probe_calibration(rng(), n_linings=40, n_draws=200)
    assert out["family_alpha_holm_0.05"] <= out["family_alpha_uncorrected_0.05"] + 1e-9
    for entry in out["per_test"].values():
        assert 0.0 <= entry["fraction_runnable"] <= 1.0


def test_the_panel_probe_covers_both_a_face_and_a_random_lining():
    out = script.probe_panel(rng(), n_draws=300)
    by_name = {r["lining"]: r for r in out["rows"]}
    assert by_name["barrel + 10-residue linkers"]["architecture_class"] == "mosaic"
    assert by_name["random 14 in 200"]["architecture_class"] == "distributed"
    assert by_name["contiguous stretch of 14"]["architecture_class"] == "unresolved"


def test_the_power_probe_shows_dilution_and_a_shape_contrast():
    out = script.probe_power(rng(), n_draws=200)
    full = [r for r in out["dilution"] if r["fraction_of_face"] == 1.0]
    thin = [r for r in out["dilution"] if r["fraction_of_face"] == 0.33]
    assert full and thin
    assert min(r["detected_rate"] for r in full) > max(r["detected_rate"] for r in thin)
    shapes = {r["shape"]: r for r in out["shapes"]}
    assert len(shapes) == 2
    ring = next(v for k, v in shapes.items() if "ring" in k)
    column = next(v for k, v in shapes.items() if "column" in k)
    assert column["lining_size"] == ring["lining_size"]
    assert column["detected_rate"] > ring["detected_rate"]


# ---------------------------------------------------------- the markdown

def test_the_markdown_reports_every_section():
    report = script.run_all(quick=True)
    text = script.markdown(report)
    for heading in ("Two identities", "does not survive the lattice",
                    "analytic threshold", "peak-width law", "has a floor",
                    "Calibration", "Recovery of known architectures", "Power"):
        assert heading in text
    assert "cannot be used for" in text


def test_the_markdown_states_what_it_is_not_evidence_for():
    text = script.markdown(script.run_all(quick=True))
    assert "no annotated sequence has been read" in text
    assert "upper bound" in text


def test_every_probe_is_listed():
    assert set(script.PROBES) == {
        "identities", "axial", "scan_inflation", "peak_width",
        "orderings_floor", "calibration", "panel", "power",
    }


# ------------------------------------------------------- the committed run

def test_the_committed_artefact_holds_every_probe():
    payload = json.loads(ARTEFACT.read_text())
    assert set(script.PROBES) <= set(payload)
    assert payload["seed"] == script.SEED


def test_the_committed_artefact_agrees_with_the_plan_s_headline_numbers():
    """The plan quotes these; a drift between them is a silent contradiction."""
    payload = json.loads(ARTEFACT.read_text())
    plan = (ROOT / "docs" / "ARCHITECTURE_PLAN.md").read_text()
    assert payload["identities"]["barcode_mismatches"] == 0
    assert "zero disagreements" in plan
    assert payload["axial"]["max_abs_deviation_at_period_2"] == 0.0
    shapes = {r["shape"]: r for r in payload["power"]["shapes"]}
    ring = next(v for k, v in shapes.items() if "ring" in k)
    column = next(v for k, v in shapes.items() if "column" in k)
    assert column["detected_rate"] > ring["detected_rate"]
    assert "Q7-P1" in plan


def test_the_cli_runs_a_quick_pass(tmp_path):
    out = tmp_path / "architecture.json"
    markdown = tmp_path / "ARCHITECTURE.md"
    assert script.main(["run", "--quick", "--out", str(out),
                        "--markdown", str(markdown)]) == 0
    assert json.loads(out.read_text())["quick"] is True
    assert "Lining architecture" in markdown.read_text()


def test_the_cli_lists_the_probes(capsys):
    assert script.main(["list"]) == 0
    printed = capsys.readouterr().out
    for name in script.PROBES:
        assert name in printed
