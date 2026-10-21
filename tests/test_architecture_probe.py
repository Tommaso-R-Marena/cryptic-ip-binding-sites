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
import pytest

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


def test_the_power_probe_reports_both_halves_and_their_better():
    """Reading one half as the detection rate is the mistake this pins shut."""
    out = script.probe_power(rng(), n_draws=200)
    for row in out["dilution"] + out["shapes"] + [out["false_positive"]]:
        for key in ("detected_topology", "detected_register", "detected_either"):
            assert key in row, key
        # 'either' is the better of the two, so it can never be below one.
        assert row["detected_either"] >= row["detected_topology"] - 1e-9
        assert row["detected_either"] >= row["detected_register"] - 1e-9


def test_the_power_probe_shows_dilution_and_a_shape_contrast():
    out = script.probe_power(rng(), n_draws=200)
    full = [r for r in out["dilution"] if r["fraction_of_face"] == 1.0]
    thin = [r for r in out["dilution"] if r["fraction_of_face"] == 0.33]
    assert full and thin
    assert min(r["detected_either"] for r in full) > max(r["detected_either"] for r in thin)
    shapes = {r["shape"]: r for r in out["shapes"]}
    assert len(shapes) == 2
    ring = next(v for k, v in shapes.items() if "ring" in k)
    column = next(v for k, v in shapes.items() if "column" in k)
    assert column["lining_size"] == ring["lining_size"]
    assert column["detected_either"] >= ring["detected_either"]


def test_a_detection_rate_never_ships_without_its_false_positive_rate():
    out = script.probe_power(rng(), n_draws=200)
    decoy = out["false_positive"]
    assert decoy["lining_size"] == out["shapes"][0]["lining_size"]
    assert decoy["detected_either"] < 0.25, "decoys should rarely be called structured"


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
        "orderings_floor", "calibration", "panel", "power", "transfer", "local",
    }


# ------------------------------------------------------- the committed run

def test_the_committed_artefact_records_the_half_swap():
    """Columns are a topology signal, rings a register one. The sharper claim."""
    payload = json.loads(ARTEFACT.read_text())
    shapes = {str(r["shape"]): r for r in payload["power"]["shapes"]}
    ring = next(v for k, v in shapes.items() if "ring" in k)
    column = next(v for k, v in shapes.items() if "column" in k)
    assert column["detected_topology"] > ring["detected_topology"]
    assert ring["detected_register"] > ring["detected_topology"]


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
    assert column["detected_either"] >= ring["detected_either"]
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


# ------------------------------------------------- transfer to irregularity

def test_the_transfer_probe_covers_jitter_and_bulges():
    out = script.probe_transfer(rng(), n_draws=150, reps=4)
    labels = [r["condition"] for r in out["rows"]]
    assert "ideal" in labels
    assert any("jitter" in label for label in labels)
    assert any("bulge" in label for label in labels)
    for row in out["rows"]:
        for key in ("detected_topology", "detected_register", "detected_either"):
            assert 0.0 <= row[key] <= 1.0


def test_jitter_leaves_the_covalent_geometry_exactly_ideal():
    """Only torsions move, so an irregular fixture is still a valid structure."""
    for jitter in (0.0, script.gf.REALISTIC_JITTER, 20.0):
        report = script.gf.bond_length_report(
            script.gf.helix_barrel(4, 16, 9.0, jitter=jitter, seed=1))
        assert report["CA_CB_mean"] == pytest.approx(1.530, abs=1e-3)
        assert report["CA_CB_max_deviation"] < 1e-9


def test_a_bulge_leaves_the_covalent_geometry_exactly_ideal():
    report = script.gf.bond_length_report(
        script.gf.helix_barrel(4, 16, 9.0, bulges={0: [8]}, seed=1))
    assert report["CA_CB_max_deviation"] < 1e-9


def test_the_committed_transfer_result_is_a_fall_not_a_rise():
    """The headline negative, pinned so it cannot be quietly reversed."""
    payload = json.loads(ARTEFACT.read_text())
    transfer = payload["transfer"]
    assert transfer["ideal_rate"] > transfer["realistic_rate"]
    assert transfer["worst_rate"] <= transfer["realistic_rate"]
    worst = min(transfer["rows"], key=lambda r: r["detected_either"])
    assert "bulge" in worst["condition"] or "jitter" in worst["condition"]


def test_the_report_says_the_transfer_result_is_a_hard_negative():
    text = script.markdown(json.loads(ARTEFACT.read_text()))
    assert "hard negative" in text
    assert "upper bound" in text


# ------------------------------------------- the local statistic, sections 10

def test_the_slab_helper_gives_the_transfer_and_local_probes_the_same_fixture():
    """The paired comparison is only meaningful if both see one lining per seed."""
    settings = {"jitter": script.gf.REALISTIC_JITTER, "bulges": None}
    first, size = script._inward_slab(settings, 3, 14, np.random.default_rng(11))
    second, again = script._inward_slab(settings, 3, 14, np.random.default_rng(11))
    assert first == second
    assert size == again
    assert len(first) == 14


def test_the_slab_helper_reports_a_face_too_small_rather_than_clipping():
    settings = {"jitter": 0.0, "bulges": None}
    pick, face = script._inward_slab(settings, 0, 10_000, np.random.default_rng(1))
    assert pick is None
    assert face > 0


def test_mcnemar_is_exact_and_two_sided():
    assert script._mcnemar(0, 0) == 1.0
    assert script._mcnemar(5, 5) == pytest.approx(1.0)
    assert script._mcnemar(6, 0) == pytest.approx(2 * 0.5 ** 6)
    assert script._mcnemar(0, 6) == pytest.approx(2 * 0.5 ** 6)


def test_the_local_probe_pairs_every_arm_and_scores_every_variant():
    out = script.probe_local(rng(), n_draws=120, reps=3)
    assert out["window_sizes"] == list(script.ar.WINDOW_SIZES)
    assert set(out["adopted"]) == set(script.VARIANTS)
    for row in out["rows"]:
        assert 0.0 <= row["pooled_detected"] <= 1.0
        for name in script.VARIANTS:
            assert 0.0 <= row[f"{name}_detected"] <= 1.0
            # A paired table has to add up: concordant plus discordant is every
            # fixture, never more.
            assert (row[f"{name}_both"] + row[f"{name}_pooled_only"]
                    + row[f"{name}_only"]) <= row["reps"]
            assert row[f"{name}_difference"] == pytest.approx(
                row[f"{name}_detected"] - row["pooled_detected"])
            assert set(row[f"{name}_fraction_by_size"]) == set(row["null_mean_by_size"])


def test_every_variant_selects_a_subset_of_the_fixed_window_sizes():
    out = script.probe_local(rng(), n_draws=120, reps=3)
    for row in out["rows"]:
        # The rank keeps exactly one size per lining, so its fractions sum to 1.
        assert sum(row["best_fraction_by_size"].values()) == pytest.approx(1.0)
        for name in script.VARIANTS:
            assert all(0.0 <= rate <= 1.0
                       for rate in row[f"{name}_fraction_by_size"].values())


def test_the_decision_rule_is_applied_and_not_asserted():
    """`adopted` must follow from the measured numbers, in both directions."""
    out = script.probe_local(rng(), n_draws=120, reps=3)
    primary = next(r for r in out["rows"] if "realistic" in r["arm"])
    for name in script.VARIANTS:
        expected = (primary[f"{name}_detected"] > primary["pooled_detected"]
                    and out["decoys"][f"{name}_false_positive"] <= 0.05)
        assert out["adopted"][name] is bool(expected)


def test_the_selection_note_says_whether_the_threshold_actually_fired():
    quiet = [{"pruned_fraction_by_size": {"6": 1.0, "8": 1.0, "10": 1.0},
              "best_fraction_by_size": {"6": 0.0, "8": 0.0, "10": 1.0}}]
    fired = [{"pruned_fraction_by_size": {"6": 0.5, "8": 1.0, "10": 1.0},
              "best_fraction_by_size": {"6": 0.0, "8": 0.2, "10": 0.8}}]
    local = {"saturated_null": 0.9}
    assert "never fired" in " ".join(script._selection_note(quiet, local))
    assert "never fired" not in " ".join(script._selection_note(fired, local))
    assert "[6]" in " ".join(script._selection_note(fired, local))


def test_the_conclusion_follows_the_measurement_in_both_directions():
    """The report must not assert the negative it was written expecting."""
    names = list(script.VARIANTS)
    lost = [{"arm": "jitter 8 deg (realistic)", "pooled_detected": 0.30,
             **{f"{n}_detected": 0.0 for n in names}}]
    won = [{"arm": "jitter 8 deg (realistic)", "pooled_detected": 0.30,
            **{f"{n}_detected": 0.0 for n in names}, "best_detected": 0.60}]
    assert "Locality is not" in script._conclusion(lost, names)
    assert "beat the pooled statistic" in script._conclusion(won, names)
    assert "held-out" in script._conclusion(won, names)
