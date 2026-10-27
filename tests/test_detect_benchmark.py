"""The benchmark harness itself, and the fixture-validity rules it rests on."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from cryptic_ip.detect import benchmark as BM
from cryptic_ip.detect import geometry_fixtures as GF

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("detect_benchmark", ROOT / "scripts" / "detect_benchmark.py")
cli = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(cli)


@pytest.fixture(scope="module")
def small_fixtures():
    """Smaller barrels, so the suite tests run in seconds."""
    return {
        "basic_lined": GF.helix_barrel(
            n_helices=6, residues_per_helix=12, barrel_radius=8.0,
            lining={0: "ARG", 2: "LYS", 4: "ARG"},
        ),
        "acidic_lined": GF.helix_barrel(
            n_helices=6, residues_per_helix=12, barrel_radius=8.0,
            lining={0: "ASP", 2: "GLU", 4: "ASP"},
        ),
        "narrow": GF.helix_barrel(n_helices=6, residues_per_helix=12, barrel_radius=7.0),
        "wide_open": GF.helix_barrel(n_helices=6, residues_per_helix=12, barrel_radius=10.0),
        "trp_plugged": GF.helix_barrel(
            n_helices=6, residues_per_helix=12, barrel_radius=7.0, lining={0: "TRP", 2: "TRP"},
        ),
    }


# -------------------------------------------------------------- the structure

def test_every_suite_declares_what_its_inputs_support():
    # A suite claiming more than its fixtures can show is the failure mode this
    # whole benchmark is designed against.
    for code, function, _ in BM.SUITES:
        del function
        assert code.startswith("B")
    report = BM.run_benchmark(only=["B6"])
    for suite in report.suites:
        assert suite.supports in BM.ALL_SUPPORTS


def test_an_unknown_suite_is_refused():
    with pytest.raises(BM.BenchmarkError, match="unknown suite"):
        BM.run_benchmark(only=["B99"])


def test_selecting_one_suite_runs_only_that_one():
    report = BM.run_benchmark(only=["B6"])
    assert [s.code for s in report.suites] == ["B6"]


def test_the_report_knows_whether_anything_failed():
    report = BM.run_benchmark(only=["B6"])
    assert report.ok
    assert report.as_dict()["n_failures"] == 0


def test_a_failing_measurement_makes_the_report_not_ok():
    suite = BM.Suite(code="BX", title="t", supports=BM.SUPPORTS_STATISTICS)
    suite.add(BM.Measurement(name="broken", value=0.0, passed=False))
    report = BM.BenchmarkReport(suites=[suite])
    assert not report.ok
    assert suite.failures


def test_measurements_without_a_verdict_do_not_fail_the_report():
    # Robustness and cost report curves; there is no principled threshold for
    # them and inventing one would be the same mistake as a tuned cutoff.
    suite = BM.Suite(code="BX", title="t", supports=BM.SUPPORTS_POSITIONS)
    suite.add(BM.Measurement(name="a curve", value=1.23, passed=None))
    assert BM.BenchmarkReport(suites=[suite]).ok


def test_the_markdown_states_what_the_benchmark_is_not():
    text = BM.run_benchmark(only=["B6"]).markdown()
    assert "say nothing" in text and "real binding sites" in text
    assert "DETECTOR_PLAN.md" in text


def test_the_report_round_trips_to_json():
    payload = BM.run_benchmark(only=["B6"]).as_dict()
    assert json.loads(json.dumps(payload))["n_suites"] == 1


# ----------------------------------------------------------------- the suites

def test_determinism_suite_passes_on_valid_fixtures(small_fixtures):
    suite = BM.suite_determinism(small_fixtures, spacing=1.2)
    assert suite.ok and suite.measurements


def test_invariance_suite_measures_rotation_and_translation(small_fixtures):
    suite = BM.suite_invariance(small_fixtures, points=192)
    names = {m.name for m in suite.measurements}
    assert "rotation discrepancy" in names
    assert "translation discrepancy" in names
    assert suite.ok


def test_robustness_suite_reports_a_curve_without_a_verdict(small_fixtures):
    suite = BM.suite_robustness(small_fixtures, noise_levels=(0.0, 0.5), points=192)
    assert suite.measurements
    assert all(m.passed is None for m in suite.measurements)


def test_discrimination_suite_separates_basic_from_acidic(small_fixtures):
    suite = BM.suite_discrimination(small_fixtures, points=192)
    gap = next(m for m in suite.measurements if m.name == "basic minus acidic well depth")
    assert gap.value > 0.0 and gap.passed


def test_cost_suite_fits_a_scaling_exponent():
    suite = BM.suite_cost(sizes=(4, 6, 8), points=192, repeats=1)
    exponents = [m for m in suite.measurements if "scaling exponent" in m.name]
    assert exponents, "the exponent is the useful number, not the absolute timings"
    assert all(m.value < 2.0 for m in exponents)


def test_recovery_suite_checks_the_estimator_against_known_truth():
    suite = BM.suite_recovery(seed=7)
    names = {m.name for m in suite.measurements}
    assert any("label frequency recovered" in n for n in names)
    assert "conformal coverage at alpha=0.1" in names
    assert "elastic network rigid-body modes" in names
    assert suite.ok


def test_head_to_head_feeds_both_scorers_the_same_measurements(small_fixtures):
    suite = BM.suite_head_to_head(small_fixtures, points=192)
    names = {m.name for m in suite.measurements}
    assert "deployed composite: worst positive minus best negative" in names
    assert "monotone score: worst positive minus best negative" in names


def test_head_to_head_admits_no_constructed_negative(small_fixtures):
    suite = BM.suite_head_to_head(small_fixtures, points=192)
    admitted = next(
        m for m in suite.measurements
        if m.name == "constructed negatives the monotone score admits"
    )
    assert admitted.passed, admitted.detail


# ------------------------------------------------------------- ground truth

def test_only_the_fixture_with_both_geometry_and_chemistry_is_a_positive():
    # The labels were wrong at first: two poly-alanine channels were marked
    # positive because they are buried, when burial is not site-ness. Both
    # scorers duly "failed" and the benchmark had found a flaw in itself.
    assert BM.GROUND_TRUTH["basic_lined"] == 1
    assert sum(BM.GROUND_TRUTH.values()) == 1


def test_every_label_publishes_its_reason():
    assert set(BM.GROUND_TRUTH) == set(BM.GROUND_TRUTH_REASON)
    for name, reason in BM.GROUND_TRUTH_REASON.items():
        assert len(reason) > 20, f"{name} has no usable reason"


def test_the_neutral_channels_are_negatives_and_say_why():
    for name in ("narrow", "trp_plugged"):
        assert BM.GROUND_TRUTH[name] == 0
        assert "neutral" in BM.GROUND_TRUTH_REASON[name]


def test_measure_fixture_feeds_both_scorers_from_one_measurement(small_fixtures):
    measured = BM.measure_fixture(small_fixtures["basic_lined"], points=192)
    assert measured is not None
    # Everything both scorers need comes from this one dictionary.
    for key in ("volume", "depth", "sasa", "basic_count", "potential", "enclosure",
                "openness", "shape", "sealed_fraction", "plddt"):
        assert key in measured


def test_a_fixture_with_no_cavity_returns_nothing_rather_than_zeroes():
    tiny = GF.helix_barrel(n_helices=3, residues_per_helix=8, barrel_radius=22.0)
    assert BM.measure_fixture(tiny, points=192) is None


# -------------------------------------------------------------------- the CLI

def test_the_cli_lists_the_suites(capsys):
    assert cli.main(["list"]) == 0
    out = capsys.readouterr().out
    for code in ("B1", "B4", "B6", "B7"):
        assert code in out


def test_the_cli_writes_both_report_formats(tmp_path, capsys):
    out_json = tmp_path / "bench.json"
    out_md = tmp_path / "bench.md"
    code = cli.main([
        "run", "--only", "B6", "--out", str(out_json), "--markdown", str(out_md),
    ])
    assert code == 0
    assert json.loads(out_json.read_text())["ok"] is True
    assert "Detector benchmark" in out_md.read_text()
    capsys.readouterr()


def test_the_cli_refuses_an_unknown_suite(capsys):
    assert cli.main(["run", "--only", "NOPE"]) == 1
    assert "unknown suite" in capsys.readouterr().out
