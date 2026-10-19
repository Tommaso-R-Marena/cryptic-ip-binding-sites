"""Tests for the figure builder.

A figure that disagrees with the report it illustrates is worse than no figure,
so the checks here are that every figure is built from the committed artefacts
and that it fails loudly rather than silently when a key it needs is missing.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "architecture_figures", ROOT / "scripts" / "architecture_figures.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(script)

ARCHITECTURE = ROOT / "results" / "architecture" / "architecture.json"
HOMOLOGY = ROOT / "results" / "architecture" / "homology.json"
FIGURES = ROOT / "results" / "architecture" / "figures"


def test_every_figure_is_declared():
    assert set(script.FIGURES) == {"power", "shape", "homology", "transfer", "local"}


def test_only_validated_palette_slots_are_used():
    """The first three categorical slots are the ones that clear all-pairs."""
    assert script.SERIES == ("#2a78d6", "#eb6834", "#1baf7a")


def test_no_figure_uses_a_fourth_series():
    text = (ROOT / "scripts" / "architecture_figures.py").read_text()
    assert "SERIES[3]" not in text


def test_the_figures_build_from_the_committed_artefacts(tmp_path):
    produced = script.build(ARCHITECTURE, HOMOLOGY, tmp_path,
                            ROOT / "results" / "architecture" / "real" / "real.json")
    assert len(produced) == len(script.FIGURES)
    for path in produced:
        assert path.exists() and path.stat().st_size > 0
        assert path.with_suffix(".pdf").exists()


def test_a_missing_key_fails_loudly(tmp_path):
    """Silence here would mean a figure drawn from nothing."""
    broken = tmp_path / "architecture.json"
    payload = json.loads(ARCHITECTURE.read_text())
    del payload["power"]["shapes"]
    broken.write_text(json.dumps(payload))
    script.style()
    with pytest.raises(KeyError):
        script.figure_shape(json.loads(broken.read_text()), tmp_path / "F2")


def test_the_committed_figures_exist_in_both_formats():
    for stem in ("F1_power", "F2_shape", "F3_homology", "F4_transfer"):
        assert (FIGURES / f"{stem}.png").exists()
        assert (FIGURES / f"{stem}.pdf").exists()


def test_detection_degrades_gracefully_across_the_real_lining_range():
    """The band is marked because the method's limit sits in it, not below it.

    Reading the better of the two halves, detection is complete at the top of
    the range and still well above the false-positive rate at the bottom. An
    earlier version of this test asserted a collapse, which was true only of
    the topology half read alone.
    """
    low, high = script.REAL_LINING_RANGE
    assert low < high
    payload = json.loads(ARCHITECTURE.read_text())
    inside = [r for r in payload["power"]["dilution"]
              if low <= r["n_positions"] <= high and r["mode"] == "spatial"]
    assert len(inside) >= 2
    decoy = payload["power"]["false_positive"]["detected_either"]
    assert max(r["detected_either"] for r in inside) > 0.9
    assert min(r["detected_either"] for r in inside) > decoy + 0.2


def test_the_figure_reads_the_composite_rate_not_one_half():
    text = (ROOT / "scripts" / "architecture_figures.py").read_text()
    assert "detected_either" in text


def test_the_shape_figure_has_exactly_two_shapes_of_one_size():
    payload = json.loads(ARCHITECTURE.read_text())
    shapes = payload["power"]["shapes"]
    assert len(shapes) == 2
    assert len({r["lining_size"] for r in shapes}) == 1


def test_the_shape_figure_shows_both_halves_because_they_swap():
    payload = json.loads(ARCHITECTURE.read_text())
    shapes = {str(r["shape"]): r for r in payload["power"]["shapes"]}
    ring = next(v for k, v in shapes.items() if "ring" in k)
    column = next(v for k, v in shapes.items() if "column" in k)
    assert ring["detected_register"] > ring["detected_topology"]
    assert column["detected_topology"] >= column["detected_register"] - 1e-9


def test_the_homology_figure_has_its_controls_and_a_square_matrix():
    payload = json.loads(HOMOLOGY.read_text())
    rows = payload["retrieval"]
    assert sum(1 for r in rows if r["is_control"]) == 2
    matrix = payload["distance_matrix"]
    assert len(matrix) == len(matrix[0]) == payload["n_linings"]
    assert len(payload["labels"]) == payload["n_linings"]
