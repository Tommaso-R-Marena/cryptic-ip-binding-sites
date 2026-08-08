"""Tests for the population-level burial survey.

The survey exists to put the burial boundary on a distribution rather than on a
handful of controls, so what matters is that its threshold estimators report
honestly: they must find a split when the data is genuinely split, and must not
manufacture one when it is not.
"""

from __future__ import annotations

import numpy as np
import pytest

from scripts.burial_survey import (
    EntryMeasurement,
    density_minimum,
    otsu_threshold,
    resolve_structures,
    summarise,
)


class TestThresholdEstimators:
    def test_density_minimum_finds_the_trough_of_a_split_distribution(self):
        rng = np.random.default_rng(0)
        buried = rng.normal(0.09, 0.02, 200)
        exposed = rng.normal(0.40, 0.05, 200)
        values = np.clip(np.concatenate([buried, exposed]), 0.0, 1.0)

        trough = density_minimum(values)
        assert trough is not None
        # The trough must lie between the modes, not on top of one.
        assert 0.12 < trough < 0.35

    def test_density_minimum_reports_none_for_a_continuous_distribution(self):
        """No trough is the honest answer when the data is unimodal.

        Returning a number here would let a boundary be justified by a statistic
        that had no evidence behind it.
        """
        rng = np.random.default_rng(1)
        values = np.clip(rng.normal(0.30, 0.06, 400), 0.0, 1.0)
        assert density_minimum(values) is None

    def test_otsu_splits_a_bimodal_distribution_between_the_modes(self):
        rng = np.random.default_rng(2)
        values = np.clip(
            np.concatenate([rng.normal(0.08, 0.02, 200), rng.normal(0.45, 0.05, 200)]),
            0.0,
            1.0,
        )
        threshold = otsu_threshold(values)
        assert threshold is not None
        assert 0.10 < threshold < 0.40

    def test_estimators_need_data(self):
        assert density_minimum([]) is None
        assert otsu_threshold([]) is None
        assert density_minimum([0.1, 0.2]) is None


class TestSummary:
    def _measurement(self, pdb_id, relative_sasa, **kwargs):
        return EntryMeasurement(
            pdb_id=pdb_id,
            relative_sasa=relative_sasa,
            burial_class=kwargs.get("burial_class", "surface"),
            series=kwargs.get("series", "InsP6"),
            comp_id=kwargs.get("comp_id", "IHP"),
            n_basic_residues=kwargs.get("n_basic_residues", 4),
        )

    def test_summary_counts_measured_and_failed_entries(self):
        measurements = [
            self._measurement("1AAA", 0.09, burial_class="cryptic"),
            self._measurement("2BBB", 0.40),
            EntryMeasurement(pdb_id="3CCC", error="no phosphorylated inositol ligand detected"),
        ]
        summary = summarise(measurements)

        assert summary["n_entries"] == 3
        assert summary["n_measured"] == 2
        assert summary["n_failed"] == 1
        assert summary["failure_reasons"] == {
            "no phosphorylated inositol ligand detected": 1
        }
        assert summary["class_counts"] == {"cryptic": 1, "surface": 1}
        assert summary["series_counts"] == {"InsP6": 2}

    def test_summary_is_strictly_serialisable(self):
        """The survey output is an artifact, so it must be valid JSON."""
        import json

        measurements = [self._measurement(f"{i:04d}", 0.1 + i * 0.01) for i in range(20)]
        summary = summarise(measurements)
        # allow_nan=False raises rather than emitting an invalid token.
        json.dumps(summary, allow_nan=False)

    def test_summary_survives_a_set_with_no_usable_measurements(self):
        summary = summarise([EntryMeasurement(pdb_id="1AAA", error="boom")])
        assert summary["n_measured"] == 0
        assert "relative_sasa" not in summary

    def test_configured_boundary_is_reported_alongside_the_estimates(self):
        """The point is comparison: what the data says vs what is configured."""
        rng = np.random.default_rng(3)
        values = np.clip(
            np.concatenate([rng.normal(0.08, 0.02, 60), rng.normal(0.45, 0.05, 60)]),
            0.0,
            1.0,
        )
        measurements = [self._measurement(f"{i:04d}", float(v)) for i, v in enumerate(values)]
        summary = summarise(measurements)

        candidates = summary["candidate_boundaries"]
        assert "density_minimum" in candidates
        assert "otsu" in candidates
        assert candidates["configured_cryptic_max"] == pytest.approx(0.12)
        assert summary["n_below_configured_boundary"] > 0


class TestStructureResolution:
    def test_missing_structures_are_skipped_not_fabricated(self, tmp_path):
        csv = tmp_path / "entries.csv"
        csv.write_text("pdb_id\n1AAA\n2BBB\n", encoding="utf-8")
        structures = tmp_path / "structures"
        structures.mkdir()
        (structures / "1AAA.pdb").write_text("", encoding="utf-8")

        resolved = resolve_structures(structures, csv)
        assert [entry["pdb_id"] for entry in resolved] == ["1AAA"]

    def test_limit_caps_the_survey(self, tmp_path):
        csv = tmp_path / "entries.csv"
        csv.write_text("pdb_id\n1AAA\n2BBB\n3CCC\n", encoding="utf-8")
        structures = tmp_path / "structures"
        structures.mkdir()
        for pdb_id in ("1AAA", "2BBB", "3CCC"):
            (structures / f"{pdb_id}.pdb").write_text("", encoding="utf-8")

        assert len(resolve_structures(structures, csv, limit=2)) == 2


def test_all_digit_identifiers_are_not_mangled_into_numbers(tmp_path):
    """Identifiers are text: 0012 must not become 12 and lose its file.

    Pandas infers an all-digit column as integers, which strips leading zeros,
    so the structure the entry names would silently never be found and the entry
    would be dropped from the survey without an error.
    """
    csv = tmp_path / "entries.csv"
    csv.write_text("pdb_id\n0012\n1ZY7\n", encoding="utf-8")
    structures = tmp_path / "structures"
    structures.mkdir()
    (structures / "0012.pdb").write_text("", encoding="utf-8")
    (structures / "1ZY7.pdb").write_text("", encoding="utf-8")

    resolved = resolve_structures(structures, csv)
    assert sorted(entry["pdb_id"] for entry in resolved) == ["0012", "1ZY7"]
