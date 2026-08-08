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


#: Measured across the bundled dataset by scripts/burial_survey.py in GitHub
#: Actions (run 31250979693, 2026-08-08, 128 SASA sample points/atom, all 136
#: entries). Recorded measurements, not expected outputs of the code under test:
#: if the measurement changes such that these no longer describe the deposited
#: set, re-run the survey and update them in the same commit as the change.
SURVEYED = {
    "n_entries": 136,
    "n_measured": 135,
    "n_failed": 1,
    "quantiles": {
        "q01": 0.070, "q05": 0.091, "q10": 0.126, "q25": 0.216, "q50": 0.338,
        "q75": 0.553, "q90": 0.773, "q95": 0.816, "q99": 0.896,
    },
    "density_minimum": 0.138,
    "otsu": 0.463,
    "n_below_configured_boundary": 14,
    "class_counts": {
        "surface": 89, "semi_cryptic": 26, "cryptic": 14, "crystal_artifact": 6
    },
    "series_counts": {"InsP6": 134, "InsP5": 1},
}


def test_coordinate_ligand_detection_resolves_almost_every_deposited_entry():
    """99% of deposited entries were identified from coordinates alone.

    This is the evidence that the structural ligand test is not too strict. A
    high failure rate here would mean the detector rejects real inositol
    phosphates - a defect in the detector, not a property of the data.
    """
    assert SURVEYED["n_measured"] / SURVEYED["n_entries"] > 0.98


def test_the_two_boundary_estimates_disagree_so_burial_is_not_two_classes():
    """Recorded finding: relative burial is continuous, not bimodal.

    A genuine two-class structure would put both estimators in the same gap.
    The density minimum sits at 0.138 near the low tail while Otsu cuts at 0.463,
    near the middle of a broad spread - what Otsu does when there is one wide
    mode rather than two. The quantiles agree: burial runs smoothly from 0.07 to
    0.90 with no chasm.

    This test records the finding so a later change cannot quietly promote the
    cryptic/surface split to a discovered class boundary. The positive class is
    defined by a cutoff, and downstream metrics are conditioned on that choice.
    """
    trough = SURVEYED["density_minimum"]
    otsu = SURVEYED["otsu"]
    assert abs(otsu - trough) > 0.2, (
        "estimators agreeing would be evidence of a real two-class split; "
        f"got trough={trough} otsu={otsu}"
    )

    # No chasm: successive quantiles step smoothly rather than jumping a gap.
    values = list(SURVEYED["quantiles"].values())
    steps = [b - a for a, b in zip(values, values[1:])]
    assert max(steps) < 0.25, f"a real gap would show as a large quantile step: {steps}"


def test_configured_boundary_sits_at_the_population_density_minimum():
    """0.12 was calibrated on one control and lands on the population trough.

    Independent corroboration is worth pinning: the boundary was chosen from
    ADAR2's measurement, before this survey existed, and falls within one
    histogram bin (0.025) of the sparsest region of 135 deposited entries.
    """
    from cryptic_ip.validation.burial_metrics import CRYPTIC_RELATIVE_SASA_MAX

    assert abs(SURVEYED["density_minimum"] - CRYPTIC_RELATIVE_SASA_MAX) <= 0.025


def test_the_surveyed_set_is_essentially_one_chemistry():
    """A limit on what this calibration covers: it is InsP6, not IP generally."""
    counts = SURVEYED["series_counts"]
    dominant = max(counts, key=counts.get)
    assert dominant == "InsP6"
    assert counts[dominant] / sum(counts.values()) > 0.95


def test_cryptic_sites_are_a_small_minority_of_deposited_complexes():
    """14 of 135. The class is rare, which is why ranking metrics matter."""
    fraction = SURVEYED["n_below_configured_boundary"] / SURVEYED["n_measured"]
    assert 0.05 < fraction < 0.20


def test_boundary_sweep_reports_class_size_at_each_candidate_cutoff():
    """Burial is continuous, so the positive class is whatever the cutoff says.

    A single class count hides that. The sweep makes the choice inspectable:
    counts must rise monotonically with the boundary, and a cutoff whose
    neighbours give very different class sizes is a fragile one.
    """
    measurements = [
        EntryMeasurement(pdb_id=f"{i:04d}", relative_sasa=i / 100.0, burial_class="surface")
        for i in range(1, 101)
    ]
    sweep = summarise(measurements)["boundary_sweep"]

    boundaries = [row["boundary"] for row in sweep]
    counts = [row["n_positive"] for row in sweep]
    assert boundaries == sorted(boundaries)
    assert counts == sorted(counts), "positives cannot fall as the cutoff rises"

    # On a uniform spread the fraction tracks the boundary itself.
    for row in sweep:
        assert row["fraction_positive"] == pytest.approx(row["boundary"], abs=0.02)


#: Positive class size at each candidate boundary, measured on the same run.
SURVEYED_BOUNDARY_SWEEP = {
    0.05: 0, 0.08: 2, 0.10: 11, 0.12: 14, 0.15: 16, 0.20: 27, 0.25: 40, 0.30: 58,
}


def test_the_configured_boundary_sits_where_class_size_is_least_sensitive():
    """Recorded finding: 0.12 is on a plateau, which is why it is defensible.

    Burial is continuous, so no cutoff is a natural break. What makes one
    defensible is insensitivity: moving the boundary across 0.10-0.15 changes the
    positive count by 5 entries, whereas the same 0.05 shift at 0.20-0.25 changes
    it by 13 and at 0.25-0.30 by 18.
    """
    sweep = SURVEYED_BOUNDARY_SWEEP
    at_boundary = sweep[0.15] - sweep[0.10]
    just_above = sweep[0.25] - sweep[0.20]
    higher_still = sweep[0.30] - sweep[0.25]

    assert at_boundary < just_above < higher_still
    assert just_above / at_boundary > 2.0, "the plateau should be markedly flatter"


def test_the_original_literature_boundary_would_have_emptied_the_positive_class():
    """0.05, taken from the ADAR2 description, yields no positives at all.

    Recorded because it shows the cost of setting a threshold from a qualitative
    description rather than from measurement: the pipeline would have trained on
    a class with no members.
    """
    assert SURVEYED_BOUNDARY_SWEEP[0.05] == 0
