"""Regression tests pinning the control calibration to measured values.

The thresholds in :mod:`cryptic_ip.validation.burial_metrics` and
:mod:`cryptic_ip.validation.control_scoring` are calibrated against measurements
of the deposited control structures, produced by ``scripts/calibrate_controls.py``
and recorded below. Pinning them here means a later change to the burial
definition, the scorer, or a threshold cannot silently invalidate the
calibration: the arithmetic that the tier-1 gate performs is reproduced from
fixed numbers, so a regression shows up without needing database access.

These are *recorded measurements*, not expected outputs of the code under test.
If the measurement pipeline changes such that these values no longer describe the
controls, re-run the calibration script and update them here in the same commit
as the change.
"""

import pytest

from cryptic_ip.validation.burial_metrics import (
    CRYPTIC_RELATIVE_SASA_MAX,
    SEMI_CRYPTIC_RELATIVE_SASA_MAX,
    classify_burial_relative,
)
from cryptic_ip.validation.control_scoring import (
    RELATIVE_BURIAL_FULL_CREDIT,
    RELATIVE_BURIAL_ZERO_CREDIT,
    burial_component_relative,
    negative_passed,
    positive_passed,
    validation_score,
)

#: Measured on the deposited structures by scripts/calibrate_controls.py
#: (GitHub Actions run 31241170491, 2026-08-08, 256 SASA sample points/atom).
#: Five controls; ``site_composite`` is the composite at the pocket with the
#: greatest ligand-atom overlap.
MEASURED = {
    "ADAR2": {
        "pdb": "1ZY7",
        "role": "positive",
        "expected_class": "cryptic",
        "relative_sasa": 0.093,
        "relative_phosphate_sasa": 0.089,
        "burial_depth": 5.76,
        "enclosure": 0.941,
        "n_basic_residues": 8,
        "pocket_volume": 1524.71,
        "site_composite": 0.643,
        # Composite of the pocket the old centre-distance rule selected.
        "nearest_centre_composite": 0.432,
    },
    "Pds5B": {
        "pdb": "5HDT",
        "role": "positive",
        "expected_class": "surface",
        "relative_sasa": 0.466,
        "relative_phosphate_sasa": 0.460,
        "burial_depth": 5.63,
        "enclosure": 0.738,
        "n_basic_residues": 8,
        "pocket_volume": 894.59,
        "site_composite": 0.561,
    },
    "HDAC1": {
        "pdb": "5ICN",
        "role": "positive",
        "expected_class": "cryptic",
        "relative_sasa": 0.089,
        "relative_phosphate_sasa": None,
        "burial_depth": 4.10,
        "enclosure": 0.906,
        "n_basic_residues": 3,
        "pocket_volume": 593.04,
        "site_composite": 0.530,
    },
    "PLCd1_PH": {
        "pdb": "1MAI",
        "role": "negative",
        "expected_class": "surface",
        "relative_sasa": 0.373,
        "relative_phosphate_sasa": 0.348,
        "burial_depth": 4.68,
        "enclosure": 0.598,
        "n_basic_residues": 5,
        "pocket_volume": 1532.16,
        "site_composite": 0.415,
    },
    "Btk_PH": {
        "pdb": "1BWN",
        "role": "negative",
        "expected_class": "surface",
        "relative_sasa": 0.253,
        "relative_phosphate_sasa": 0.266,
        "burial_depth": 5.98,
        "enclosure": 0.668,
        "n_basic_residues": 3,
        "pocket_volume": 491.44,
        "site_composite": 0.498,
    },
}

#: Controls whose ligand is genuinely sequestered, and those that are not.
BURIED_CONTROLS = ("ADAR2", "HDAC1")
EXPOSED_CONTROLS = ("Pds5B", "PLCd1_PH", "Btk_PH")


def test_cryptic_boundary_admits_the_paradigm_buried_site():
    """ADAR2 at 0.093 must fall inside the cryptic boundary, with margin."""
    adar2 = MEASURED["ADAR2"]["relative_sasa"]
    assert adar2 <= CRYPTIC_RELATIVE_SASA_MAX
    # At least 20 % headroom, so ordinary measurement variation cannot flip it.
    assert adar2 <= 0.8 * CRYPTIC_RELATIVE_SASA_MAX


def test_cryptic_boundary_excludes_the_paradigm_surface_site():
    """PLC-delta-1 at 0.373 must be well outside both buried classes."""
    plc = MEASURED["PLCd1_PH"]["relative_sasa"]
    assert plc > SEMI_CRYPTIC_RELATIVE_SASA_MAX
    assert classify_burial_relative(plc, relative_phosphate_sasa=0.348) == "surface"


def test_controls_are_classified_as_measured():
    assert (
        classify_burial_relative(
            MEASURED["ADAR2"]["relative_sasa"],
            relative_phosphate_sasa=MEASURED["ADAR2"]["relative_phosphate_sasa"],
        )
        == "cryptic"
    )


def test_burial_boundary_sits_in_the_gap_between_buried_and_exposed_controls():
    """The boundary must fall in the gap the panel leaves, not outside the range.

    Across five controls the buried ligands measure 0.089-0.093 and the exposed
    ones 0.253-0.466, leaving a clear gap. The boundary sits inside it.
    """
    buried = [MEASURED[name]["relative_sasa"] for name in BURIED_CONTROLS]
    exposed = [MEASURED[name]["relative_sasa"] for name in EXPOSED_CONTROLS]
    assert max(buried) < CRYPTIC_RELATIVE_SASA_MAX < min(exposed)


def test_every_control_classifies_as_the_panel_measured_it():
    for name, control in MEASURED.items():
        assert (
            classify_burial_relative(
                control["relative_sasa"],
                relative_phosphate_sasa=control["relative_phosphate_sasa"],
            )
            == control["expected_class"]
        ), name


def test_burial_credit_separates_buried_from_exposed_across_the_panel():
    for name in BURIED_CONTROLS:
        assert burial_component_relative(
            MEASURED[name]["relative_sasa"]
        ) == pytest.approx(1.0), name
    for name in EXPOSED_CONTROLS:
        assert burial_component_relative(MEASURED[name]["relative_sasa"]) < 0.6, name
    assert RELATIVE_BURIAL_FULL_CREDIT < RELATIVE_BURIAL_ZERO_CREDIT


def test_volume_window_covers_every_real_ligand_site():
    """The window applies to cavity volume, so it must span the observed sites.

    Under the original 300-800 A^3 window - the volume of the *ligand* rather
    than the cavity - ADAR2's site scored 0.16 while the Btk surface negative
    scored 1.00, so the component penalised the paradigm positive and rewarded a
    negative.
    """
    from cryptic_ip.analysis.scorer import PocketScorer

    scorer = PocketScorer()
    for name, control in MEASURED.items():
        assert scorer.score_volume(control["pocket_volume"]) > 0.9, name

    # The defect it replaces: the old window inverted the ranking.
    old_window = PocketScorer(
        parameters=type(scorer.parameters)(volume_optimum_high=800.0)
    )
    assert old_window.score_volume(MEASURED["ADAR2"]["pocket_volume"]) < old_window.score_volume(
        MEASURED["Btk_PH"]["pocket_volume"]
    )


def test_tier1_gate_passes_on_the_measured_values():
    """Reproduce the gate's arithmetic from the recorded measurements."""
    adar2, plc = MEASURED["ADAR2"], MEASURED["PLCd1_PH"]

    positive = validation_score(
        "positive",
        adar2["site_composite"],
        None,
        adar2["n_basic_residues"],
        relative_sasa=adar2["relative_sasa"],
    )
    negative = validation_score(
        "negative",
        plc["site_composite"],
        None,
        plc["n_basic_residues"],
        relative_sasa=plc["relative_sasa"],
    )

    assert positive_passed(
        positive,
        None,
        adar2["n_basic_residues"],
        adar2["site_composite"],
        burial_class="cryptic",
        relative_sasa=adar2["relative_sasa"],
    )
    assert negative_passed(
        negative, None, plc["site_composite"], relative_sasa=plc["relative_sasa"]
    )
    # The tier-1 gate requires separation above 0.50.
    assert positive - negative > 0.50


def test_gate_would_fail_on_the_pocket_the_old_rule_selected():
    """Selecting by centre distance picks a pocket that cannot clear the gate.

    This is why the control validators now select by ligand-atom overlap: the
    nearest-centre pocket in ADAR2 scores 0.432, below the 0.50 composite the
    positive-control criterion requires, even though the pocket that actually
    holds the InsP6 scores 0.644.
    """
    adar2 = MEASURED["ADAR2"]
    wrong = validation_score(
        "positive",
        adar2["nearest_centre_composite"],
        None,
        adar2["n_basic_residues"],
        relative_sasa=adar2["relative_sasa"],
    )
    assert not positive_passed(
        wrong,
        None,
        adar2["n_basic_residues"],
        adar2["nearest_centre_composite"],
        burial_class="cryptic",
        relative_sasa=adar2["relative_sasa"],
    )


def test_depth_does_not_separate_the_control_panel():
    """Recorded finding: burial depth is uninformative on real structures.

    On idealised synthetic spheres depth separates buried from surface sites
    cleanly (22 A against 5 A). On the five deposited controls the values are
    completely interleaved - and the *largest* depth in the panel belongs to a
    surface negative:

        ADAR2 5.76 (buried)   Pds5B 5.63   HDAC1 4.10 (buried)
        PLCd1 4.68 (surface)  Btk   5.98 (surface)

    Depth is a minimum over all solvent-exposed atoms, and a real protein
    surface is irregular enough that some exposed atom lies within a few
    Angstrom of almost any interior point, so the minimum saturates. An
    idealised sphere has no such irregularity, which is exactly why the
    synthetic benchmark could not reveal this.

    This test records the finding so a future change cannot quietly promote
    depth to a primary burial criterion.
    """
    buried_depths = [MEASURED[name]["burial_depth"] for name in BURIED_CONTROLS]
    exposed_depths = [MEASURED[name]["burial_depth"] for name in EXPOSED_CONTROLS]
    # No threshold on depth can separate the two groups.
    assert max(exposed_depths) > max(buried_depths)
    assert min(exposed_depths) > min(buried_depths)


def test_enclosure_separates_the_control_panel():
    """Enclosure, unlike depth, does separate buried from exposed controls."""
    buried = [MEASURED[name]["enclosure"] for name in BURIED_CONTROLS]
    exposed = [MEASURED[name]["enclosure"] for name in EXPOSED_CONTROLS]
    assert min(buried) > max(exposed), (
        "enclosure should separate the panel: "
        f"buried {sorted(buried)} vs exposed {sorted(exposed)}"
    )


def test_calibration_measurements_are_strictly_valid_json():
    """Undefined ratios must serialise as null, not as a bare NaN token.

    A phosphate SASA ratio is undefined when the matched ligand carries no
    phosphate, and Python represents that as NaN. ``json.dumps`` emits NaN as a
    bare ``NaN`` token by default, which is not part of the JSON grammar: the
    measurements file round-tripped through Python but was rejected by jq,
    JavaScript and R. Since that file is a published artifact of the calibration
    run, it has to be readable by a strict parser.
    """
    import json

    from scripts.calibrate_controls import _json_safe

    payload = _json_safe(
        {
            "relative_sasa": 0.093,
            "relative_phosphate_sasa": float("nan"),
            "depth": float("inf"),
            "instances": [{"enclosure": float("nan")}],
            "n": 3,
        }
    )
    # allow_nan=False makes the encoder raise rather than emit an invalid token.
    text = json.dumps(payload, allow_nan=False)
    reloaded = json.loads(text)

    assert reloaded["relative_phosphate_sasa"] is None
    assert reloaded["depth"] is None
    assert reloaded["instances"][0]["enclosure"] is None
    assert reloaded["relative_sasa"] == pytest.approx(0.093)
    assert reloaded["n"] == 3
