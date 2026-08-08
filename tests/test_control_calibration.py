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
#: (GitHub Actions run 31252064686, 2026-08-08, 256 SASA sample points/atom).
#: Composites and site volumes were refreshed after the volume window was
#: corrected, so every value here comes from one run of one scorer.
#: Five controls; ``site_composite`` is the composite at the pocket with the
#: greatest ligand-atom overlap.
MEASURED = {
    "ADAR2": {
        "comp_id": "IHP",
        "pdb": "1ZY7",
        "role": "positive",
        "expected_class": "cryptic",
        "relative_sasa": 0.093,
        "relative_phosphate_sasa": 0.089,
        "burial_depth": 5.76,
        "enclosure": 0.941,
        "n_basic_residues": 8,
        "pocket_volume": 1487.53,
        "site_composite": 0.727,
        # Composite of the pocket the old centre-distance rule selected,
        # measured under the scorer of the time (run 31241170491). Retained as
        # a historical demonstration that centre-distance selection grades a
        # pocket that does not hold the ligand, not as a current measurement.
        "nearest_centre_composite": 0.432,
    },
    "Pds5B": {
        "comp_id": "IHP",
        "pdb": "5HDT",
        "role": "positive",
        "expected_class": "surface",
        "relative_sasa": 0.466,
        "relative_phosphate_sasa": 0.460,
        "burial_depth": 5.631,
        "enclosure": 0.738,
        "n_basic_residues": 8,
        "pocket_volume": 908.96,
        "site_composite": 0.582,
    },
    # Corrected. The earlier entry recorded 6A0 at relative SASA 0.089 and
    # classified HDAC1 as cryptic, but 6A0 carries no phosphate: it is not an
    # inositol phosphate, and burial had been measured on it because the old
    # identifier whitelist admitted it and it happened to be the most buried
    # matching copy. Identifying ligands from coordinates excludes it, and the
    # most buried *phosphorylated* copy in 5ICN is a solvent-exposed InsP6.
    # HDAC1 is therefore an exposed control, not a buried one.
    #
    # Its depth, basic count, site volume and composite are measured at the
    # InsP6 site, not carried over from the 6A0 measurement they replaced.
    "HDAC1": {
        "comp_id": "IHP",
        "pdb": "5ICN",
        "role": "positive",
        "expected_class": "surface",
        "relative_sasa": 0.436,
        "relative_phosphate_sasa": 0.444,
        "burial_depth": 4.392,
        "enclosure": 0.629,
        "n_basic_residues": 4,
        "pocket_volume": 939.34,
        "site_composite": 0.677,
    },
    "PLCd1_PH": {
        "comp_id": "I3P",
        "pdb": "1MAI",
        "role": "negative",
        "expected_class": "surface",
        "relative_sasa": 0.373,
        "relative_phosphate_sasa": 0.348,
        "burial_depth": 4.68,
        "enclosure": 0.598,
        "n_basic_residues": 5,
        "pocket_volume": 1533.9,
        "site_composite": 0.499,
    },
    "Btk_PH": {
        "comp_id": "4IP",
        "pdb": "1BWN",
        "role": "negative",
        "expected_class": "surface",
        "relative_sasa": 0.253,
        "relative_phosphate_sasa": 0.266,
        "burial_depth": 5.977,
        "enclosure": 0.668,
        "n_basic_residues": 3,
        "pocket_volume": 499.97,
        "site_composite": 0.498,
    },
}

#: Which chemical component each control's burial was measured on. Recorded
#: because the panel is only meaningful if it measured the right molecule: the
#: identifier alone does not say whether a component is an inositol *phosphate*,
#: which is why the pipeline now identifies ligands from coordinates and reports
#: the phosphate count. HDAC1 (5ICN) matched 6A0 and reports an undefined
#: phosphate SASA, meaning no phosphate group was resolved on it; the calibration
#: digest now prints the detected series so this can be read off directly rather
#: than inferred.
#:
#: Controls whose ligand is genuinely sequestered, and those that are not.
BURIED_CONTROLS = ("ADAR2",)
EXPOSED_CONTROLS = ("Pds5B", "HDAC1", "PLCd1_PH", "Btk_PH")


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

    The gap is real but the panel is now thin on the buried side. Correcting
    HDAC1 - whose burial had been measured on a component carrying no phosphate -
    moved it from the buried group to the exposed one, leaving ADAR2 at 0.093 as
    the *only* sequestered control against four exposed ones at 0.253-0.466.

    The boundary of 0.12 still sits inside the gap, but it now rests on a single
    structure. Widening the buried side is the most valuable thing that can be
    done to this calibration; scripts/burial_survey.py measures the whole
    deposited set for that purpose.
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
    with_volume = {
        name: control
        for name, control in MEASURED.items()
        if control.get("pocket_volume") is not None
    }
    # Guard against the test quietly emptying itself if entries lose the field.
    assert len(with_volume) >= 3, f"too few recorded volumes: {sorted(with_volume)}"
    for name, control in with_volume.items():
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

        ADAR2 5.76 (buried)   Pds5B 5.63 (surface)
        PLCd1 4.68 (surface)  Btk   5.98 (surface)

    Depth is a minimum over all solvent-exposed atoms, and a real protein
    surface is irregular enough that some exposed atom lies within a few
    Angstrom of almost any interior point, so the minimum saturates. An
    idealised sphere has no such irregularity, which is exactly why the
    synthetic benchmark could not reveal this.

    This test records the finding so a future change cannot quietly promote
    depth to a primary burial criterion.
    """
    buried_depths = [
        MEASURED[name]["burial_depth"]
        for name in BURIED_CONTROLS
        if MEASURED[name].get("burial_depth") is not None
    ]
    exposed_depths = [
        MEASURED[name]["burial_depth"]
        for name in EXPOSED_CONTROLS
        if MEASURED[name].get("burial_depth") is not None
    ]
    assert buried_depths and len(exposed_depths) >= 2

    # No threshold on depth can isolate the buried control: its depth falls
    # strictly inside the spread of the exposed ones, so any cutoff that admits
    # it also admits a surface site.
    assert min(exposed_depths) < min(buried_depths)
    assert max(exposed_depths) > max(buried_depths)


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
