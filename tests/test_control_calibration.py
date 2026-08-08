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
#: (GitHub Actions run 31238630563, 2026-08-08, 256 SASA sample points/atom).
MEASURED = {
    "ADAR2": {
        "pdb": "1ZY7",
        "role": "positive",
        "relative_sasa": 0.093,
        "relative_phosphate_sasa": 0.089,
        "burial_depth": 5.76,
        "enclosure": 0.941,
        "n_basic_residues": 8,
        # Composite at the pocket with 100 % ligand-atom overlap.
        "site_composite": 0.644,
        # Composite of the pocket the old centre-distance rule selected.
        "nearest_centre_composite": 0.432,
    },
    "PLCd1_PH": {
        "pdb": "1MAI",
        "role": "negative",
        "relative_sasa": 0.373,
        "relative_phosphate_sasa": 0.348,
        "burial_depth": 4.68,
        "enclosure": 0.598,
        "n_basic_residues": 5,
        "site_composite": 0.415,
    },
}


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


def test_burial_boundary_sits_between_the_two_controls():
    """The boundary must separate the controls, not sit outside their range."""
    assert (
        MEASURED["ADAR2"]["relative_sasa"]
        < CRYPTIC_RELATIVE_SASA_MAX
        < MEASURED["PLCd1_PH"]["relative_sasa"]
    )


def test_burial_credit_separates_the_controls():
    assert burial_component_relative(MEASURED["ADAR2"]["relative_sasa"]) == pytest.approx(1.0)
    assert burial_component_relative(MEASURED["PLCd1_PH"]["relative_sasa"]) < 0.2
    assert RELATIVE_BURIAL_FULL_CREDIT < RELATIVE_BURIAL_ZERO_CREDIT


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


def test_enclosure_separates_the_controls_more_cleanly_than_depth():
    """Recorded finding: burial depth barely separates real controls.

    On idealised synthetic spheres, depth cleanly separates buried from surface
    sites. On the deposited structures it does not - 5.76 A versus 4.68 A -
    because a real protein surface is irregular enough that some exposed atom is
    almost always within a few Angstrom of any interior point, and depth is a
    minimum over all exposed atoms. Enclosure, which integrates over directions
    rather than taking a minimum, separates the same two structures cleanly.
    This test records that finding so a future change cannot quietly rely on
    depth as the primary burial discriminator.
    """
    depth_gap = MEASURED["ADAR2"]["burial_depth"] - MEASURED["PLCd1_PH"]["burial_depth"]
    enclosure_gap = MEASURED["ADAR2"]["enclosure"] - MEASURED["PLCd1_PH"]["enclosure"]

    # Relative separation, since the two quantities have different units.
    depth_relative = depth_gap / MEASURED["PLCd1_PH"]["burial_depth"]
    enclosure_relative = enclosure_gap / MEASURED["PLCd1_PH"]["enclosure"]
    assert enclosure_relative > depth_relative
    assert depth_relative < 0.30
