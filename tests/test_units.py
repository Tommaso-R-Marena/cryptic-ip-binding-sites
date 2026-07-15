"""Tests for the pocket-analysis units registry."""

import pandas as pd

from cryptic_ip.analysis import units


def test_known_feature_units():
    assert units.unit_for("burial_depth") == units.ANGSTROM
    assert units.unit_for("volume") == units.ANGSTROM_CU
    assert units.unit_for("sasa") == units.ANGSTROM_SQ
    assert units.unit_for("electrostatic_potential") == "kT/e"
    assert units.unit_for("basic_residues") == "count"


def test_unknown_column_has_no_unit():
    assert units.unit_for("does_not_exist") == ""


def test_label_with_unit_appends_unit_only_when_present():
    assert units.label_with_unit("burial_depth") == f"burial_depth ({units.ANGSTROM})"
    assert units.label_with_unit("volume") == f"volume ({units.ANGSTROM_CU})"
    # Identifier-style columns are returned unchanged.
    assert units.label_with_unit("pocket_id") == "pocket_id"
    assert units.label_with_unit("classification") == "classification"


def test_label_columns_maps_each_column():
    cols = ["pocket_id", "burial_depth", "sasa"]
    mapping = units.label_columns(cols)
    assert mapping == {
        "pocket_id": "pocket_id",
        "burial_depth": f"burial_depth ({units.ANGSTROM})",
        "sasa": f"sasa ({units.ANGSTROM_SQ})",
    }


def test_label_columns_usable_for_dataframe_rename():
    frame = pd.DataFrame({"pocket_id": [1], "burial_depth": [18.0], "sasa": [3.0]})
    renamed = frame.rename(columns=units.label_columns(frame.columns))
    assert f"burial_depth ({units.ANGSTROM})" in renamed.columns
    assert f"sasa ({units.ANGSTROM_SQ})" in renamed.columns
    # Values are untouched by relabeling.
    assert renamed.iloc[0][f"burial_depth ({units.ANGSTROM})"] == 18.0


def test_units_mapping_excludes_unitless_columns():
    mapping = units.units_mapping(["pocket_id", "burial_depth", "electrostatic_potential"])
    assert mapping == {
        "burial_depth": units.ANGSTROM,
        "electrostatic_potential": "kT/e",
    }


def test_units_mapping_defaults_to_full_registry():
    mapping = units.units_mapping()
    assert "burial_depth" in mapping and "volume" in mapping
    # Unitless identifiers are excluded.
    assert "pocket_id" not in mapping
