"""Tests for MD validation workflow utilities."""

from cryptic_ip.validation.md_validation import (
    OpenMMMDValidationPipeline,
    PocketStabilityThresholds,
)


def test_parse_helpers():
    pipeline = OpenMMMDValidationPipeline(output_dir="/tmp/md_validation_test")

    residues = pipeline._parse_residues("10, 12,15")
    center = pipeline._parse_center("1.0,2.5,3.25")

    assert residues == [10, 12, 15]
    assert center == (1.0, 2.5, 3.25)


def test_classify_stably_buried():
    pipeline = OpenMMMDValidationPipeline(
        thresholds=PocketStabilityThresholds(
            sasa_stably_buried=1.0,
            sasa_exposed=5.0,
            rmsf_stably_buried=0.2,
            rmsf_exposed=0.4,
            waters_stably_buried=1.0,
            waters_exposed=4.0,
        )
    )

    cls = pipeline.classify_pocket_stability(
        {
            "avg_pocket_sasa_nm2": 0.5,
            "avg_pocket_rmsf_nm": 0.15,
            "avg_waters_in_pocket": 0.4,
        }
    )

    assert cls == "stably buried"


def test_classify_exposed():
    pipeline = OpenMMMDValidationPipeline()

    cls = pipeline.classify_pocket_stability(
        {
            "avg_pocket_sasa_nm2": 6.0,
            "avg_pocket_rmsf_nm": 0.6,
            "avg_waters_in_pocket": 7.0,
        }
    )

    assert cls == "exposed"


def test_classify_transiently_accessible():
    pipeline = OpenMMMDValidationPipeline()

    cls = pipeline.classify_pocket_stability(
        {
            "avg_pocket_sasa_nm2": 2.0,
            "avg_pocket_rmsf_nm": 0.3,
            "avg_waters_in_pocket": 2.0,
        }
    )

    assert cls == "transiently accessible"
