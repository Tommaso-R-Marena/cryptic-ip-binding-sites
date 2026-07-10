"""Tests for advanced burial metrics."""

from pathlib import Path

import pytest

from cryptic_ip.validation.burial_metrics import classify_burial, compute_burial_metrics
from cryptic_ip.validation.control_scoring import (
    positive_passed,
    validation_score,
)


@pytest.fixture
def adar2_pdb(adar2_structure_path):
    return adar2_structure_path


def test_classify_burial_cryptic():
    assert classify_burial(0.0) == "cryptic"
    assert classify_burial(30.0) == "semi_cryptic"
    assert classify_burial(80.0) == "surface"


def test_adar2_burial_metrics(adar2_pdb):
    metrics = compute_burial_metrics(adar2_pdb)
    assert metrics.burial_class == "cryptic"
    assert metrics.ligand_sasa is not None
    assert metrics.ligand_sasa < 10.0
    assert metrics.phosphate_sasa is not None
    assert metrics.delta_sasa is not None


def test_semi_cryptic_positive_criteria():
    val = validation_score("positive", 0.75, 36.0, 8)
    assert positive_passed(val, 36.0, 8, 0.75, burial_class="semi_cryptic")
