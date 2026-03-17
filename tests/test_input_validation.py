from pathlib import Path

import pytest

from cryptic_ip.errors import UnsupportedFormatError, ValidationError
from cryptic_ip.utils.input_validation import parse_or_raise, validate_structure_file, ExportConfig


def test_validate_structure_file_missing(tmp_path):
    with pytest.raises(ValidationError):
        validate_structure_file(tmp_path / "missing.pdb")


def test_validate_structure_file_format(tmp_path):
    bad = tmp_path / "input.txt"
    bad.write_text("x")
    with pytest.raises(UnsupportedFormatError):
        validate_structure_file(bad)


def test_parse_or_raise_export_config():
    model = parse_or_raise(ExportConfig, {"export_format": "csv"}, "bad")
    assert model.export_format == "csv"
    with pytest.raises(ValidationError):
        parse_or_raise(ExportConfig, {"export_format": "xml"}, "bad")
