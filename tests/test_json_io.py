"""Tests for strict JSON output of published artifacts."""

import json
from pathlib import Path

import numpy as np
import pytest

from cryptic_ip.utils.json_io import dumps_strict, json_safe, write_json_strict


def test_undefined_values_become_null_not_bare_nan():
    text = dumps_strict({"auroc": float("nan"), "depth": float("inf"), "ok": 0.5})
    assert "NaN" not in text and "Infinity" not in text
    assert json.loads(text) == {"auroc": None, "depth": None, "ok": 0.5}


def test_numpy_types_and_paths_are_converted():
    payload = json_safe(
        {"n": np.int64(3), "x": np.float32(0.25), "arr": np.array([1.0, np.nan]),
         "path": Path("models/m.pkl"), "flag": np.bool_(True)}
    )
    assert payload == {
        "n": 3, "x": 0.25, "arr": [1.0, None], "path": "models/m.pkl", "flag": True
    }


def test_nesting_is_handled_throughout():
    payload = json_safe({"a": [{"b": (float("nan"), 1.0)}]})
    assert payload == {"a": [{"b": [None, 1.0]}]}


def test_written_file_is_readable_by_a_strict_parser(tmp_path):
    path = write_json_strict(tmp_path / "meta.json", {"v": float("nan")})
    # parse_constant is only invoked for NaN/Infinity tokens; it must never fire.
    json.loads(path.read_text(), parse_constant=pytest.fail)
