"""Strictly valid JSON output for published artifacts.

Undefined quantities are routine in this pipeline - a phosphate SASA ratio for
a ligand with no phosphate, an AUROC for a fold containing one class - and
Python evaluates them to ``NaN``. ``json.dumps`` writes ``NaN`` as a bare token
by default, which is not part of the JSON grammar: Python reads its own output
back happily, so nothing complains, but jq, JavaScript and R all reject the
file. Every artifact other tools are meant to read goes through here instead.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Union

import numpy as np


def json_safe(value: Any) -> Any:
    """Convert a nested structure into one that serialises as strict JSON.

    Non-finite floats become ``None`` (``null``), which is both valid and the
    correct meaning for an undefined value. NumPy scalars and arrays become
    built-in types, and paths become strings.

    Args:
        value: Arbitrary nested structure.

    Returns:
        An equivalent structure containing only JSON-native values.
    """
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def dumps_strict(value: Any, **kwargs: Any) -> str:
    """Serialise to JSON, refusing to emit anything outside the grammar.

    Args:
        value: Structure to serialise.
        **kwargs: Passed through to :func:`json.dumps`.

    Returns:
        The JSON text.
    """
    kwargs.setdefault("indent", 2)
    return json.dumps(json_safe(value), allow_nan=False, **kwargs)


def write_json_strict(path: Union[str, Path], value: Any, **kwargs: Any) -> Path:
    """Write ``value`` to ``path`` as strictly valid JSON.

    Args:
        path: Destination file.
        value: Structure to serialise.
        **kwargs: Passed through to :func:`json.dumps`.

    Returns:
        The path written.
    """
    target = Path(path)
    target.write_text(dumps_strict(value, **kwargs), encoding="utf-8")
    return target
