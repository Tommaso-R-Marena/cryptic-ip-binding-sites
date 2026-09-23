"""Physical units for pocket-analysis features.

Central registry so that units are attached consistently wherever pocket
descriptors are shown to a user (CLI tables, the Streamlit app) or exported as
metadata. Machine-readable column names in persisted CSV/JSON outputs are left
unchanged so downstream parsing and schema validation keep working; units are
provided alongside as labels or as a separate mapping.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional

# Unicode symbols so the analysis reports proper structural-biology units.
ANGSTROM = "\u00c5"  # Å
ANGSTROM_SQ = "\u00c5\u00b2"  # Å²
ANGSTROM_CU = "\u00c5\u00b3"  # Å³

#: Canonical unit label for each pocket-analysis feature/column. Empty string
#: means the column is an identifier/label with no physical unit.
FEATURE_UNITS: Dict[str, str] = {
    "pocket_id": "",
    "rank": "",
    "classification": "",
    "uniprot_id": "",
    "protein_file": "",
    "composite_score": "0-1",
    "volume": ANGSTROM_CU,
    # fpocket local hydrophobic density proxy (dimensionless), distinct from the
    # geometric burial depth below.
    "depth": "fpocket density",
    "burial_depth": ANGSTROM,
    "sasa": ANGSTROM_SQ,
    "basic_residues": "count",
    "residue_count": "count",
    "electrostatic_potential": "kT/e",
    "plddt_confidence": "pLDDT 0-100",
}


def unit_for(column: str) -> str:
    """Return the unit label for a column ("" when it has no physical unit)."""
    return FEATURE_UNITS.get(column, "")


def label_with_unit(column: str) -> str:
    """Return a display header like ``"burial_depth (Å)"`` (unitless columns unchanged)."""
    unit = unit_for(column)
    return f"{column} ({unit})" if unit else column


def label_columns(columns: Iterable[str]) -> Dict[str, str]:
    """Map each column name to its unit-annotated display header.

    Suitable for ``DataFrame.rename(columns=...)`` when rendering tables for
    users without mutating the underlying machine-readable column names.
    """
    return {column: label_with_unit(column) for column in columns}


def units_mapping(columns: Optional[Iterable[str]] = None) -> Dict[str, str]:
    """Return a ``{column: unit}`` mapping, excluding columns without a unit.

    When ``columns`` is omitted the full registry is described. Intended for
    embedding as machine-readable metadata in exported results.
    """
    selected: List[str] = list(columns) if columns is not None else list(FEATURE_UNITS)
    return {column: unit_for(column) for column in selected if unit_for(column)}
