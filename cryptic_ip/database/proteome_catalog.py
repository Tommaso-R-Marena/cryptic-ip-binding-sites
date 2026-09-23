"""Phase 2: catalogue an AlphaFold proteome and check it is screen-ready.

A corrupted file or a naming inconsistency causes *silent* failures at proteome
scale - a structure that cannot be parsed simply contributes no pockets, and the
hit rate is quietly computed over fewer proteins than it claims. This module
turns the Phase 2 quality-control checklist into measurements:

1. File count against the AlphaFold DB's expected count.
2. No zero-byte or truncated files.
3. Consistent naming: ``AF-{UniProtID}-F{n}-model_v{k}``.
4. Each model is geometrically a protein (replaces the manual PyMOL spot-check;
   see :func:`ca_geometry_fraction`).
5. pLDDT extracted and stored per model.
6. A master CSV linking UniProt ID, file, organism, length and mean pLDDT.

Parsing is plain text over the fixed-column PDB format rather than a full
structure parse, because a catalogue of ~40,000 models has to be cheap to build.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

#: Proteomes screened, with the protein counts the project plan expects from the
#: AlphaFold DB. The counts are for AlphaFold DB v4; a later release may differ,
#: which the QC report states rather than treating as an error.
PROTEOMES: Dict[str, Dict[str, object]] = {
    "yeast": {
        "proteome_id": "UP000002311",
        "taxon": "559292",
        "organism": "Saccharomyces cerevisiae",
        "expected_proteins": 6049,
        "ip6_uM": 20.0,
    },
    "human": {
        "proteome_id": "UP000005640",
        "taxon": "9606",
        "organism": "Homo sapiens",
        "expected_proteins": 23391,
        "ip6_uM": 25.0,
    },
    "dictyostelium": {
        "proteome_id": "UP000002195",
        "taxon": "44689",
        "organism": "Dictyostelium discoideum",
        "expected_proteins": 12622,
        "ip6_uM": 520.0,
    },
}

AF_NAME = re.compile(
    r"^AF-(?P<uniprot_id>[A-Z0-9]+(?:-\d+)?)-F(?P<fragment>\d+)-model_v(?P<version>\d+)"
    r"\.(?P<fmt>pdb|cif)(?P<gz>\.gz)?$"
)

#: Smallest plausible model file. The project plan's "> 1 KB" rule.
MIN_FILE_BYTES = 1024

#: Consecutive C-alpha atoms in a real polypeptide sit 3.8 A apart (trans
#: peptide). A model whose consecutive C-alpha distances fall outside this band
#: is not a protein chain - a corrupt coordinate block, a unit error, or a
#: file that is not what its name says.
CA_CA_RANGE = (3.6, 4.0)
#: Fraction of consecutive pairs that must fall in :data:`CA_CA_RANGE`. Cis
#: prolines (~2.9 A) are the only legitimate exception and are rare.
MIN_CA_GEOMETRY_FRACTION = 0.95

PLDDT_CUTOFF = 70.0


@dataclass
class ModelRecord:
    """One AlphaFold model file, as measured."""

    uniprot_id: str
    fragment: int
    version: int
    filename: str
    path: str
    n_bytes: int
    name_ok: bool
    length: int = 0
    mean_plddt: float = float("nan")
    fraction_plddt_70: float = float("nan")
    ca_geometry_fraction: float = float("nan")
    complete: bool = False
    error: str = ""

    @property
    def qc_pass(self) -> bool:
        return (
            self.name_ok
            and self.complete
            and self.n_bytes >= MIN_FILE_BYTES
            and self.length > 0
            and self.ca_geometry_fraction >= MIN_CA_GEOMETRY_FRACTION
        )


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def ca_geometry_fraction(ca_coords: np.ndarray, resseqs: np.ndarray, chains: np.ndarray) -> float:
    """Fraction of sequence-adjacent C-alpha pairs at a peptide-bond distance.

    Only pairs adjacent in sequence *and* chain are compared, so a chain break
    or a numbering gap is not mistaken for bad geometry.
    """
    if len(ca_coords) < 2:
        return float("nan")
    adjacent = (np.diff(resseqs) == 1) & (chains[1:] == chains[:-1])
    if not adjacent.any():
        return float("nan")
    distances = np.linalg.norm(np.diff(ca_coords, axis=0), axis=1)[adjacent]
    low, high = CA_CA_RANGE
    return float(np.mean((distances >= low) & (distances <= high)))


def measure_model(path: Path) -> ModelRecord:
    """Measure one model file without a full structure parse."""
    path = Path(path)
    match = AF_NAME.match(path.name)
    record = ModelRecord(
        uniprot_id=match.group("uniprot_id") if match else path.name.split("-")[1] if path.name.count("-") >= 2 else path.stem,
        fragment=int(match.group("fragment")) if match else 0,
        version=int(match.group("version")) if match else 0,
        filename=path.name,
        path=str(path),
        n_bytes=path.stat().st_size if path.exists() else 0,
        name_ok=bool(match),
    )
    if record.n_bytes == 0:
        record.error = "zero-byte file"
        return record
    if match and match.group("fmt") != "pdb":
        record.error = "not a PDB-format model"
        return record

    coords: List[List[float]] = []
    resseqs: List[int] = []
    chains: List[str] = []
    plddt: List[float] = []
    ended = False
    try:
        with _open_text(path) as handle:
            for line in handle:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    resseqs.append(int(line[22:26]))
                    chains.append(line[21])
                    plddt.append(float(line[60:66]))
                elif line.startswith("END"):
                    ended = True
    except (OSError, EOFError, ValueError) as exc:
        # A truncated gzip stream raises EOFError; a cut-off coordinate line
        # raises ValueError. Either way the file is not usable as it stands.
        record.error = f"unreadable: {type(exc).__name__}: {exc}"
        return record

    record.complete = ended
    if not ended:
        record.error = "truncated: no END record"
    record.length = len(coords)
    if coords:
        arr = np.asarray(plddt, dtype=float)
        record.mean_plddt = float(arr.mean())
        record.fraction_plddt_70 = float(np.mean(arr >= PLDDT_CUTOFF))
        record.ca_geometry_fraction = ca_geometry_fraction(
            np.asarray(coords, dtype=float), np.asarray(resseqs), np.asarray(chains)
        )
    elif not record.error:
        record.error = "no C-alpha atoms"
    return record


def build_catalog(
    paths: Iterable[Path],
    organism_key: str,
) -> pd.DataFrame:
    """Build the master catalogue for one proteome."""
    info = PROTEOMES[organism_key]
    rows = []
    for path in sorted(paths):
        record = measure_model(Path(path))
        row = asdict(record)
        row["qc_pass"] = record.qc_pass
        row["organism_key"] = organism_key
        row["organism"] = info["organism"]
        row["proteome_id"] = info["proteome_id"]
        rows.append(row)
    columns = [
        "uniprot_id", "fragment", "version", "organism_key", "organism", "proteome_id",
        "filename", "path", "n_bytes", "length", "mean_plddt", "fraction_plddt_70",
        "ca_geometry_fraction", "name_ok", "complete", "qc_pass", "error",
    ]
    return pd.DataFrame(rows, columns=columns)


def qc_report(catalog: pd.DataFrame, organism_key: str) -> Dict[str, object]:
    """Evaluate the Phase 2 checklist on a catalogue.

    Each item records the measurement and whether it passes. A count mismatch
    against the plan's expected number is reported, not raised: the expected
    counts are for AlphaFold DB v4, and a later release legitimately differs.
    """
    info = PROTEOMES[organism_key]
    expected = int(info["expected_proteins"])
    proteins = catalog["uniprot_id"].nunique() if not catalog.empty else 0
    f1 = catalog[catalog["fragment"] == 1] if not catalog.empty else catalog
    multi = (
        catalog.groupby("uniprot_id")["fragment"].max().gt(1).sum() if not catalog.empty else 0
    )
    zero = int((catalog["n_bytes"] == 0).sum()) if not catalog.empty else 0
    small = int((catalog["n_bytes"] < MIN_FILE_BYTES).sum()) if not catalog.empty else 0
    truncated = int((~catalog["complete"]).sum()) if not catalog.empty else 0
    bad_names = int((~catalog["name_ok"]).sum()) if not catalog.empty else 0
    geometry_fail = (
        int((catalog["ca_geometry_fraction"].fillna(0) < MIN_CA_GEOMETRY_FRACTION).sum())
        if not catalog.empty
        else 0
    )
    plddt_missing = int(catalog["mean_plddt"].isna().sum()) if not catalog.empty else 0
    versions = sorted(int(v) for v in catalog["version"].unique()) if not catalog.empty else []

    def ratio(a: float, b: float) -> Optional[float]:
        return float(a / b) if b else None

    return {
        "organism_key": organism_key,
        "organism": info["organism"],
        "proteome_id": info["proteome_id"],
        "alphafold_versions": versions,
        "checklist": {
            "file_count": {
                "expected_proteins": expected,
                "observed_proteins": int(proteins),
                "observed_files": int(len(catalog)),
                "f1_models": int(len(f1)),
                "multi_fragment_proteins": int(multi),
                "ratio_to_expected": ratio(proteins, expected),
                # Within 2 %: releases add and retire a small number of entries.
                "passes": bool(expected and abs(proteins - expected) / expected <= 0.02),
            },
            "no_zero_byte_or_truncated": {
                "zero_byte": zero,
                "under_1kb": small,
                "truncated": truncated,
                "passes": zero == 0 and small == 0 and truncated == 0,
            },
            "consistent_naming": {"nonconforming": bad_names, "passes": bad_names == 0},
            "geometry_is_protein": {
                "failing_models": geometry_fail,
                "criterion": (
                    f">= {MIN_CA_GEOMETRY_FRACTION:.0%} of sequence-adjacent C-alpha pairs "
                    f"within {CA_CA_RANGE[0]}-{CA_CA_RANGE[1]} A"
                ),
                "passes": geometry_fail == 0,
            },
            "plddt_extracted": {"missing": plddt_missing, "passes": plddt_missing == 0},
        },
        "screenable_f1_models": int(f1["qc_pass"].sum()) if len(f1) else 0,
        "length": _describe(catalog["length"]) if not catalog.empty else {},
        "mean_plddt": _describe(catalog["mean_plddt"]) if not catalog.empty else {},
        "fraction_models_mean_plddt_ge_70": (
            float((catalog["mean_plddt"] >= PLDDT_CUTOFF).mean()) if not catalog.empty else None
        ),
    }


def _describe(series: pd.Series) -> Dict[str, float]:
    values = series.dropna().astype(float)
    if values.empty:
        return {}
    return {
        "median": float(values.median()),
        "mean": float(values.mean()),
        "q10": float(values.quantile(0.10)),
        "q90": float(values.quantile(0.90)),
        "max": float(values.max()),
    }


def shard(catalog: pd.DataFrame, index: int, count: int) -> pd.DataFrame:
    """Return the screenable F1 models assigned to shard ``index`` of ``count``.

    Models are dealt round-robin in descending length, so every shard receives
    a similar mix of long and short proteins and the shards finish together -
    runtime grows faster than linearly with length, so contiguous slices of a
    length-sorted list would leave one shard with all the giants.
    """
    if not 0 <= index < count:
        raise ValueError(f"shard index {index} outside 0..{count - 1}")
    eligible = catalog[(catalog["fragment"] == 1) & catalog["qc_pass"].astype(bool)]
    ordered = eligible.sort_values(["length", "uniprot_id"], ascending=[False, True])
    return ordered.iloc[index::count].reset_index(drop=True)
