#!/usr/bin/env python3
"""Measure per-copy relative burial across every deposited IP complex available.

Why this exists
---------------
The burial boundary that separates a cryptic site from a surface site is the
single most consequential constant in this pipeline: it defines the positive
class. It is currently calibrated on five hand-picked control structures, which
is enough to show that a gap exists but not enough to say where in that gap the
boundary belongs, nor whether the distribution is even bimodal.

The bundled validation set holds 136 deposited entries containing inositol
phosphates. Their ``classification`` column, however, was produced by the earlier
burial definition - solvent-accessible surface summed over *all* copies of a
ligand and compared against an absolute cutoff. That measure grows with the
number of copies in the asymmetric unit and with ligand size, so it cannot be
compared across entries, and it labelled 133 of 136 entries "Surface". Those
labels are therefore not usable as evidence.

This script re-measures the whole set with the corrected definition - burial
computed per ligand copy and normalised by that copy's own isolated surface -
and reports the resulting distribution rather than a verdict. It answers three
questions a five-structure panel cannot:

* Is relative burial bimodal, as a cryptic/surface dichotomy would require, or
  continuous?
* Where does the density minimum between the modes actually fall?
* How many entries would change class under the calibrated boundary?

The output is a distribution and a set of candidate thresholds with the evidence
for each, not a recommendation. Choosing the boundary remains a judgement, but
this makes it a judgement about data.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.analysis.inositol_detection import (  # noqa: E402
    detect_inositol_residues,
)
from cryptic_ip.analysis.structure_arrays import load_structure_arrays  # noqa: E402
from cryptic_ip.validation.burial_metrics import (  # noqa: E402
    CRYPTIC_RELATIVE_SASA_MAX,
    SEMI_CRYPTIC_RELATIVE_SASA_MAX,
    compute_ligand_burial,
)

LOGGER = logging.getLogger("burial_survey")

#: Bin width for the density estimate used to locate a minimum between modes.
#: Relative SASA lies in [0, 1]; 40 bins gives 0.025 resolution, fine enough to
#: place a boundary to two decimals without chasing single-entry noise.
N_DENSITY_BINS = 40

#: Candidate boundaries are searched only in this interior range. A "boundary"
#: at 0.02 or 0.95 would be an artefact of the tails rather than a class split.
SEARCH_LOW = 0.05
SEARCH_HIGH = 0.60

#: Each side of a candidate trough must hold at least this share of the data.
#: Below this a "mode" is a handful of entries, and a boundary drawn against it
#: would be fitting noise.
MIN_MODE_MASS_FRACTION = 0.10

#: A trough must hold at most this fraction of the smaller flanking peak.
#: Without it, the shoulder of a single broad mode qualifies as a gap.
MAX_DIP_RATIO = 0.5

#: Boundaries at which class membership is reported. Because burial turns out to
#: be continuous rather than two-class, the positive class is *defined* by
#: whichever cutoff is chosen, and a single count hides how sensitive that
#: definition is. Reporting the sweep makes the choice inspectable: a boundary
#: whose neighbours give wildly different class sizes is a fragile one.
BOUNDARY_SWEEP = (0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30)


@dataclass
class EntryMeasurement:
    """Per-entry burial measurement.

    Attributes:
        pdb_id: Entry identifier.
        comp_id: Component the most-buried copy belongs to.
        series: Inositol phosphate series of that copy, e.g. ``"InsP6"``.
        n_copies: Ligand copies measured in the entry.
        relative_sasa: Relative SASA of the most buried copy.
        relative_phosphate_sasa: The same restricted to phosphate groups.
        enclosure: Enclosure fraction of the most buried copy.
        burial_depth: Burial depth of the most buried copy (Å).
        n_basic_residues: Basic residues coordinating that copy.
        burial_class: Class assigned by the current thresholds.
        error: Why the entry could not be measured, when applicable.
    """

    pdb_id: str
    comp_id: Optional[str] = None
    series: Optional[str] = None
    n_copies: int = 0
    relative_sasa: Optional[float] = None
    relative_phosphate_sasa: Optional[float] = None
    enclosure: Optional[float] = None
    burial_depth: Optional[float] = None
    n_basic_residues: Optional[int] = None
    burial_class: Optional[str] = None
    error: Optional[str] = None


def _finite(value: Any) -> Optional[float]:
    """Return a float when finite, otherwise ``None``."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def measure_entry(path: Path, pdb_id: str, *, sasa_points: int) -> EntryMeasurement:
    """Measure the most buried inositol phosphate copy in one entry.

    Args:
        path: Structure file.
        pdb_id: Entry identifier, for reporting.
        sasa_points: SASA sample points per atom.

    Returns:
        The measurement, carrying ``error`` when the entry could not be used.
    """
    try:
        instances = compute_ligand_burial(path, n_points=sasa_points)
    except Exception as exc:  # noqa: BLE001 - one bad entry must not stop a survey
        return EntryMeasurement(pdb_id=pdb_id, error=f"{type(exc).__name__}: {exc}")

    if not instances:
        return EntryMeasurement(
            pdb_id=pdb_id, error="no phosphorylated inositol ligand detected"
        )

    best = instances[0]
    series = None
    try:
        arrays = load_structure_arrays(path)
        for residue in detect_inositol_residues(arrays, require_phosphate=False):
            if residue.comp_id == best.comp_id:
                series = residue.series
                break
    except Exception:  # noqa: BLE001 - series is provenance, not the measurement
        series = None

    return EntryMeasurement(
        pdb_id=pdb_id,
        comp_id=best.comp_id,
        series=series,
        n_copies=len(instances),
        relative_sasa=_finite(best.relative_sasa),
        relative_phosphate_sasa=_finite(best.relative_phosphate_sasa),
        enclosure=_finite(best.enclosure),
        burial_depth=_finite(best.burial_depth),
        n_basic_residues=int(best.n_basic_residues),
        burial_class=best.burial_class,
    )


def _worker(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Process-pool entry point."""
    return asdict(
        measure_entry(
            Path(payload["path"]), payload["pdb_id"], sasa_points=payload["sasa_points"]
        )
    )


def density_minimum(
    values: Sequence[float],
    *,
    n_bins: int = N_DENSITY_BINS,
    low: float = SEARCH_LOW,
    high: float = SEARCH_HIGH,
    min_mass_fraction: float = MIN_MODE_MASS_FRACTION,
    dip_ratio: float = MAX_DIP_RATIO,
) -> Optional[float]:
    """Locate the trough between two modes, or report that there is none.

    A cryptic/surface dichotomy implies two populations separated by a sparse
    region, and the natural boundary is that sparse region. A minimum on its own
    is not evidence of one - every distribution has a smallest interior bin - so
    a candidate must satisfy both conditions that make a trough a trough:

    * **mass on both sides**, each at least ``min_mass_fraction`` of the total,
      so a thin tail cannot pass as a population;
    * **an actual dip**, the candidate bin holding at most ``dip_ratio`` of the
      smaller flanking peak, so the shoulder of a single mode cannot pass as a
      gap between two.

    Returning ``None`` is a real result: it says the measurements give no
    evidence of a split, and a boundary drawn on them would be arbitrary.

    Args:
        values: Relative SASA measurements.
        n_bins: Histogram bins across [0, 1].
        low: Lower end of the interior search range.
        high: Upper end of the interior search range.
        min_mass_fraction: Minimum share of the data required on each side.
        dip_ratio: Maximum trough height as a fraction of the smaller peak.

    Returns:
        Midpoint of the emptiest qualifying bin, or ``None`` when none qualifies.
    """
    finite = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    if finite.size < 3:
        return None

    counts, edges = np.histogram(finite, bins=n_bins, range=(0.0, 1.0))
    centres = 0.5 * (edges[:-1] + edges[1:])
    total = int(counts.sum())
    if total == 0:
        return None

    best_index: Optional[int] = None
    best_count = np.inf
    for index in np.flatnonzero((centres >= low) & (centres <= high)):
        below = counts[:index]
        above = counts[index + 1 :]
        if below.size == 0 or above.size == 0:
            continue
        if below.sum() < min_mass_fraction * total:
            continue
        if above.sum() < min_mass_fraction * total:
            continue
        smaller_peak = min(int(below.max()), int(above.max()))
        if smaller_peak <= 0:
            continue
        if counts[index] > dip_ratio * smaller_peak:
            continue
        if counts[index] < best_count:
            best_count = counts[index]
            best_index = int(index)

    if best_index is None:
        return None
    return float(centres[best_index])


def otsu_threshold(values: Sequence[float], *, n_bins: int = N_DENSITY_BINS) -> Optional[float]:
    """Return the threshold maximising between-class variance.

    Otsu's criterion picks the split that makes the two resulting groups as
    internally homogeneous as possible. It always returns a value, so it is
    reported *alongside* :func:`density_minimum` rather than instead of it: the
    two agreeing is evidence of a real split, the two disagreeing is evidence
    that the distribution is closer to continuous.

    Args:
        values: Relative SASA measurements.
        n_bins: Histogram bins across [0, 1].

    Returns:
        Threshold value, or ``None`` when there is too little data.
    """
    finite = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    if finite.size < 3:
        return None

    counts, edges = np.histogram(finite, bins=n_bins, range=(0.0, 1.0))
    centres = 0.5 * (edges[:-1] + edges[1:])
    total = counts.sum()
    if total == 0:
        return None

    weights = counts / total
    cumulative = np.cumsum(weights)
    means = np.cumsum(weights * centres)
    grand_mean = means[-1]

    denominator = cumulative * (1.0 - cumulative)
    with np.errstate(divide="ignore", invalid="ignore"):
        between = (grand_mean * cumulative - means) ** 2 / denominator
    between[~np.isfinite(between)] = -np.inf
    if not np.any(np.isfinite(between)):
        return None
    return float(centres[int(np.argmax(between))])


def summarise(measurements: Sequence[EntryMeasurement]) -> Dict[str, Any]:
    """Summarise the survey into a distribution and candidate boundaries.

    Args:
        measurements: Per-entry measurements.

    Returns:
        A JSON-serialisable summary.
    """
    usable = [m for m in measurements if m.relative_sasa is not None]
    values = [m.relative_sasa for m in usable]

    summary: Dict[str, Any] = {
        "n_entries": len(measurements),
        "n_measured": len(usable),
        "n_failed": len(measurements) - len(usable),
        "failure_reasons": {},
    }
    for measurement in measurements:
        if measurement.relative_sasa is None:
            reason = measurement.error or "unknown"
            summary["failure_reasons"][reason] = (
                summary["failure_reasons"].get(reason, 0) + 1
            )

    if not values:
        return summary

    array = np.asarray(values, dtype=float)
    summary["relative_sasa"] = {
        "min": float(array.min()),
        "max": float(array.max()),
        "mean": float(array.mean()),
        "median": float(np.median(array)),
        "quantiles": {
            f"q{int(q * 100):02d}": float(np.quantile(array, q))
            for q in (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
        },
    }

    counts, edges = np.histogram(array, bins=N_DENSITY_BINS, range=(0.0, 1.0))
    summary["histogram"] = {
        "bin_edges": [float(e) for e in edges],
        "counts": [int(c) for c in counts],
    }

    trough = density_minimum(values)
    otsu = otsu_threshold(values)
    summary["candidate_boundaries"] = {
        "density_minimum": trough,
        "otsu": otsu,
        "configured_cryptic_max": float(CRYPTIC_RELATIVE_SASA_MAX),
        "configured_semi_cryptic_max": float(SEMI_CRYPTIC_RELATIVE_SASA_MAX),
    }
    summary["is_bimodal_by_density"] = trough is not None

    by_class: Dict[str, int] = {}
    for measurement in usable:
        key = measurement.burial_class or "unknown"
        by_class[key] = by_class.get(key, 0) + 1
    summary["class_counts"] = by_class

    by_series: Dict[str, int] = {}
    for measurement in usable:
        key = measurement.series or "unclassified"
        by_series[key] = by_series.get(key, 0) + 1
    summary["series_counts"] = by_series

    summary["n_below_configured_boundary"] = int(
        (array <= CRYPTIC_RELATIVE_SASA_MAX).sum()
    )

    # How many entries the positive class would contain at each candidate
    # boundary, and how fast that count moves. On a continuous distribution the
    # class size is a function of the cutoff, so this is the honest way to show
    # what the choice costs.
    sweep = []
    for boundary in BOUNDARY_SWEEP:
        n_positive = int((array <= boundary).sum())
        sweep.append(
            {
                "boundary": float(boundary),
                "n_positive": n_positive,
                "fraction_positive": float(n_positive / array.size),
            }
        )
    summary["boundary_sweep"] = sweep
    return summary


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--structures-dir", type=Path, default=Path("data/validation"))
    parser.add_argument(
        "--entry-csv",
        type=Path,
        default=Path("data/validation/ip_binding_validation_dataset.csv"),
        help="Dataset listing entries to survey; missing structures are skipped",
    )
    parser.add_argument(
        "--output-json", type=Path, default=Path("results/calibration/burial_survey.json")
    )
    parser.add_argument(
        "--output-csv", type=Path, default=Path("results/calibration/burial_survey.csv")
    )
    parser.add_argument("--sasa-points", type=int, default=256)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument(
        "--limit", type=int, default=0, help="Survey at most this many entries (0 = all)"
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    return parser.parse_args(argv)


def resolve_structures(
    structures_dir: Path, entry_csv: Path, *, limit: int = 0
) -> List[Dict[str, str]]:
    """Match dataset entries to structure files present on disk.

    Args:
        structures_dir: Directory searched for coordinate files.
        entry_csv: Dataset listing ``pdb_id`` values.
        limit: Maximum entries to return (0 for all).

    Returns:
        ``{"pdb_id", "path"}`` for each entry whose structure is available.
    """
    import pandas as pd

    if not entry_csv.exists():
        raise FileNotFoundError(f"entry CSV not found: {entry_csv}")

    # Read identifiers as text: pandas would infer an all-digit identifier as an
    # integer and strip its leading zeros, so the file it names would never be
    # found. PDB identifiers are alphanumeric tokens, not numbers.
    frame = pd.read_csv(entry_csv, dtype={"pdb_id": str})
    pdb_ids = [str(v).strip().upper() for v in frame["pdb_id"].dropna().unique()]

    resolved: List[Dict[str, str]] = []
    for pdb_id in pdb_ids:
        for suffix in (".pdb", ".cif", ".ent", ".mmcif"):
            for candidate in (
                structures_dir / f"{pdb_id}{suffix}",
                structures_dir / "raw" / "structures" / f"{pdb_id}{suffix}",
            ):
                if candidate.exists():
                    resolved.append({"pdb_id": pdb_id, "path": str(candidate)})
                    break
            else:
                continue
            break
        if limit and len(resolved) >= limit:
            break
    return resolved


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the survey."""
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    targets = resolve_structures(args.structures_dir, args.entry_csv, limit=args.limit)
    if not targets:
        LOGGER.error(
            "No structures found under %s for entries in %s. Download them first.",
            args.structures_dir,
            args.entry_csv,
        )
        return 1

    LOGGER.info("Surveying %d entries with structures on disk", len(targets))
    payloads = [dict(target, sasa_points=args.sasa_points) for target in targets]

    records: List[Dict[str, Any]] = []
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(_worker, payload): payload for payload in payloads}
            for future in as_completed(futures):
                records.append(future.result())
    else:
        records = [_worker(payload) for payload in payloads]

    records.sort(key=lambda record: str(record["pdb_id"]))
    measurements = [EntryMeasurement(**record) for record in records]
    summary = summarise(measurements)

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps({"summary": summary, "entries": records}, indent=2, allow_nan=False),
        encoding="utf-8",
    )

    import pandas as pd

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(args.output_csv, index=False)

    print("\n" + "=" * 78)
    print("PER-COPY RELATIVE BURIAL SURVEY")
    print("=" * 78)
    print(f"  entries with structures : {summary['n_entries']}")
    print(f"  measured                : {summary['n_measured']}")
    print(f"  failed                  : {summary['n_failed']}")
    for reason, count in sorted(summary["failure_reasons"].items()):
        print(f"      {count:4d}  {reason}")

    if summary.get("relative_sasa"):
        stats = summary["relative_sasa"]
        print("\n  relative SASA distribution")
        print(f"      min {stats['min']:.3f}   median {stats['median']:.3f}   max {stats['max']:.3f}")
        quantiles = stats["quantiles"]
        print("      " + "  ".join(f"{k} {v:.3f}" for k, v in quantiles.items()))

        candidates = summary["candidate_boundaries"]
        print("\n  candidate boundaries")
        trough = candidates["density_minimum"]
        print(f"      density minimum : {'none - no trough found' if trough is None else f'{trough:.3f}'}")
        print(f"      Otsu            : {candidates['otsu']:.3f}")
        print(f"      configured       : {candidates['configured_cryptic_max']:.3f}")
        print(f"\n  bimodal by density: {summary['is_bimodal_by_density']}")
        print(f"  entries at or below the configured boundary: {summary['n_below_configured_boundary']}")
        print(f"  class counts: {summary['class_counts']}")
        print(f"  series counts: {summary['series_counts']}")

    print(f"\nWrote {args.output_json} and {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
