#!/usr/bin/env python3
"""Extract pocket descriptors and training labels from a structure collection.

This is the bridge between the ground-truth dataset and model training. For each
structure it runs pocket detection, computes the full descriptor suite, and
labels each pocket by **atom-level overlap** with the observed ligand copies -
not by centroid distance, and not by a structure-level burial class.

Design points
-------------
* **Ambiguous pockets are excluded, not called negative.** Pockets that partially
  overlap the ligand are written with label ``-1`` and dropped at training time.
  Forcing them to 0 injects label noise exactly where the decision boundary lies.
* **Detector recall is reported.** If a structure contains a ligand but no pocket
  matches it, that is a pocket-detection failure and it bounds every downstream
  result. It is counted and printed rather than being absorbed into the negatives.
* **Grouping keys are emitted.** ``group_key`` (UniProt accession when known,
  otherwise the entry id) travels with every row so training can split by protein
  rather than by pocket.
* **Resumable.** Per-structure results are cached as JSON, so an interrupted run
  restarts where it stopped - which matters when the collection is thousands of
  structures.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.analysis.analyzer import ProteinAnalyzer  # noqa: E402
from cryptic_ip.analysis.features import FEATURE_NAMES  # noqa: E402
from cryptic_ip.analysis.labeling import (  # noqa: E402
    LigandSite,
    assign_pocket_labels,
    summarise_labels,
)
from cryptic_ip.analysis.structure_arrays import load_structure_arrays  # noqa: E402
from cryptic_ip.validation.burial_metrics import (  # noqa: E402
    compute_ligand_burial,
    find_ligand_instances,
)
from cryptic_ip.validation.structure_context import LIGAND_RESNAMES  # noqa: E402

LOGGER = logging.getLogger("extract_pocket_features")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--structures-dir",
        type=Path,
        default=Path("data/validation/raw/structures"),
        help="Directory of PDB/mmCIF structure files",
    )
    parser.add_argument(
        "--entry-csv",
        type=Path,
        default=Path("data/validation/ip_binding_validation_dataset.csv"),
        help="Per-entry dataset CSV supplying metadata and grouping keys",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/ml_training/pocket_features.csv"),
        help="Output feature/label table",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("results/ml_training/feature_cache"),
        help="Per-structure JSON cache enabling resume",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=Path("results/ml_training/labeling_summary.json"),
        help="Label and detector-recall summary output",
    )
    parser.add_argument(
        "--comp-ids",
        nargs="*",
        default=None,
        help="Ligand component identifiers; defaults to the built-in registry",
    )
    parser.add_argument(
        "--require-cryptic",
        action="store_true",
        help=(
            "Only cryptic ligand copies produce positives; surface copies become "
            "ambiguous and are excluded rather than treated as negatives."
        ),
    )
    parser.add_argument("--jobs", type=int, default=4, help="Parallel worker processes")
    parser.add_argument(
        "--sasa-points", type=int, default=256, help="SASA sample points per atom"
    )
    parser.add_argument(
        "--max-structures", type=int, default=None, help="Cap structures processed"
    )
    parser.add_argument("--min-alpha-spheres", type=int, default=3, help="fpocket -m parameter")
    parser.add_argument(
        "--skip-electrostatics",
        action="store_true",
        default=True,
        help="Skip APBS; the screened Coulomb surrogate is always computed",
    )
    parser.add_argument(
        "--with-electrostatics",
        dest="skip_electrostatics",
        action="store_false",
        help="Run APBS for the Poisson-Boltzmann potential feature",
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    return parser.parse_args(argv)


def _ligand_sites(path: Path, comp_ids: Sequence[str], sasa_points: int) -> List[LigandSite]:
    """Locate ligand copies and attach their burial class.

    Args:
        path: Structure file.
        comp_ids: Ligand component identifiers.
        sasa_points: SASA sample points per atom.

    Returns:
        Ligand copies, most buried first.
    """
    arrays = load_structure_arrays(path)
    found = find_ligand_instances(arrays, comp_ids)
    if not found:
        return []

    burial_by_id = {
        instance.instance_id: instance.burial_class
        for instance in compute_ligand_burial(path, comp_ids=list(comp_ids), n_points=sasa_points)
    }

    sites: List[LigandSite] = []
    for key, comp_id, atom_indices in found:
        suffix = f"{key[2]}{key[3]}".strip()
        instance_id = f"{comp_id}_{key[1]}_{suffix}"
        sites.append(
            LigandSite(
                instance_id=instance_id,
                comp_id=comp_id,
                coords=arrays.coords[atom_indices],
                burial_class=burial_by_id.get(instance_id, "unknown"),
            )
        )
    return sites


def process_structure(
    task: Tuple[str, str, Sequence[str], Dict[str, Any]]
) -> Tuple[str, List[Dict[str, Any]], Dict[str, Any]]:
    """Worker: detect, describe and label every pocket of one structure.

    Args:
        task: ``(structure_id, path, comp_ids, options)``.

    Returns:
        ``(structure_id, rows, diagnostics)``.
    """
    structure_id, path_str, comp_ids, options = task
    path = Path(path_str)
    diagnostics: Dict[str, Any] = {"structure_id": structure_id, "error": None}

    try:
        sites = _ligand_sites(path, comp_ids, int(options.get("sasa_points", 256)))
        diagnostics["n_ligand_instances"] = len(sites)

        analyzer = ProteinAnalyzer(
            str(path), skip_electrostatics=bool(options.get("skip_electrostatics", True))
        )
        analyzer.detect_pockets(min_alpha_sphere=int(options.get("min_alpha_spheres", 3)))
        if not options.get("skip_electrostatics", True):
            analyzer.calculate_electrostatics()
        analyzer.calculate_sasa()

        pocket_records: List[Dict[str, Any]] = []
        pocket_geometry: List[Tuple[int, np.ndarray, np.ndarray]] = []
        for pocket_id in analyzer.pockets["pocket_id"]:
            record = analyzer.analyze_pocket(int(pocket_id))
            pocket_records.append(record)
            centre = np.asarray(record["center"], dtype=float)
            spheres = analyzer._pocket_alpha_spheres(int(pocket_id))
            points = spheres if spheres is not None and len(spheres) else centre.reshape(1, 3)
            pocket_geometry.append((int(pocket_id), centre, points))

        assignments = assign_pocket_labels(
            pocket_geometry,
            sites,
            require_cryptic=bool(options.get("require_cryptic", False)),
        )
        assignment_by_id = {a.pocket_id: a for a in assignments}

        rows: List[Dict[str, Any]] = []
        for record in pocket_records:
            pocket_id = int(record["pocket_id"])
            row: Dict[str, Any] = {
                "structure_id": structure_id,
                "pocket_id": pocket_id,
            }
            row.update({name: record.get(name, np.nan) for name in FEATURE_NAMES})
            row["mean_local_hydrophobic_density"] = record.get(
                "mean_local_hydrophobic_density", np.nan
            )
            assignment = assignment_by_id.get(pocket_id)
            if assignment is not None:
                row.update(assignment.to_row())
                row["pocket_id"] = pocket_id
            rows.append(row)

        diagnostics["n_pockets"] = len(rows)
        diagnostics["n_positive"] = sum(1 for row in rows if row.get("label") == 1)
        analyzer.cleanup()
        return structure_id, rows, diagnostics

    except Exception as exc:  # noqa: BLE001 - one structure must not abort the run
        diagnostics["error"] = f"{type(exc).__name__}: {exc}"
        return structure_id, [], diagnostics


def load_entry_metadata(entry_csv: Path) -> Dict[str, Dict[str, Any]]:
    """Load per-entry metadata used for grouping and annotation.

    Args:
        entry_csv: Per-entry dataset CSV.

    Returns:
        Metadata keyed by upper-case structure identifier; empty when absent.
    """
    if not entry_csv.exists():
        LOGGER.warning("Entry metadata %s not found; grouping falls back to file stem", entry_csv)
        return {}
    frame = pd.read_csv(entry_csv)
    out: Dict[str, Dict[str, Any]] = {}
    for _, row in frame.iterrows():
        identifier = str(row.get("pdb_id", "")).upper()
        if identifier:
            out[identifier] = row.to_dict()
    return out


def group_key_for(structure_id: str, metadata: Dict[str, Any]) -> str:
    """Choose the grouping key that prevents leakage between related structures.

    The UniProt accession is preferred over the entry identifier: the PDB holds
    many entries of the same protein, and splitting by entry would place the same
    protein on both sides of a fold. Falls back to the entry identifier when no
    accession is annotated.

    Args:
        structure_id: Structure identifier.
        metadata: Entry metadata row.

    Returns:
        The group key.
    """
    accessions = str(metadata.get("uniprot_ids", "") or "").strip()
    if accessions:
        # Sorting makes the key deterministic for multi-chain complexes.
        return "|".join(sorted(accessions.split(";")))
    return structure_id.upper()


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    comp_ids = list(args.comp_ids) if args.comp_ids else sorted(LIGAND_RESNAMES)
    structures = sorted(
        [
            path
            for path in args.structures_dir.glob("*")
            if path.suffix.lower() in {".pdb", ".cif", ".ent", ".mmcif"}
        ]
    )
    if args.max_structures:
        structures = structures[: args.max_structures]
    if not structures:
        LOGGER.error("No structure files found in %s", args.structures_dir)
        return 2

    metadata = load_entry_metadata(args.entry_csv)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    options = {
        "sasa_points": args.sasa_points,
        "skip_electrostatics": args.skip_electrostatics,
        "min_alpha_spheres": args.min_alpha_spheres,
        "require_cryptic": args.require_cryptic,
    }

    tasks: List[Tuple[str, str, Sequence[str], Dict[str, Any]]] = []
    cached_rows: Dict[str, List[Dict[str, Any]]] = {}
    diagnostics: List[Dict[str, Any]] = []

    for path in structures:
        structure_id = path.stem.upper()
        cache_path = args.cache_dir / f"{structure_id}.json"
        if cache_path.exists():
            try:
                payload = json.loads(cache_path.read_text(encoding="utf-8"))
                cached_rows[structure_id] = payload["rows"]
                diagnostics.append(payload["diagnostics"])
                continue
            except (json.JSONDecodeError, KeyError, OSError):
                cache_path.unlink(missing_ok=True)
        tasks.append((structure_id, str(path), comp_ids, options))

    LOGGER.info(
        "Processing %d structures (%d already cached)", len(tasks), len(cached_rows)
    )

    def store(structure_id: str, rows: List[Dict[str, Any]], diag: Dict[str, Any]) -> None:
        cached_rows[structure_id] = rows
        diagnostics.append(diag)
        (args.cache_dir / f"{structure_id}.json").write_text(
            json.dumps({"rows": rows, "diagnostics": diag}, default=float), encoding="utf-8"
        )

    if args.jobs > 1 and tasks:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(process_structure, task): task[0] for task in tasks}
            for index, future in enumerate(as_completed(futures), start=1):
                structure_id, rows, diag = future.result()
                store(structure_id, rows, diag)
                if index % 25 == 0:
                    LOGGER.info("  %d/%d structures processed", index, len(tasks))
    else:
        for index, task in enumerate(tasks, start=1):
            structure_id, rows, diag = process_structure(task)
            store(structure_id, rows, diag)
            if index % 25 == 0:
                LOGGER.info("  %d/%d structures processed", index, len(tasks))

    all_rows: List[Dict[str, Any]] = []
    per_structure_assignments: Dict[str, List[Any]] = {}
    for structure_id, rows in cached_rows.items():
        entry_meta = metadata.get(structure_id, {})
        for row in rows:
            row = dict(row)
            row["structure_id"] = structure_id
            row["group_key"] = group_key_for(structure_id, entry_meta)
            row["uniprot_ids"] = entry_meta.get("uniprot_ids", "")
            row["organism"] = entry_meta.get("organism", "")
            row["resolution"] = entry_meta.get("resolution", np.nan)
            row["structure_classification"] = entry_meta.get("classification", "")
            all_rows.append(row)

    if not all_rows:
        LOGGER.error("No pockets were extracted; check that fpocket is installed and working.")
        return 3

    frame = pd.DataFrame(all_rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_csv, index=False)

    # Detector-recall accounting: which ligand-bearing structures yielded a site?
    structures_with_ligand = [
        diag["structure_id"]
        for diag in diagnostics
        if diag.get("n_ligand_instances", 0) > 0 and not diag.get("error")
    ]
    from cryptic_ip.analysis.labeling import PocketAssignment, PocketLabel

    for structure_id, rows in cached_rows.items():
        per_structure_assignments[structure_id] = [
            PocketAssignment(
                pocket_id=int(row.get("pocket_id", 0)),
                label=PocketLabel(int(row.get("label", 0))),
                overlap_fraction=float(row.get("overlap_fraction", 0.0) or 0.0),
            )
            for row in rows
        ]
    summary = summarise_labels(per_structure_assignments, structures_with_ligand)

    payload = {
        **summary.to_dict(),
        "n_structures_processed": len(cached_rows),
        "n_failures": sum(1 for diag in diagnostics if diag.get("error")),
        "failures": {
            diag["structure_id"]: diag["error"] for diag in diagnostics if diag.get("error")
        },
        "n_groups": int(frame["group_key"].nunique()),
        "feature_columns": list(FEATURE_NAMES),
        "output_csv": str(args.output_csv),
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    LOGGER.info(
        "Wrote %d pockets from %d structures (%d groups): %d positive, %d negative, %d ambiguous",
        len(frame),
        summary.n_structures,
        payload["n_groups"],
        summary.n_positive,
        summary.n_negative,
        summary.n_ambiguous,
    )
    if np.isfinite(summary.site_recall):
        LOGGER.info(
            "Pocket-detector recall on known ligand sites: %.1f%% (%d/%d structures)",
            100.0 * summary.site_recall,
            summary.n_structures_with_site,
            summary.n_structures_with_site + summary.n_structures_with_ligand_but_no_site,
        )
        if summary.site_recall < 0.8:
            LOGGER.warning(
                "Detector recall is below 80%%. Every missed site is unreachable by the "
                "classifier, so this bounds the whole pipeline; consider lowering "
                "--min-alpha-spheres."
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
