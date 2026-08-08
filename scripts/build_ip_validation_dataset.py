#!/usr/bin/env python3
"""Build the inositol phosphate ground-truth dataset from the RCSB PDB.

What this collects
------------------
1. **Every** inositol phosphate chemical component in the PDB, discovered by
   full-text search of the chemical component dictionary and validated against
   the component API, rather than a hand-written list of eight identifiers.
2. Every structure containing at least one of those components, across **all**
   experimental methods. Restricting to X-ray - as the previous version did -
   discarded the cryo-EM structures of large assemblies where inositol phosphates
   most often act as structural cofactors.
3. Optionally, a matched set of high-resolution structures containing **no**
   inositol phosphate. These supply protein-level negatives: without them every
   negative pocket comes from an IP-binding protein, which makes the benchmark
   easier than a real proteome screen.
4. Per **ligand copy** burial measurements - relative SASA, relative phosphate
   SASA, burial depth, enclosure, coordination counts - rather than a single
   number summed over copies.

Why per-copy matters
--------------------
The previous builder summed ligand SASA over every copy of a residue name in an
entry and classified the entry from that sum. Copy number is a crystallisation
artefact, so a structure with six InsP6 copies scored six times the SASA of one
with a single copy, and nearly every entry was labelled ``Surface``. That single
defect propagated into the training labels and left the shipped model at chance.

Outputs
-------
``--output-csv``
    One row per ligand copy, with burial measurements and entry metadata.
``--entry-csv``
    One row per entry, summarising its most buried copy.
``--manifest``
    JSON provenance: ligand registry, queries issued, file checksums, API
    statistics, and software versions.

Network access to ``search.rcsb.org``, ``data.rcsb.org`` and ``files.rcsb.org``
is required. The script fails with an actionable message if they are unreachable
rather than silently writing an empty dataset.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import platform
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

# Allow running as a plain script from the repository root.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.database.ip_ligands import (  # noqa: E402
    IPLigand,
    comp_ids as ligand_comp_ids,
    discover_ip_ligands,
    phosphorylated_comp_ids,
)
from cryptic_ip.database.rcsb_client import (  # noqa: E402
    RcsbClient,
    RcsbUnavailableError,
    sha256_file,
)
from cryptic_ip.validation.burial_metrics import compute_ligand_burial  # noqa: E402

LOGGER = logging.getLogger("build_ip_validation_dataset")

LIGAND_FIELDS: Tuple[str, ...] = (
    "pdb_id",
    "instance_id",
    "comp_id",
    "chain_id",
    "resseq",
    "icode",
    "model_id",
    "series",
    "n_atoms",
    "n_phosphorus",
    "sasa_complex",
    "sasa_isolated",
    "relative_sasa",
    "phosphate_sasa_complex",
    "phosphate_sasa_isolated",
    "relative_phosphate_sasa",
    "burial_depth",
    "enclosure",
    "n_protein_contacts",
    "n_contact_chains",
    "n_basic_residues",
    "n_basic_nitrogens",
    "n_hydroxyl_residues",
    "mean_occupancy",
    "mean_bfactor",
    "radius_of_gyration",
    "burial_class",
    "is_probable_artifact",
    "uniprot_ids",
    "organism",
    "taxonomy_id",
    "resolution",
    "experimental_method",
    "r_free",
    "release_date",
    "title",
    "structure_path",
)

ENTRY_FIELDS: Tuple[str, ...] = (
    "pdb_id",
    "n_ligand_instances",
    "n_cryptic_instances",
    "best_instance_id",
    "best_comp_id",
    "best_relative_sasa",
    "best_relative_phosphate_sasa",
    "best_burial_depth",
    "best_enclosure",
    "best_n_basic_residues",
    "classification",
    "uniprot_ids",
    "organism",
    "taxonomy_id",
    "resolution",
    "experimental_method",
    "r_free",
    "release_date",
    "has_ip_ligand",
    "structure_path",
)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/validation/ip_ligand_instances.csv"),
        help="Per-ligand-copy output CSV",
    )
    parser.add_argument(
        "--entry-csv",
        type=Path,
        default=Path("data/validation/ip_binding_validation_dataset.csv"),
        help="Per-entry summary CSV",
    )
    parser.add_argument(
        "--download-dir",
        type=Path,
        default=Path("data/validation/raw"),
        help="Directory for downloaded structures and cached API responses",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/validation/dataset_manifest.json"),
        help="Provenance manifest output path",
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=None,
        help="Cap the number of ligand-bearing entries (for smoke runs)",
    )
    parser.add_argument(
        "--n-decoys",
        type=int,
        default=0,
        help=(
            "Number of inositol-phosphate-free entries to collect as protein-level "
            "negatives. 0 disables decoy collection."
        ),
    )
    parser.add_argument(
        "--decoy-max-resolution",
        type=float,
        default=2.0,
        help="Resolution ceiling for decoy entries (Angstrom)",
    )
    parser.add_argument(
        "--experimental-methods",
        nargs="*",
        default=None,
        help=(
            "Restrict to these exptl.method values. Default accepts every method, "
            "which keeps cryo-EM structures of large IP-dependent assemblies."
        ),
    )
    parser.add_argument(
        "--max-resolution",
        type=float,
        default=None,
        help="Optional resolution ceiling for ligand-bearing entries (Angstrom)",
    )
    parser.add_argument(
        "--include-unphosphorylated",
        action="store_true",
        help="Include free inositols (InsP0) as reference ligands",
    )
    parser.add_argument("--jobs", type=int, default=4, help="Parallel worker processes")
    parser.add_argument(
        "--download-workers", type=int, default=4, help="Concurrent download threads"
    )
    parser.add_argument(
        "--sasa-points",
        type=int,
        default=512,
        help="SASA sample points per atom; higher is more precise and slower",
    )
    parser.add_argument(
        "--min-proteins",
        type=int,
        default=20,
        help="Warn when fewer unique UniProt accessions are collected",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level",
    )
    return parser.parse_args(argv)


def _first_resolution(entry: Dict[str, Any]) -> Optional[float]:
    """Extract the reported resolution from an entry metadata document."""
    info = entry.get("rcsb_entry_info") or {}
    values = info.get("resolution_combined") or []
    if values:
        try:
            return float(values[0])
        except (TypeError, ValueError):
            return None
    return None


def _entry_annotations(entry: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Flatten an entry metadata document into scalar annotation fields.

    Metadata is best-effort: coordinates drive every measurement, so a metadata
    gap reduces annotation richness without invalidating any number.

    Args:
        entry: GraphQL entry document, or ``None``.

    Returns:
        Annotation fields with empty strings for anything unavailable.
    """
    if not entry:
        return {
            "uniprot_ids": "",
            "organism": "",
            "taxonomy_id": "",
            "resolution": "",
            "experimental_method": "",
            "r_free": "",
            "release_date": "",
            "title": "",
        }

    accessions: List[str] = []
    organisms: List[str] = []
    taxa: List[str] = []
    for polymer in entry.get("polymer_entities") or []:
        identifiers = polymer.get("rcsb_polymer_entity_container_identifiers") or {}
        for reference in identifiers.get("reference_sequence_identifiers") or []:
            if (reference or {}).get("database_name") == "UniProt":
                accession = (reference or {}).get("database_accession")
                if accession:
                    accessions.append(str(accession))
        for source in polymer.get("rcsb_entity_source_organism") or []:
            name = (source or {}).get("ncbi_scientific_name")
            if name:
                organisms.append(str(name))
            taxon = (source or {}).get("ncbi_taxonomy_id")
            if taxon:
                taxa.append(str(taxon))

    methods = [m.get("method") for m in (entry.get("exptl") or []) if m and m.get("method")]
    refine = (entry.get("refine") or [{}])
    r_free = None
    if refine and isinstance(refine, list) and refine[0]:
        r_free = refine[0].get("ls_R_factor_R_free")
    accession_info = entry.get("rcsb_accession_info") or {}
    resolution = _first_resolution(entry)

    return {
        "uniprot_ids": ";".join(sorted(set(accessions))),
        "organism": ";".join(sorted(set(organisms))),
        "taxonomy_id": ";".join(sorted(set(taxa))),
        "resolution": "" if resolution is None else f"{resolution:.2f}",
        "experimental_method": ";".join(sorted(set(m for m in methods if m))),
        "r_free": "" if r_free in (None, "") else f"{float(r_free):.4f}",
        "release_date": str(accession_info.get("initial_release_date") or ""),
        "title": str((entry.get("struct") or {}).get("title") or ""),
    }


def _measure_structure(
    args: Tuple[str, str, Sequence[str], int]
) -> Tuple[str, List[Dict[str, Any]], Optional[str]]:
    """Worker: measure every ligand copy in one structure file.

    Defined at module scope so it can be pickled for a process pool.

    Args:
        args: ``(pdb_id, structure_path, comp_ids, sasa_points)``.

    Returns:
        ``(pdb_id, instance_dicts, error_message)``.
    """
    pdb_id, structure_path, comp_ids, sasa_points = args
    try:
        instances = compute_ligand_burial(
            Path(structure_path), comp_ids=list(comp_ids), n_points=int(sasa_points)
        )
        return pdb_id, [instance.to_dict() for instance in instances], None
    except Exception as exc:  # noqa: BLE001 - one bad entry must not abort the run
        return pdb_id, [], f"{type(exc).__name__}: {exc}"


def resolve_ligand_registry(
    client: RcsbClient, *, include_unphosphorylated: bool
) -> Tuple[List[IPLigand], Dict[str, Any]]:
    """Discover and validate the inositol phosphate component vocabulary."""
    ligands, provenance = discover_ip_ligands(
        client, include_unphosphorylated=include_unphosphorylated
    )
    if not ligands:
        raise RcsbUnavailableError(
            "No inositol phosphate components could be resolved. This usually means "
            "the RCSB chemical component API is unreachable; check network access "
            "and any outbound proxy policy."
        )
    LOGGER.info(
        "Ligand vocabulary: %s",
        ", ".join(f"{lig.comp_id}({lig.series})" for lig in ligands),
    )
    return ligands, provenance


def collect_entries(
    client: RcsbClient,
    ligands: Sequence[IPLigand],
    args: argparse.Namespace,
) -> Tuple[List[str], List[str]]:
    """Search for ligand-bearing entries and, optionally, decoy entries."""
    target_ids = list(phosphorylated_comp_ids(ligands) or ligand_comp_ids(ligands))
    entries = client.search_entries_with_components(
        target_ids,
        experimental_methods=args.experimental_methods,
        max_resolution=args.max_resolution,
        max_results=args.max_entries,
    )
    LOGGER.info("Found %d entries containing an inositol phosphate", len(entries))

    decoys: List[str] = []
    if args.n_decoys > 0:
        decoys = client.search_decoy_entries(
            exclude_comp_ids=list(ligand_comp_ids(ligands)),
            max_resolution=args.decoy_max_resolution,
            max_results=args.n_decoys,
        )
        decoys = [pdb_id for pdb_id in decoys if pdb_id not in set(entries)][: args.n_decoys]
        LOGGER.info("Selected %d inositol-phosphate-free decoy entries", len(decoys))
    return entries, decoys


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, Any]]) -> int:
    """Write rows to CSV, restricted to the declared field names.

    Args:
        path: Output path.
        fieldnames: Column order.
        rows: Row mappings.

    Returns:
        Number of rows written.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
            count += 1
    return count


def build_dataset(args: argparse.Namespace) -> Dict[str, Any]:
    """Run the full collection and write every output.

    Args:
        args: Parsed command-line arguments.

    Returns:
        The provenance manifest.
    """
    started = datetime.now(timezone.utc)
    structures_dir = args.download_dir / "structures"
    cache_dir = args.download_dir / "api_cache"
    client = RcsbClient(cache_dir=cache_dir)

    ligands, ligand_provenance = resolve_ligand_registry(
        client, include_unphosphorylated=args.include_unphosphorylated
    )
    target_comp_ids = list(ligand_comp_ids(ligands))

    entries, decoys = collect_entries(client, ligands, args)
    all_ids = list(entries) + list(decoys)
    if not all_ids:
        raise RcsbUnavailableError("Search returned no entries; nothing to build.")

    LOGGER.info("Downloading %d structures", len(all_ids))
    downloads = client.download_structures(
        all_ids, structures_dir, max_workers=args.download_workers
    )
    download_by_id = {record.identifier: record for record in downloads}

    LOGGER.info("Fetching entry metadata")
    metadata = client.fetch_entry_metadata(all_ids)

    ligand_rows: List[Dict[str, Any]] = []
    entry_rows: List[Dict[str, Any]] = []
    failures: Dict[str, str] = {}
    unique_proteins: set = set()

    tasks = [
        (pdb_id, download_by_id[pdb_id].path, target_comp_ids, args.sasa_points)
        for pdb_id in all_ids
        if pdb_id in download_by_id
    ]
    LOGGER.info("Measuring ligand burial for %d structures on %d workers", len(tasks), args.jobs)

    measured: Dict[str, List[Dict[str, Any]]] = {}
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(_measure_structure, task): task[0] for task in tasks}
            for index, future in enumerate(as_completed(futures), start=1):
                pdb_id, instances, error = future.result()
                if error:
                    failures[pdb_id] = error
                measured[pdb_id] = instances
                if index % 50 == 0:
                    LOGGER.info("  measured %d/%d", index, len(tasks))
    else:
        for index, task in enumerate(tasks, start=1):
            pdb_id, instances, error = _measure_structure(task)
            if error:
                failures[pdb_id] = error
            measured[pdb_id] = instances
            if index % 50 == 0:
                LOGGER.info("  measured %d/%d", index, len(tasks))

    decoy_set = set(decoys)
    series_by_comp = {lig.comp_id: lig.series for lig in ligands}

    for pdb_id in all_ids:
        record = download_by_id.get(pdb_id)
        if record is None:
            failures.setdefault(pdb_id, "structure file unavailable")
            continue
        annotations = _entry_annotations(metadata.get(pdb_id))
        instances = measured.get(pdb_id, [])
        if annotations["uniprot_ids"]:
            unique_proteins.update(annotations["uniprot_ids"].split(";"))

        for instance in instances:
            row = dict(instance)
            row.update(annotations)
            row["pdb_id"] = pdb_id
            row["series"] = series_by_comp.get(str(instance.get("comp_id")), "")
            row["structure_path"] = record.path
            ligand_rows.append(row)

        if instances:
            best = instances[0]
            entry_rows.append(
                {
                    "pdb_id": pdb_id,
                    "n_ligand_instances": len(instances),
                    "n_cryptic_instances": sum(
                        1 for inst in instances if inst.get("burial_class") == "cryptic"
                    ),
                    "best_instance_id": best.get("instance_id"),
                    "best_comp_id": best.get("comp_id"),
                    "best_relative_sasa": best.get("relative_sasa"),
                    "best_relative_phosphate_sasa": best.get("relative_phosphate_sasa"),
                    "best_burial_depth": best.get("burial_depth"),
                    "best_enclosure": best.get("enclosure"),
                    "best_n_basic_residues": best.get("n_basic_residues"),
                    "classification": best.get("burial_class"),
                    "has_ip_ligand": True,
                    "structure_path": record.path,
                    **annotations,
                }
            )
        else:
            entry_rows.append(
                {
                    "pdb_id": pdb_id,
                    "n_ligand_instances": 0,
                    "n_cryptic_instances": 0,
                    "classification": "decoy" if pdb_id in decoy_set else "no_ligand_parsed",
                    "has_ip_ligand": False,
                    "structure_path": record.path,
                    **annotations,
                }
            )

    n_ligand_rows = write_csv(args.output_csv, LIGAND_FIELDS, ligand_rows)
    n_entry_rows = write_csv(args.entry_csv, ENTRY_FIELDS, entry_rows)

    class_counts: Dict[str, int] = {}
    for row in ligand_rows:
        key = str(row.get("burial_class", "unknown"))
        class_counts[key] = class_counts.get(key, 0) + 1

    manifest = {
        "generated_at_utc": started.isoformat(),
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "ligand_registry": [lig.to_dict() for lig in ligands],
        "ligand_discovery": ligand_provenance,
        "search": {
            "n_ligand_entries": len(entries),
            "n_decoy_entries": len(decoys),
            "experimental_methods": args.experimental_methods,
            "max_resolution": args.max_resolution,
            "decoy_max_resolution": args.decoy_max_resolution,
        },
        "downloads": {
            "n_requested": len(all_ids),
            "n_retrieved": len(downloads),
            "files": [record.to_dict() for record in downloads],
        },
        "measurements": {
            "n_ligand_instances": n_ligand_rows,
            "n_entries": n_entry_rows,
            "burial_class_counts": class_counts,
            "n_unique_uniprot": len(unique_proteins),
            "sasa_points_per_atom": args.sasa_points,
            "failures": failures,
        },
        "api_statistics": client.stats.to_dict(),
        "outputs": {
            "ligand_csv": str(args.output_csv),
            "entry_csv": str(args.entry_csv),
        },
    }
    for key, path in (("ligand_csv", args.output_csv), ("entry_csv", args.entry_csv)):
        if path.exists():
            manifest["outputs"][f"{key}_sha256"] = sha256_file(path)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    LOGGER.info(
        "Wrote %d ligand copies across %d entries (%s)",
        n_ligand_rows,
        n_entry_rows,
        ", ".join(f"{name}={count}" for name, count in sorted(class_counts.items())),
    )
    LOGGER.info("Unique UniProt accessions: %d", len(unique_proteins))
    if len(unique_proteins) < args.min_proteins:
        LOGGER.warning(
            "Collected %d proteins, below the requested minimum of %d.",
            len(unique_proteins),
            args.min_proteins,
        )
    if failures:
        LOGGER.warning("%d structures failed measurement; see the manifest.", len(failures))
    return manifest


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point.

    Returns:
        ``0`` on success, ``2`` when the RCSB services are unreachable.
    """
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    try:
        build_dataset(args)
    except RcsbUnavailableError as exc:
        LOGGER.error("Dataset build failed: %s", exc)
        LOGGER.error(
            "This step requires outbound access to search.rcsb.org, data.rcsb.org "
            "and files.rcsb.org. In a restricted environment, run it where those "
            "hosts are reachable and copy data/validation/ across."
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
