#!/usr/bin/env python3
"""Build a validation dataset for inositol phosphate (IP) binding structures.

This script:
1. Queries the RCSB PDB Search API for X-ray structures containing IP3/IP4/IP5/IP6
2. Downloads structure files and per-ligand metadata
3. Computes ligand SASA using FreeSASA
4. Classifies ligands as cryptic / semi-cryptic / surface
5. Writes a reproducible CSV dataset suitable for methods reporting
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set

import requests

LOGGER = logging.getLogger("build_ip_validation_dataset")

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_CORE_URL = "https://data.rcsb.org/rest/v1/core"
RCSB_FILES_URL = "https://files.rcsb.org/download"
TARGET_LIGANDS = ("IP3", "IP4", "IP5", "IP6")


@dataclass
class DatasetRow:
    """One dataset row for output CSV."""

    pdb_id: str
    uniprot_id: str
    ligand_type: str
    sasa: float
    classification: str
    resolution: Optional[float]
    organism: str


class RcsbClient:
    """Small RCSB API client with retry/rate-limit handling."""

    def __init__(self, delay_s: float = 0.1, max_retries: int = 5, timeout_s: int = 30):
        self.session = requests.Session()
        self.delay_s = delay_s
        self.max_retries = max_retries
        self.timeout_s = timeout_s

    def _request(self, method: str, url: str, **kwargs) -> requests.Response:
        last_error: Optional[Exception] = None
        for attempt in range(1, self.max_retries + 1):
            try:
                response = self.session.request(method, url, timeout=self.timeout_s, **kwargs)
                if response.status_code == 429:
                    retry_after = response.headers.get("Retry-After")
                    wait = float(retry_after) if retry_after else min(2**attempt, 30)
                    LOGGER.warning("Rate limited on %s; sleeping %.1fs", url, wait)
                    time.sleep(wait)
                    continue
                response.raise_for_status()
                time.sleep(self.delay_s)
                return response
            except requests.RequestException as err:
                last_error = err
                wait = min(2**attempt, 30)
                LOGGER.warning("Request error on %s (attempt %d/%d): %s", url, attempt, self.max_retries, err)
                time.sleep(wait)
        raise RuntimeError(f"Request failed after retries: {url}") from last_error

    def get_json(self, url: str) -> Dict:
        return self._request("GET", url).json()

    def post_json(self, url: str, payload: Dict) -> Dict:
        return self._request("POST", url, json=payload).json()

    def download_file(self, url: str, out_path: Path) -> Path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if out_path.exists():
            return out_path
        response = self._request("GET", url)
        out_path.write_bytes(response.content)
        return out_path


def build_search_query(ligands: Sequence[str]) -> Dict:
    """Build RCSB search query for X-ray entries with requested ligands."""
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION",
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_container_identifiers.nonpolymer_comp_ids",
                        "operator": "in",
                        "value": list(ligands),
                    },
                },
            ],
        },
        "request_options": {
            "results_content_type": ["experimental"],
            "return_all_hits": True,
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
        "return_type": "entry",
    }


def query_candidate_pdb_ids(client: RcsbClient, ligands: Sequence[str]) -> List[str]:
    payload = build_search_query(ligands)
    result = client.post_json(RCSB_SEARCH_URL, payload)
    pdb_ids = [x["identifier"].upper() for x in result.get("result_set", [])]
    unique_ids = sorted(set(pdb_ids))
    LOGGER.info("Found %d candidate PDB entries", len(unique_ids))
    return unique_ids


def parse_target_ligands_from_pdb(pdb_path: Path, target_ligands: Set[str]) -> Set[str]:
    found: Set[str] = set()
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("HETATM", "ATOM  ")):
                continue
            resname = line[17:20].strip().upper()
            if resname in target_ligands:
                found.add(resname)
    return found


def ligand_sasa(pdb_path: Path, ligand_id: str) -> float:
    """Compute total SASA for all residues matching ligand_id."""
    try:
        import freesasa
    except ImportError as exc:  # pragma: no cover - import guard
        raise RuntimeError("freesasa is required. Install via conda/pip before running this script.") from exc

    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)
    selection = ("lig", f"resn {ligand_id}")
    selected = freesasa.selectArea((selection,), structure, result)
    return float(selected.get("lig", 0.0))


def classify_sasa(value: float) -> str:
    if value < 5.0:
        return "Cryptic"
    if value <= 20.0:
        return "Semi-cryptic"
    return "Surface"


def extract_uniprot_ids(entry_data: Dict, polymer_entities: Iterable[Dict]) -> str:
    ids: Set[str] = set()
    for entity in polymer_entities:
        refs = entity.get("rcsb_polymer_entity_container_identifiers", {}).get(
            "reference_sequence_identifiers", []
        )
        for ref in refs:
            if ref.get("database_name") == "UniProt":
                accession = ref.get("database_accession")
                if accession:
                    ids.add(accession)
    if ids:
        return ";".join(sorted(ids))

    # fallback to entry-level external references if polymer annotations are absent
    for ref in entry_data.get("struct_ref", []):
        if ref.get("db_name", "").lower() == "uniprot":
            ref_id = ref.get("pdbx_db_accession")
            if ref_id:
                ids.add(ref_id)
    return ";".join(sorted(ids)) if ids else "NA"


def fetch_polymer_entities(client: RcsbClient, pdb_id: str, entry_data: Dict) -> List[Dict]:
    entities: List[Dict] = []
    entity_ids = entry_data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
    for entity_id in entity_ids:
        try:
            url = f"{RCSB_CORE_URL}/polymer_entity/{pdb_id}/{entity_id}"
            entities.append(client.get_json(url))
        except Exception as err:  # noqa: BLE001 - continue on missing entity
            LOGGER.warning("Unable to fetch polymer entity %s/%s: %s", pdb_id, entity_id, err)
    return entities


def fetch_ligand_metadata(client: RcsbClient, ligand_id: str, metadata_dir: Path) -> None:
    out_path = metadata_dir / f"{ligand_id}.json"
    if out_path.exists():
        return
    url = f"{RCSB_CORE_URL}/chemcomp/{ligand_id}"
    try:
        data = client.get_json(url)
        out_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    except Exception as err:  # noqa: BLE001
        LOGGER.warning("Unable to fetch ligand metadata for %s: %s", ligand_id, err)


def build_dataset(output_csv: Path, download_dir: Path, min_proteins: int) -> None:
    client = RcsbClient()
    structures_dir = download_dir / "structures"
    metadata_dir = download_dir / "ligand_metadata"
    structures_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    rows: List[DatasetRow] = []
    unique_proteins: Set[str] = set()

    for pdb_id in query_candidate_pdb_ids(client, TARGET_LIGANDS):
        try:
            structure_path = client.download_file(
                f"{RCSB_FILES_URL}/{pdb_id}.pdb", structures_dir / f"{pdb_id}.pdb"
            )
            entry_data = client.get_json(f"{RCSB_CORE_URL}/entry/{pdb_id}")
            polymer_entities = fetch_polymer_entities(client, pdb_id, entry_data)

            uniprot = extract_uniprot_ids(entry_data, polymer_entities)
            target_ligands = parse_target_ligands_from_pdb(structure_path, set(TARGET_LIGANDS))
            if not target_ligands:
                continue

            resolution = None
            res_list = entry_data.get("rcsb_entry_info", {}).get("resolution_combined", [])
            if res_list:
                resolution = float(res_list[0])
            organism = ""
            src_org = entry_data.get("rcsb_entry_container_identifiers", {}).get(
                "source_organism_scientific_name", []
            )
            if src_org:
                organism = src_org[0]

            for ligand_id in sorted(target_ligands):
                fetch_ligand_metadata(client, ligand_id, metadata_dir)
                sasa = ligand_sasa(structure_path, ligand_id)
                rows.append(
                    DatasetRow(
                        pdb_id=pdb_id,
                        uniprot_id=uniprot,
                        ligand_type=ligand_id,
                        sasa=sasa,
                        classification=classify_sasa(sasa),
                        resolution=resolution,
                        organism=organism or "NA",
                    )
                )

            if uniprot != "NA":
                unique_proteins.update(uniprot.split(";"))
        except Exception as err:  # noqa: BLE001
            LOGGER.warning("Skipping %s due to error: %s", pdb_id, err)

    rows.sort(key=lambda row: (row.pdb_id, row.ligand_type))
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pdb_id",
                "uniprot_id",
                "ligand_type",
                "sasa",
                "classification",
                "resolution",
                "organism",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "pdb_id": row.pdb_id,
                    "uniprot_id": row.uniprot_id,
                    "ligand_type": row.ligand_type,
                    "sasa": f"{row.sasa:.3f}",
                    "classification": row.classification,
                    "resolution": "" if row.resolution is None else f"{row.resolution:.2f}",
                    "organism": row.organism,
                }
            )

    LOGGER.info("Wrote %d rows to %s", len(rows), output_csv)
    LOGGER.info("Unique UniProt proteins: %d", len(unique_proteins))
    if len(unique_proteins) < min_proteins:
        LOGGER.warning(
            "Collected %d proteins, below requested minimum of %d. Consider broadening filters.",
            len(unique_proteins),
            min_proteins,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/validation/ip_binding_validation_dataset.csv"),
        help="Output CSV path",
    )
    parser.add_argument(
        "--download-dir",
        type=Path,
        default=Path("data/validation/raw"),
        help="Directory for downloaded structures and metadata",
    )
    parser.add_argument(
        "--min-proteins",
        type=int,
        default=20,
        help="Warn if fewer unique proteins are collected",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    build_dataset(args.output_csv, args.download_dir, args.min_proteins)


if __name__ == "__main__":
    main()
