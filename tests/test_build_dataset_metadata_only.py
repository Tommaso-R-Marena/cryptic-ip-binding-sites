"""The dataset builder's metadata-only mode (entry table without downloads)."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone

import pandas as pd

from scripts import build_ip_validation_dataset as builder


class _Stats:
    def to_dict(self):
        return {"requests": 1}


class _Client:
    stats = _Stats()

    def fetch_entry_metadata(self, ids):
        return {
            "1ZY7": {
                "rcsb_accession_info": {"initial_release_date": "2005-09-06T00:00:00Z"},
                "rcsb_entry_info": {"resolution_combined": [1.7]},
                "exptl": [{"method": "X-RAY DIFFRACTION"}],
                "polymer_entities": [{
                    "rcsb_polymer_entity_container_identifiers": {
                        "reference_sequence_identifiers": [{"database_name": "UniProt", "database_accession": "P78563"}]
                    },
                    "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Homo sapiens", "ncbi_taxonomy_id": 9606}],
                }],
            }
        }


def test_writes_every_entry_and_reports_missing_metadata(tmp_path):
    args = argparse.Namespace(
        entry_csv=tmp_path / "entries.csv", manifest=tmp_path / "manifest.json",
        experimental_methods=["X-RAY DIFFRACTION"], max_resolution=3.5,
    )
    manifest = builder._write_metadata_only(
        args, _Client(), [], {}, ["1ZY7", "9ZZZ"], [], datetime.now(timezone.utc)
    )
    table = pd.read_csv(args.entry_csv).set_index("pdb_id")
    assert table.loc["1ZY7", "uniprot_ids"] == "P78563"
    assert str(table.loc["1ZY7", "release_date"]).startswith("2005-09-06")
    assert pd.isna(table.loc["9ZZZ", "release_date"])
    assert manifest["entries_without_metadata"] == ["9ZZZ"]
