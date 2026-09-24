#!/usr/bin/env python3
"""Assign every benchmark entry to its sequence and strict homology groups.

Reads the entry table and the structures, searches every protein chain against
every other (MMseqs2; Foldseek for the strict grouping), and writes the entry
table back with ``homology_group`` (sequence) and ``homology_group_strict``
columns, plus a JSON report of the group size distributions. See
:mod:`cryptic_ip.benchmark.homology` and docs/ANALYSIS_PLAN.md, section 3.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path
from typing import List, Optional, Sequence

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.benchmark import homology  # noqa: E402

LOGGER = logging.getLogger("homology_groups")
STRUCTURE_SUFFIXES = (".pdb", ".cif", ".ent", ".mmcif")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entry-csv", type=Path, required=True)
    parser.add_argument("--structures-dir", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--report-json", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, default=Path("results/homology"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--mmseqs-hits", type=Path, default=None, help="Use this hit table instead of running MMseqs2")
    parser.add_argument("--foldseek-hits", type=Path, default=None, help="Use this hit table instead of running Foldseek")
    return parser.parse_args(argv)


def _structure_files(directory: Path, entries: Sequence[str]) -> List[Path]:
    wanted = {entry.upper() for entry in entries}
    return sorted(
        path
        for path in directory.iterdir()
        if path.suffix.lower() in STRUCTURE_SUFFIXES and path.stem.upper() in wanted
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    table = pd.read_csv(args.entry_csv, dtype={"pdb_id": str})
    entries = [str(e).upper() for e in table["pdb_id"]]
    files = _structure_files(args.structures_dir, entries)
    have = {path.stem.upper() for path in files}
    missing = sorted(set(entries) - have)
    if missing:
        LOGGER.warning("%d entries have no structure file and stay singletons: %s", len(missing), missing[:20])

    args.work_dir.mkdir(parents=True, exist_ok=True)
    mmseqs_hits = args.mmseqs_hits
    if mmseqs_hits is None:
        records = {}
        for path in files:
            try:
                records[path.stem.upper()] = homology.chain_sequences(path)
            except Exception as exc:  # recorded: a parse failure must not hide an entry
                LOGGER.warning("%s: no sequences (%s)", path.name, exc)
        fasta = args.work_dir / "chains.fasta"
        n_chains = homology.write_fasta(records, fasta)
        LOGGER.info("Searching %d protein chains from %d entries", n_chains, len(records))
        mmseqs_hits = homology.run_mmseqs(fasta, args.work_dir / "mmseqs", threads=args.threads)
    seq_links = homology.sequence_links(mmseqs_hits)

    foldseek_hits = args.foldseek_hits
    if foldseek_hits is None:
        foldseek_hits = homology.run_foldseek(files, args.work_dir / "foldseek", threads=args.threads)
    struct_links = homology.structure_links(foldseek_hits)

    sequence_groups = homology.connected_groups(entries, seq_links)
    strict_groups = homology.connected_groups(entries, seq_links | struct_links)
    table["homology_group"] = [sequence_groups[e] for e in entries]
    table["homology_group_strict"] = [strict_groups[e] for e in entries]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_csv, index=False)

    report = {
        "thresholds": {
            "min_identity": homology.MIN_IDENTITY,
            "min_shorter_coverage": homology.MIN_SHORTER_COVERAGE,
            "max_evalue": homology.MAX_EVALUE,
            "min_tm_score": homology.MIN_TM_SCORE,
            "min_chain_length": homology.MIN_CHAIN_LENGTH,
        },
        "n_sequence_links": len(seq_links),
        "n_structure_links": len(struct_links),
        "entries_without_structure": missing,
        "sequence": asdict(homology.GroupReport.of(sequence_groups)),
        "strict": asdict(homology.GroupReport.of(strict_groups)),
    }
    args.report_json.parent.mkdir(parents=True, exist_ok=True)
    args.report_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    LOGGER.info("Groups: %s", json.dumps({k: report[k] for k in ("sequence", "strict")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
