#!/usr/bin/env python3
"""Turn the project's committed result tables into study inputs.

Two products, and which one you get depends on what the source actually
carries:

``architecture``
    One record per pocket with its lining **positions**, the chain they sit on,
    a grouping key and any covariates. Every statistic in
    `cryptic_ip/detect/architecture.py` is a function of positions alone, so
    this file is sufficient for study Q7 and is emitted whenever a table has a
    residue list.

``sequence``
    The same records with a chain **sequence** attached, which study Q6 needs
    because its statistics are functions of the letters. Emitted only when
    sequences are supplied. There is no fallback and no placeholder: a
    sequence-inputs file with invented letters would produce a study result
    that means nothing, so without sequences the file is simply not written.

Two refusals are enforced here rather than written in a document.

**Chain-blind positions are refused for multi-chain sources.** The repo's
pocket residue numbers come from `PocketAnalyzer.get_pocket_residues`, whose
own docstring says they are chain-blind. On a single-chain AlphaFold model that
costs nothing because there is one chain to be blind about; on a PDB entry with
several chains the same number names several residues and a window built from
it silently mixes them. The builder therefore needs either a chain column or an
explicit declaration that the source is single-chain, and refuses otherwise.

**A missing column stops the run.** Column names are passed in and verified
against the file, never guessed, because a builder that quietly skips rows it
cannot parse produces a study that goes green on less data than it claims.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

#: Sources whose structures carry exactly one polymer chain, so a chain-blind
#: residue number is unambiguous. AlphaFold models are one chain by
#: construction; anything else has to say so or supply a chain column.
SINGLE_CHAIN_SOURCES = ("alphafold",)

#: The chain label used when a source is declared single-chain. It is recorded
#: in every record so a downstream reader never has to assume it.
IMPLICIT_CHAIN = "A"


class InputError(ValueError):
    """A source table from which honest inputs cannot be built."""


@dataclass(frozen=True)
class Columns:
    """Which column holds what, named explicitly rather than inferred."""

    identifier: str
    residues: str
    group: str
    chain: Optional[str] = None
    label: Optional[str] = None
    covariates: Sequence[str] = ()

    def required(self) -> List[str]:
        names = [self.identifier, self.residues, self.group]
        for optional in (self.chain, self.label):
            if optional:
                names.append(optional)
        names.extend(self.covariates)
        return names


def parse_residues(value: object) -> List[int]:
    """Residue numbers from a delimited cell, in the format the repo writes.

    `proteome_screen.py` joins them with commas; other tables use spaces or
    semicolons. Negative and insertion-coded numbers are not parsed: a residue
    number this cannot read is an error rather than a silently dropped residue.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    tokens = [t for t in re.split(r"[,;\s]+", text) if t]
    out: List[int] = []
    for token in tokens:
        if not re.fullmatch(r"-?\d+", token):
            raise InputError(
                f"cannot read {token!r} as a residue number; the cell was {text[:60]!r}"
            )
        out.append(int(token))
    return sorted(set(out))


def verify_columns(frame: pd.DataFrame, columns: Columns) -> None:
    missing = [name for name in columns.required() if name not in frame.columns]
    if missing:
        raise InputError(
            f"the table is missing {missing}; it has {sorted(frame.columns)[:12]}..."
        )


def build_records(
    frame: pd.DataFrame,
    columns: Columns,
    *,
    source: str,
    single_chain: bool,
    min_positions: int = 2,
) -> Dict[str, object]:
    """One record per pocket, with the reason for every row that is dropped."""
    verify_columns(frame, columns)
    if not single_chain and not columns.chain:
        raise InputError(
            "this source is not declared single-chain and no chain column was "
            "given; the repo's pocket residue numbers are chain-blind, so a "
            "window built from them would mix chains"
        )
    records: List[Dict[str, object]] = []
    dropped: Dict[str, int] = {}

    def drop(reason: str) -> None:
        dropped[reason] = dropped.get(reason, 0) + 1

    for _index, row in frame.iterrows():
        positions = parse_residues(row[columns.residues])
        if len(positions) < min_positions:
            drop(f"fewer than {min_positions} lining positions")
            continue
        chain = str(row[columns.chain]) if columns.chain else IMPLICIT_CHAIN
        group = row[columns.group]
        if group is None or (isinstance(group, float) and pd.isna(group)):
            drop("no grouping key")
            continue
        record: Dict[str, object] = {
            "id": str(row[columns.identifier]),
            "group": str(group),
            "chain": chain,
            "chain_is_implicit": bool(columns.chain is None),
            "source": source,
            "residues": [[chain, int(p)] for p in positions],
            "n_positions": len(positions),
            "span": int(positions[-1] - positions[0]),
        }
        if columns.label:
            record["label"] = str(row[columns.label])
        for name in columns.covariates:
            value = row[name]
            record[name] = None if pd.isna(value) else (
                float(value) if isinstance(value, (int, float)) else str(value))
        records.append(record)
    return {
        "source": source,
        "single_chain": single_chain,
        "n_rows": int(len(frame)),
        "n_records": len(records),
        "dropped": dict(sorted(dropped.items())),
        "records": records,
    }


def attach_sequences(payload: Dict[str, object], sequences: Dict[str, str]) -> Dict[str, object]:
    """Add a chain sequence to each record, dropping any record without one.

    A record whose sequence is missing is dropped with a reason rather than
    carried with a placeholder, and a record whose positions fall outside its
    sequence is dropped too - that is the signature of a numbering offset, and
    guessing the offset is how a whole study ends up shifted by a constant.
    """
    kept, dropped = [], dict(payload.get("dropped") or {})

    def drop(reason: str) -> None:
        dropped[reason] = dropped.get(reason, 0) + 1

    for record in payload["records"]:  # type: ignore[union-attr]
        sequence = sequences.get(str(record["id"]))
        if not sequence:
            drop("no sequence for this identifier")
            continue
        positions = [int(p) for _chain, p in record["residues"]]  # type: ignore[union-attr]
        if min(positions) < 1 or max(positions) > len(sequence):
            drop("lining positions fall outside the sequence")
            continue
        kept.append({**record, "sequence": sequence, "offset": 1})
    return {**payload, "records": kept, "n_records": len(kept), "dropped": dropped}


def read_fasta(path: Path) -> Dict[str, str]:
    """Identifier to sequence, keyed on the first whitespace-delimited token."""
    sequences: Dict[str, str] = {}
    name, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        elif name:
            chunks.append(line.strip())
    if name:
        sequences[name] = "".join(chunks)
    return sequences


def write_jsonl(records: Iterable[Dict[str, object]], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
            count += 1
    return count


#: The committed tables this builder knows, with their real column names.
#: Verified against the files in the repository rather than assumed, and a
#: table whose columns have moved will fail `verify_columns` loudly.
KNOWN = {
    "specificity-candidates": {
        "path": "results/specificity/candidates.csv",
        "columns": Columns(
            identifier="uniprot_id", residues="top_pocket_residues",
            group="cluster", label="explained",
            covariates=("plddt_mean", "annotated", "learned_score")),
        "single_chain": True,
        "note": "AlphaFold proteome-screen candidates: one chain per model",
    },
}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--table", type=Path, default=None,
                        help="a CSV to read; omit to use --known")
    parser.add_argument("--known", choices=sorted(KNOWN), default=None,
                        help="one of the committed tables, with its verified columns")
    parser.add_argument("--identifier", default=None)
    parser.add_argument("--residues", default=None)
    parser.add_argument("--group", default=None)
    parser.add_argument("--chain", default=None)
    parser.add_argument("--label", default=None)
    parser.add_argument("--covariate", action="append", default=[])
    parser.add_argument("--single-chain", action="store_true",
                        help="declare that every structure has one polymer chain")
    parser.add_argument("--sequences", type=Path, default=None,
                        help="FASTA; without it no sequence-inputs file is written")
    parser.add_argument("--architecture-out", type=Path, required=True)
    parser.add_argument("--sequence-out", type=Path, default=None)
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args(argv)

    if args.known:
        spec = KNOWN[args.known]
        table = Path(spec["path"])  # type: ignore[arg-type]
        columns = spec["columns"]
        single_chain = bool(spec["single_chain"])
        source = args.known
    else:
        if not args.table:
            raise SystemExit("give --table or --known")
        missing = [n for n in ("identifier", "residues", "group") if not getattr(args, n)]
        if missing:
            raise SystemExit(f"--{' and --'.join(missing)} are required with --table")
        table = args.table
        columns = Columns(identifier=args.identifier, residues=args.residues,
                          group=args.group, chain=args.chain, label=args.label,
                          covariates=tuple(args.covariate))
        single_chain = args.single_chain
        source = table.stem

    frame = pd.read_csv(table)
    try:
        payload = build_records(frame, columns, source=source, single_chain=single_chain)
    except InputError as error:
        raise SystemExit(f"{table}: {error}") from error

    written = write_jsonl(payload["records"], args.architecture_out)  # type: ignore[arg-type]
    print(f"architecture inputs: {written} records -> {args.architecture_out}")

    if args.sequences:
        with_sequences = attach_sequences(payload, read_fasta(args.sequences))
        target = args.sequence_out or args.architecture_out.with_name("sequence_inputs.jsonl")
        count = write_jsonl(with_sequences["records"], target)  # type: ignore[arg-type]
        print(f"sequence inputs: {count} records -> {target}")
        payload = {**payload, "sequence_records": count,
                   "dropped": with_sequences["dropped"]}
    else:
        print("no --sequences given, so no sequence-inputs file was written; "
              "study Q6 needs letters and none are invented here")

    summary = {k: v for k, v in payload.items() if k != "records"}
    summary["n_architecture_records"] = written
    if args.report:
        args.report.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
