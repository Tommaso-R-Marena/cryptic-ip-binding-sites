"""Homology groups: which entries must never be split across train and test.

Two entries are linked when any protein chain of one is homologous to any
protein chain of the other; groups are the connected components of that graph.
Every chain takes part, not only those touching the ligand - a negative pocket
memorised in one fold is as much a leak as a positive one.

Two groupings are built (docs/ANALYSIS_PLAN.md, section 3):

* **sequence** - MMseqs2, at least 30 % identity over at least 50 % of the
  shorter chain, E <= 1e-3;
* **strict** - the sequence links plus Foldseek structural links, alignment
  TM-score >= 0.5 over at least 50 % of the shorter chain, E <= 1e-3, which
  joins homologues too remote for sequence search.

The search tools are run by :func:`run_mmseqs` and :func:`run_foldseek`; the
decisions about what counts as a link are made here, on their tabular output,
so they are testable without the tools.
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np

LOGGER = logging.getLogger(__name__)

MIN_IDENTITY = 0.30
MIN_SHORTER_COVERAGE = 0.50
MAX_EVALUE = 1e-3
MIN_TM_SCORE = 0.50
#: Chains shorter than this carry too little sequence to align meaningfully.
MIN_CHAIN_LENGTH = 20

MMSEQS_COLUMNS = ("query", "target", "fident", "qstart", "qend", "tstart", "tend", "qlen", "tlen", "evalue")
FOLDSEEK_COLUMNS = ("query", "target", "alntmscore", "qstart", "qend", "tstart", "tend", "qlen", "tlen", "evalue")

THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common modified residues, mapped to their parent so homology is not hidden.
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C", "CME": "C", "HYP": "P",
    "MLY": "K", "M3L": "K", "KCX": "K", "LLP": "K", "SEC": "U", "PYL": "O",
}


def entry_of(name: str) -> str:
    """Entry identifier from a chain name written by this module or by Foldseek.

    Names are ``ENTRY|CHAIN`` (FASTA written here) or ``<file>_<chain>`` with
    the file's extensions (Foldseek). PDB identifiers contain neither ``|``,
    ``_`` nor ``.``.
    """
    head = str(name).split("|", 1)[0].split("_", 1)[0].split(".", 1)[0]
    return head.strip().upper()


def chain_sequences(path: Path) -> Dict[str, str]:
    """One-letter sequences of the protein chains of a structure, from C-alpha atoms."""
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    arrays = load_structure_arrays(path)
    mask = arrays.is_polymer & (arrays.atom_names == "CA")
    sequences: Dict[str, List[str]] = {}
    seen: Set[Tuple[str, int, str]] = set()
    for index in np.flatnonzero(mask):
        chain = str(arrays.chain_ids[index])
        key = (chain, int(arrays.resseqs[index]), str(arrays.icodes[index]))
        if key in seen:  # alternate locations
            continue
        seen.add(key)
        sequences.setdefault(chain, []).append(THREE_TO_ONE.get(str(arrays.resnames[index]), "X"))
    return {
        chain: "".join(residues)
        for chain, residues in sequences.items()
        if len(residues) >= MIN_CHAIN_LENGTH and set(residues) != {"X"}
    }


def write_fasta(records: Mapping[str, Mapping[str, str]], path: Path) -> int:
    """Write ``{entry: {chain: sequence}}`` as FASTA named ``ENTRY|CHAIN``; returns the count."""
    n = 0
    with open(path, "w", encoding="utf-8") as handle:
        for entry in sorted(records):
            for chain in sorted(records[entry]):
                handle.write(f">{entry.upper()}|{chain}\n{records[entry][chain]}\n")
                n += 1
    return n


def _shorter_coverage(row: Mapping[str, str]) -> float:
    qlen, tlen = float(row["qlen"]), float(row["tlen"])
    if qlen <= tlen:
        return (float(row["qend"]) - float(row["qstart"]) + 1.0) / qlen
    return (float(row["tend"]) - float(row["tstart"]) + 1.0) / tlen


def _read_hits(path: Path, columns: Sequence[str]) -> Iterable[Dict[str, str]]:
    with open(path, encoding="utf-8") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if not fields:
                continue
            if len(fields) != len(columns):
                raise ValueError(f"{path}: expected {len(columns)} columns, got {len(fields)}")
            yield dict(zip(columns, fields))


def sequence_links(
    hits_path: Path,
    *,
    min_identity: float = MIN_IDENTITY,
    min_coverage: float = MIN_SHORTER_COVERAGE,
    max_evalue: float = MAX_EVALUE,
) -> Set[Tuple[str, str]]:
    """Entry pairs linked by an MMseqs2 hit meeting every threshold."""
    links: Set[Tuple[str, str]] = set()
    for row in _read_hits(hits_path, MMSEQS_COLUMNS):
        a, b = entry_of(row["query"]), entry_of(row["target"])
        if a == b:
            continue
        identity = float(row["fident"])
        identity = identity / 100.0 if identity > 1.0 else identity  # percent or fraction
        if (
            identity >= min_identity
            and float(row["evalue"]) <= max_evalue
            and _shorter_coverage(row) >= min_coverage
        ):
            links.add(tuple(sorted((a, b))))
    return links


def structure_links(
    hits_path: Path,
    *,
    min_tm_score: float = MIN_TM_SCORE,
    min_coverage: float = MIN_SHORTER_COVERAGE,
    max_evalue: float = MAX_EVALUE,
) -> Set[Tuple[str, str]]:
    """Entry pairs linked by a Foldseek hit meeting every threshold."""
    links: Set[Tuple[str, str]] = set()
    for row in _read_hits(hits_path, FOLDSEEK_COLUMNS):
        a, b = entry_of(row["query"]), entry_of(row["target"])
        if a == b:
            continue
        if (
            float(row["alntmscore"]) >= min_tm_score
            and float(row["evalue"]) <= max_evalue
            and _shorter_coverage(row) >= min_coverage
        ):
            links.add(tuple(sorted((a, b))))
    return links


class _UnionFind:
    def __init__(self, items: Iterable[str]) -> None:
        self.parent = {item: item for item in items}

    def find(self, item: str) -> str:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != root:  # path compression
            self.parent[item], item = root, self.parent[item]
        return root

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            # The smaller identifier becomes the root, so group names are stable.
            if rb < ra:
                ra, rb = rb, ra
            self.parent[rb] = ra


def connected_groups(entries: Iterable[str], links: Iterable[Tuple[str, str]]) -> Dict[str, str]:
    """Group name for every entry: ``G:<smallest entry id in its component>``.

    Links naming an entry outside ``entries`` are ignored (it was not kept).
    """
    uf = _UnionFind(sorted({str(e).upper() for e in entries}))
    for a, b in links:
        if a in uf.parent and b in uf.parent:
            uf.union(a, b)
    return {entry: f"G:{uf.find(entry)}" for entry in uf.parent}


@dataclass
class GroupReport:
    """Size distribution of one grouping."""

    n_entries: int
    n_groups: int
    largest_group: int
    largest_group_fraction: float
    n_singletons: int

    @classmethod
    def of(cls, groups: Mapping[str, str]) -> "GroupReport":
        sizes: Dict[str, int] = {}
        for group in groups.values():
            sizes[group] = sizes.get(group, 0) + 1
        largest = max(sizes.values()) if sizes else 0
        return cls(
            n_entries=len(groups),
            n_groups=len(sizes),
            largest_group=largest,
            largest_group_fraction=largest / len(groups) if groups else float("nan"),
            n_singletons=sum(1 for size in sizes.values() if size == 1),
        )


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise RuntimeError(f"{tool} is required for homology grouping and was not found on PATH")
    return path


def run_mmseqs(fasta: Path, work_dir: Path, *, threads: int = 4) -> Path:
    """All-against-all MMseqs2 search; returns the hit table path."""
    work_dir.mkdir(parents=True, exist_ok=True)
    out = work_dir / "mmseqs_hits.tsv"
    subprocess.run(
        [
            _require("mmseqs"), "easy-search", str(fasta), str(fasta), str(out), str(work_dir / "tmp"),
            "-e", str(MAX_EVALUE),
            "-s", "7.5",  # most sensitive preset: remote homologues must not be missed
            "--max-seqs", "100000",  # keep every hit, not the top 300
            "--threads", str(threads),
            "--format-output", ",".join(MMSEQS_COLUMNS),
        ],
        check=True,
    )
    return out


def run_foldseek(structures: Sequence[Path], work_dir: Path, *, threads: int = 4) -> Path:
    """All-against-all Foldseek search over every chain; returns the hit table path."""
    work_dir.mkdir(parents=True, exist_ok=True)
    listing = work_dir / "structures.txt"
    listing.write_text("".join(f"{path}\n" for path in structures), encoding="utf-8")
    db = work_dir / "db"
    foldseek = _require("foldseek")
    subprocess.run(
        [foldseek, "createdb", str(listing), str(db), "--chain-name-mode", "1", "--threads", str(threads)],
        check=True,
    )
    aln, out = work_dir / "aln", work_dir / "foldseek_hits.tsv"
    subprocess.run(
        [
            foldseek, "search", str(db), str(db), str(aln), str(work_dir / "tmp"),
            "-e", str(MAX_EVALUE),
            "-s", "9.5",  # most sensitive preset
            "--max-seqs", "100000",
            "--threads", str(threads),
        ],
        check=True,
    )
    subprocess.run(
        [
            foldseek, "convertalis", str(db), str(db), str(aln), str(out),
            "--format-output", ",".join(FOLDSEEK_COLUMNS),
        ],
        check=True,
    )
    return out

