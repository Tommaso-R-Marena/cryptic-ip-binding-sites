#!/usr/bin/env python3
"""Generate curated candidate dossiers from a proteome-screen hit table.

Turns raw screen hits into reviewable, manuscript-ready candidate reports by
combining the pocket metrics with UniProt functional annotation and the
pocket-lining basic (Arg/Lys/His) residues that would coordinate an inositol
phosphate.

Usage:
    python scripts/characterize_candidates.py \
        --hits-csv results/yeast_pilot/yeast_pilot_hits.csv \
        --output-dir results/candidates
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.analysis.units import ANGSTROM, ANGSTROM_CU, ANGSTROM_SQ  # noqa: E402
from cryptic_ip.validation.structure_context import BASIC_RESNAMES  # noqa: E402

BASIC_ONE_LETTER = {"ARG": "R", "LYS": "K", "HIS": "H"}


def _to_float(value: Any) -> Optional[float]:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(result) else result


def build_candidate_record(
    hit: Dict[str, Any],
    *,
    rank: int,
    uniprot_info: Optional[Dict[str, Any]] = None,
    coordinating_residues: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Assemble a candidate record from a hit row plus optional lookups (pure)."""
    record: Dict[str, Any] = {
        "rank": rank,
        "uniprot_id": hit.get("uniprot_id"),
        "composite_score": _to_float(hit.get("composite_score")),
        "classification": hit.get("classification", ""),
        "pocket_id": int(hit["pocket_id"]) if _to_float(hit.get("pocket_id")) is not None else None,
        "volume": _to_float(hit.get("volume")),
        "burial_depth": _to_float(hit.get("burial_depth")),
        "sasa": _to_float(hit.get("sasa")),
        "basic_residues": (
            int(hit["basic_residues"]) if _to_float(hit.get("basic_residues")) is not None else None
        ),
        "plddt_confidence": _to_float(hit.get("plddt_confidence")),
        "electrostatic_potential": _to_float(hit.get("electrostatic_potential")),
        "structure_path": hit.get("structure_path"),
    }
    if uniprot_info:
        record.update(
            {
                "gene_name": uniprot_info.get("gene_name", ""),
                "protein_name": uniprot_info.get("protein_name", ""),
                "organism": uniprot_info.get("organism", ""),
                "sequence_length": uniprot_info.get("sequence_length", 0),
                "function": uniprot_info.get("function", ""),
                "subcellular_location": uniprot_info.get("subcellular_location", []),
                "go_molecular_function": (uniprot_info.get("go_terms", {}) or {}).get(
                    "molecular_function", []
                ),
            }
        )
    if coordinating_residues is not None:
        record["coordinating_residues"] = coordinating_residues
    return record


def _fmt(value: Optional[float], digits: int = 2) -> str:
    return f"{value:.{digits}f}" if isinstance(value, (int, float)) else "n/a"


def format_candidate_markdown(record: Dict[str, Any]) -> str:
    """Render a single candidate dossier as Markdown (pure)."""
    name = record.get("protein_name") or record.get("uniprot_id", "Unknown")
    gene = record.get("gene_name") or "?"
    header = f"### {record.get('rank', '?')}. {record.get('uniprot_id', '?')} — {name}"

    lines = [
        header,
        "",
        f"- **Gene / organism**: {gene} / {record.get('organism', 'n/a')}",
        f"- **Composite score**: {_fmt(record.get('composite_score'), 3)} "
        f"({record.get('classification', '')})",
        f"- **Pocket volume**: {_fmt(record.get('volume'), 1)} {ANGSTROM_CU}",
        f"- **Pocket SASA**: {_fmt(record.get('sasa'), 2)} {ANGSTROM_SQ}  "
        f"| **Burial depth**: {_fmt(record.get('burial_depth'), 2)} {ANGSTROM}",
        f"- **Basic residues near pocket**: {record.get('basic_residues', 'n/a')}",
        f"- **AlphaFold pLDDT (pocket)**: {_fmt(record.get('plddt_confidence'), 1)}",
    ]
    if record.get("coordinating_residues"):
        lines.append(
            f"- **Candidate coordinating residues**: {', '.join(record['coordinating_residues'])}"
        )
    if record.get("subcellular_location"):
        lines.append(f"- **Subcellular location**: {', '.join(record['subcellular_location'])}")
    if record.get("function"):
        function = record["function"]
        if len(function) > 400:
            function = function[:400].rstrip() + "…"
        lines += ["", f"> {function}"]
    lines.append("")
    return "\n".join(lines)


def get_coordinating_residues(structure_path: str, pocket_id: int) -> List[str]:
    """Return pocket-lining basic residues (e.g. ['ARG123', 'LYS130']) via re-analysis."""
    from cryptic_ip.analysis import ProteinAnalyzer

    analyzer = ProteinAnalyzer(structure_path, skip_electrostatics=True)
    analyzer.detect_pockets()
    pocket_residue_numbers = set(analyzer.get_pocket_residues(int(pocket_id), distance_cutoff=5.0))

    coordinating: List[tuple[int, str]] = []
    for model in analyzer.structure:
        for chain in model:
            for residue in chain:
                if (
                    residue.id[1] in pocket_residue_numbers
                    and residue.get_resname() in BASIC_RESNAMES
                ):
                    letter = BASIC_ONE_LETTER.get(residue.get_resname(), "?")
                    coordinating.append((int(residue.id[1]), f"{letter}{residue.id[1]}"))
    return [label for _, label in sorted(coordinating)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/candidates"))
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--no-uniprot", action="store_true", help="Skip UniProt annotation lookup")
    parser.add_argument(
        "--no-structure",
        action="store_true",
        help="Skip fpocket re-analysis for coordinating residues",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    hits = pd.read_csv(args.hits_csv)
    if "composite_score" in hits.columns:
        hits = hits.sort_values("composite_score", ascending=False)
    hits = hits.head(args.top_n)

    uniprot_client = None
    if not args.no_uniprot:
        from cryptic_ip.database.uniprot_client import UniProtClient

        uniprot_client = UniProtClient()

    records: List[Dict[str, Any]] = []
    for rank, (_, hit) in enumerate(hits.iterrows(), start=1):
        hit_dict = hit.to_dict()
        uniprot_info = None
        if uniprot_client is not None and hit_dict.get("uniprot_id"):
            try:
                uniprot_info = uniprot_client.get_protein_info(str(hit_dict["uniprot_id"]))
            except Exception as exc:  # noqa: BLE001 - annotation is best-effort
                print(f"UniProt lookup failed for {hit_dict.get('uniprot_id')}: {exc}")

        coordinating = None
        if (
            not args.no_structure
            and hit_dict.get("structure_path")
            and pd.notna(hit_dict.get("pocket_id"))
        ):
            path = Path(str(hit_dict["structure_path"]))
            if path.exists():
                try:
                    coordinating = get_coordinating_residues(str(path), int(hit_dict["pocket_id"]))
                except Exception as exc:  # noqa: BLE001 - structural re-analysis is best-effort
                    print(f"Coordinating-residue analysis failed for {path.name}: {exc}")

        records.append(
            build_candidate_record(
                hit_dict, rank=rank, uniprot_info=uniprot_info, coordinating_residues=coordinating
            )
        )

    summary_df = pd.DataFrame(records)
    summary_csv = args.output_dir / "candidates_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    (args.output_dir / "candidates.json").write_text(
        json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    md_lines = [
        "# Cryptic IP-binding site candidate dossiers",
        "",
        f"Curated from `{args.hits_csv}` ({len(records)} candidate(s)).",
        "",
        "Each candidate is a high-scoring buried pocket flagged for experimental "
        "review (DSF thermal stability ± IP6, site-directed mutagenesis of the "
        "coordinating residues, mass spec for IP occupancy).",
        "",
    ]
    md_lines.extend(format_candidate_markdown(record) for record in records)
    (args.output_dir / "candidate_dossiers.md").write_text("\n".join(md_lines), encoding="utf-8")

    print(f"Wrote {len(records)} candidate dossier(s) to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
