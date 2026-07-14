#!/usr/bin/env python3
"""Yeast proteome pilot screen (default: 500 AlphaFold structures)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.analysis import ProteinAnalyzer
from cryptic_ip.analysis.filters import CandidateFilter
from cryptic_ip.database.batch_processing import AlphaFoldBatchDownloader, AnalysisCache, ParallelProcessor

YEAST_PROTEOME_ID = "UP000002311"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("results/yeast_pilot"))
    parser.add_argument("--structures-dir", type=Path, default=Path("data/structures/yeast_pilot"))
    parser.add_argument("--n-proteins", type=int, default=500)
    parser.add_argument("--score-threshold", type=float, default=0.75)
    parser.add_argument("--min-plddt", type=float, default=70.0)
    parser.add_argument("--min-basic", type=int, default=4)
    parser.add_argument("--max-sasa", type=float, default=10.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-electrostatics", action="store_true", default=True)
    parser.add_argument("--with-electrostatics", action="store_true", help="Run APBS per structure (slow)")
    return parser.parse_args()


def _analyze_structure(item: Dict[str, Any]) -> Dict[str, Any]:
    pdb_path = Path(item["filepath"])
    uniprot_id = item["uniprot_id"]
    work_dir = Path(item["work_dir"]) / uniprot_id
    analyzer = ProteinAnalyzer(
        str(pdb_path),
        work_dir=str(work_dir),
        skip_electrostatics=item.get("skip_electrostatics", True),
    )
    scored = analyzer.run_pipeline(include_electrostatics=not item.get("skip_electrostatics", True))
    if scored.empty:
        return {"uniprot_id": uniprot_id, "hits": []}

    filt = CandidateFilter(min_score=item["score_threshold"], min_plddt=item["min_plddt"])
    ranked = filt.filter_cryptic_candidates(
        scored,
        structure_path=str(pdb_path),
        min_basic=item.get("min_basic", 4),
        max_sasa=item.get("max_sasa", 10.0),
    )
    hits = ranked.head(3).to_dict(orient="records")
    for hit in hits:
        hit["uniprot_id"] = uniprot_id
        hit["structure_path"] = str(pdb_path)
    return {"uniprot_id": uniprot_id, "hits": hits}


def main() -> int:
    args = parse_args()
    skip_electrostatics = args.skip_electrostatics and not args.with_electrostatics
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.structures_dir.mkdir(parents=True, exist_ok=True)

    downloader = AlphaFoldBatchDownloader(output_dir=args.structures_dir)
    if not args.skip_download:
        uniprot_ids = downloader.fetch_proteome_uniprot_ids(YEAST_PROTEOME_ID)[: args.n_proteins]
        (args.output_dir / "pilot_uniprot_ids.txt").write_text("\n".join(uniprot_ids), encoding="utf-8")
        for uniprot_id in uniprot_ids:
            try:
                downloader.af_client.fetch_structure(uniprot_id)
            except Exception as exc:
                print(f"Download failed for {uniprot_id}: {exc}")

    pdb_files = sorted(args.structures_dir.glob("AF-*-F1-model_v*.pdb"))[: args.n_proteins]
    if not pdb_files:
        raise RuntimeError(f"No AlphaFold structures found in {args.structures_dir}")

    work_root = args.output_dir / "work"
    items = [
        {
            "uniprot_id": path.name.split("-")[1],
            "filepath": str(path),
            "work_dir": str(work_root),
            "score_threshold": args.score_threshold,
            "min_plddt": args.min_plddt,
            "min_basic": args.min_basic,
            "max_sasa": args.max_sasa,
            "skip_electrostatics": skip_electrostatics,
        }
        for path in pdb_files
    ]

    processor = ParallelProcessor(
        analyze_function=_analyze_structure,
        workers=args.workers,
        checkpoint_path=args.output_dir / "screen_checkpoint.json",
    )
    results = processor.run(items)

    hit_rows: List[Dict[str, Any]] = []
    for result in results:
        hit_rows.extend(result.get("hits", []))

    hits_df = pd.DataFrame(hit_rows)
    if hits_df.empty:
        hits_df = pd.DataFrame(
            columns=[
                "pocket_id",
                "composite_score",
                "basic_residues",
                "sasa",
                "volume",
                "rank",
                "classification",
                "uniprot_id",
                "structure_path",
            ]
        )
    hits_path = args.output_dir / "yeast_pilot_hits.csv"
    hits_df.to_csv(hits_path, index=False)

    proteins_with_hits = int(hits_df["uniprot_id"].nunique()) if not hits_df.empty else 0
    summary = {
        "proteome_id": YEAST_PROTEOME_ID,
        "structures_screened": len(pdb_files),
        "proteins_with_hits": proteins_with_hits,
        "total_hits": int(len(hits_df)),
        "hit_rate": float(proteins_with_hits / len(pdb_files)) if len(pdb_files) else 0.0,
        "score_threshold": args.score_threshold,
        "min_plddt": args.min_plddt,
        "min_basic": args.min_basic,
        "max_sasa": args.max_sasa,
    }
    (args.output_dir / "yeast_pilot_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    cache = AnalysisCache(args.output_dir / "screen_cache.sqlite", pipeline_version="yeast_pilot_v1")
    for result in results:
        if result.get("hits"):
            cache.set_cached_result(result["uniprot_id"], "yeast", result)
    cache.close()

    print(json.dumps(summary, indent=2))
    print(f"Hits written to {hits_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
