#!/usr/bin/env python3
"""Phase 2 and Phase 3: catalogue, screen and summarise AlphaFold proteomes.

Subcommands, in pipeline order:

``catalog``
    Measure every model in a proteome directory, write the master CSV and the
    Phase 2 QC report.
``screen``
    Screen one shard of a catalogue. Every pocket is recorded with its full
    descriptor set - not only those that pass - so hit calling can be varied
    afterwards without re-running the structures.
``annotate``
    Fetch UniProt annotations (gene, name, keywords, location, binding sites)
    for a proteome, for the enrichment test and the known-binder check.
``aggregate``
    Combine shards and proteomes: hit rates with exact intervals, a threshold
    sweep, the organism comparison, keyword enrichment, recovery of known
    inositol phosphate binders, and the ranked lists for manual inspection.

Electrostatics use the screened-Coulomb surrogate rather than APBS: an APBS
grid per structure is not affordable at proteome scale. The Phase 1 report
records APBS on the controls.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.analysis.proteome_stats import (  # noqa: E402
    HitCriteria,
    adjusted_comparison,
    hit_rates,
    keyword_enrichment,
    known_ip_annotation,
    pairwise_fisher,
    protein_table,
    rank_percentile,
    threshold_sweep,
)
from cryptic_ip.database.proteome_catalog import (  # noqa: E402
    PROTEOMES,
    build_catalog,
    qc_report,
    shard,
)
from cryptic_ip.utils.json_io import write_json_strict  # noqa: E402

#: Proteins known to bind an inositol phosphate, used to test whether a
#: proteome-wide ranking puts them near the top. Only the ADAR2 site is a
#: buried monomeric site; the HDAC sites sit at the co-repressor interface and
#: are not expected to look buried in a monomer model, and are included so that
#: expectation is tested rather than assumed.
KNOWN_BINDERS: Dict[str, Dict[str, str]] = {
    "P78563": {"gene": "ADARB1 (ADAR2)", "site": "buried InsP6, folding cofactor (1ZY7)"},
    "P55265": {"gene": "ADAR (ADAR1)", "site": "InsP6 pocket conserved with ADAR2"},
    "Q9BUB4": {"gene": "ADAT1", "site": "InsP6 required for activity (no structure)"},
    "Q9NTI5": {"gene": "PDS5B", "site": "InsP6 in insect-cell crystal structure (5HDT)"},
    "Q13547": {"gene": "HDAC1", "site": "InsP4 at MTA1 co-repressor interface (5ICN)"},
    "O15379": {"gene": "HDAC3", "site": "InsP4 at SMRT co-repressor interface (4A69)"},
}

ENRICHMENT_KEYWORDS = (
    "Nucleus",
    "RNA-binding",
    "DNA-binding",
    "Chromatin regulator",
    "Transcription regulation",
    "Mitochondrion",
    "Membrane",
    "Metal-binding",
)

#: Columns kept per pocket. The full descriptor set, so any rescoring is
#: possible later, plus the identity and location of the pocket.
POCKET_COLUMNS_DROP = ("center", "pocket_residue_numbers")


# --------------------------------------------------------------------- catalog
def cmd_catalog(args: argparse.Namespace) -> int:
    directory = Path(args.structures_dir)
    paths = [
        p
        for pattern in ("AF-*.pdb.gz", "AF-*.pdb")
        for p in directory.rglob(pattern)
    ]
    catalog = build_catalog(paths, args.organism)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(out / f"{args.organism}_catalog.csv", index=False)
    report = qc_report(catalog, args.organism)
    if args.source:
        report["source"] = args.source
    write_json_strict(out / f"{args.organism}_qc.json", report, indent=2)
    print(json.dumps(report, indent=2, default=str))
    return 0


# ---------------------------------------------------------------------- screen
def _screen_one(item: Dict[str, str]) -> Dict[str, object]:
    """Screen one model in an isolated working directory."""
    from cryptic_ip.analysis import ProteinAnalyzer

    started = time.time()
    source = Path(item["path"])
    work = Path(tempfile.mkdtemp(prefix=f"screen_{item['uniprot_id']}_"))
    status = {"uniprot_id": item["uniprot_id"], "n_pockets": 0, "seconds": 0.0, "error": ""}
    try:
        pdb_path = work / source.name.replace(".gz", "")
        if source.suffix == ".gz":
            with gzip.open(source, "rb") as src, open(pdb_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
        else:
            shutil.copy(source, pdb_path)
        analyzer = ProteinAnalyzer(str(pdb_path), work_dir=str(work / "work"), skip_electrostatics=True)
        scored = analyzer.run_pipeline(include_electrostatics=False)
        rows = []
        for record in scored.to_dict(orient="records"):
            center = record.get("center") or (np.nan, np.nan, np.nan)
            residues = record.get("pocket_residue_numbers") or []
            for key in POCKET_COLUMNS_DROP:
                record.pop(key, None)
            record["center_x"], record["center_y"], record["center_z"] = (float(c) for c in center)
            record["pocket_residues"] = ",".join(str(r) for r in residues)
            record["uniprot_id"] = item["uniprot_id"]
            rows.append(record)
        status["n_pockets"] = len(rows)
        return {"status": status, "pockets": rows}
    except Exception as exc:  # recorded, never fatal: a proteome has edge cases
        status["error"] = f"{type(exc).__name__}: {exc}"[:500]
        return {"status": status, "pockets": []}
    finally:
        status["seconds"] = round(time.time() - started, 2)
        shutil.rmtree(work, ignore_errors=True)


def cmd_screen(args: argparse.Namespace) -> int:
    catalog = pd.read_csv(args.catalog)
    assigned = shard(catalog, args.shard_index, args.shard_count)
    if args.limit:
        assigned = assigned.head(args.limit)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{args.organism}_shard{args.shard_index:03d}of{args.shard_count:03d}"
    print(f"{tag}: {len(assigned)} models, {args.workers} workers", flush=True)

    items = [
        {"uniprot_id": row.uniprot_id, "path": row.path}
        for row in assigned.itertuples(index=False)
    ]
    statuses: List[Dict[str, object]] = []
    part: List[Dict[str, object]] = []
    part_index = 0
    started = time.time()

    def flush() -> None:
        nonlocal part, part_index
        if part:
            frame = pd.DataFrame(part)
            frame.insert(0, "organism_key", args.organism)
            frame.to_csv(out / f"{tag}_pockets_part{part_index:03d}.csv.gz", index=False)
            part_index += 1
            part = []
        pd.DataFrame(statuses).to_csv(out / f"{tag}_status.csv", index=False)

    def record(done: int, result: Dict[str, object]) -> None:
        statuses.append(result["status"])
        part.extend(result["pockets"])
        if done % 100 == 0 or done == len(items):
            flush()
            rate = done / max(time.time() - started, 1e-9)
            errors = sum(1 for s in statuses if s["error"])
            print(f"{tag}: {done}/{len(items)} ({rate:.2f}/s, {errors} errors)", flush=True)

    if args.workers <= 1:
        for done, item in enumerate(items, start=1):
            record(done, _screen_one(item))
    else:
        # Fresh worker processes periodically, so memory from one giant model
        # is returned to the system rather than carried through the shard.
        with ProcessPoolExecutor(max_workers=args.workers, max_tasks_per_child=25) as pool:
            futures = [pool.submit(_screen_one, item) for item in items]
            for done, future in enumerate(as_completed(futures), start=1):
                record(done, future.result())
    flush()
    return 0


# -------------------------------------------------------------------- annotate
UNIPROT_FIELDS = (
    "accession",
    "gene_primary",
    "protein_name",
    "keyword",
    "cc_subcellular_location",
    "ft_binding",
    "cc_function",
)
ANNOTATION_COLUMNS = (
    "uniprot_id",
    "gene",
    "protein_name",
    "keywords",
    "subcellular_location",
    "binding_site",
    "function",
)


def cmd_annotate(args: argparse.Namespace) -> int:
    info = PROTEOMES[args.organism]
    query = urllib.parse.urlencode(
        {
            "query": f"proteome:{info['proteome_id']}",
            "fields": ",".join(UNIPROT_FIELDS),
            "format": "tsv",
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/stream?{query}"
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    target = out / f"{args.organism}_uniprot.tsv"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=600) as response:
                target.write_bytes(response.read())
            break
        except Exception as exc:
            print(f"UniProt attempt {attempt + 1} failed: {exc}")
            time.sleep(2 ** (attempt + 1))
    else:
        return 1
    frame = pd.read_csv(target, sep="\t", dtype=str)
    frame.columns = list(ANNOTATION_COLUMNS[: len(frame.columns)])
    frame.to_csv(target, sep="\t", index=False)
    print(f"{args.organism}: {len(frame)} UniProt entries -> {target}")
    return 0


# ------------------------------------------------------------------- aggregate
def _read_many(paths: List[Path], **kwargs) -> pd.DataFrame:
    frames = [pd.read_csv(p, **kwargs) for p in paths if p.stat().st_size > 0]
    frames = [f for f in frames if not f.empty]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def cmd_aggregate(args: argparse.Namespace) -> int:
    shards = Path(args.shards_dir)
    catalogs = Path(args.catalog_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    criteria = HitCriteria(min_score=args.min_score)

    pockets = _read_many(sorted(shards.rglob("*_pockets_part*.csv.gz")), low_memory=False)
    status = _read_many(sorted(shards.rglob("*_status.csv")))
    catalog = _read_many(sorted(catalogs.rglob("*_catalog.csv")))
    if status.empty:
        print("no shard status files found")
        return 1
    status["error"] = status["error"].fillna("")
    f1 = catalog[catalog["fragment"] == 1].drop_duplicates("uniprot_id")
    screened = status[status["error"] == ""].merge(
        f1[["uniprot_id", "organism_key", "length", "mean_plddt", "fraction_plddt_70"]],
        on="uniprot_id",
        how="left",
    )
    failed = status[status["error"] != ""]

    proteins = protein_table(pockets, screened, criteria)
    annotations = _read_many(sorted(catalogs.rglob("*_uniprot.tsv")), sep="\t", dtype=str)
    if not annotations.empty:
        annotations = annotations.drop_duplicates("uniprot_id")
        proteins = proteins.merge(
            annotations[["uniprot_id", "gene", "protein_name", "keywords"]], on="uniprot_id", how="left"
        )
        known = known_ip_annotation(annotations)
        proteins["annotated_ip_binder"] = proteins["uniprot_id"].map(known).fillna(False).astype(bool)

    rates = hit_rates(proteins)
    sweep = threshold_sweep(pockets, screened, np.round(np.arange(0.50, 0.951, 0.05), 2), criteria)
    fisher = pairwise_fisher(proteins) if proteins["organism_key"].nunique() > 1 else pd.DataFrame()
    adjusted = adjusted_comparison(proteins) if proteins["organism_key"].nunique() > 1 else {}

    # Candidates: every passing pocket, best first.
    passing = pockets[criteria.passes(pockets)] if not pockets.empty else pockets
    keep = [c for c in (
        "organism_key", "uniprot_id", "pocket_id", "composite_score", "volume", "sasa",
        "basic_residues", "burial_depth", "enclosure", "mean_relative_sasa",
        "coulomb_potential_kt", "plddt_mean", "plddt_min", "center_x", "center_y",
        "center_z", "pocket_residues",
    ) if c in pockets.columns]
    candidates = passing[keep].sort_values("composite_score", ascending=False)
    if "gene" in proteins.columns and not candidates.empty:
        candidates = candidates.merge(
            proteins[["uniprot_id", "gene", "protein_name", "annotated_ip_binder"]],
            on="uniprot_id",
            how="left",
        )
    candidates.to_csv(out / "candidates.csv", index=False)

    # Top 50 per proteome by best confident pocket, for manual inspection
    # whether or not anything clears the strict filter.
    top = (
        proteins.dropna(subset=["best_confident_score"])
        .sort_values("best_confident_score", ascending=False)
        .groupby("organism_key")
        .head(50)
    )
    top.to_csv(out / "top50_per_proteome.csv", index=False)

    enrichment = []
    if not annotations.empty:
        for organism, group in proteins.groupby("organism_key"):
            group = group.reset_index(drop=True)
            eligible = group[group["eligible"]].reset_index(drop=True)
            if eligible.empty:
                continue
            # Hits are likely too few to test; the top 1 % of eligible proteins
            # by best confident pocket score is a ranking-based set that always
            # exists, and is reported alongside.
            cutoff = eligible["best_confident_score"].quantile(0.99)
            for label, mask, base in (
                ("strict hits", eligible["is_hit"], eligible),
                ("top 1% by score", eligible["best_confident_score"] >= cutoff, eligible),
            ):
                if int(mask.sum()) == 0:
                    continue
                table = keyword_enrichment(base, annotations, mask, ENRICHMENT_KEYWORDS)
                table.insert(0, "set", label)
                table.insert(0, "organism_key", organism)
                enrichment.append(table)
    enrichment_df = pd.concat(enrichment, ignore_index=True) if enrichment else pd.DataFrame()

    known_rows = []
    for uniprot_id, meta in KNOWN_BINDERS.items():
        row = proteins[proteins["uniprot_id"] == uniprot_id]
        if row.empty:
            known_rows.append({"uniprot_id": uniprot_id, **meta, "screened": False})
            continue
        r = row.iloc[0]
        known_rows.append(
            {
                "uniprot_id": uniprot_id,
                **meta,
                "screened": True,
                "best_confident_score": r["best_confident_score"],
                "rank_percentile": rank_percentile(proteins, uniprot_id),
                "is_hit": bool(r["is_hit"]),
                "blocking_gate": r["blocking_gate"],
            }
        )
    known_df = pd.DataFrame(known_rows)

    proteins.to_csv(out / "proteins.csv.gz", index=False)
    rates.to_csv(out / "hit_rates.csv", index=False)
    sweep.to_csv(out / "threshold_sweep.csv", index=False)
    fisher.to_csv(out / "pairwise_fisher.csv", index=False)
    enrichment_df.to_csv(out / "keyword_enrichment.csv", index=False)
    known_df.to_csv(out / "known_binders.csv", index=False)
    failed.to_csv(out / "failed_structures.csv", index=False)

    ip6 = {k: v["ip6_uM"] for k, v in PROTEOMES.items()}
    summary = {
        "criteria": criteria.__dict__,
        "structures": {
            "screened_ok": int(len(screened)),
            "failed": int(len(failed)),
            "pockets": int(len(pockets)),
            "median_seconds_per_structure": float(status["seconds"].median()),
        },
        "hit_rates": rates.to_dict(orient="records"),
        "ip6_uM": ip6,
        "pairwise_fisher": fisher.to_dict(orient="records"),
        "adjusted_comparison": adjusted,
        "blocking_gate_counts": {
            organism: group["blocking_gate"].value_counts().to_dict()
            for organism, group in proteins.groupby("organism_key")
        },
        "known_binders": known_df.to_dict(orient="records"),
    }
    write_json_strict(out / "screen_summary.json", summary, indent=2)
    _digest(summary, rates, sweep, fisher, adjusted, candidates, top, known_df, enrichment_df, out)
    return 0


def _pct(x) -> str:
    return "n/a" if x is None or pd.isna(x) else f"{100 * float(x):.2f}%"


def _digest(summary, rates, sweep, fisher, adjusted, candidates, top, known_df, enrichment_df, out) -> None:
    lines = ["## Proteome screen\n"]
    s = summary["structures"]
    lines.append(
        f"{s['screened_ok']} structures screened ({s['failed']} failed), {s['pockets']} pockets, "
        f"median {s['median_seconds_per_structure']:.1f} s per structure.\n"
    )
    lines.append("### Hit rates (strict criteria)\n")
    lines.append("| proteome | IP6 (uM) | screened | hits | rate [95% CI] | eligible | rate among eligible [95% CI] |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for r in rates.to_dict(orient="records"):
        lines.append(
            f"| {r['organism_key']} | {summary['ip6_uM'].get(r['organism_key'], '?')} | {r['screened']} | "
            f"{r['hits']} | {_pct(r['hit_rate'])} [{_pct(r['ci_low'])}-{_pct(r['ci_high'])}] | "
            f"{r['eligible']} | {_pct(r['hit_rate_eligible'])} "
            f"[{_pct(r['ci_low_eligible'])}-{_pct(r['ci_high_eligible'])}] |"
        )
    if not sweep.empty:
        lines.append("\n### Hits across score thresholds (other gates fixed)\n")
        pivot = sweep.pivot(index="min_score", columns="organism_key", values="hits")
        lines.append(pivot.to_markdown())
    if not fisher.empty:
        lines.append("\n### Dictyostelium vs the others (Fisher exact)\n")
        lines.append(fisher.to_markdown(index=False, floatfmt=".3g"))
    if adjusted:
        lines.append("\n### Adjusted for length and model confidence\n")
        lines.append("```\n" + json.dumps(adjusted, indent=2, default=str)[:4000] + "\n```")
    lines.append("\n### Why proteins were not hits (closest pocket's first failing gate)\n")
    lines.append("```\n" + json.dumps(summary["blocking_gate_counts"], indent=2) + "\n```")
    lines.append("\n### Known inositol phosphate binders\n")
    lines.append(known_df.to_markdown(index=False, floatfmt=".3f") if not known_df.empty else "none")
    lines.append("\n### Strict candidates\n")
    lines.append(candidates.head(40).to_markdown(index=False, floatfmt=".3f") if not candidates.empty else "none")
    lines.append("\n### Top 10 per proteome by best confident pocket (for inspection)\n")
    cols = [c for c in ("organism_key", "uniprot_id", "gene", "best_confident_score", "is_hit",
                        "blocking_gate", "length", "annotated_ip_binder") if c in top.columns]
    lines.append(top.groupby("organism_key").head(10)[cols].to_markdown(index=False, floatfmt=".3f"))
    if not enrichment_df.empty:
        lines.append("\n### Keyword enrichment\n")
        lines.append(enrichment_df.to_markdown(index=False, floatfmt=".3g"))
    text = "\n".join(lines) + "\n"
    (out / "DIGEST.md").write_text(text, encoding="utf-8")
    print(text)


# ------------------------------------------------------------------------ main
def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("catalog")
    p.add_argument("--organism", required=True, choices=sorted(PROTEOMES))
    p.add_argument("--structures-dir", required=True)
    p.add_argument("--output-dir", default="results/proteome_screen/catalog")
    p.add_argument("--source", default="", help="Where the models came from, recorded in the QC report")
    p.set_defaults(func=cmd_catalog)

    p = sub.add_parser("screen")
    p.add_argument("--organism", required=True, choices=sorted(PROTEOMES))
    p.add_argument("--catalog", required=True)
    p.add_argument("--shard-index", type=int, default=0)
    p.add_argument("--shard-count", type=int, default=1)
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--limit", type=int, default=0, help="Screen only the first N (testing)")
    p.add_argument("--output-dir", default="results/proteome_screen/shards")
    p.set_defaults(func=cmd_screen)

    p = sub.add_parser("annotate")
    p.add_argument("--organism", required=True, choices=sorted(PROTEOMES))
    p.add_argument("--output-dir", default="results/proteome_screen/catalog")
    p.set_defaults(func=cmd_annotate)

    p = sub.add_parser("aggregate")
    p.add_argument("--shards-dir", default="results/proteome_screen/shards")
    p.add_argument("--catalog-dir", default="results/proteome_screen/catalog")
    p.add_argument("--output-dir", default="results/proteome_screen/summary")
    p.add_argument("--min-score", type=float, default=0.75)
    p.set_defaults(func=cmd_aggregate)
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
