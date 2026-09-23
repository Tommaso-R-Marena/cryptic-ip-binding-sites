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
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.analysis.proteome_stats import (  # noqa: E402
    CALIBRATED_CRITERIA,
    PLAN_CRITERIA,
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
from cryptic_ip.database.async_fetch import FetchJob, fetch_all  # noqa: E402
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
    from cryptic_ip.analysis.geometry import hull_depths
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    started = time.time()
    source = Path(item["path"])
    work = Path(tempfile.mkdtemp(prefix=f"screen_{item['uniprot_id']}_"))
    status = {
        "uniprot_id": item["uniprot_id"],
        "n_pockets": 0,
        "seconds": 0.0,
        "fpocket_seconds": float("nan"),
        "error": "",
    }
    try:
        pdb_path = work / source.name.replace(".gz", "")
        if source.suffix == ".gz":
            with gzip.open(source, "rb") as src, open(pdb_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
        else:
            shutil.copy(source, pdb_path)
        analyzer = ProteinAnalyzer(str(pdb_path), work_dir=str(work / "work"), skip_electrostatics=True)
        # run_pipeline without electrostatics, with fpocket timed on its own:
        # the split says whether a slow shard is pocket detection or descriptors.
        detect_started = time.time()
        analyzer.detect_pockets()
        status["fpocket_seconds"] = round(time.time() - detect_started, 2)
        scored = analyzer.score_all_pockets()
        # Each pocket's hull depth is one of its descriptors. The deepest
        # atom's hull depth is kept as the protein's own scale, so a pocket's
        # depth can be read relative to the size of the protein it sits in.
        arrays = load_structure_arrays(pdb_path)
        protein = arrays.coords[arrays.is_polymer]
        try:
            protein_scale = float(hull_depths(protein, protein).max())
        except Exception:  # degenerate (e.g. planar) coordinates have no 3-D hull
            protein_scale = float("nan")
        rows = []
        for record in scored.to_dict(orient="records"):
            center = record.get("center") or (np.nan, np.nan, np.nan)
            residues = record.get("pocket_residue_numbers") or []
            for key in POCKET_COLUMNS_DROP:
                record.pop(key, None)
            record["center_x"], record["center_y"], record["center_z"] = (float(c) for c in center)
            record["pocket_residues"] = ",".join(str(r) for r in residues)
            record["uniprot_id"] = item["uniprot_id"]
            record["protein_max_hull_depth"] = protein_scale
            rows.append(record)
        status["n_pockets"] = len(rows)
        return {"status": status, "pockets": rows}
    except Exception as exc:  # recorded, never fatal: a proteome has edge cases
        status["error"] = f"{type(exc).__name__}: {exc}"[:500]
        return {"status": status, "pockets": []}
    finally:
        status["seconds"] = round(time.time() - started, 2)
        shutil.rmtree(work, ignore_errors=True)


#: Report when no model has finished for this long, naming those in flight.
HEARTBEAT_S = 600.0
#: Log every model that takes at least this long.
SLOW_S = 120.0


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
        status = result["status"]
        statuses.append(status)
        part.extend(result["pockets"])
        if status["seconds"] >= SLOW_S:
            print(f"{tag}: slow model {status['uniprot_id']}: {status['seconds']:.0f} s "
                  f"(fpocket {status['fpocket_seconds']:.0f} s), {status['n_pockets']} pockets "
                  f"{status['error']}", flush=True)
        if done % 50 == 0 or done == len(items):
            flush()
            rate = done / max(time.time() - started, 1e-9)
            errors = sum(1 for s in statuses if s["error"])
            print(f"{tag}: {done}/{len(items)} ({rate:.2f}/s, {errors} errors)", flush=True)

    if args.workers <= 1:
        for done, item in enumerate(items, start=1):
            record(done, _screen_one(item))
    else:
        # Fresh worker processes periodically, so memory from one giant model
        # is returned to the system rather than carried through the shard
        # (Python 3.11+; earlier versions keep their workers).
        recycle = {"max_tasks_per_child": 25} if sys.version_info >= (3, 11) else {}
        with ProcessPoolExecutor(max_workers=args.workers, **recycle) as pool:
            futures = {pool.submit(_screen_one, item): item for item in items}
            pending = set(futures)
            done = 0
            while pending:
                finished, pending = wait(pending, timeout=HEARTBEAT_S, return_when=FIRST_COMPLETED)
                if not finished:
                    # Name what is running, so a stall is attributable to models.
                    running = [futures[f]["uniprot_id"] for f in pending if f.running()]
                    print(
                        f"{tag}: nothing finished in {HEARTBEAT_S:.0f} s; "
                        f"{done}/{len(items)} done; running {', '.join(running)}",
                        flush=True,
                    )
                for future in finished:
                    done += 1
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
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
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
            # Gzip in transit: the human annotation table is ~100 MB as text.
            "compressed": "true",
        }
    )
    url = f"{UNIPROT_STREAM}?{query}"
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    target = out / f"{args.organism}_uniprot.tsv"
    # Retried with backoff, decompressed and checked before being written.
    (result,) = fetch_all([FetchJob(args.organism, url, target)], concurrency=1, timeout=900)
    if not result.ok:
        print(f"UniProt download failed: {result.error}")
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


def _analyse(
    label: str,
    criteria: HitCriteria,
    pockets: pd.DataFrame,
    screened: pd.DataFrame,
    annotations: pd.DataFrame,
    out: Path,
) -> Dict[str, object]:
    """Hit calling and every downstream comparison, under one hit definition."""
    proteins = protein_table(pockets, screened, criteria)
    if not annotations.empty:
        proteins = proteins.merge(
            annotations[["uniprot_id", "gene", "protein_name", "keywords"]], on="uniprot_id", how="left"
        )
        known = known_ip_annotation(annotations)
        proteins["annotated_ip_binder"] = proteins["uniprot_id"].map(known).fillna(False).astype(bool)
    if "hull_depth" in pockets.columns and not pockets.empty:
        confident = pockets[pockets["plddt_mean"].astype(float).fillna(0) >= criteria.min_plddt]
        proteins["max_confident_hull_depth"] = (
            proteins["uniprot_id"].map(confident.groupby("uniprot_id")["hull_depth"].max())
        )

    multi = proteins["organism_key"].nunique() > 1
    rates = hit_rates(proteins)
    sweep = threshold_sweep(pockets, screened, np.round(np.arange(0.50, 0.951, 0.05), 2), criteria)
    fisher = pairwise_fisher(proteins) if multi else pd.DataFrame()
    adjusted = adjusted_comparison(proteins) if multi else {}

    passing = pockets[criteria.passes(pockets)] if not pockets.empty else pockets
    keep = [c for c in (
        "organism_key", "uniprot_id", "pocket_id", "composite_score", "volume", "sasa",
        "basic_residues", "burial_depth", "hull_depth", "enclosure", "mean_relative_sasa",
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

    # Top 50 per proteome by best confident pocket, for manual inspection
    # whether or not anything clears the filter.
    top = (
        proteins.dropna(subset=["best_confident_score"])
        .sort_values("best_confident_score", ascending=False)
        .groupby("organism_key")
        .head(50)
    )

    enrichment = []
    if not annotations.empty:
        for organism, group in proteins.groupby("organism_key"):
            eligible = group[group["eligible"]].reset_index(drop=True)
            if eligible.empty:
                continue
            # Hits may be too few to test; the top 1 % of eligible proteins by
            # best confident pocket score always exists and is reported beside.
            cutoff = eligible["best_confident_score"].quantile(0.99)
            for set_label, mask in (
                ("hits", eligible["is_hit"]),
                ("top 1% by score", eligible["best_confident_score"] >= cutoff),
            ):
                if int(mask.sum()) == 0:
                    continue
                table = keyword_enrichment(eligible, annotations, mask, ENRICHMENT_KEYWORDS)
                table.insert(0, "set", set_label)
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
                "max_confident_hull_depth": r.get("max_confident_hull_depth"),
                "is_hit": bool(r["is_hit"]),
                "blocking_gate": r["blocking_gate"],
            }
        )
    known_df = pd.DataFrame(known_rows)

    target = out / label
    target.mkdir(parents=True, exist_ok=True)
    proteins.to_csv(target / "proteins.csv.gz", index=False)
    rates.to_csv(target / "hit_rates.csv", index=False)
    sweep.to_csv(target / "threshold_sweep.csv", index=False)
    fisher.to_csv(target / "pairwise_fisher.csv", index=False)
    enrichment_df.to_csv(target / "keyword_enrichment.csv", index=False)
    known_df.to_csv(target / "known_binders.csv", index=False)
    candidates.to_csv(target / "candidates.csv", index=False)
    top.to_csv(target / "top50_per_proteome.csv", index=False)

    return {
        "label": label,
        "criteria": dict(criteria.__dict__),
        "proteins": proteins,
        "rates": rates,
        "sweep": sweep,
        "fisher": fisher,
        "adjusted": adjusted,
        "candidates": candidates,
        "top": top,
        "known": known_df,
        "enrichment": enrichment_df,
    }


def cmd_aggregate(args: argparse.Namespace) -> int:
    shards = Path(args.shards_dir)
    catalogs = Path(args.catalog_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    pockets = _read_many(sorted(shards.rglob("*_pockets_part*.csv.gz")), low_memory=False)
    scorer_parameters = None
    if not args.keep_screen_scores and not pockets.empty:
        # Every pocket is stored with its full descriptor set, so it is scored
        # here with the current scorer rather than trusted from the screen: a
        # scorer change then needs a re-aggregation, never a re-screen. The
        # screen's own score is kept alongside.
        from cryptic_ip.analysis.scorer import PocketScorer

        scorer = PocketScorer()
        pockets["composite_score_screen"] = pockets.get("composite_score")
        pockets["composite_score"] = scorer.score_frame(pockets)
        scorer_parameters = dict(scorer.parameters.__dict__)
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
    failed.to_csv(out / "failed_structures.csv", index=False)
    # Models a shard never reached (it timed out, or its job failed) are
    # neither screened nor failed: they are reported, and left out of every
    # denominator, rather than disappearing without a count.
    eligible = f1[f1["qc_pass"].astype(bool)] if "qc_pass" in f1.columns else f1
    unscreened = eligible[~eligible["uniprot_id"].isin(status["uniprot_id"])]
    unscreened[["uniprot_id", "organism_key", "length"]].to_csv(out / "unscreened_structures.csv", index=False)
    annotations = _read_many(sorted(catalogs.rglob("*_uniprot.tsv")), sep="\t", dtype=str)
    if not annotations.empty:
        annotations = annotations.drop_duplicates("uniprot_id")

    analyses = [_analyse("plan", PLAN_CRITERIA.replace(min_score=args.min_score), pockets, screened, annotations, out)]
    if "hull_depth" in pockets.columns:
        analyses.append(_analyse("calibrated", CALIBRATED_CRITERIA, pockets, screened, annotations, out))

    ip6 = {k: v["ip6_uM"] for k, v in PROTEOMES.items()}
    summary = {
        "structures": {
            "screened_ok": int(len(screened)),
            "failed": int(len(failed)),
            "unscreened": {k: int(v) for k, v in unscreened["organism_key"].value_counts().items()},
            "pockets": int(len(pockets)),
            "median_seconds_per_structure": float(status["seconds"].median()),
        },
        "scoring": (
            {"rescored": True, "parameters": scorer_parameters}
            if scorer_parameters is not None
            else {"rescored": False}
        ),
        "ip6_uM": ip6,
    }
    for analysis in analyses:
        summary[analysis["label"]] = {
            "criteria": analysis["criteria"],
            "hit_rates": analysis["rates"].to_dict(orient="records"),
            "pairwise_fisher": analysis["fisher"].to_dict(orient="records"),
            "adjusted_comparison": analysis["adjusted"],
            "blocking_gate_counts": {
                organism: group["blocking_gate"].value_counts().to_dict()
                for organism, group in analysis["proteins"].groupby("organism_key")
            },
            "known_binders": analysis["known"].to_dict(orient="records"),
        }
    write_json_strict(out / "screen_summary.json", summary, indent=2)
    _digest(summary, analyses, out)
    return 0


def _pct(x) -> str:
    return "n/a" if x is None or pd.isna(x) else f"{100 * float(x):.2f}%"


def _digest(summary, analyses, out) -> None:
    lines = ["## Proteome screen\n"]
    s = summary["structures"]
    lines.append(
        f"{s['screened_ok']} structures screened ({s['failed']} failed), {s['pockets']} pockets, "
        f"median {s['median_seconds_per_structure']:.1f} s per structure.\n"
    )
    scoring = summary.get("scoring", {})
    if scoring.get("rescored"):
        lines.append(
            "Pocket scores recomputed at aggregation from the stored descriptors "
            f"(depth measure: {scoring['parameters'].get('depth_measure')}).\n"
        )
    if s["unscreened"]:
        missing = ", ".join(f"{k} {v}" for k, v in sorted(s["unscreened"].items()))
        lines.append(
            f"**Incomplete:** QC-passing models never screened ({missing}); "
            "rates below are over the models that were.\n"
        )
    for analysis in analyses:
        c = analysis["criteria"]
        lines.append(f"## Hit definition: {analysis['label']}\n")
        lines.append(
            f"score >= {c['min_score']}, lining SASA <= {c['max_sasa']}, basic >= {c['min_basic']}, "
            f"volume {c['min_volume']:.0f}-{c['max_volume']:.0f}, pLDDT >= {c['min_plddt']}, "
            f"hull depth >= {c['min_hull_depth']}\n"
        )
        lines.append("| proteome | IP6 (uM) | screened | hits | rate [95% CI] | eligible | rate among eligible [95% CI] |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for r in analysis["rates"].to_dict(orient="records"):
            lines.append(
                f"| {r['organism_key']} | {summary['ip6_uM'].get(r['organism_key'], '?')} | {r['screened']} | "
                f"{r['hits']} | {_pct(r['hit_rate'])} [{_pct(r['ci_low'])}-{_pct(r['ci_high'])}] | "
                f"{r['eligible']} | {_pct(r['hit_rate_eligible'])} "
                f"[{_pct(r['ci_low_eligible'])}-{_pct(r['ci_high_eligible'])}] |"
            )
        sweep = analysis["sweep"]
        if not sweep.empty:
            lines.append("\n#### Hits across score thresholds (other gates fixed)\n")
            lines.append(sweep.pivot(index="min_score", columns="organism_key", values="hits").to_markdown())
        if not analysis["fisher"].empty:
            lines.append("\n#### Dictyostelium vs the others (Fisher exact)\n")
            lines.append(analysis["fisher"].to_markdown(index=False, floatfmt=".3g"))
        if analysis["adjusted"]:
            lines.append("\n#### Adjusted for length and model confidence\n")
            lines.append("```\n" + json.dumps(analysis["adjusted"], indent=2, default=str)[:4000] + "\n```")
        gates = summary[analysis["label"]]["blocking_gate_counts"]
        lines.append("\n#### Why proteins were not hits (closest pocket's first failing gate)\n")
        lines.append("```\n" + json.dumps(gates, indent=2) + "\n```")
        lines.append("\n#### Known inositol phosphate binders\n")
        known = analysis["known"]
        lines.append(known.to_markdown(index=False, floatfmt=".3f") if not known.empty else "none")
        lines.append("\n#### Candidates (best 40)\n")
        candidates = analysis["candidates"]
        lines.append(candidates.head(40).to_markdown(index=False, floatfmt=".3f") if not candidates.empty else "none")
        top = analysis["top"]
        cols = [col for col in ("organism_key", "uniprot_id", "gene", "best_confident_score",
                                "max_confident_hull_depth", "is_hit", "blocking_gate", "length",
                                "annotated_ip_binder") if col in top.columns]
        lines.append("\n#### Top 10 per proteome by best confident pocket\n")
        lines.append(top.groupby("organism_key").head(10)[cols].to_markdown(index=False, floatfmt=".3f"))
        if not analysis["enrichment"].empty:
            lines.append("\n#### Keyword enrichment\n")
            lines.append(analysis["enrichment"].to_markdown(index=False, floatfmt=".3g"))
        lines.append("")
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
    p.add_argument(
        "--keep-screen-scores",
        action="store_true",
        help="use the composite score computed during the screen instead of rescoring",
    )
    p.set_defaults(func=cmd_aggregate)
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
