#!/usr/bin/env python3
"""Entry table for the transfer dataset of phosphate-dense ligands (docs/TRANSFER_PLAN.md).

Resolves the ligand class by its formula rule, finds X-ray entries at the
plan's resolution that hold a class ligand, removes every entry that holds any
inositol phosphate (so the transfer set and the inositol phosphate benchmark
share no entry), draws the plan's seeded sample, and writes the entry table
with metadata, the class component list (for extraction's ``--comp-ids``) and
a manifest recording every decision.

    python scripts/build_polyanion_dataset.py --entry-csv data/transfer/entries.csv \\
        --comp-ids-out data/transfer/comp_ids.txt --manifest data/transfer/manifest.json
"""

from __future__ import annotations

import argparse
import json
import logging
import random
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.database.ip_ligands import comp_ids as ip_comp_ids  # noqa: E402
from cryptic_ip.database.ip_ligands import discover_ip_ligands  # noqa: E402
from cryptic_ip.database.polyanion_ligands import resolve_polyanion_ligands  # noqa: E402
from cryptic_ip.database.rcsb_client import RcsbClient, RcsbUnavailableError  # noqa: E402

LOGGER = logging.getLogger("build_polyanion_dataset")

# Pre-registered (docs/TRANSFER_PLAN.md).
SAMPLE_ENTRIES = 1200
SAMPLE_SEED = 20260925
MAX_RESOLUTION = 2.5
METHODS = ("X-RAY DIFFRACTION",)


def sample(entries: Sequence[str], n: int, seed: int) -> List[str]:
    """A seeded sample that depends only on the set of entries, not their order."""
    pool = sorted(set(entries))
    return sorted(random.Random(seed).sample(pool, min(n, len(pool))))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--entry-csv", type=Path, required=True)
    parser.add_argument("--comp-ids-out", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, default=Path("data/transfer/api_cache"))
    parser.add_argument("--sample-entries", type=int, default=SAMPLE_ENTRIES)
    parser.add_argument("--sample-seed", type=int, default=SAMPLE_SEED)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")

    started = datetime.now(timezone.utc)
    client = RcsbClient(cache_dir=args.cache_dir)
    try:
        ligands, provenance = resolve_polyanion_ligands(client)
        if not ligands:
            raise RcsbUnavailableError("no class ligand resolved; is the chemical component API reachable?")
        class_ids = [lig.comp_id for lig in ligands]
        found = client.search_entries_with_components(
            class_ids, experimental_methods=list(METHODS), max_resolution=MAX_RESOLUTION)
        ip_ligands, _ = discover_ip_ligands(client, include_unphosphorylated=True)
        with_ip = set(client.search_entries_with_components(list(ip_comp_ids(ip_ligands))))
        eligible = sorted(set(found) - with_ip)
        chosen = sample(eligible, args.sample_entries, args.sample_seed)
        LOGGER.info("%d entries hold a class ligand; %d also hold an inositol phosphate and are removed; "
                    "%d sampled", len(set(found)), len(set(found) & with_ip), len(chosen))

        # Reuse the benchmark builder's metadata writer, so the two tables share one schema.
        from scripts.build_ip_validation_dataset import _write_metadata_only

        writer_args = argparse.Namespace(
            entry_csv=args.entry_csv, manifest=args.manifest,
            experimental_methods=list(METHODS), max_resolution=MAX_RESOLUTION,
        )
        manifest = _write_metadata_only(writer_args, client, ligands, provenance, chosen, [], started)
    except RcsbUnavailableError as exc:
        LOGGER.error("transfer dataset build failed: %s", exc)
        return 2

    manifest["transfer"] = {
        "plan": "docs/TRANSFER_PLAN.md",
        "class_comp_ids": class_ids,
        "entries_with_class_ligand": len(set(found)),
        "removed_for_inositol_phosphate": len(set(found) & with_ip),
        "eligible": len(eligible),
        "sampled": len(chosen),
        "sample_seed": args.sample_seed,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    args.comp_ids_out.parent.mkdir(parents=True, exist_ok=True)
    args.comp_ids_out.write_text(" ".join(class_ids) + "\n")
    print(json.dumps(manifest["transfer"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
