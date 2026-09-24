#!/usr/bin/env python3
"""Fetch structures concurrently, verified, with a manifest.

One command for every workflow that downloads structures, replacing
sequential ``wget``/``urllib`` loops that neither retried, nor checked what
they received, nor wrote atomically - so a cut-off download could be kept and
reused as a cached structure.

Examples::

    # RCSB entries, legacy PDB preferred, mmCIF when an entry has no PDB file
    python scripts/fetch_structures.py rcsb --out data/validation 1ZY7 1MAI
    python scripts/fetch_structures.py rcsb --out data/validation/raw/structures \\
        --ids-csv data/validation/ip_binding_validation_dataset.csv --column pdb_id

    # current AlphaFold models by UniProt accession
    python scripts/fetch_structures.py alphafold --out models --ids-file ids.txt

Exits non-zero when any download *failed* (retries exhausted); identifiers the
source does not have are reported but are not failures, unless
``--require-all`` is given.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.database.async_fetch import (  # noqa: E402
    fetch_alphafold_models,
    fetch_rcsb_structures,
    write_manifest,
)


def _identifiers(args: argparse.Namespace) -> List[str]:
    ids = list(args.ids or [])
    if args.ids_file:
        ids += [line.strip() for line in Path(args.ids_file).read_text().splitlines() if line.strip()]
    if args.ids_csv:
        import pandas as pd

        frame = pd.read_csv(args.ids_csv, dtype=str)
        ids += frame[args.column].dropna().astype(str).tolist()
    if args.ids_json:
        data = json.loads(Path(args.ids_json).read_text())
        for key in (args.json_key or "").split("."):
            if key:
                data = data[key]
        ids += [str(x) for x in data]
    return list(dict.fromkeys(i.strip().upper() for i in ids if i.strip()))


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", choices=["rcsb", "alphafold"])
    parser.add_argument("ids", nargs="*")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--ids-file")
    parser.add_argument("--ids-csv")
    parser.add_argument("--column", default="pdb_id")
    parser.add_argument("--ids-json", help="JSON file holding a list of identifiers")
    parser.add_argument("--json-key", help="Dotted key of the list inside --ids-json")
    parser.add_argument(
        "--prefer", default="pdb,cif",
        help="RCSB formats in order of preference (default pdb,cif: legacy PDB, "
        "falling back to mmCIF for entries too large to have one)",
    )
    parser.add_argument("--format", default="pdb", choices=["pdb", "cif"], help="AlphaFold model format")
    parser.add_argument("--concurrency", type=int, default=0, help="0 = source default")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--require-all", action="store_true")
    # Intermixed: identifiers may follow the options ("rcsb --out d 1ABC").
    args = parser.parse_intermixed_args(argv)

    ids = _identifiers(args)
    if not ids:
        print("no identifiers given", file=sys.stderr)
        return 2
    manifest = args.manifest or args.out / f"{args.source}_manifest.json"
    if args.source == "rcsb":
        results = fetch_rcsb_structures(
            ids, args.out, prefer=[f.strip() for f in args.prefer.split(",") if f.strip()],
            concurrency=args.concurrency or 8,
        )
    else:
        results = fetch_alphafold_models(ids, args.out, fmt=args.format, concurrency=args.concurrency or 16)

    summary = write_manifest(results.values(), manifest)
    print(json.dumps({"source": args.source, **summary, "manifest": str(manifest)}))
    for result in results.values():
        if not result.ok:
            kind = "absent" if result.not_found else "FAILED"
            print(f"  {kind} {result.key}: {result.error} ({result.attempts} attempts)")
    if summary["failed"] or (args.require_all and summary["not_found"]):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
