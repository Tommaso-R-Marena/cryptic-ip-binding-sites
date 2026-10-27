#!/usr/bin/env python3
"""Fetch a structure set across several attempts, waiting between them.

The strict yeast re-screen needs all 499 of its AlphaFold models before it
means anything, and the endpoint rate-limits by address. One pass gets a
fraction of them, so the set is assembled over several attempts with the
partial download kept between each - which the workflow's cache already does.
What was missing was anything to make the next attempt happen: a human
re-ran the workflow, and the set sat incomplete for days between nudges.

This supplies the retry policy, so one job assembles what it can:

    fetch_until_complete.py alphafold --out DIR --ids-json F --json-key processed \
        --attempts 5 --sleep 1800 --progress progress.json

Three rules decide when to stop, and the middle one is the point:

* **complete** - nothing is missing, so stop and succeed;
* **stalled** - two attempts in a row added nothing, so the endpoint's window
  is longer than this job can usefully wait. Stop and fail, rather than
  burning the remaining budget on sleeps that cannot help. The partial set is
  kept, so a later run resumes from it;
* **exhausted** - the attempt budget ran out with progress still being made,
  which is a reason to run again rather than to wait longer here.

A deficit of `not_found` is never retried. That is the endpoint saying an
accession has no model, which no amount of waiting changes, so it is reported
as a separate and permanent shortfall.

The progress file is rewritten after EVERY attempt, so a job killed at its
budget still records how much arrived and when - which is the only measurement
anyone has of how fast the quota recovers.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

#: Stop reasons, in the order the loop can reach them.
COMPLETE = "complete"
NOT_FOUND_ONLY = "only absent models remain"
STALLED = "stalled: two attempts added nothing"
EXHAUSTED = "attempt budget exhausted"


def _fetcher() -> Callable[[List[str]], int]:
    """`fetch_structures.main`, loaded by path like the other sibling scripts."""
    spec = importlib.util.spec_from_file_location(
        "fetch_structures", ROOT / "scripts" / "fetch_structures.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules["fetch_structures"] = module
    spec.loader.exec_module(module)
    return module.main


def read_manifest(path: Path) -> Dict[str, int]:
    """How many models are present, absent upstream, and still owed.

    `ok` counts a model on disk whether this attempt fetched it or a previous
    one did, because that is what the screen needs. `failed` is what another
    attempt could still get; `not_found` is what no attempt can.
    """
    if not path.exists() or path.stat().st_size == 0:
        return {"total": 0, "ok": 0, "not_found": 0, "failed": 0}
    try:
        records = json.loads(path.read_text())
    except json.JSONDecodeError:
        return {"total": 0, "ok": 0, "not_found": 0, "failed": 0}
    rows = records if isinstance(records, list) else records.get("results", [])
    ok = sum(1 for r in rows if r.get("ok"))
    not_found = sum(1 for r in rows if r.get("not_found"))
    return {"total": len(rows), "ok": ok, "not_found": not_found,
            "failed": len(rows) - ok - not_found}


def assemble(argv: Sequence[str], *, manifest: Path, attempts: int, sleep_s: float,
             progress: Optional[Path] = None,
             fetch: Optional[Callable[[List[str]], int]] = None,
             sleeper: Callable[[float], None] = time.sleep) -> Dict[str, object]:
    """Call the fetcher until the set is complete, stalled, or out of attempts."""
    call = fetch if fetch is not None else _fetcher()
    rows: List[Dict[str, object]] = []
    before = read_manifest(manifest)["ok"]
    reason, barren = EXHAUSTED, 0

    for attempt in range(1, int(attempts) + 1):
        code = call(list(argv))
        counts = read_manifest(manifest)
        gained = counts["ok"] - before
        before = counts["ok"]
        rows.append({"attempt": attempt, "exit_code": int(code), "gained": gained,
                     **counts})
        if progress:
            write_progress(progress, rows, reason="in progress")
        print(f"attempt {attempt}/{attempts}: +{gained} -> {counts['ok']} of "
              f"{counts['total']} present, {counts['failed']} owed, "
              f"{counts['not_found']} absent upstream", flush=True)

        if counts["total"] and counts["failed"] == 0 and counts["not_found"] == 0:
            reason = COMPLETE
            break
        if counts["total"] and counts["failed"] == 0:
            # Waiting cannot conjure a model the endpoint does not have.
            reason = NOT_FOUND_ONLY
            break
        barren = barren + 1 if gained <= 0 else 0
        if barren >= 2:
            reason = STALLED
            break
        if attempt < int(attempts):
            print(f"waiting {sleep_s:.0f}s before the next attempt", flush=True)
            sleeper(float(sleep_s))

    final = read_manifest(manifest)
    report = {"reason": reason, "complete": reason == COMPLETE,
              "attempts_used": len(rows), "attempts_allowed": int(attempts),
              "sleep_seconds": float(sleep_s), "attempts": rows, **final}
    if progress:
        write_progress(progress, rows, reason=reason, final=final)
    return report


def write_progress(path: Path, rows: List[Dict[str, object]], *, reason: str,
                   final: Optional[Dict[str, int]] = None) -> None:
    """The attempt-by-attempt record, rewritten after every attempt.

    Written as the loop goes rather than at the end, so a job killed at its
    budget still leaves the only measurement anyone has of the endpoint's
    recovery rate.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: Dict[str, object] = {"reason": reason, "attempts": rows}
    if final:
        payload.update(final)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", choices=["rcsb", "alphafold"])
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--ids-json")
    parser.add_argument("--json-key")
    parser.add_argument("--ids-csv")
    parser.add_argument("--column")
    parser.add_argument("--format")
    parser.add_argument("--concurrency", type=int, default=0)
    parser.add_argument("--attempts", type=int, default=5)
    parser.add_argument("--sleep", type=float, default=1800.0)
    parser.add_argument("--progress", type=Path, default=None)
    args = parser.parse_args(list(argv) if argv is not None else None)

    passed: List[str] = [args.source, "--out", str(args.out)]
    for flag, value in (("--ids-json", args.ids_json), ("--json-key", args.json_key),
                        ("--ids-csv", args.ids_csv), ("--column", args.column),
                        ("--format", args.format)):
        if value:
            passed += [flag, str(value)]
    if args.concurrency:
        passed += ["--concurrency", str(args.concurrency)]

    report = assemble(passed, manifest=args.out / f"{args.source}_manifest.json",
                      attempts=args.attempts, sleep_s=args.sleep, progress=args.progress)
    print(json.dumps({k: v for k, v in report.items() if k != "attempts"}, indent=2))
    return 0 if report["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
