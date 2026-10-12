#!/usr/bin/env python3
"""Run the detector benchmark and write its report.

    python scripts/detect_benchmark.py run
    python scripts/detect_benchmark.py run --fast --only B1 B2
    python scripts/detect_benchmark.py run --out results/detector/benchmark.json \
                                           --markdown results/detector/BENCHMARK.md

These are engineering measurements on synthetic fixtures with valid covalent
geometry. They show that the detector is deterministic, that its stated
invariances hold, how it degrades under coordinate noise, what it costs and
how that cost scales, and that its statistical machinery recovers parameters
whose true value is known by construction.

They do not show that it finds real binding sites. That is Q1-Q5 in
docs/DETECTOR_PLAN.md, asked of real structures, and none of those has run.
A benchmark that blurred the two would be the most dangerous artefact in this
repository, because it would look exactly like validation.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.detect import benchmark as BM  # noqa: E402


def cmd_run(args: argparse.Namespace) -> int:
    report = BM.run_benchmark(seed=args.seed, only=args.only, fast=args.fast)
    text = report.markdown()
    print(text)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(report.as_dict(), indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.out}")
    if args.markdown:
        Path(args.markdown).parent.mkdir(parents=True, exist_ok=True)
        Path(args.markdown).write_text(text)
        print(f"wrote {args.markdown}")
    if not report.ok:
        failures = [m.name for suite in report.suites for m in suite.failures]
        print(f"\n{len(failures)} measurement(s) failed: {', '.join(failures)}")
        return 1
    return 0


def cmd_list(args: argparse.Namespace) -> int:
    del args
    for code, function, _ in BM.SUITES:
        doc = (function.__doc__ or "").strip().splitlines()[0]
        print(f"{code}  {doc}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("run", help="run the suites")
    p.add_argument("--seed", type=int, default=20261012)
    p.add_argument("--only", nargs="+", help="suite codes, e.g. B1 B4")
    p.add_argument("--fast", action="store_true", help="coarser settings, for a smoke run")
    p.add_argument("--out", help="write the report as JSON")
    p.add_argument("--markdown", help="write the report as markdown")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("list", help="list the suites")
    p.set_defaults(func=cmd_list)

    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except BM.BenchmarkError as exc:
        print(f"error: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
