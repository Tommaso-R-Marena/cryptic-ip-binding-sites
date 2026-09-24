#!/usr/bin/env python3
"""Extract a result block from a GitHub Actions job log.

Workflows print their results between ``BEGIN_<NAME>_JSON`` / ``END_<NAME>_JSON``
(a JSON document) or ``BEGIN_<NAME>_B64`` / ``END_<NAME>_B64`` (a gzip file,
base64-encoded, wrapped over several lines). This reads a log - plain text, or
the JSON a log-fetching tool saved (``{"logs_content": "..."}``) - strips the
runner's timestamps and writes the block to a file, so recorded results are
copied from the run, never retyped.

    python scripts/extract_log_block.py read saved_log.txt REDOCKING_JSON -o results/redocking/redocking.json
    python scripts/extract_log_block.py read saved_log.txt REDOCKING_COPIES_B64 -o results/redocking/copies.csv

A workflow prints a block with the same tool, so both ends agree on the format:

    python scripts/extract_log_block.py emit REDOCKING_JSON results/redocking/redocking.json
"""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
import sys
from pathlib import Path
from typing import List, Optional, Sequence

TIMESTAMP = re.compile(r"^﻿?\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z ?")


def log_text(raw: str) -> str:
    """The log's text, whether ``raw`` is the log itself or a tool's JSON wrapper."""
    stripped = raw.lstrip("﻿").lstrip()
    if stripped.startswith("{"):
        try:
            data = json.loads(stripped)
        except json.JSONDecodeError:
            return raw
        for key in ("logs_content", "content", "log"):
            if isinstance(data.get(key), str):
                return data[key]
    return raw


def block_lines(text: str, name: str) -> List[str]:
    """Lines strictly between ``BEGIN_<name>`` and ``END_<name>``, timestamps removed (last block wins)."""
    begin, end = f"BEGIN_{name}", f"END_{name}"
    blocks: List[List[str]] = []
    current: Optional[List[str]] = None
    for line in text.splitlines():
        line = TIMESTAMP.sub("", line).rstrip("\r")
        if line.strip() == begin:
            current = []
        elif line.strip() == end and current is not None:
            blocks.append(current)
            current = None
        elif current is not None:
            current.append(line)
    if not blocks:
        raise SystemExit(f"no complete {begin} ... {end} block in the log")
    return blocks[-1]


def extract(text: str, name: str) -> bytes:
    lines = block_lines(text, name)
    if name.endswith("_B64"):
        return gzip.decompress(base64.b64decode("".join(line.strip() for line in lines)))
    payload = "\n".join(lines).strip()
    return (json.dumps(json.loads(payload), indent=2) + "\n").encode()


def emit(name: str, path: Path, width: int = 100) -> str:
    """The block for ``path``: compact JSON, or gzip + base64 wrapped at ``width``."""
    raw = path.read_bytes()
    if name.endswith("_B64"):
        encoded = base64.b64encode(gzip.compress(raw, mtime=0)).decode()
        body = "\n".join(encoded[i:i + width] for i in range(0, len(encoded), width))
    else:
        body = json.dumps(json.loads(raw), separators=(",", ":"), allow_nan=True)
    return f"BEGIN_{name}\n{body}\nEND_{name}"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    r = sub.add_parser("read")
    r.add_argument("log", type=Path)
    r.add_argument("name", help="block name without BEGIN_/END_, e.g. REDOCKING_JSON")
    r.add_argument("-o", "--output", type=Path, required=True)
    e = sub.add_parser("emit")
    e.add_argument("name")
    e.add_argument("path", type=Path)
    args = parser.parse_args(argv)
    if args.command == "emit":
        print(emit(args.name, args.path))
        return 0
    data = extract(log_text(args.log.read_text(encoding="utf-8", errors="replace")), args.name)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(data)
    print(f"{args.name}: {len(data)} bytes -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
