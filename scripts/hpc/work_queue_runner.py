"""Simple filesystem-backed work queue runner for cluster nodes."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", required=True, help="JSONL queue file")
    parser.add_argument("--results", required=True, help="JSONL output file")
    args = parser.parse_args()

    queue_path = Path(args.queue)
    result_path = Path(args.results)
    result_path.parent.mkdir(parents=True, exist_ok=True)

    with queue_path.open() as qh, result_path.open("a") as rh:
        for line in qh:
            if not line.strip():
                continue
            item = json.loads(line)
            started = time.time()
            rh.write(json.dumps({"uniprot_id": item["uniprot_id"], "status": "queued", "started": started}) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
