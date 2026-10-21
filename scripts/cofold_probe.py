#!/usr/bin/env python3
"""Feasibility probe for Boltz-2 co-folding on CPU (no study, no decisions).

Times one IP6 co-fold at several protein lengths so study M can be scoped to what
the available compute actually supports. Measuring runtime is not peeking at a
result: nothing here scores a candidate or compares an arm.

    python scripts/cofold_probe.py --out probe.json [--lengths 100 300 600]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence

#: A real, well-behaved fold (E. coli homoserine dehydrogenase N-terminal region)
#: truncated to the requested length. Identity is irrelevant: this times the model,
#: it does not ask a biological question.
BASE = ("MKVLAAGIVGLNLGGSLAKELVKRGHEVTVYDVNQEAVDHLVAQGATAVASPAEAAKDADLVILAVPAEAVEAVL"
        "FGENGLLEGLRPGSLLIDMSTIAPLASREISQALAEKGIHMLDAPVSGGVGGAEAGTLTFMVGGDAAVFERVKPL"
        "FEALGKNITLVGGNGDGQTAKVANQIIVALNIAAVSEALTLATKAGVDPARVREALMGGFASSKILEVHGERMIK"
        "RTFNPGFRIDLHIKDLANALDTARGVGAQLPITAAVMEMMQAAHADGLGDQDHSAVACVYEKLAGVQVKRNGDQR")


def sequence_of(length: int) -> str:
    """A sequence of exactly ``length`` residues, tiling BASE if needed."""
    return (BASE * (length // len(BASE) + 1))[:length]


def write_input(path: Path, sequence: str, ccd: str, affinity: bool) -> Path:
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f"      sequence: {sequence}",
             "      msa: empty",
             "  - ligand:", "      id: B", f"      ccd: {ccd}"]
    if affinity:
        lines += ["properties:", "  - affinity:", "      binder: B"]
    path.write_text("\n".join(lines) + "\n")
    return path


def run_one(length: int, work: Path, ccd: str, timeout: int, affinity: bool) -> Dict[str, object]:
    work.mkdir(parents=True, exist_ok=True)
    yaml_path = write_input(work / f"probe_{length}.yaml", sequence_of(length), ccd, affinity)
    cmd = ["boltz", "predict", str(yaml_path), "--accelerator", "cpu", "--out_dir", str(work / "out"),
           "--diffusion_samples", "1", "--recycling_steps", "0", "--output_format", "pdb"]
    record: Dict[str, object] = {"length": length, "ccd": ccd, "affinity_requested": affinity,
                                 "command": " ".join(cmd)}
    start = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        record["returncode"] = proc.returncode
        record["ok"] = proc.returncode == 0
        if proc.returncode != 0:
            record["stderr_tail"] = (proc.stderr or proc.stdout)[-1500:]
    except subprocess.TimeoutExpired:
        record["ok"] = False
        record["timed_out_after_s"] = timeout
    record["seconds"] = round(time.time() - start, 1)
    hits = sorted((work / "out").rglob("*.pdb")) + sorted((work / "out").rglob("affinity*.json"))
    record["outputs"] = [h.name for h in hits][:8]
    for affinity_json in sorted((work / "out").rglob("affinity*.json")):
        try:
            record["affinity"] = json.loads(affinity_json.read_text())
        except Exception as exc:  # noqa: BLE001 - probe only
            record["affinity_error"] = str(exc)[:200]
        break
    return record


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--lengths", type=int, nargs="+", default=[100, 300, 600])
    parser.add_argument("--ccd", default="IHP", help="IP6 in the PDB chemical component dictionary")
    parser.add_argument("--timeout", type=int, default=5400)
    parser.add_argument("--work", type=Path, default=Path("probe_work"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--no-affinity", action="store_true")
    args = parser.parse_args(argv)

    results: List[Dict[str, object]] = []
    for length in args.lengths:
        record = run_one(length, args.work / str(length), args.ccd, args.timeout, not args.no_affinity)
        results.append(record)
        print(json.dumps(record, default=str)[:900], flush=True)
        if not record.get("ok"):
            print(f"stopping: length {length} did not complete", flush=True)
            break
    payload = {"probe": "boltz-2 cpu feasibility", "results": results,
               "note": "runtime measurement only; no candidate is scored and no arm is compared"}
    args.out.write_text(json.dumps(payload, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
