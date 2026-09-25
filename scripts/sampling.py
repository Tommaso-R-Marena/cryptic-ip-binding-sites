#!/usr/bin/env python3
"""Study G, docking (docs/SAMPLING_PLAN.md): the same sites at a larger search budget.

    python scripts/sampling.py dock --arm e128 --census census.csv --structures-dir s \
        --ccd-dir ccd --shard-index 0 --shard-count 30 --out-dir arms

Receptor, ligand, box, starting pose, pose list and E_el are exactly study F's
(``scripts/rerank.py``); only Vina's exhaustiveness changes. Arm ``e128`` runs
seeds 1-3 at exhaustiveness 128 on every copy; arm ``e512`` runs seed 1 at 512 on
the pre-registered stratified subset. The exhaustiveness-32 arm is not re-docked:
the report reads it from study F's pose lists.
"""

from __future__ import annotations

import argparse
import json
import logging
import multiprocessing as mp
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

LOGGER = logging.getLogger("sampling")

N_POSES = 40
ENERGY_RANGE = 10.0

#: arm -> (exhaustiveness, seeds, subset only)
ARMS: Dict[str, tuple] = {"e128": (128, (1, 2, 3), False), "e512": (512, (1,), True)}
SUBSET_SIZE = 80


def subset_keys(frame: pd.DataFrame, size: int = SUBSET_SIZE) -> List[str]:
    """The plan's E512 subset: within each burial class, every k-th copy by ``copy_key``.

    Deterministic, uses no docking outcome, and keeps each class's share of the whole.
    """
    total = len(frame)
    out: List[str] = []
    for cls, group in sorted(frame.groupby("burial_class"), key=lambda kv: str(kv[0])):
        keys = sorted(group["copy_key"].astype(str))
        want = max(1, round(size * len(keys) / total)) if total else 0
        step = max(1, len(keys) / want) if want else 1
        out += [keys[min(len(keys) - 1, int(i * step))] for i in range(min(want, len(keys)))]
    return sorted(dict.fromkeys(out))


def dock_copy(ctx, exhaustiveness: int, seeds: Sequence[int]) -> List[Dict[str, object]]:
    """One copy at one exhaustiveness: per seed the pose list with Vina score, RMSD and E_el.

    This mirrors ``scripts/rerank.py``'s routine deliberately rather than importing it:
    study F is finished and recorded, and editing its script would re-run its workflow.
    ``tests/test_sampling.py`` pins the two against each other on the same settings.
    """
    from cryptic_ip.docking import engine
    from cryptic_ip.docking.ligand import start_pose, to_pdbqt
    from cryptic_ip.rescoring.electrostatics import ReceptorField, read_pdbqt, read_pdbqt_models

    receptor, _ = ctx.receptor()
    field = ReceptorField(read_pdbqt(Path(receptor.pdbqt).read_text()))
    size = [ctx.side] * 3
    runs: List[Dict[str, object]] = []
    for seed in seeds:
        pose = start_pose(ctx.ligands["primary"], ctx.centre, seed, crystal=ctx.crystal)
        result = engine.dock(receptor, to_pdbqt(pose.mol), ctx.centre, size, seed=seed,
                             exhaustiveness=exhaustiveness, n_poses=N_POSES, energy_range=ENERGY_RANGE)
        table = engine.pose_table(result, ctx.crystal, site_centroid=ctx.centre)
        models = read_pdbqt_models(result.poses_pdbqt)
        if len(models) != len(table.scores):
            raise RuntimeError(f"{len(models)} PDBQT models but {len(table.scores)} scored poses")
        runs.append({"seed": seed, "vina": table.scores, "rmsd": table.rmsd,
                     "eel": [field.energy(m) for m in models], "seconds": result.seconds})
    return runs


def cmd_dock_one(args: argparse.Namespace) -> int:
    """Child process: one copy, one arm."""
    from redocking import Context, _json_default

    exhaustiveness, seeds, _ = ARMS[args.arm]
    row = json.loads(args.row_json.read_text())
    record: Dict[str, object] = {"copy_key": row["copy_key"], "pdb_id": row["pdb_id"], "arm": args.arm,
                                 "exhaustiveness": exhaustiveness}
    start = time.time()
    try:
        ctx = Context(row, args.structures_dir, args.ccd_dir, args.work_dir)
        record["complete"] = bool(ctx.complete)
        record["runs"] = dock_copy(ctx, exhaustiveness, seeds)
    except Exception as exc:  # noqa: BLE001 - recorded as the copy's failure reason
        record["error"] = f"{type(exc).__name__}: {exc}"[:500]
        record["traceback"] = traceback.format_exc()[-1500:]
    record["seconds"] = time.time() - start
    args.out_json.write_text(json.dumps(record, default=_json_default))
    return 0


def _child(argv: List[str]) -> None:
    main(argv)


def cmd_dock(args: argparse.Namespace) -> int:
    from redocking import _isnan, _json_default, shard_entries

    census = pd.read_csv(args.census, dtype={"icode": str}, keep_default_na=False, na_values=[""])
    census["icode"] = census["icode"].fillna("")
    flags = census[["selected", "primary_set"]].astype(str).apply(lambda c: c.str.lower() == "true")
    selected = census[flags.all(axis=1)].copy()
    if ARMS[args.arm][2]:
        keys = set(subset_keys(selected))
        selected = selected[selected["copy_key"].astype(str).isin(keys)]
        LOGGER.info("%s: the pre-registered subset holds %d copies", args.arm, len(selected))
    mine = shard_entries(selected["pdb_id"].unique().tolist(), args.shard_index, args.shard_count)
    rows = selected[selected["pdb_id"].isin(mine)].sort_values(["pdb_id", "chain", "resseq"])
    if args.limit:
        rows = rows.head(args.limit)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    deadline = time.time() + args.budget_seconds
    done = not_reached = 0
    with open(args.out_dir / f"{args.arm}_{args.shard_index}.jsonl", "w") as fh:
        for _, row in rows.iterrows():
            data = {k: (None if _isnan(v) else v) for k, v in row.to_dict().items()}
            if time.time() > deadline:
                fh.write(json.dumps({"copy_key": data["copy_key"], "pdb_id": data["pdb_id"], "arm": args.arm,
                                     "error": "not reached: shard time budget"}) + "\n")
                not_reached += 1
                continue
            key = f"{str(data['copy_key']).replace(':', '_')}_{args.arm}"
            row_json, out_json = args.work_dir / f"{key}.row.json", args.work_dir / f"{key}.out.json"
            row_json.write_text(json.dumps(data, default=_json_default))
            child = ["dock-one", "--arm", args.arm, "--row-json", str(row_json), "--out-json", str(out_json),
                     "--structures-dir", str(args.structures_dir), "--ccd-dir", str(args.ccd_dir),
                     "--work-dir", str(args.work_dir)]
            proc = mp.get_context("fork").Process(target=_child, args=(child,))
            proc.start()
            limit = min(args.copy_timeout, max(60.0, deadline + args.grace_seconds - time.time()))
            proc.join(limit)
            if proc.is_alive():
                proc.terminate()
                proc.join(10)
                record = {"copy_key": data["copy_key"], "pdb_id": data["pdb_id"], "arm": args.arm,
                          "error": f"timeout after {limit:.0f} s"}
            elif out_json.exists():
                record = json.loads(out_json.read_text())
            else:
                record = {"copy_key": data["copy_key"], "pdb_id": data["pdb_id"], "arm": args.arm,
                          "error": f"child exited with code {proc.exitcode} and no record"}
            fh.write(json.dumps(record, default=_json_default) + "\n")
            fh.flush()
            done += 1
            LOGGER.info("%s %s: %s (%.0f s)", args.arm, data["copy_key"], record.get("error", "ok"),
                        record.get("seconds", 0.0))
    print(json.dumps({"arm": args.arm, "shard": args.shard_index, "copies": int(len(rows)), "run": done,
                      "not_reached": not_reached}))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("dock", "dock-one"):
        d = sub.add_parser(name)
        d.add_argument("--arm", choices=sorted(ARMS), required=True)
        d.add_argument("--structures-dir", type=Path, required=True)
        d.add_argument("--ccd-dir", type=Path, required=True)
        d.add_argument("--work-dir", type=Path, default=Path("/tmp/sampling_work"))
        if name == "dock":
            d.add_argument("--census", type=Path, required=True)
            d.add_argument("--shard-index", type=int, default=0)
            d.add_argument("--shard-count", type=int, default=1)
            d.add_argument("--out-dir", type=Path, required=True)
            d.add_argument("--budget-seconds", type=float, default=5 * 3600 + 20 * 60 - 1800)
            d.add_argument("--grace-seconds", type=float, default=1800)
            d.add_argument("--copy-timeout", type=float, default=5400)
            d.add_argument("--limit", type=int, default=0)
        else:
            d.add_argument("--row-json", type=Path, required=True)
            d.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args(argv)
    return cmd_dock(args) if args.command == "dock" else cmd_dock_one(args)


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
