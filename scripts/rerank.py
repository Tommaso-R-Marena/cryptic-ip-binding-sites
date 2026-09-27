#!/usr/bin/env python3
"""Study F, docking (docs/RERANK_PLAN.md): Vina pose lists with an electrostatic energy per pose.

    python scripts/rerank.py dock --census census.csv --structures-dir s --ccd-dir ccd \
        --shard-index 0 --shard-count 30 --out-dir arms

Receptor, ligand, box and starting pose come from scripts/redocking.py's
``Context`` (docs/REDOCKING_PLAN.md). Per copy and seed (1-3): Vina at
exhaustiveness 32, at most 40 poses within 10 kcal/mol; per pose the Vina score,
the symmetric RMSD to the crystal copy and E_el (cryptic_ip.rescoring.electrostatics).
The minimised crystal pose is scored the same way. Every copy runs in a child
process with a time limit; a shard stops starting copies when its budget is spent.
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

LOGGER = logging.getLogger("rerank")
N_POSES = 40
ENERGY_RANGE = 10.0
EXHAUSTIVENESS = [32]


def dock_copy(ctx) -> List[Dict[str, object]]:
    from cryptic_ip.docking import engine
    from cryptic_ip.docking.ligand import pose_on_crystal, start_pose, to_pdbqt
    from cryptic_ip.rescoring.electrostatics import ReceptorField, read_pdbqt, read_pdbqt_models

    receptor, _ = ctx.receptor()
    field = ReceptorField(read_pdbqt(Path(receptor.pdbqt).read_text()))
    size = [ctx.side] * 3
    runs: List[Dict[str, object]] = []
    for seed in engine.SEEDS:
        pose = start_pose(ctx.ligands["primary"], ctx.centre, seed, crystal=ctx.crystal)
        ligand_text = to_pdbqt(pose.mol)
        result = engine.dock(receptor, ligand_text, ctx.centre, size, seed=seed, exhaustiveness=EXHAUSTIVENESS[0],
                             n_poses=N_POSES, energy_range=ENERGY_RANGE)
        table = engine.pose_table(result, ctx.crystal, site_centroid=ctx.centre)
        models = read_pdbqt_models(result.poses_pdbqt)
        if len(models) != len(table.scores):
            raise RuntimeError(f"{len(models)} PDBQT models but {len(table.scores)} scored poses")
        runs.append({"seed": seed, "vina": table.scores, "rmsd": table.rmsd,
                     "eel": [field.energy(m) for m in models], "ligand_net_charge": read_pdbqt(ligand_text).net,
                     "seconds": result.seconds})
    if ctx.complete:
        placed = pose_on_crystal(ctx.ligands["primary"], ctx.crystal)
        _, after, text = engine.score_and_minimise(receptor, to_pdbqt(placed), ctx.centre, size)
        runs.append({"seed": 0, "crystal_minimised_vina": after,
                     "crystal_minimised_eel": field.energy(read_pdbqt_models(text)[0])})
    return runs


def cmd_dock_one(args: argparse.Namespace) -> int:
    from redocking import Context, _json_default

    row = json.loads(args.row_json.read_text())
    EXHAUSTIVENESS[0] = args.exhaustiveness
    record: Dict[str, object] = {"copy_key": row["copy_key"], "pdb_id": row["pdb_id"]}
    start = time.time()
    try:
        ctx = Context(row, args.structures_dir, args.ccd_dir, args.work_dir)
        record["complete"] = bool(ctx.complete)
        record["runs"] = dock_copy(ctx)
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
    truthy = census[["selected", "primary_set"]].astype(str).apply(lambda c: c.str.lower() == "true")
    selected = census[truthy.all(axis=1)]
    mine = shard_entries(selected["pdb_id"].unique().tolist(), args.shard_index, args.shard_count)
    rows = selected[selected["pdb_id"].isin(mine)].sort_values(["pdb_id", "chain", "resseq"])
    if args.limit:
        rows = rows.head(args.limit)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    deadline = time.time() + args.budget_seconds
    done = not_reached = 0
    with open(args.out_dir / f"rerank_{args.shard_index}.jsonl", "w") as fh:
        for _, row in rows.iterrows():
            data = {k: (None if _isnan(v) else v) for k, v in row.to_dict().items()}
            if time.time() > deadline:
                fh.write(json.dumps({"copy_key": data["copy_key"], "pdb_id": data["pdb_id"],
                                     "error": "not reached: shard time budget"}) + "\n")
                not_reached += 1
                continue
            key = str(data["copy_key"]).replace(":", "_")
            row_json, out_json = args.work_dir / f"{key}.row.json", args.work_dir / f"{key}.out.json"
            row_json.write_text(json.dumps(data, default=_json_default))
            child_args = ["dock-one", "--row-json", str(row_json), "--out-json", str(out_json),
                          "--structures-dir", str(args.structures_dir), "--ccd-dir", str(args.ccd_dir),
                          "--work-dir", str(args.work_dir), "--exhaustiveness", str(args.exhaustiveness)]
            proc = mp.get_context("fork").Process(target=_child, args=(child_args,))
            proc.start()
            limit = min(args.copy_timeout, max(60.0, deadline + args.grace_seconds - time.time()))
            proc.join(limit)
            if proc.is_alive():
                proc.terminate()
                proc.join(10)
                record = {"copy_key": data["copy_key"], "pdb_id": data["pdb_id"],
                          "error": f"timeout after {limit:.0f} s"}
            elif out_json.exists():
                record = json.loads(out_json.read_text())
            else:
                record = {"copy_key": data["copy_key"], "pdb_id": data["pdb_id"],
                          "error": f"child exited with code {proc.exitcode} and no record"}
            fh.write(json.dumps(record, default=_json_default) + "\n")
            fh.flush()
            done += 1
            LOGGER.info("rerank %s: %s (%.0f s)", data["copy_key"], record.get("error", "ok"),
                        record.get("seconds", 0.0))
    print(json.dumps({"shard": args.shard_index, "copies": int(len(rows)), "run": done, "not_reached": not_reached}))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("dock", "dock-one"):
        d = sub.add_parser(name)
        d.add_argument("--structures-dir", type=Path, required=True)
        d.add_argument("--ccd-dir", type=Path, required=True)
        d.add_argument("--work-dir", type=Path, default=Path("/tmp/rerank_work"))
        d.add_argument("--exhaustiveness", type=int, default=32, help="the plan fixes 32; tests only")
        if name == "dock":
            d.add_argument("--census", type=Path, required=True)
            d.add_argument("--shard-index", type=int, default=0)
            d.add_argument("--shard-count", type=int, default=1)
            d.add_argument("--out-dir", type=Path, required=True)
            d.add_argument("--budget-seconds", type=float, default=5 * 3600 + 20 * 60 - 1800)
            d.add_argument("--grace-seconds", type=float, default=1200)
            d.add_argument("--copy-timeout", type=float, default=3600)
            d.add_argument("--limit", type=int, default=0)
        else:
            d.add_argument("--row-json", type=Path, required=True)
            d.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args(argv)
    return cmd_dock(args) if args.command == "dock" else cmd_dock_one(args)


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
