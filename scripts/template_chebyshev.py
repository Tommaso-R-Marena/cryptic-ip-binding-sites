#!/usr/bin/env python3
"""Study O (docs/TEMPLATE2_PLAN.md): the corrected IP6 template transplant.

Study N passed both its primary tests and told us nothing, for three reasons recorded in
``results/template/PROVENANCE.md``: its shortlist threshold collapsed onto the censoring
floor, its placements were not buried, and its search was not exhaustive because distance
pruning is unsound under an RMSD acceptance criterion. This study fixes each.

* **Chebyshev acceptance.** A placement is admissible when *every* matched anchor lies
  within :data:`EPS` of its partner, and pruning runs at ``2 * EPS``. ``pruning_complete``
  in ``formal/RequestProject/Pruning.lean`` proves that discards no admissible
  correspondence, so the search is exhaustive by construction rather than by assertion.
* **Burial of the placed ligand**, measured by the same function that produced the
  census's own ``relative_sasa`` column, against a threshold read off the crystal IP6
  burial classes before any candidate was scored.
* **Controls matched on basic-residue count**, which the fit score structurally requires
  and which study N's controls were never matched on.

The geometric primitives are imported from ``scripts/template_fit.py`` rather than
reimplemented, so study N and study O share one Kabsch and one correspondence search.

Subcommands::

    templates  anchors, ligand coordinates and ligand atom identities per IHP copy
    pool       each arm's queries plus a pool of candidate controls (no structures)
    counts     basic-residue count of each pocket, from the AlphaFold models (sharded)
    select     keep the pool member whose basic count is closest to its candidate's
    fit        score one arm's queries against the template library (sharded)
    guard      O1 alone, so the workflow can gate the candidate arm on it
    report     O1, O2 and O3
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import template_fit as tf  # noqa: E402  the shared geometric primitives
from cryptic_ip.docking.stats import MIN_GROUPS, auc_estimate, paired_difference  # noqa: E402

LOGGER = logging.getLogger("template_chebyshev")

#: Fixed by docs/TEMPLATE2_PLAN.md. None may be tuned on results.
SEED = 20261007
N_BOOTSTRAP = 2000
#: Chebyshev acceptance: every matched anchor within EPS of its partner.
EPS = 2.5
#: Pruning tolerance. ``pruning_complete`` proves 2*EPS discards nothing admissible.
PRUNE_TOLERANCE = 2.0 * EPS
K_ANCHORS = 4
#: Burial of the *placed* ligand. The largest relative SASA among the 64 crystal IHP
#: copies this project calls cryptic or semi-cryptic; the stricter cryptic-only cut is
#: reported as a sensitivity and cannot change a decision.
BURIED_MAX_REL_SASA = 0.2363
CRYPTIC_MAX_REL_SASA = 0.1044
#: Controls drawn per candidate before selecting on basic-residue count.
CONTROL_POOL = 5
#: O2's no-difference margin, in angstroms.
NO_DIFFERENCE = 0.25
TEMPLATE_COMP_ID = "IHP"
ARMS = ("annotated", "annotated_controls", "candidates", "controls")
QUERY_ARMS = {"annotated": "annotated_controls", "candidates": "controls"}


# ------------------------------------------------------------------ templates
def ligand_identity(arrays, ligand_atoms: np.ndarray) -> Tuple[List[str], List[str]]:
    """Atom names and elements of a copy's heavy atoms, in the same order as its coords.

    Needed because burial is measured by writing the placed ligand back out as a PDB and
    calling the project's own SASA routine, which identifies the copy from its atoms.
    """
    heavy = ligand_atoms[arrays.elements[ligand_atoms] != "H"]
    names = [str(n).strip() for n in arrays.atom_names[heavy]]
    elements = [str(e).strip().upper() for e in arrays.elements[heavy]]
    return names, elements


def cmd_templates(args: argparse.Namespace) -> int:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.validation.burial_metrics import find_ligand_instances

    census = pd.read_csv(args.census)
    # Stated in words that match this predicate, per amendment 2: the flag is empty for
    # cryo-EM copies, which have no unit cell, so the symmetry test is inapplicable to
    # them rather than failed, and excluding them would discard exactly the large
    # assemblies where IP6 acts as a structural cofactor.
    keep = census[(census["comp_id"] == TEMPLATE_COMP_ID) & (census["status"] == "eligible")
                  & (census["complete"].astype(str).str.lower() == "true")
                  & (census["symmetry_contact"].astype(str).str.lower() != "true")]
    LOGGER.info("%d eligible %s copies in %d entries, %d strict homology groups",
                len(keep), TEMPLATE_COMP_ID, keep["pdb_id"].nunique(),
                keep["homology_group_strict"].nunique())
    out, dropped = [], []
    for pdb_id, block in keep.groupby("pdb_id"):
        path = next(iter(sorted(Path(args.structures).glob(f"{pdb_id}.*"))), None)
        if path is None:
            dropped.append({"copy_key": f"{pdb_id}:*", "reason": "structure missing"})
            continue
        arrays = load_structure_arrays(path)
        found = {(key[1], int(key[2]), key[3] or ""): atoms
                 for key, _comp, atoms in find_ligand_instances(arrays, comp_ids=(TEMPLATE_COMP_ID,))}
        for _, row in block.iterrows():
            key = (tf._text(row["chain"]), int(row["resseq"]), tf._text(row.get("icode")))
            atoms = found.get(key)
            if atoms is None:
                dropped.append({"copy_key": row["copy_key"], "reason": "copy not found in structure"})
                continue
            points, labels, ligand_xyz = tf.template_anchors(arrays, atoms)
            if len(points) < K_ANCHORS:
                dropped.append({"copy_key": row["copy_key"], "reason": f"{len(points)} anchors"})
                continue
            names, elements = ligand_identity(arrays, atoms)
            out.append({"copy_key": row["copy_key"], "pdb_id": pdb_id,
                        "homology_group_strict": row.get("homology_group_strict"),
                        "uniprot_ids": tf._text(row.get("uniprot_ids")),
                        "burial_class": row.get("burial_class"),
                        "relative_sasa": row.get("relative_sasa"),
                        "anchors": points, "anchor_labels": labels,
                        "ligand": ligand_xyz, "ligand_names": names, "ligand_elements": elements})
    LOGGER.info("%d templates kept, %d dropped", len(out), len(dropped))
    tf.write_json(Path(args.out), {"templates": out, "dropped": dropped,
                                   "n_kept": len(out), "n_dropped": len(dropped),
                                   "n_entries": int(keep["pdb_id"].nunique()),
                                   "n_strict_groups": int(keep["homology_group_strict"].nunique())})
    return 0


# ------------------------------------------------------------------ arms
def arm_rows(proteins: pd.DataFrame, candidates: Optional[pd.DataFrame], arm: str) -> pd.DataFrame:
    """One arm's queries, before control selection."""
    import triage

    if arm in ("candidates", "controls"):
        queries = candidates[["uniprot_id", "organism_key", "top_pocket_residues", "plddt_mean",
                              "hull_depth", "combined", "cluster"]].copy()
    else:
        queries = proteins[proteins["annotated"].astype(str).str.lower() == "true"].copy()
    pool = triage.matched_controls(queries, proteins, n_per_candidate=CONTROL_POOL)
    # As in studies I and N: a protein drawn as a control that is also a query counts as a
    # control only.
    clash = set(pool["uniprot_id"]) & set(queries["uniprot_id"])
    if clash:
        LOGGER.info("%d queries are also drawn as controls and count as controls only: %s",
                    len(clash), sorted(clash))
        queries = queries[~queries["uniprot_id"].isin(clash)]
    if arm in ("controls", "annotated_controls"):
        frame = pool
    else:
        frame = queries.assign(matched_to=queries["uniprot_id"], pool_rank=-1)
    keep = ["uniprot_id", "top_pocket_residues", "cluster", "matched_to", "pool_rank"]
    frame = frame[[c for c in keep if c in frame.columns]].copy()
    frame["arm"] = arm
    return frame


def cmd_pool(args: argparse.Namespace) -> int:
    proteins = pd.read_csv(args.proteins)
    candidates = pd.read_csv(args.candidates) if args.candidates else None
    frames = [arm_rows(proteins, candidates, arm) for arm in ARMS]
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(args.out, index=False)
    LOGGER.info("pool: %s", out.groupby("arm").size().to_dict())
    return 0


def basic_count(arrays, residues: Sequence[int]) -> int:
    """How many K/R/H the pocket holds - the quantity the fit score requires."""
    wanted = {(str(c), int(r)) for c in set(arrays.chain_ids.astype(str)) for r in residues}
    points, _labels = tf.residue_points(arrays, wanted)
    return int(len(points))


def cmd_counts(args: argparse.Namespace) -> int:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    rows = pd.read_csv(args.rows).drop_duplicates("uniprot_id")
    shard = tf.shard(list(rows.to_dict("records")), args.shard, args.shards)
    if args.list_only:
        # The workflow must fetch this shard's models before it can count them, and the
        # shard split has to be the same one the counting uses.
        accessions = sorted({str(r["uniprot_id"]) for r in shard})
        Path(args.out).write_text("\n".join(accessions) + "\n")
        LOGGER.info("%d accessions for shard %d/%d", len(accessions), args.shard, args.shards)
        return 0
    out = []
    for row in shard:
        acc = str(row["uniprot_id"])
        record = {"uniprot_id": acc}
        path = next(iter(sorted(Path(args.af_dir).glob(f"AF-{acc}-F1-*.pdb"))), None)
        if path is None:
            record["error"] = "no AlphaFold model"
        else:
            try:
                arrays = load_structure_arrays(path)
                record["n_basic"] = basic_count(arrays, tf.parse_positions(row.get("top_pocket_residues")))
            except Exception as exc:  # a structure we cannot read is recorded, not fatal
                record["error"] = f"{type(exc).__name__}: {exc}"
        out.append(record)
    LOGGER.info("counted %d of %d accessions", len(out), len(rows))
    tf.write_json(Path(args.out), out)
    return 0


def select_controls(pool: pd.DataFrame, counts: Dict[str, int]) -> pd.DataFrame:
    """Keep, per query, the pool member whose basic count is closest to the query's.

    Ties break on accession, so the choice does not depend on row order.
    """
    pool = pool.copy()
    pool["n_basic"] = pool["uniprot_id"].map(counts)
    kept = [pool[pool["arm"].isin(QUERY_ARMS.keys())]]
    for query_arm, control_arm in QUERY_ARMS.items():
        queries = pool[pool["arm"] == query_arm].set_index("uniprot_id")
        controls = pool[pool["arm"] == control_arm]
        picks = []
        for target, block in controls.groupby("matched_to"):
            if target not in queries.index:
                continue
            want = queries.loc[target, "n_basic"]
            block = block.dropna(subset=["n_basic"])
            if block.empty or pd.isna(want):
                continue
            block = block.assign(gap=(block["n_basic"] - want).abs())
            picks.append(block.sort_values(["gap", "uniprot_id"]).iloc[0])
        if picks:
            kept.append(pd.DataFrame(picks))
    return pd.concat(kept, ignore_index=True)


def cmd_select(args: argparse.Namespace) -> int:
    pool = pd.read_csv(args.pool)
    counts: Dict[str, int] = {}
    for part in sorted(Path(args.counts).rglob("*.json*")):
        for record in tf.read_json(part):
            if "n_basic" in record:
                counts[str(record["uniprot_id"])] = int(record["n_basic"])
    chosen = select_controls(pool, counts)
    chosen.to_csv(args.out, index=False)
    for query_arm, control_arm in QUERY_ARMS.items():
        q = chosen[chosen["arm"] == query_arm]["n_basic"].dropna()
        c = chosen[chosen["arm"] == control_arm]["n_basic"].dropna()
        LOGGER.info("%s: n=%d mean basic %.2f | %s: n=%d mean basic %.2f",
                    query_arm, len(q), q.mean() if len(q) else float("nan"),
                    control_arm, len(c), c.mean() if len(c) else float("nan"))
    return 0


# ------------------------------------------------------------------ fitting
def placed_relative_sasa(arrays, placed: np.ndarray, names: Sequence[str],
                         elements: Sequence[str]) -> Optional[float]:
    """Relative SASA of the transplanted ligand, by the project's own burial routine.

    The receptor and the placed ligand are written out as one PDB and measured with
    :func:`compute_ligand_burial`, the same function that produced the census's
    ``relative_sasa`` column - so the threshold in the plan and the measurement here are
    on the same footing.
    """
    from cryptic_ip.validation.burial_metrics import compute_ligand_burial

    lines, serial = [], 1
    polymer = arrays.is_polymer & (arrays.elements != "H")
    for idx in np.flatnonzero(polymer):
        name = str(arrays.atom_names[idx]).strip()
        x, y, z = arrays.coords[idx]
        lines.append(f"ATOM  {serial:5d} {name:<4s} {str(arrays.resnames[idx]).strip():>3s} "
                     f"A{int(arrays.resseqs[idx]):4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00"
                     f"          {str(arrays.elements[idx]).strip():>2s}")
        serial += 1
    for name, element, (x, y, z) in zip(names, elements, placed):
        lines.append(f"HETATM{serial:5d} {name:<4s} {TEMPLATE_COMP_ID} X9999    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}")
        serial += 1
    with tempfile.TemporaryDirectory() as work:
        path = Path(work) / "placed.pdb"
        path.write_text("\n".join(lines) + "\nEND\n")
        try:
            records = compute_ligand_burial(path, comp_ids=(TEMPLATE_COMP_ID,))
        except Exception as exc:  # pragma: no cover - depends on the SASA backend
            LOGGER.debug("burial failed: %s", exc)
            return None
    return float(records[0].relative_sasa) if records else None


def best_chebyshev(template_xyz: np.ndarray, ligand_xyz: np.ndarray, query_xyz: np.ndarray,
                   tree, *, k: int = K_ANCHORS, eps: float = EPS) -> Optional[Dict[str, object]]:
    """Lowest achievable worst-anchor deviation among clash-free placements.

    Pruning runs at ``2 * eps``, which ``pruning_complete`` proves loses no placement whose
    worst anchor is within ``eps``. Candidates are then ranked by that worst deviation and
    clash-tested in order, so the first admissible one is the minimum.
    """
    maps = tf.correspondences(template_xyz, query_xyz, k, tolerance=2.0 * eps)
    if not maps:
        return None
    mobile = np.stack([template_xyz[list(t)] for t, _ in maps])
    target = np.stack([query_xyz[list(q)] for _, q in maps])
    rotation, translation, _rmsd = tf.kabsch_batch(mobile, target)
    fitted = np.einsum("nki,nij->nkj", mobile, rotation) + translation[:, None, :]
    worst = np.linalg.norm(fitted - target, axis=2).max(axis=1)
    order = np.argsort(worst)
    order = order[worst[order] <= eps]
    for start in range(0, len(order), tf.CLASH_BATCH):
        batch = order[start:start + tf.CLASH_BATCH]
        placed = np.einsum("ai,nij->naj", ligand_xyz, rotation[batch]) + translation[batch][:, None, :]
        hits = tree.query_ball_point(placed.reshape(-1, 3), tf.CLASH_DISTANCE, return_length=True)
        ok = ~np.any(hits.reshape(len(batch), -1) > 0, axis=1)
        if not ok.any():
            continue
        idx = int(batch[np.flatnonzero(ok)[0]])
        return {"deviation": float(worst[idx]), "n_anchors": int(k),
                "template_anchors": [int(i) for i in maps[idx][0]],
                "query_anchors": [int(i) for i in maps[idx][1]],
                "ligand_xyz": ligand_xyz @ rotation[idx] + translation[idx]}
    return None


def score_query(templates: Sequence[dict], arrays, residues: Sequence[int],
                *, eps: float = EPS) -> Dict[str, object]:
    from scipy.spatial import cKDTree

    wanted = {(str(c), int(r)) for c in set(arrays.chain_ids.astype(str)) for r in residues}
    points, labels = tf.residue_points(arrays, wanted)
    n_basic = int(len(points))
    if n_basic < K_ANCHORS:
        return {"score": -eps, "deviation": eps, "n_basic": n_basic, "censored": True,
                "reason": "fewer basic residues in the pocket than matched anchors"}
    centre = points.mean(axis=0)
    points, labels = tf.nearest(points, labels, centre, tf.MAX_QUERY_ANCHORS)
    heavy = arrays.is_polymer & (arrays.elements != "H")
    tree = cKDTree(arrays.coords[heavy])
    best: Optional[Dict[str, object]] = None
    for t in templates:
        fit = best_chebyshev(np.asarray(t["anchors"], dtype=float),
                             np.asarray(t["ligand"], dtype=float), points, tree, eps=eps)
        if fit is None:
            continue
        if best is None or fit["deviation"] < best["deviation"]:
            placed = fit.pop("ligand_xyz")
            fit["template"] = t["copy_key"]
            fit["template_burial"] = t.get("burial_class")
            fit["query_anchor_labels"] = [labels[i] for i in fit["query_anchors"]]
            fit["relative_sasa"] = placed_relative_sasa(
                arrays, placed, t.get("ligand_names") or [], t.get("ligand_elements") or [])
            best = fit
    if best is None:
        return {"score": -eps, "deviation": eps, "n_basic": n_basic, "censored": True,
                "reason": "no clash-free placement within eps"}
    rel = best.get("relative_sasa")
    best.update({"score": -best["deviation"], "n_basic": n_basic, "censored": False,
                 "buried": bool(rel is not None and rel <= BURIED_MAX_REL_SASA),
                 "cryptic": bool(rel is not None and rel <= CRYPTIC_MAX_REL_SASA)})
    return best


def cmd_fit(args: argparse.Namespace) -> int:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    templates = tf.read_json(Path(args.templates))["templates"]
    rows = pd.read_csv(args.arms)
    rows = rows[rows["arm"] == args.arm]
    shard = tf.shard(list(rows.to_dict("records")), args.shard, args.shards)
    LOGGER.info("arm %s: %d of %d queries in shard %d/%d", args.arm, len(shard), len(rows),
                args.shard, args.shards)
    out = []
    for row in shard:
        acc = str(row["uniprot_id"])
        record = {"uniprot_id": acc, "arm": args.arm, "cluster": row.get("cluster"),
                  "matched_to": row.get("matched_to")}
        path = next(iter(sorted(Path(args.af_dir).glob(f"AF-{acc}-F1-*.pdb"))), None)
        if path is None:
            record.update({"score": None, "error": "no AlphaFold model"})
            out.append(record)
            continue
        usable, removed = tf.eligible_templates(templates, acc, row.get("cluster"))
        record.update({"templates_removed": removed, "n_templates": len(usable)})
        try:
            arrays = load_structure_arrays(path)
            record.update(score_query(usable, arrays, tf.parse_positions(row.get("top_pocket_residues"))))
        except Exception as exc:
            record.update({"score": None, "error": f"{type(exc).__name__}: {exc}"})
        out.append(record)
        LOGGER.info("%s %s score=%s buried=%s", args.arm, acc, record.get("score"), record.get("buried"))
    tf.write_json(Path(args.out), out)
    return 0


def cmd_merge(args: argparse.Namespace) -> int:
    records: List[dict] = []
    for part in sorted(Path(args.parts).rglob("*.json*")):
        payload = tf.read_json(part)
        if isinstance(payload, list):
            records.extend(payload)
    LOGGER.info("%d fit records merged", len(records))
    tf.write_json(Path(args.out), records)
    return 0


# ------------------------------------------------------------------ report
def arm_records(records: Sequence[dict], arm: str) -> pd.DataFrame:
    frame = pd.DataFrame([r for r in records if r.get("arm") == arm])
    if frame.empty or "score" not in frame.columns:
        return pd.DataFrame()
    return frame[frame["score"].notna()].copy()


def guard(records: Sequence[dict], n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    """O1: does the fit score separate annotated binders from matched controls?"""
    binders = arm_records(records, "annotated")
    controls = arm_records(records, "annotated_controls")
    if binders.empty or controls.empty:
        return {"verdict": "not evaluable", "reason": "an arm has no scored query",
                "n_binders": int(len(binders)), "n_controls": int(len(controls))}
    frame = pd.concat([binders.assign(label=1), controls.assign(label=0)], ignore_index=True)
    auc = auc_estimate(frame["label"].tolist(), frame["score"].tolist(), tf.groups_of(frame),
                       n_bootstrap=n_bootstrap, seed=SEED)
    out = {"auc": auc, "n": int(len(frame)), "n_binders": int(len(binders)),
           "n_controls": int(len(controls))}
    if not auc.get("evidence") or "roc_auc" not in auc:
        out["verdict"] = "not evaluable"
        out["reason"] = f"{auc.get('groups')} clusters, fewer than the {MIN_GROUPS} required"
        return out
    out["verdict"] = "pass" if auc["roc_auc"]["low"] > 0.5 else "fail"
    return out


def paired(records: Sequence[dict], n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    """O2: paired candidate-minus-control difference in fit score."""
    cands = arm_records(records, "candidates")
    controls = arm_records(records, "controls")
    if cands.empty or controls.empty:
        return {"decision": "not evaluable", "reason": "an arm has no scored query"}
    cands = cands.drop_duplicates("uniprot_id").set_index("uniprot_id")
    pairs = controls[controls["matched_to"].isin(cands.index)]
    if pairs.empty:
        return {"decision": "not evaluable", "reason": "no candidate kept its matched control"}
    partner = cands.loc[pairs["matched_to"]]
    est = paired_difference(partner["score"].to_numpy(dtype=float),
                            pairs["score"].to_numpy(dtype=float),
                            tf.groups_of(partner.reset_index()),
                            n_bootstrap=n_bootstrap, seed=SEED)
    out = {"estimate": est, "n_pairs": int(len(pairs))}
    if not est.get("evidence") or "per_group" not in est:
        out["decision"] = "not evaluable"
        out["reason"] = f"{est.get('groups')} clusters, fewer than the {MIN_GROUPS} required"
        return out
    low, high = est["per_group"]["low"], est["per_group"]["high"]
    if low > 0:
        out["decision"] = "better in candidates"
    elif high < 0:
        out["decision"] = "worse in candidates"
    elif low > -NO_DIFFERENCE and high < NO_DIFFERENCE:
        out["decision"] = "no difference"
    else:
        out["decision"] = "inconclusive"
    return out


def shortlist(records: Sequence[dict]) -> Dict[str, object]:
    """O3: absolute criteria - clash-free, within eps on every anchor, and buried."""
    def passing(frame: pd.DataFrame) -> pd.DataFrame:
        if frame.empty:
            return frame
        ok = (~frame["censored"].astype(bool)) & frame["buried"].astype(bool)
        return frame[ok]

    cands, controls = arm_records(records, "candidates"), arm_records(records, "controls")
    members = []
    fields = ("uniprot_id", "score", "deviation", "n_anchors", "n_basic", "template",
              "template_burial", "relative_sasa", "cryptic", "query_anchor_labels")
    for _, row in passing(cands).sort_values("deviation").iterrows():
        members.append({k: row.get(k) for k in fields})
    return {"eps": EPS, "buried_max_relative_sasa": BURIED_MAX_REL_SASA,
            "cryptic_max_relative_sasa": CRYPTIC_MAX_REL_SASA,
            "members": members,
            "n_candidates_scored": int(len(cands)), "n_controls_scored": int(len(controls)),
            "n_candidates_passing": int(len(passing(cands))),
            "n_controls_passing": int(len(passing(controls))),
            "n_candidates_cryptic": int(passing(cands)["cryptic"].sum()) if not cands.empty else 0,
            "n_controls_cryptic": int(passing(controls)["cryptic"].sum()) if not controls.empty else 0}


def balance(records: Sequence[dict]) -> Dict[str, object]:
    """How well the basic-count matching worked - the defect this study exists to fix."""
    out: Dict[str, object] = {}
    for query_arm, control_arm in QUERY_ARMS.items():
        q, c = arm_records(records, query_arm), arm_records(records, control_arm)
        if q.empty or c.empty or "n_basic" not in q or "n_basic" not in c:
            continue
        out[query_arm] = {"query_mean_basic": float(q["n_basic"].mean()),
                          "control_mean_basic": float(c["n_basic"].mean()),
                          "query_censored": int(q["censored"].astype(bool).sum()),
                          "control_censored": int(c["censored"].astype(bool).sum()),
                          "n_query": int(len(q)), "n_control": int(len(c))}
    return out


def build(records: Sequence[dict], n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    from cryptic_ip.benchmark import protocol

    o1 = guard(records, n_bootstrap=n_bootstrap)
    result: Dict[str, object] = {"seed": SEED, "n_bootstrap": n_bootstrap, "eps": EPS,
                                 "prune_tolerance": PRUNE_TOLERANCE, "O1": o1,
                                 "n_records": len(records), "balance": balance(records),
                                 "errors": sorted({str(r["error"]) for r in records if r.get("error")})}
    if o1.get("verdict") != "pass":
        result["O2"] = {"decision": "not run",
                        "reason": "the calibration guard did not pass, so no candidate is scored"}
        result["O3"] = {"reason": "not reported: the calibration guard did not pass", "members": []}
        return result
    o2 = paired(records, n_bootstrap=n_bootstrap)
    result["O2"] = o2
    result["O3"] = shortlist(records)
    pvalues = {}
    if "roc_auc" in o1.get("auc", {}):
        pvalues["O1"] = o1["auc"]["roc_auc"]["p_value"]
    if "per_group" in o2.get("estimate", {}):
        pvalues["O2"] = o2["estimate"]["per_group"]["p_value"]
    if pvalues:
        result["holm"] = protocol.holm(pvalues)
    return result


def markdown(r: Dict[str, object]) -> str:
    lines = ["# Study O: can a real IP6 site be transplanted into a buried candidate pocket?", "",
             f"Pre-registered in `docs/TEMPLATE2_PLAN.md`. Seed {r['seed']}, {r['n_bootstrap']} "
             "resamples of MMseqs2 30 % clusters.",
             f"Chebyshev acceptance at eps = {r['eps']} A with pruning at {r['prune_tolerance']} A, "
             "which `pruning_complete` proves discards no admissible placement.", ""]
    o1 = r["O1"]
    lines += ["## O1 - calibration guard", "", f"Verdict: **{o1.get('verdict')}**."]
    roc = o1.get("auc", {}).get("roc_auc")
    if roc:
        lines.append(f"AUC {roc['point']:.3f} [{roc['low']:.3f}, {roc['high']:.3f}] over "
                     f"{o1['n_binders']} annotated binders and {o1['n_controls']} matched controls "
                     f"in {o1['auc'].get('groups')} clusters.")
    if o1.get("reason"):
        lines.append(f"Reason: {o1['reason']}.")
    o2 = r["O2"]
    lines += ["", "## O2 - candidates against matched controls", "",
              f"Decision: **{o2.get('decision')}**."]
    per_group = o2.get("estimate", {}).get("per_group")
    if per_group:
        lines.append(f"Paired difference {per_group['point']:.3f} A "
                     f"[{per_group['low']:.3f}, {per_group['high']:.3f}] over {o2['n_pairs']} pairs "
                     f"in {o2['estimate'].get('groups')} clusters.")
    if o2.get("reason"):
        lines.append(f"Reason: {o2['reason']}.")
    o3 = r["O3"]
    lines += ["", "## O3 - template-compatible candidates", ""]
    if o3.get("reason"):
        lines.append(str(o3["reason"]) + ".")
    else:
        lines.append(f"Criteria, all absolute and fixed before any candidate was scored: clash-free, "
                     f"every matched anchor within {o3['eps']} A, and the placed ligand buried at "
                     f"relative SASA <= {o3['buried_max_relative_sasa']}.")
        lines.append(f"**{o3['n_candidates_passing']} of {o3['n_candidates_scored']} candidates** meet "
                     f"them, against **{o3['n_controls_passing']} of {o3['n_controls_scored']} "
                     f"matched controls**. At the stricter cryptic cut "
                     f"(<= {o3['cryptic_max_relative_sasa']}): {o3['n_candidates_cryptic']} candidates "
                     f"and {o3['n_controls_cryptic']} controls.")
        if o3["members"]:
            lines += ["", "| accession | worst anchor (A) | basic | template | relative SASA | cryptic |",
                      "| --- | --- | --- | --- | --- | --- |"]
            for m in o3["members"]:
                rel = m.get("relative_sasa")
                lines.append(f"| {m['uniprot_id']} | {float(m['deviation']):.3f} | {m.get('n_basic')} | "
                             f"{m.get('template')} | "
                             f"{'n/a' if rel is None else format(float(rel), '.3f')} | "
                             f"{'yes' if m.get('cryptic') else 'no'} |")
        lines += ["", "These are pockets that **can** host a real IP6 site in a buried position.",
                  "The plan fixes the word: they are template-compatible, not predicted binders."]
    if r.get("balance"):
        lines += ["", "## Basic-residue balance after matching", "",
                  "| arm | mean basic (query) | mean basic (control) | censored query | censored control |",
                  "| --- | --- | --- | --- | --- |"]
        for arm, b in r["balance"].items():
            lines.append(f"| {arm} | {b['query_mean_basic']:.2f} | {b['control_mean_basic']:.2f} | "
                         f"{b['query_censored']}/{b['n_query']} | {b['control_censored']}/{b['n_control']} |")
    if r.get("errors"):
        lines += ["", "## Queries not scored", ""] + [f"- {e}" for e in r["errors"]]
    return "\n".join(lines) + "\n"


def cmd_guard(args: argparse.Namespace) -> int:
    records = tf.read_json(Path(args.fits))
    result = guard(records, n_bootstrap=args.n_bootstrap)
    tf.write_json(Path(args.out), result)
    print("BEGIN_TEMPLATE2_GUARD_JSON")
    print(json.dumps(result, sort_keys=True, default=tf._json_default))
    print("END_TEMPLATE2_GUARD_JSON")
    if args.github_output:
        with open(args.github_output, "a") as fh:
            fh.write(f"verdict={result.get('verdict')}\n")
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    records = tf.read_json(Path(args.fits))
    result = build(records, n_bootstrap=args.n_bootstrap)
    tf.write_json(Path(args.out), result)
    Path(args.markdown).write_text(markdown(result))
    LOGGER.info("O1 %s, O2 %s", result["O1"].get("verdict"), result["O2"].get("decision"))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("templates", help="anchors, ligands and ligand atom identities")
    p.add_argument("--census", required=True)
    p.add_argument("--structures", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_templates)

    p = sub.add_parser("pool", help="each arm's queries plus a pool of candidate controls")
    p.add_argument("--proteins", required=True)
    p.add_argument("--candidates")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_pool)

    p = sub.add_parser("counts", help="basic-residue count per pocket")
    p.add_argument("--rows", required=True)
    p.add_argument("--af-dir", default="")
    p.add_argument("--list-only", action="store_true",
                   help="write this shard's accessions and exit, so the models can be fetched")
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=1)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_counts)

    p = sub.add_parser("select", help="keep the basic-count-closest control per query")
    p.add_argument("--pool", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_select)

    p = sub.add_parser("fit", help="score one arm against the template library")
    p.add_argument("--templates", required=True)
    p.add_argument("--arms", required=True)
    p.add_argument("--af-dir", required=True)
    p.add_argument("--arm", required=True, choices=list(ARMS))
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=1)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_fit)

    p = sub.add_parser("merge", help="concatenate fit shards")
    p.add_argument("--parts", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser("guard", help="O1 alone, to gate the candidate arm")
    p.add_argument("--fits", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    p.add_argument("--github-output", default=None)
    p.set_defaults(func=cmd_guard)

    p = sub.add_parser("report", help="O1, O2 and O3")
    p.add_argument("--fits", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--markdown", required=True)
    p.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    p.set_defaults(func=cmd_report)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
