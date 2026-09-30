#!/usr/bin/env python3
"""Study N (docs/TEMPLATE_PLAN.md): can a real IP6 site be transplanted onto a pocket?

The question studies A-L never answered is whether a candidate pocket is
geometrically capable of holding IP6. Docking cannot answer it here: study F
attributed 230 of 242 failures to scoring, and study G showed a 16x search budget
recovers nothing. So this study drops the scoring function entirely.

For every crystal IP6 copy in the redocking census, the basic residues that
coordinate its phosphates give a small constellation of anchor points. For a query
pocket, the same points are taken from its own K/R/H residues. A rigid map between
four anchors, by Kabsch superposition, carries the crystal ligand onto the query;
the placement is admissible only when no ligand heavy atom lands inside 2.2 A of a
protein heavy atom. The query's score is the lowest anchor RMSD over all admissible
placements and all eligible templates.

This is a necessary-condition test. Passing says the pocket *can* host IP6, which is
a hypothesis for an experiment; it is not evidence of binding, and the plan fixes the
word for it: template-compatible.

Subcommands::

    templates  anchors and ligand coordinates for every eligible IHP copy
    fit        score one arm's queries against the template library (sharded)
    merge      concatenate fit shards
    report     N1 (the calibration guard), N2, N3; printed between markers
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.docking.stats import MIN_GROUPS, auc_estimate, paired_difference  # noqa: E402

LOGGER = logging.getLogger("template_fit")

#: Every constant below is fixed by docs/TEMPLATE_PLAN.md and none may be tuned on results.
SEED = 20261006
N_BOOTSTRAP = 2000
#: A side-chain nitrogen this close to a ligand phosphate oxygen makes its residue an anchor.
ANCHOR_DISTANCE = 4.0
#: Anchors kept per template and per query, nearest the ligand and pocket centroid.
MAX_TEMPLATE_ANCHORS = 6
MAX_QUERY_ANCHORS = 10
#: Matched anchors per fit: the primary, and the secondary descriptor.
K_PRIMARY = 4
K_SECONDARY = 5
#: A correspondence survives pruning only if every pairwise distance agrees this closely.
DISTANCE_TOLERANCE = 1.5
#: A ligand heavy atom closer than this to a protein heavy atom is a clash.
CLASH_DISTANCE = 2.2
#: Fits worse than this are censored; the score is reported as -RMSD, so larger is better.
RMSD_CEILING = 4.0
#: N2's no-difference margin, in angstroms.
NO_DIFFERENCE = 0.25
#: A candidate is template-compatible above this percentile of the control arm.
COMPATIBLE_PERCENTILE = 95.0

BASIC_NITROGENS = {
    "LYS": ("NZ",),
    "ARG": ("NE", "NH1", "NH2"),
    "HIS": ("ND1", "NE2"),
}
TEMPLATE_COMP_ID = "IHP"


# ------------------------------------------------------------------ geometry
def kabsch_batch(mobile: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Optimal rigid transforms for a stack of point correspondences.

    Args:
        mobile: ``(n, k, 3)`` points to move.
        target: ``(n, k, 3)`` points to move onto.

    Returns:
        ``(rotation (n, 3, 3), translation (n, 3), rmsd (n,))`` such that
        ``mobile @ rotation + translation`` is the fitted placement.
    """
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    mob_c = mobile.mean(axis=1, keepdims=True)
    tar_c = target.mean(axis=1, keepdims=True)
    p = mobile - mob_c
    q = target - tar_c
    cov = np.einsum("nki,nkj->nij", p, q)
    u, _, vt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(np.einsum("nij,njk->nik", u, vt)))
    correction = np.repeat(np.eye(3)[None, :, :], len(mobile), axis=0)
    correction[:, 2, 2] = d
    rotation = np.einsum("nij,njk,nkl->nil", u, correction, vt)
    fitted = np.einsum("nki,nij->nkj", p, rotation)
    rmsd = np.sqrt(np.mean(np.sum((fitted - q) ** 2, axis=2), axis=1))
    translation = tar_c[:, 0, :] - np.einsum("ni,nij->nj", mob_c[:, 0, :], rotation)
    return rotation, translation, rmsd


def correspondences(template_xyz: np.ndarray, query_xyz: np.ndarray, k: int,
                    tolerance: float = DISTANCE_TOLERANCE) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
    """Injective template-to-query anchor maps whose pairwise distances all agree.

    Pruning on distances is exact for a rigid map: a superposition cannot make two
    points match if the distance between them differs by more than the tolerance, so
    nothing admissible is lost and the enumeration stays small.
    """
    template_xyz = np.asarray(template_xyz, dtype=float)
    query_xyz = np.asarray(query_xyz, dtype=float)
    t, q = len(template_xyz), len(query_xyz)
    if t < k or q < k:
        return []
    td = np.linalg.norm(template_xyz[:, None, :] - template_xyz[None, :, :], axis=2)
    qd = np.linalg.norm(query_xyz[:, None, :] - query_xyz[None, :, :], axis=2)
    out: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = []
    for subset in combinations(range(t), k):
        stack: List[List[int]] = [[]]
        while stack:
            partial = stack.pop()
            depth = len(partial)
            if depth == k:
                out.append((subset, tuple(partial)))
                continue
            for cand in range(q):
                if cand in partial:
                    continue
                if all(abs(td[subset[depth], subset[i]] - qd[cand, partial[i]]) <= tolerance
                       for i in range(depth)):
                    stack.append(partial + [cand])
    return out


def clash_free(ligand_xyz: np.ndarray, tree, cutoff: float = CLASH_DISTANCE) -> bool:
    """True when no ligand heavy atom lies within ``cutoff`` of a protein heavy atom."""
    return not np.any(tree.query_ball_point(ligand_xyz, cutoff, return_length=True))


#: Placements clash-tested per batch. The batch is an implementation detail, not a cap:
#: candidates are tested in ascending RMSD, so the first admissible one found is still
#: the global minimum. Batching only keeps the coordinate array bounded.
CLASH_BATCH = 4096


def best_fit(template_xyz: np.ndarray, ligand_xyz: np.ndarray, query_xyz: np.ndarray, tree,
             *, k: int = K_PRIMARY, ceiling: float = RMSD_CEILING) -> Optional[Dict[str, object]]:
    """Lowest-RMSD admissible placement of one template on one query pocket.

    Correspondences are ranked by anchor RMSD and clash-tested in that order, so the
    first admissible one *is* the minimum and the rest need no test.
    """
    maps = correspondences(template_xyz, query_xyz, k)
    if not maps:
        return None
    mobile = np.stack([template_xyz[list(t)] for t, _ in maps])
    target = np.stack([query_xyz[list(q)] for _, q in maps])
    rotation, translation, rmsd = kabsch_batch(mobile, target)
    order = np.argsort(rmsd)
    order = order[rmsd[order] <= ceiling]
    for start in range(0, len(order), CLASH_BATCH):
        batch = order[start:start + CLASH_BATCH]
        placed = np.einsum("ai,nij->naj", ligand_xyz, rotation[batch]) + translation[batch][:, None, :]
        hits = tree.query_ball_point(placed.reshape(-1, 3), CLASH_DISTANCE, return_length=True)
        admissible = ~np.any(hits.reshape(len(batch), -1) > 0, axis=1)
        if not admissible.any():
            continue
        idx = int(batch[np.flatnonzero(admissible)[0]])
        return {"rmsd": float(rmsd[idx]), "n_anchors": int(k),
                "template_anchors": [int(i) for i in maps[idx][0]],
                "query_anchors": [int(i) for i in maps[idx][1]],
                "ligand_xyz": ligand_xyz @ rotation[idx] + translation[idx]}
    return None


# ------------------------------------------------------------------ anchors
def residue_points(arrays, residues: Iterable[Tuple[str, int]]) -> Tuple[np.ndarray, List[str]]:
    """One point per basic residue: the mean of its own side-chain nitrogens."""
    wanted = {(str(c), int(r)) for c, r in residues}
    points, labels = [], []
    for key in np.unique(arrays.residue_index):
        sel = arrays.residue_index == key
        resname = str(arrays.resnames[sel][0]).upper()
        chain = str(arrays.chain_ids[sel][0])
        resseq = int(arrays.resseqs[sel][0])
        if resname not in BASIC_NITROGENS or (wanted and (chain, resseq) not in wanted):
            continue
        names = BASIC_NITROGENS[resname]
        atoms = sel & np.isin(np.char.strip(arrays.atom_names.astype(str)), names)
        if not atoms.any():
            continue
        points.append(arrays.coords[atoms].mean(axis=0))
        labels.append(f"{chain}:{resname}{resseq}")
    if not points:
        return np.zeros((0, 3)), []
    return np.asarray(points, dtype=float), labels


def nearest(points: np.ndarray, labels: List[str], centre: np.ndarray, limit: int) -> Tuple[np.ndarray, List[str]]:
    """The ``limit`` points closest to ``centre``, which bounds the correspondence search."""
    if len(points) <= limit:
        return points, labels
    order = np.argsort(np.linalg.norm(points - centre, axis=1))[:limit]
    return points[order], [labels[i] for i in order]


def template_anchors(arrays, ligand_atoms: np.ndarray) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Anchor points, their labels, and the ligand heavy atoms, for one crystal copy.

    A residue anchors the ligand when one of its side-chain nitrogens comes within
    :data:`ANCHOR_DISTANCE` of a *phosphate oxygen* - the oxygens are what the basic
    side chains actually contact, and using them rather than the whole ligand keeps a
    residue packed against the inositol ring from counting as an anchor.
    """
    heavy = ligand_atoms[arrays.elements[ligand_atoms] != "H"]
    ligand_xyz = arrays.coords[heavy]
    phosphorus = arrays.coords[heavy[arrays.elements[heavy] == "P"]]
    oxygens = arrays.coords[heavy[arrays.elements[heavy] == "O"]]
    if len(phosphorus) and len(oxygens):
        near_p = np.linalg.norm(oxygens[:, None, :] - phosphorus[None, :, :], axis=2).min(axis=1) <= 1.8
        oxygens = oxygens[near_p]
    if not len(oxygens):
        return np.zeros((0, 3)), [], ligand_xyz
    points, labels = residue_points(arrays, ())
    if not len(points):
        return np.zeros((0, 3)), [], ligand_xyz
    close = np.linalg.norm(points[:, None, :] - oxygens[None, :, :], axis=2).min(axis=1) <= ANCHOR_DISTANCE
    points, labels = points[close], [labels[i] for i in np.flatnonzero(close)]
    points, labels = nearest(points, labels, ligand_xyz.mean(axis=0), MAX_TEMPLATE_ANCHORS)
    return points, labels, ligand_xyz


# ------------------------------------------------------------------ I/O
def write_json(path: Path, payload: object) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True, default=_json_default)
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


def read_json(path: Path) -> object:
    raw = Path(path).read_bytes()
    if str(path).endswith(".gz"):
        raw = gzip.decompress(raw)
    return json.loads(raw)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(repr(o))


def _text(value: object) -> str:
    """A CSV cell as text, with pandas' NaN read as empty rather than the string "nan"."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return str(value).strip()


def shard(items: Sequence, index: int, count: int) -> List:
    return [x for i, x in enumerate(items) if i % max(count, 1) == index]


def parse_positions(text: object) -> List[int]:
    """Pocket residue numbers, as the screen wrote them."""
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return []
    return [int(float(p)) for p in str(text).replace(";", ",").split(",") if p.strip()]


# ------------------------------------------------------------------ templates
def cmd_templates(args: argparse.Namespace) -> int:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.validation.burial_metrics import find_ligand_instances

    census = pd.read_csv(args.census)
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
        # Exact name matching, not structural re-detection: the census already chose
        # these copies by name, and re-detecting could return a differently named one.
        found = {(key[1], int(key[2]), key[3] or ""): atoms
                 for key, comp, atoms in find_ligand_instances(arrays, comp_ids=(TEMPLATE_COMP_ID,))}
        for _, row in block.iterrows():
            key = (_text(row["chain"]), int(row["resseq"]), _text(row.get("icode")))
            atoms = found.get(key)
            if atoms is None:
                dropped.append({"copy_key": row["copy_key"], "reason": "copy not found in structure"})
                continue
            points, labels, ligand_xyz = template_anchors(arrays, atoms)
            if len(points) < K_PRIMARY:
                dropped.append({"copy_key": row["copy_key"], "reason": f"{len(points)} anchors"})
                continue
            out.append({"copy_key": row["copy_key"], "pdb_id": pdb_id,
                        "homology_group_strict": row.get("homology_group_strict"),
                        "uniprot_ids": _text(row.get("uniprot_ids")),
                        "burial_class": row.get("burial_class"),
                        "anchors": points, "anchor_labels": labels, "ligand": ligand_xyz})
    LOGGER.info("%d templates kept, %d dropped", len(out), len(dropped))
    write_json(Path(args.out), {"templates": out, "dropped": dropped,
                                "n_kept": len(out), "n_dropped": len(dropped)})
    return 0


# ------------------------------------------------------------------ fitting
def eligible_templates(templates: Sequence[dict], accession: str, cluster: object) -> Tuple[List[dict], Dict[str, int]]:
    """Templates left after the plan's three leakage rules, with the counts removed by each."""
    kept, removed = [], {"same_accession": 0, "same_cluster": 0}
    cluster_key = None if cluster is None or (isinstance(cluster, float) and np.isnan(cluster)) else str(cluster)
    for t in templates:
        accessions = {a.strip() for a in _text(t.get("uniprot_ids")).replace(";", ",").split(",") if a.strip()}
        if accession in accessions:
            removed["same_accession"] += 1
            continue
        if cluster_key is not None and cluster_key in accessions:
            removed["same_cluster"] += 1
            continue
        kept.append(t)
    return kept, removed


def score_query(templates: Sequence[dict], arrays, residues: Sequence[int], *, k: int = K_PRIMARY) -> Dict[str, object]:
    """The query's fit score: the lowest admissible anchor RMSD over every template."""
    from scipy.spatial import cKDTree

    wanted = {(str(c), int(r)) for c in set(arrays.chain_ids.astype(str)) for r in residues}
    points, labels = residue_points(arrays, wanted)
    if len(points) < k:
        return {"score": -RMSD_CEILING, "rmsd": RMSD_CEILING, "n_query_anchors": int(len(points)),
                "censored": True, "reason": "fewer basic residues in the pocket than matched anchors"}
    centre = points.mean(axis=0)
    points, labels = nearest(points, labels, centre, MAX_QUERY_ANCHORS)
    heavy = arrays.is_polymer & (arrays.elements != "H")
    tree = cKDTree(arrays.coords[heavy])
    best: Optional[Dict[str, object]] = None
    for t in templates:
        fit = best_fit(np.asarray(t["anchors"], dtype=float), np.asarray(t["ligand"], dtype=float),
                       points, tree, k=k)
        if fit is None:
            continue
        if best is None or fit["rmsd"] < best["rmsd"]:
            placed = fit.pop("ligand_xyz")
            fit["template"] = t["copy_key"]
            fit["template_burial"] = t.get("burial_class")
            fit["query_anchor_labels"] = [labels[i] for i in fit["query_anchors"]]
            fit["buried_fraction"] = float(np.mean(
                tree.query_ball_point(placed, 8.0, return_length=True) >= 60))
            best = fit
    if best is None:
        return {"score": -RMSD_CEILING, "rmsd": RMSD_CEILING, "n_query_anchors": int(len(points)),
                "censored": True, "reason": "no admissible placement within the ceiling"}
    best.update({"score": -best["rmsd"], "n_query_anchors": int(len(points)), "censored": False})
    return best


def arm_frame(args: argparse.Namespace) -> pd.DataFrame:
    """The queries of one arm, as accession / cluster / pocket-residue rows."""
    sys.path.insert(0, str(ROOT / "scripts"))
    import importlib.util

    spec = importlib.util.spec_from_file_location("triage_mod", ROOT / "scripts" / "triage.py")
    triage = importlib.util.module_from_spec(spec)
    sys.modules["triage_mod"] = triage
    spec.loader.exec_module(triage)

    proteins = pd.read_csv(args.proteins)
    if args.arm in ("candidates", "controls"):
        cands = pd.read_csv(args.candidates)
        cands = cands[["uniprot_id", "organism_key", "top_pocket_residues", "plddt_mean", "hull_depth",
                       "combined", "cluster"]]
    else:
        cands = proteins[proteins["annotated"].astype(str).str.lower() == "true"].copy()
    controls = triage.matched_controls(cands, proteins)
    # A protein drawn as a control that is also an annotated binder counts as a control only.
    clash = set(controls["uniprot_id"]) & set(cands["uniprot_id"])
    if clash:
        LOGGER.info("%d queries are also drawn as controls and count as controls only: %s",
                    len(clash), sorted(clash))
        cands = cands[~cands["uniprot_id"].isin(clash)]
    if args.arm == "controls" or args.arm == "annotated_controls":
        frame = controls
    else:
        frame = cands.assign(role="query", matched_to=cands["uniprot_id"])
    keep = ["uniprot_id", "top_pocket_residues", "cluster", "matched_to"]
    frame = frame[[c for c in keep if c in frame.columns]].copy()
    frame["arm"] = args.arm
    return frame


def cmd_fit(args: argparse.Namespace) -> int:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    templates = read_json(Path(args.templates))["templates"]
    frame = arm_frame(args)
    rows = shard(list(frame.to_dict("records")), args.shard, args.shards)
    LOGGER.info("arm %s: %d of %d queries in shard %d/%d", args.arm, len(rows), len(frame),
                args.shard, args.shards)
    out = []
    for row in rows:
        acc = str(row["uniprot_id"])
        record = {"uniprot_id": acc, "arm": args.arm, "cluster": row.get("cluster"),
                  "matched_to": row.get("matched_to")}
        path = next(iter(sorted(Path(args.af_dir).glob(f"AF-{acc}-F1-*.pdb"))), None)
        if path is None:
            record.update({"score": None, "error": "no AlphaFold model"})
            out.append(record)
            continue
        usable, removed = eligible_templates(templates, acc, row.get("cluster"))
        record["templates_removed"] = removed
        record["n_templates"] = len(usable)
        try:
            arrays = load_structure_arrays(path)
            record.update(score_query(usable, arrays, parse_positions(row.get("top_pocket_residues"))))
        except Exception as exc:  # a structure this study cannot read is recorded, not fatal
            record.update({"score": None, "error": f"{type(exc).__name__}: {exc}"})
        out.append(record)
        LOGGER.info("%s %s score=%s", args.arm, acc, record.get("score"))
    write_json(Path(args.out), out)
    return 0


def cmd_merge(args: argparse.Namespace) -> int:
    records: List[dict] = []
    for part in sorted(Path(args.parts).rglob("*.json*")):
        payload = read_json(part)
        if isinstance(payload, list):
            records.extend(payload)
    LOGGER.info("%d fit records merged", len(records))
    write_json(Path(args.out), records)
    return 0


# ------------------------------------------------------------------ report
def arm_records(records, arm: str) -> pd.DataFrame:
    frame = pd.DataFrame([r for r in records if r.get("arm") == arm])
    if frame.empty or "score" not in frame.columns:
        return pd.DataFrame()
    return frame[frame["score"].notna()].copy()


def groups_of(frame: pd.DataFrame) -> List[str]:
    """Resampling groups: the MMseqs2 30 % cluster, or the accession when it has none."""
    return frame["cluster"].fillna(frame["uniprot_id"]).astype(str).tolist()


def guard(records, n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    """N1: does the fit score separate annotated binders from their matched controls?"""
    binders = arm_records(records, "annotated")
    controls = arm_records(records, "annotated_controls")
    if binders.empty or controls.empty:
        return {"verdict": "not evaluable", "reason": "an arm has no scored query",
                "n_binders": int(len(binders)), "n_controls": int(len(controls))}
    frame = pd.concat([binders.assign(label=1), controls.assign(label=0)], ignore_index=True)
    auc = auc_estimate(frame["label"].tolist(), frame["score"].tolist(), groups_of(frame),
                       n_bootstrap=n_bootstrap, seed=SEED)
    out = {"auc": auc, "n": int(len(frame)), "n_binders": int(len(binders)),
           "n_controls": int(len(controls))}
    if not auc.get("evidence") or "roc_auc" not in auc:
        out["verdict"] = "not evaluable"
        out["reason"] = (f"{auc.get('groups')} clusters, fewer than the {MIN_GROUPS} "
                         f"this project requires for evidence")
        return out
    out["verdict"] = "pass" if auc["roc_auc"]["low"] > 0.5 else "fail"
    return out


def paired(records, n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    """N2: the paired candidate-minus-control difference in fit score."""
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
                            groups_of(partner.reset_index()),
                            n_bootstrap=n_bootstrap, seed=SEED)
    out = {"estimate": est, "n_pairs": int(len(pairs))}
    if not est.get("evidence") or "per_group" not in est:
        out["decision"] = "not evaluable"
        out["reason"] = (f"{est.get('groups')} clusters, fewer than the {MIN_GROUPS} "
                         f"this project requires for evidence")
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


def shortlist(records) -> Dict[str, object]:
    """N3: the candidates above the control arm's 95th percentile, reported descriptively."""
    cands = arm_records(records, "candidates")
    controls = arm_records(records, "controls")
    if cands.empty or controls.empty:
        return {"threshold": None, "members": [], "reason": "not reported: an arm has no scored query"}
    threshold = float(np.percentile(controls["score"].to_numpy(dtype=float), COMPATIBLE_PERCENTILE))
    fields = ("uniprot_id", "score", "rmsd", "n_anchors", "template", "template_burial",
              "buried_fraction", "query_anchor_labels")
    members = [{k: row.get(k) for k in fields}
               for _, row in cands[cands["score"] > threshold].sort_values("score", ascending=False).iterrows()]
    return {"threshold": threshold, "percentile": COMPATIBLE_PERCENTILE, "members": members,
            "n_candidates_scored": int(len(cands)), "n_controls_scored": int(len(controls))}


def build(records, n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    from cryptic_ip.benchmark import protocol

    n1 = guard(records, n_bootstrap=n_bootstrap)
    result: Dict[str, object] = {"seed": SEED, "n_bootstrap": n_bootstrap, "N1": n1,
                                 "n_records": len(records),
                                 "errors": sorted({str(r["error"]) for r in records if r.get("error")})}
    if n1.get("verdict") != "pass":
        result["N2"] = {"decision": "not run",
                        "reason": "the calibration guard did not pass, so no candidate is scored"}
        result["N3"] = {"reason": "not reported: the calibration guard did not pass", "members": []}
        return result
    n2 = paired(records, n_bootstrap=n_bootstrap)
    result["N2"] = n2
    result["N3"] = shortlist(records)
    pvalues = {}
    if "roc_auc" in n1.get("auc", {}):
        pvalues["N1"] = n1["auc"]["roc_auc"]["p_value"]
    if "per_group" in n2.get("estimate", {}):
        pvalues["N2"] = n2["estimate"]["per_group"]["p_value"]
    if pvalues:
        result["holm"] = protocol.holm(pvalues)
    return result


def markdown(r: Dict[str, object]) -> str:
    lines = ["# Study N: can a real IP6 site be transplanted onto a candidate pocket?", "",
             "Pre-registered in `docs/TEMPLATE_PLAN.md`. Seed "
             f"{r['seed']}, {r['n_bootstrap']} resamples of MMseqs2 30 % clusters.", ""]
    n1 = r["N1"]
    lines += ["## N1 - calibration guard", "", f"Verdict: **{n1.get('verdict')}**."]
    roc = n1.get("auc", {}).get("roc_auc")
    if roc:
        lines.append(f"AUC {roc['point']:.3f} [{roc['low']:.3f}, {roc['high']:.3f}] over "
                     f"{n1['n_binders']} annotated binders and {n1['n_controls']} matched controls "
                     f"in {n1['auc'].get('groups')} clusters.")
    if n1.get("reason"):
        lines.append(f"Reason: {n1['reason']}.")
    n2 = r["N2"]
    lines += ["", "## N2 - candidates against matched controls", "",
              f"Decision: **{n2.get('decision')}**."]
    per_group = n2.get("estimate", {}).get("per_group")
    if per_group:
        lines.append(f"Paired difference {per_group['point']:.3f} A "
                     f"[{per_group['low']:.3f}, {per_group['high']:.3f}] over {n2['n_pairs']} pairs "
                     f"in {n2['estimate'].get('groups')} clusters.")
    if n2.get("reason"):
        lines.append(f"Reason: {n2['reason']}.")
    n3 = r["N3"]
    lines += ["", "## N3 - template-compatible candidates", ""]
    if n3.get("reason"):
        lines.append(str(n3["reason"]) + ".")
    else:
        lines.append(f"Threshold: the {n3['percentile']:.0f}th percentile of "
                     f"{n3['n_controls_scored']} scored controls, at {n3['threshold']:.3f}.")
        lines += ["", "| accession | score | RMSD (A) | template | buried fraction |",
                  "| --- | --- | --- | --- | --- |"]
        for m in n3["members"]:
            rmsd = m.get("rmsd")
            rmsd = -float(m["score"]) if rmsd is None else float(rmsd)
            lines.append(f"| {m['uniprot_id']} | {float(m['score']):.3f} | {rmsd:.3f} | "
                         f"{m.get('template')} | {float(m.get('buried_fraction') or 0.0):.2f} |")
        lines += ["", "These are pockets that **can** host IP6 under a rigid transplant of a real site.",
                  "The plan fixes the word: they are template-compatible, not predicted binders."]
    if r.get("errors"):
        lines += ["", "## Queries not scored", ""] + [f"- {e}" for e in r["errors"]]
    return "\n".join(lines) + "\n"


def cmd_accessions(args: argparse.Namespace) -> int:
    """The accessions one arm's shard needs, so the models can be fetched before fitting."""
    rows = shard(list(arm_frame(args).to_dict("records")), args.shard, args.shards)
    accessions = sorted({str(r["uniprot_id"]) for r in rows})
    Path(args.out).write_text("\n".join(accessions) + "\n")
    LOGGER.info("%d accessions for arm %s shard %d/%d", len(accessions), args.arm, args.shard, args.shards)
    return 0


def cmd_guard(args: argparse.Namespace) -> int:
    """N1 alone, so the workflow can gate the candidate arm on its verdict."""
    records = read_json(Path(args.fits))
    result = guard(records, n_bootstrap=args.n_bootstrap)
    write_json(Path(args.out), result)
    print("BEGIN_TEMPLATE_GUARD_JSON")
    print(json.dumps(result, sort_keys=True, default=_json_default))
    print("END_TEMPLATE_GUARD_JSON")
    if args.github_output:
        with open(args.github_output, "a") as fh:
            fh.write(f"verdict={result.get('verdict')}\n")
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    records = read_json(Path(args.fits))
    result = build(records, n_bootstrap=args.n_bootstrap)
    write_json(Path(args.out), result)
    Path(args.markdown).write_text(markdown(result))
    # The blocks are emitted by the workflow, markdown first and JSON last: the log
    # tool returns a bounded *tail*, so the large base64 block must not come after
    # the numbers that have to be read back.
    LOGGER.info("N1 %s, N2 %s", result["N1"].get("verdict"), result["N2"].get("decision"))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("templates", help="anchors and ligand coordinates for every eligible IHP copy")
    p.add_argument("--census", required=True)
    p.add_argument("--structures", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_templates)

    p = sub.add_parser("fit", help="score one arm's queries against the template library")
    p.add_argument("--templates", required=True)
    p.add_argument("--proteins", required=True)
    p.add_argument("--candidates")
    p.add_argument("--af-dir", required=True)
    p.add_argument("--arm", required=True,
                   choices=["annotated", "annotated_controls", "candidates", "controls"])
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=1)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_fit)

    p = sub.add_parser("merge", help="concatenate fit shards")
    p.add_argument("--parts", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser("accessions", help="the accessions one arm's shard needs")
    p.add_argument("--proteins", required=True)
    p.add_argument("--candidates")
    p.add_argument("--arm", required=True,
                   choices=["annotated", "annotated_controls", "candidates", "controls"])
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=1)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_accessions)

    p = sub.add_parser("guard", help="N1 alone, to gate the candidate arm")
    p.add_argument("--fits", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    p.add_argument("--github-output", default=None)
    p.set_defaults(func=cmd_guard)

    p = sub.add_parser("report", help="N1, N2 and N3")
    p.add_argument("--fits", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--markdown", required=True)
    p.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    p.set_defaults(func=cmd_report)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
