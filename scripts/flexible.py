#!/usr/bin/env python3
"""Study J (docs/FLEXIBLE_PLAN.md): does a flexible receptor recover study A's failures?

Every docking run in studies A, F and G held the receptor rigid, so each pocket's side
chains sat in whatever rotamer the crystal caught them in. Study F's 0.110 top-pose success
is therefore optimistic for prospective use and may understate what the protocol could do
with rotamer freedom. This study lets the pocket's side chains move and measures the paired
difference.

Two pieces of machinery need care, and both are tested:

* **Torsion trees.** A flexible Vina run needs a second PDBQT holding each movable side
  chain as a ``ROOT``/``BRANCH`` tree. Meeko builds those; see
  ``docs/FLEXIBLE_PLAN_AMENDMENT_1.md`` for why both arms are therefore docked through the
  Meeko receptor rather than reusing study G's rigid arm.
* **Pose parsing.** With flexible residues, Vina returns the moved side chains *inside each
  pose model*, between ``BEGIN_RES`` and ``END_RES``. Those atoms are not ligand atoms and
  would corrupt every RMSD if they reached the pose reader, so they are stripped first.

``cryptic_ip/docking/`` is deliberately not modified: its paths trigger the full redocking
benchmark, and studies A, F and G must stay reproducible. The Vina call here mirrors
``engine.dock`` and hands its result to ``engine.pose_table``, so scoring and RMSD are the
shared code.

Subcommands::

    residues  which side chains each copy makes flexible, and why
    dock      both arms for a shard of copies
    report    J1, J2 and J3; printed between markers
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

LOGGER = logging.getLogger("flexible")

#: Fixed by docs/FLEXIBLE_PLAN.md; none may be tuned on results.
SEED = 20261003
N_BOOTSTRAP = 2000
EXHAUSTIVENESS = 32
N_POSES = 40
ENERGY_RANGE = 10.0
SUCCESS_RMSD = 2.0
NO_MATERIAL_CHANGE = 0.05
#: A residue is flexible when a side-chain heavy atom comes within this of the ligand.
FLEX_RADIUS = 4.0
#: Vina's flexible sampling degrades beyond roughly this many movable side chains.
MAX_FLEX = 8
#: Amendment 2: a residue Meeko cannot type is deleted only this far beyond the box
#: edge, which is Vina's own interaction cutoff, so no pose can feel it.
DROP_RADIUS = 8.0
#: No rotatable side chain to speak of, so never made flexible.
NO_SIDECHAIN = ("GLY", "ALA", "PRO")
#: A CYS whose SG is this close to another SG is in a disulphide and stays fixed.
DISULPHIDE = 2.5
BACKBONE = ("N", "CA", "C", "O", "OXT")
ARMS = ("rigid", "flex")


# ------------------------------------------------------------------ residue choice
def _text(value: object) -> str:
    """A CSV cell as text, with pandas' NaN read as empty rather than the string "nan"."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return str(value).strip()


def flexible_residues(arrays, ligand_xyz: np.ndarray, *, radius: float = FLEX_RADIUS,
                      cap: int = MAX_FLEX) -> List[Dict[str, object]]:
    """The plan's rule, in order: within ``radius``, side-chain-bearing, not a disulphide CYS.

    Ordered by minimum heavy-atom distance to the ligand and capped at ``cap`` residues
    nearest the ligand centroid, both fixed by the plan.
    """
    from scipy.spatial import cKDTree

    ligand_xyz = np.asarray(ligand_xyz, dtype=float)
    heavy = arrays.is_polymer & (arrays.elements != "H")
    names = np.char.strip(arrays.atom_names.astype(str))
    sidechain = heavy & ~np.isin(names, BACKBONE)
    if not sidechain.any():
        return []
    tree = cKDTree(ligand_xyz)
    idx = np.flatnonzero(sidechain)
    distances, _ = tree.query(arrays.coords[idx])

    sulphurs = np.flatnonzero(heavy & (names == "SG"))
    bonded = set()
    if len(sulphurs) > 1:
        sg_tree = cKDTree(arrays.coords[sulphurs])
        for a, b in sg_tree.query_pairs(DISULPHIDE):
            bonded.add(int(arrays.residue_index[sulphurs[a]]))
            bonded.add(int(arrays.residue_index[sulphurs[b]]))

    best: Dict[int, Dict[str, object]] = {}
    for position, atom in enumerate(idx):
        if distances[position] > radius:
            continue
        key = int(arrays.residue_index[atom])
        resname = str(arrays.resnames[atom]).upper()
        if resname in NO_SIDECHAIN:
            continue
        if resname == "CYS" and key in bonded:
            continue
        record = best.get(key)
        if record is None or distances[position] < record["min_distance"]:
            best[key] = {"residue_index": key, "resname": resname,
                         "chain": str(arrays.chain_ids[atom]),
                         "resseq": int(arrays.resseqs[atom]),
                         "icode": str(arrays.icodes[atom]).strip(),
                         "min_distance": float(distances[position])}
    chosen = sorted(best.values(), key=lambda r: (r["min_distance"], r["chain"], r["resseq"]))
    if len(chosen) <= cap:
        return chosen
    # The cap keeps the residues nearest the ligand *centroid*, as the plan states, while
    # the returned order stays by minimum heavy-atom distance.
    centre = ligand_xyz.mean(axis=0)
    by_centroid = sorted(chosen, key=lambda r: (
        float(np.linalg.norm(arrays.coords[arrays.residue_index == r["residue_index"]].mean(axis=0) - centre)),
        r["chain"], r["resseq"]))
    keep = {r["residue_index"] for r in by_centroid[:cap]}
    return [r for r in chosen if r["residue_index"] in keep]


def meeko_residue_id(record: Dict[str, object]) -> str:
    """Meeko addresses monomers as ``chain:resseq`` (with an insertion code if present)."""
    icode = str(record.get("icode") or "")
    return f"{record['chain']}:{record['resseq']}{icode}"


# ------------------------------------------------------------------ receptor
def canonical_pqr(text: str) -> str:
    """Re-space a PDB2PQR file so Meeko's PQR reader cannot mis-split its fields.

    PDB2PQR writes the PQR on PDB columns, where the chain occupies column 22 and the
    residue number columns 23 to 26. A four-digit residue number therefore fills its field
    and the two run together - ``GLY A2401`` - while Meeko reads the PQR by splitting on
    whitespace::

        token = items.pop(0)
        try: resnum = int(token)          # "A2401" is not an int
        except ValueError:
            chainid = token               # so the chain becomes "A2401"
            resnum = int(items.pop(0))    # and the x coordinate becomes the residue number

    which raises ``ValueError: invalid literal for int() with base 10: '0.000'``. That is
    four of the eight failures in shard 1 of run 36828409167, and exactly the four copies
    whose residue numbers reach four digits. Studies A, F and G never saw it because they
    read the same file by column.

    Had the int() happened to succeed the damage would have been worse than a crash: the
    chain would have been "A2401" and every monomer key wrong, so no flexible residue would
    have been found and both arms would have docked rigidly under a "flex" label. The
    fields go out separated by single spaces, read off the columns PDB2PQR wrote them in.
    """
    out: List[str] = []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            out.append(line)
            continue
        fields = _pqr_fields(line)
        if fields is None:  # not the layout we expect; leave it for Meeko to read as before
            out.append(line)
            continue
        record, serial, name, resname, chain, resseq, icode, tail = fields
        head = [record, serial, name, resname]
        if chain:
            head.append(chain)
        head.append(resseq)
        if icode:
            head.append(icode)
        out.append(" ".join(head + tail))
    return "\n".join(out) + ("\n" if out else "")


def _pqr_fields(line: str):
    """``(record, serial, name, resname, chain, resseq, icode, [x, y, z, charge, radius])``.

    ``None`` when the line does not read as PDB columns, so a PQR written some other way is
    passed through untouched rather than mangled.
    """
    if len(line) < 54:
        return None
    serial, name, resname = line[6:11].strip(), line[12:16].strip(), line[17:20].strip()
    chain, resseq, icode = line[21:22].strip(), line[22:26].strip(), line[26:27].strip()
    try:
        int(resseq)
        xyz = [f"{float(line[30 + 8 * i:38 + 8 * i]):.3f}" for i in range(3)]
    except ValueError:
        return None
    rest = line[54:].split()
    if len(rest) != 2 or not name or not resname:
        return None
    try:
        charge, radius = (f"{float(value):.4f}" for value in rest)
    except ValueError:
        return None
    return (line[:6].strip(), serial, name, resname, chain, resseq, icode,
            xyz + [charge, radius])


H_ANOMALY = re.compile(r"Residue (\S+) matched with template")


def residue_box_distance(pqr_text: str, box_centre: Sequence[float],
                         box_size: Sequence[float]) -> Dict[str, float]:
    """Each residue's smallest distance from the docking box, by ``chain:resseq`` key.

    Zero for a residue with an atom inside the box. Read off the canonical PQR, so it is
    the same coordinates Meeko is given.
    """
    low = np.asarray(box_centre, dtype=float) - np.asarray(box_size, dtype=float) / 2.0
    high = np.asarray(box_centre, dtype=float) + np.asarray(box_size, dtype=float) / 2.0
    best: Dict[str, float] = {}
    for line in pqr_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        fields = line.split()
        if len(fields) < 10:
            continue
        chain, resseq = fields[4], fields[5]
        try:
            xyz = np.asarray([float(v) for v in fields[-5:-2]], dtype=float)
        except ValueError:
            continue
        gap = float(np.linalg.norm(np.maximum(np.maximum(low - xyz, xyz - high), 0.0)))
        key = f"{chain}:{resseq}"
        if gap < best.get(key, np.inf):
            best[key] = gap
    return best


def prepare_pair(protonated: Tuple[Path, Path], workdir: Path,
                 flexres: Sequence[Dict[str, object]], *,
                 box_centre: Optional[Sequence[float]] = None,
                 box_size: Optional[Sequence[float]] = None,
                 drop_radius: float = DROP_RADIUS
                 ) -> Tuple[Path, Optional[Path], List[str], List[str]]:
    """Rigid and flexible PDBQT from the project's own PDB2QR output, via Meeko.

    With ``flexres`` empty the writer emits an empty flex string and ``None`` is returned in
    its place, so the rigid arm takes the single-receptor path - the pinning the plan asks
    for, so that "rigid" means rigid.

    Meeko refuses a receptor outright over residues elsewhere in the protein: ones it has
    no template for, and ones whose PQR hydrogens disagree with the template it matched,
    which makes the PQR charges inapplicable. Three of shard 1's ten copies died that way
    in run 36828409167. Given the box, such a residue is deleted when it lies more than
    ``drop_radius`` beyond the box edge, and the copy still fails when one lies closer; see
    docs/FLEXIBLE_PLAN_AMENDMENT_2.md. The deleted ids are returned, not swallowed.
    """
    from meeko import MoleculePreparation, PDBQTWriterLegacy, Polymer

    workdir.mkdir(parents=True, exist_ok=True)
    _pdb, pqr = protonated
    mk = MoleculePreparation()
    text = canonical_pqr(Path(pqr).read_text())
    boxed = box_centre is not None and box_size is not None
    options: Dict[str, object] = {}
    if boxed:
        options = {"box_center": [float(c) for c in box_centre],
                   "box_size": [float(s) for s in box_size],
                   "delete_bad_res_from_box_radius": float(drop_radius)}
    dropped: List[str] = []
    try:
        polymer = Polymer.from_pqr_string(text, mk_prep=mk, **options)
    except Exception as exc:  # noqa: BLE001 - one recoverable shape, re-raised otherwise
        anomalous = sorted(set(H_ANOMALY.findall(str(exc))))
        if not boxed or not anomalous:
            raise
        gaps = residue_box_distance(text, box_centre, box_size)
        near = [rid for rid in anomalous if gaps.get(rid, 0.0) <= drop_radius]
        if near:
            raise RuntimeError(
                f"hydrogens disagree with Meeko's template within {drop_radius} A of the "
                f"box, at {', '.join(near)}") from exc
        polymer = Polymer.from_pqr_string(text, mk_prep=mk, residues_to_delete=anomalous,
                                          **options)
        dropped = anomalous
        LOGGER.warning("deleted %d residue(s) with PQR/template hydrogen discrepancies, all "
                       "more than %.1f A beyond the box: %s", len(dropped), drop_radius,
                       ", ".join(dropped))
    applied: List[str] = []
    for record in flexres:
        rid = meeko_residue_id(record)
        if rid not in polymer.monomers:
            LOGGER.warning("flexible residue %s is not a monomer Meeko built; skipped", rid)
            continue
        polymer.flexibilize_sidechain(rid, mk)
        applied.append(rid)
    rigid_text, flex_text = (str(part) for part in PDBQTWriterLegacy.write_string_from_polymer(polymer)[:2])
    rigid = workdir / "rigid.pdbqt"
    rigid.write_text(rigid_text)
    if not flex_text.strip():
        return rigid, None, applied, dropped
    flex = workdir / "flex.pdbqt"
    flex.write_text(flex_text)
    return rigid, flex, applied, dropped


def strip_flex_residues(pdbqt: str) -> str:
    """Drop ``BEGIN_RES``/``END_RES`` blocks from a pose file.

    Vina writes each moved side chain inside the pose model. Those atoms belong to the
    receptor, not the ligand, and would corrupt every RMSD if the pose reader saw them.
    """
    out, skipping = [], False
    for line in pdbqt.splitlines():
        token = line.strip()
        if token.startswith("BEGIN_RES"):
            skipping = True
            continue
        if token.startswith("END_RES"):
            skipping = False
            continue
        if not skipping:
            out.append(line)
    return "\n".join(out) + ("\n" if out else "")


def dock_arm(rigid: Path, flex: Optional[Path], ligand_pdbqt: str, centre: Sequence[float],
             size: Sequence[float], *, seed: int) -> object:
    """One Vina run, mirroring ``engine.dock`` but able to pass a flexible receptor."""
    import time

    from vina import Vina

    from cryptic_ip.docking import engine

    start = time.time()
    v = Vina(sf_name="vina", seed=seed, cpu=0, verbosity=0)
    if flex is None:
        v.set_receptor(rigid_pdbqt_filename=str(rigid))
    else:
        v.set_receptor(rigid_pdbqt_filename=str(rigid), flex_pdbqt_filename=str(flex))
    v.set_ligand_from_string(ligand_pdbqt)
    v.compute_vina_maps(center=[float(c) for c in centre], box_size=[float(s) for s in size])
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
    text = v.poses(n_poses=N_POSES, energy_range=ENERGY_RANGE)
    energies = np.asarray(v.energies(n_poses=N_POSES, energy_range=ENERGY_RANGE), dtype=float)
    return engine.DockResult("vina", seed, energies, strip_flex_residues(text),
                             [float(c) for c in centre], [float(s) for s in size],
                             time.time() - start)


# ------------------------------------------------------------------ outcomes
def success(rmsd: Sequence[float]) -> Optional[bool]:
    """Top-pose success at 2 Å: the first pose in Vina's own ranking."""
    values = [r for r in rmsd]
    if not values or values[0] is None or not np.isfinite(values[0]):
        return None
    return bool(values[0] <= SUCCESS_RMSD)


def best_of_list(rmsd: Sequence[float]) -> Optional[bool]:
    """J2's ceiling: does *any* retained pose reach 2 Å?"""
    values = [r for r in rmsd if r is not None and np.isfinite(r)]
    if not values:
        return None
    return bool(min(values) <= SUCCESS_RMSD)


def per_copy(records: Sequence[dict]) -> pd.DataFrame:
    """One row per copy and arm: success and ceiling, averaged over the seeds."""
    rows = []
    for record in records:
        if record.get("error") or not record.get("runs"):
            continue
        for arm, runs in record["runs"].items():
            tops = [success(r.get("rmsd") or []) for r in runs]
            ceilings = [best_of_list(r.get("rmsd") or []) for r in runs]
            tops = [t for t in tops if t is not None]
            ceilings = [c for c in ceilings if c is not None]
            if not tops:
                continue
            rows.append({"copy_key": record["copy_key"], "arm": arm,
                         "homology_group_strict": record.get("homology_group_strict"),
                         "burial_class": record.get("burial_class"),
                         "n_flex": int(record.get("n_flex") or 0),
                         "success": float(np.mean(tops)), "ceiling": float(np.mean(ceilings))})
    return pd.DataFrame(rows)


def paired(frame: pd.DataFrame, column: str, n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    """Paired flex - rigid difference in ``column``, resampling whole strict groups."""
    from cryptic_ip.docking.stats import MIN_GROUPS, paired_difference

    # Guard before pivoting, not after: pivot_table raises KeyError on an empty frame, and
    # the report job runs on always() precisely so that a run whose shards all died still
    # produces a "not evaluable" report rather than a second failure.
    needed = {"copy_key", "homology_group_strict", "arm", column}
    if frame.empty or not needed <= set(frame.columns):
        return {"decision": "not evaluable", "reason": "no copy was scored"}
    wide = frame.pivot_table(index=["copy_key", "homology_group_strict"], columns="arm",
                             values=column, aggfunc="first").dropna()
    if wide.empty or not {"flex", "rigid"} <= set(wide.columns):
        return {"decision": "not evaluable", "reason": "an arm has no scored copy"}
    groups = [g for _c, g in wide.index]
    est = paired_difference(wide["flex"].to_numpy(dtype=float), wide["rigid"].to_numpy(dtype=float),
                            groups, n_bootstrap=n_bootstrap, seed=SEED)
    out = {"estimate": est, "n_copies": int(len(wide)),
           "flex_mean": float(wide["flex"].mean()), "rigid_mean": float(wide["rigid"].mean())}
    if not est.get("evidence") or "per_group" not in est:
        out["decision"] = "not evaluable"
        out["reason"] = f"{est.get('groups')} strict groups, fewer than the {MIN_GROUPS} required"
        return out
    low, high = est["per_group"]["low"], est["per_group"]["high"]
    if low > 0:
        out["decision"] = "better"
    elif high < 0:
        out["decision"] = "worse"
    elif low > -NO_MATERIAL_CHANGE and high < NO_MATERIAL_CHANGE:
        out["decision"] = "no material change"
    else:
        out["decision"] = "inconclusive"
    return out


def error_classes(records: Sequence[dict]) -> Dict[str, object]:
    """One entry per exception type, with a count and one worked example.

    The per-copy traceback is the only thing that says *where* a copy died, and it reached
    nobody: it sat in the shard artifact while the report printed the message alone. Run
    36828409167 cost a second run to recover what this would have printed the first time.
    """
    out: Dict[str, object] = {}
    for record in records:
        error = record.get("error")
        if not error:
            continue
        name = str(error).split(":", 1)[0].strip() or "unknown"
        entry = out.setdefault(name, {"n": 0, "copies": []})
        entry["n"] = int(entry["n"]) + 1
        copies = entry["copies"]
        if len(copies) < 3:
            copies.append(str(record.get("copy_key")))
        if "traceback" not in entry and record.get("traceback"):
            entry["example_copy"] = str(record.get("copy_key"))
            entry["traceback"] = str(record["traceback"])[-1200:]
    return out


def build(records: Sequence[dict], n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    from cryptic_ip.benchmark import protocol

    frame = per_copy(records)
    result: Dict[str, object] = {
        "seed": SEED, "n_bootstrap": n_bootstrap, "exhaustiveness": EXHAUSTIVENESS,
        "flex_radius": FLEX_RADIUS, "max_flex": MAX_FLEX,
        "n_records": len(records), "n_copies_scored": int(frame["copy_key"].nunique()) if not frame.empty else 0,
        "errors": sorted({str(r["error"]) for r in records if r.get("error")}),
        "error_classes": error_classes(records),
    }
    if not frame.empty:
        result["flex_residues_per_copy"] = {
            "mean": float(frame[frame.arm == "flex"]["n_flex"].mean()),
            "max": int(frame[frame.arm == "flex"]["n_flex"].max()),
            "zero": int((frame[frame.arm == "flex"]["n_flex"] == 0).sum())}
    # Amendment 2: how much receptor was deleted to make Meeko accept it, never implicit.
    deleted = {str(r["copy_key"]): list(r["residues_dropped"]) for r in records
               if r.get("residues_dropped")}
    result["residues_dropped"] = {
        "copies": len(deleted), "residues": int(sum(len(v) for v in deleted.values())),
        "by_copy": dict(sorted(deleted.items())[:20]), "drop_radius": DROP_RADIUS}
    result["J1"] = paired(frame, "success", n_bootstrap=n_bootstrap)
    result["J2"] = paired(frame, "ceiling", n_bootstrap=n_bootstrap)
    # Amendment 1: the rigid arm under this receptor pipeline, against study F's 0.110
    # under the project-prepared one, so the pipeline change is auditable.
    if not frame.empty:
        rigid = frame[frame.arm == "rigid"]
        result["rigid_arm_check"] = {"meeko_receptor_success": float(rigid["success"].mean()),
                                     "study_f_success": 0.110, "n_copies": int(len(rigid))}
    result["J3"] = {}
    if not frame.empty:
        for burial, block in frame.groupby("burial_class"):
            result["J3"][str(burial)] = paired(block, "success", n_bootstrap=n_bootstrap)
    pvalues = {}
    for name in ("J1", "J2"):
        est = result[name].get("estimate", {})
        if "per_group" in est:
            pvalues[name] = est["per_group"]["p_value"]
    if pvalues:
        result["holm"] = protocol.holm(pvalues)
    return result


def markdown(r: Dict[str, object]) -> str:
    lines = ["# Study J: does a flexible receptor recover study A's failures?", "",
             f"Pre-registered in `docs/FLEXIBLE_PLAN.md`, amended in "
             f"`docs/FLEXIBLE_PLAN_AMENDMENT_1.md` and "
             f"`docs/FLEXIBLE_PLAN_AMENDMENT_2.md`. Seed {r['seed']}, {r['n_bootstrap']} "
             f"resamples of strict homology groups, exhaustiveness {r['exhaustiveness']}.",
             f"Flexible side chains: within {r['flex_radius']} A of the ligand, at most "
             f"{r['max_flex']}.", ""]
    flex = r.get("flex_residues_per_copy")
    if flex:
        lines.append(f"Movable side chains per copy: mean {flex['mean']:.2f}, max {flex['max']}, "
                     f"and {flex['zero']} copies had none.")
    for name, title in (("J1", "J1 - top-pose success at 2 A (primary)"),
                        ("J2", "J2 - best-of-list ceiling (secondary)")):
        block = r[name]
        lines += ["", f"## {title}", "", f"Decision: **{block.get('decision')}**."]
        per_group = block.get("estimate", {}).get("per_group")
        if per_group:
            lines.append(f"flex {block['flex_mean']:.3f} against rigid {block['rigid_mean']:.3f}; "
                         f"paired difference {per_group['point']:+.3f} "
                         f"[{per_group['low']:+.3f}, {per_group['high']:+.3f}] over "
                         f"{block['n_copies']} copies in {block['estimate'].get('groups')} groups.")
        if block.get("reason"):
            lines.append(f"Reason: {block['reason']}.")
    check = r.get("rigid_arm_check")
    if check:
        lines += ["", "## The receptor pipeline, audited", "",
                  f"The rigid arm scores {check['meeko_receptor_success']:.3f} under this study's "
                  f"receptor, against study F's {check['study_f_success']:.3f} under the "
                  f"project-prepared one, over {check['n_copies']} copies. A material gap is a "
                  "finding about the pipeline, not about side-chain freedom, and is not folded "
                  "into J1."]
    if r.get("J3"):
        lines += ["", "## J3 - by burial class (exploratory, not Holm-corrected)", "",
                  "| burial class | decision | difference | copies | groups |", "| --- | --- | --- | --- | --- |"]
        for burial, block in sorted(r["J3"].items()):
            per_group = block.get("estimate", {}).get("per_group")
            shown = (f"{per_group['point']:+.3f} [{per_group['low']:+.3f}, {per_group['high']:+.3f}]"
                     if per_group else "-")
            lines.append(f"| {burial} | {block.get('decision')} | {shown} | "
                         f"{block.get('n_copies', '-')} | {block.get('estimate', {}).get('groups', '-')} |")
    dropped = r.get("residues_dropped") or {}
    if dropped.get("copies"):
        lines += ["", "## Receptor residues deleted to satisfy Meeko", "",
                  f"{dropped['residues']} residue(s) across {dropped['copies']} copies, each "
                  f"more than {dropped['drop_radius']} A beyond the docking box, where Vina's "
                  "own interaction cutoff puts them out of reach of any pose. Both arms of a "
                  "copy see the same deletions or the copy fails. Per "
                  "docs/FLEXIBLE_PLAN_AMENDMENT_2.md."]
    if r.get("error_classes"):
        lines += ["", "## Copies not docked, by cause", "",
                  "| cause | copies | examples |", "| --- | --- | --- |"]
        ordered = sorted(r["error_classes"].items(), key=lambda kv: (-int(kv[1]["n"]), kv[0]))
        for name, entry in ordered:
            lines.append(f"| {name} | {entry['n']} | {', '.join(entry['copies'])} |")
        for name, entry in ordered:
            if entry.get("traceback"):
                lines += ["", f"### {name}, at {entry.get('example_copy')}", "", "```",
                          *str(entry["traceback"]).splitlines(), "```"]
    if r.get("errors"):
        lines += ["", "## Copies not docked", ""] + [f"- {e}" for e in r["errors"]]
    return "\n".join(lines) + "\n"


def dock_copy(ctx, flexres: Sequence[Dict[str, object]]
              ) -> Tuple[Dict[str, List[dict]], List[str], List[str]]:
    """Both arms for one copy: the same receptor, differing only in movable side chains."""
    from cryptic_ip.docking import engine
    from cryptic_ip.docking.ligand import start_pose, to_pdbqt

    # Runs the project's PDB2PQR path and caches receptor_h.pdb / receptor.pqr, which is
    # the protonation both arms are then built from.
    ctx.receptor()
    protonated = (ctx.work / "receptor_h.pdb", ctx.work / "receptor.pqr")
    if not all(path.exists() for path in protonated):
        raise RuntimeError("the protonated receptor pair was not cached")

    size = [ctx.side] * 3
    box = {"box_centre": [float(c) for c in ctx.centre], "box_size": size}

    # Prepared twice into separate directories: flexibilising mutates the polymer, and the
    # rigid arm must see a receptor with no flexible residues at all.
    rigid_only, none_flex, _, dropped_rigid = prepare_pair(protonated, ctx.work / "arm_rigid",
                                                           (), **box)
    if none_flex is not None:
        raise RuntimeError("the rigid arm was given a flexible receptor")
    rigid_part, flex_part, applied, dropped_flex = prepare_pair(protonated,
                                                                ctx.work / "arm_flex",
                                                                flexres, **box)
    # The arms must differ in side-chain freedom and nothing else, so a residue deleted
    # from one receptor and not the other is a defect, not a difference to report.
    if sorted(dropped_rigid) != sorted(dropped_flex):
        raise RuntimeError(f"the arms dropped different residues: {sorted(dropped_rigid)} "
                           f"against {sorted(dropped_flex)}")

    runs: Dict[str, List[dict]] = {arm: [] for arm in ARMS}
    for arm, (rigid, flex) in (("rigid", (rigid_only, None)), ("flex", (rigid_part, flex_part))):
        for seed in engine.SEEDS:
            pose = start_pose(ctx.ligands["primary"], ctx.centre, seed, crystal=ctx.crystal)
            result = dock_arm(rigid, flex, to_pdbqt(pose.mol), ctx.centre, size, seed=seed)
            table = engine.pose_table(result, ctx.crystal, site_centroid=ctx.centre)
            runs[arm].append({"seed": seed, "vina": table.scores, "rmsd": table.rmsd,
                              "seconds": result.seconds})
    return runs, applied, dropped_flex


def cmd_dock_one(args: argparse.Namespace) -> int:
    import time
    import traceback

    from redocking import Context, _json_default

    row = json.loads(args.row_json.read_text())
    record: Dict[str, object] = {"copy_key": row["copy_key"], "pdb_id": row["pdb_id"],
                                 "homology_group_strict": row.get("homology_group_strict"),
                                 "burial_class": row.get("burial_class")}
    start = time.time()
    try:
        ctx = Context(row, args.structures_dir, args.ccd_dir, args.work_dir)
        flexres = flexible_residues(ctx.arrays, ctx.xyz)
        record["flex_residues"] = flexres
        runs, applied, dropped = dock_copy(ctx, flexres)
        record["n_flex"] = len(applied)
        record["flex_applied"] = applied
        record["residues_dropped"] = dropped
        record["runs"] = runs
    except Exception as exc:  # noqa: BLE001 - recorded as this copy's failure reason
        record["error"] = f"{type(exc).__name__}: {exc}"[:500]
        record["traceback"] = traceback.format_exc()[-1500:]
    record["seconds"] = time.time() - start
    args.out_json.write_text(json.dumps(record, default=_json_default))
    return 0


def cmd_dock(args: argparse.Namespace) -> int:
    """A shard of copies, each in its own process so one crash cannot take the shard.

    Three things protect the shard's work, all of them learned from run 36828409167, where
    thirteen shards each threw away every copy they had already docked:

    * a copy that overruns its timeout is this copy's failure, not the shard's;
    * the output file is rewritten after every copy, so a shard the runner kills still
      uploads what it had;
    * the shard stops on its own wall-clock budget rather than being killed by
      ``timeout-minutes``, and the copies it never reached say so.
    """
    import subprocess
    import tempfile
    import time

    from redocking import _json_default, shard_entries

    census = pd.read_csv(args.census)
    rows = census[census["selected"].astype(str).str.lower() == "true"]
    keys = shard_entries(sorted(rows["copy_key"].astype(str)), args.shard, args.shards)
    LOGGER.info("shard %d/%d: %d copies", args.shard, args.shards, len(keys))
    deadline = time.time() + args.shard_budget
    out: List[dict] = []
    target = Path(args.out)

    def flush() -> None:
        target.write_text(json.dumps(out, default=_json_default))

    flush()
    with tempfile.TemporaryDirectory() as tmp:
        for index, key in enumerate(keys):
            left = deadline - time.time()
            if left < args.min_copy_seconds:
                for unreached in keys[index:]:
                    out.append({"copy_key": unreached,
                                "error": "not reached: the shard's budget ran out"})
                LOGGER.info("budget spent with %d copies unreached", len(keys) - index)
                break
            row = rows[rows["copy_key"].astype(str) == key].iloc[0].to_dict()
            row_json = Path(tmp) / "row.json"
            out_json = Path(tmp) / "out.json"
            row_json.write_text(json.dumps(row, default=_json_default))
            timeout = min(args.per_copy_timeout, left)
            try:
                proc = subprocess.run(
                    [sys.executable, str(Path(__file__).resolve()), "dock-one",
                     "--row-json", str(row_json), "--out-json", str(out_json),
                     "--structures-dir", str(args.structures_dir), "--ccd-dir", str(args.ccd_dir),
                     "--work-dir", str(args.work_dir)],
                    capture_output=True, text=True, timeout=timeout)
            except subprocess.TimeoutExpired:
                # One copy overrunning is this copy's failure. Before this was caught, it
                # raised out of the loop and discarded every copy the shard had docked.
                out.append({"copy_key": key, "error": f"timed out after {timeout:.0f} s"})
            else:
                if out_json.exists():
                    out.append(json.loads(out_json.read_text()))
                    out_json.unlink()
                else:
                    out.append({"copy_key": key, "error": f"child exited {proc.returncode}",
                                "stderr": (proc.stderr or "")[-800:]})
            flush()
            # The reason, not just the fact: an error visible only inside the artifact
            # cannot be diagnosed from the logs.
            LOGGER.info("%s -> %s", key, out[-1].get("error") or "ok")
    flush()
    return 0


# ------------------------------------------------------------------ commands
def cmd_residues(args: argparse.Namespace) -> int:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    import redocking

    census = pd.read_csv(args.census)
    rows = census[census["selected"].astype(str).str.lower() == "true"]
    out = []
    for _, row in rows.iterrows():
        path = redocking.structure_path(Path(args.structures), str(row["pdb_id"]))
        if path is None:
            out.append({"copy_key": row["copy_key"], "error": "structure missing"})
            continue
        try:
            arrays = load_structure_arrays(path)
            # find_copy *raises* when the copy is absent, and the icode cell is NaN for
            # every copy without an insertion code - `NaN or ""` is NaN, so formatting it
            # produces the string "nan" and matches nothing. Both cost a whole shard once.
            _key, _comp_id, atoms = redocking.find_copy(
                arrays, _text(row["chain"]), int(row["resseq"]), _text(row.get("icode")))
            heavy = atoms[arrays.elements[atoms] != "H"]
            chosen = flexible_residues(arrays, arrays.coords[heavy])
        except Exception as exc:  # noqa: BLE001 - one copy's failure is not the shard's
            out.append({"copy_key": row["copy_key"], "error": f"{type(exc).__name__}: {exc}"[:300]})
            continue
        out.append({"copy_key": row["copy_key"], "n_flex": len(chosen), "residues": chosen})
    Path(args.out).write_text(json.dumps(out, indent=2, default=str))
    counts = [r.get("n_flex", 0) for r in out]
    LOGGER.info("%d copies, mean %.2f flexible side chains, max %d",
                len(out), float(np.mean(counts)) if counts else 0.0, max(counts) if counts else 0)
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    import template_fit as tf

    records = tf.read_json(Path(args.records))
    result = build(records, n_bootstrap=args.n_bootstrap)
    tf.write_json(Path(args.out), result)
    Path(args.markdown).write_text(markdown(result))
    LOGGER.info("J1 %s, J2 %s", result["J1"].get("decision"), result["J2"].get("decision"))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("residues", help="which side chains each copy makes flexible")
    p.add_argument("--census", required=True)
    p.add_argument("--structures", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_residues)

    p = sub.add_parser("dock", help="both arms for a shard of copies")
    p.add_argument("--census", required=True)
    p.add_argument("--structures-dir", required=True, type=Path)
    p.add_argument("--ccd-dir", required=True, type=Path)
    p.add_argument("--work-dir", required=True, type=Path)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=1)
    p.add_argument("--per-copy-timeout", type=int, default=5400)
    # Under the job's own timeout-minutes, so the shard stops itself and writes its output
    # instead of being killed with the file unwritten.
    p.add_argument("--shard-budget", type=int, default=19800)
    p.add_argument("--min-copy-seconds", type=int, default=600)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_dock)

    p = sub.add_parser("dock-one", help="one copy, run as a child process")
    p.add_argument("--row-json", required=True, type=Path)
    p.add_argument("--out-json", required=True, type=Path)
    p.add_argument("--structures-dir", required=True, type=Path)
    p.add_argument("--ccd-dir", required=True, type=Path)
    p.add_argument("--work-dir", required=True, type=Path)
    p.set_defaults(func=cmd_dock_one)

    p = sub.add_parser("report", help="J1, J2 and J3")
    p.add_argument("--records", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--markdown", required=True)
    p.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    p.set_defaults(func=cmd_report)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
