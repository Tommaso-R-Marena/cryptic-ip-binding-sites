#!/usr/bin/env python3
"""The redocking benchmark (docs/REDOCKING_PLAN.md).

``census``   every IP copy of a shard's entries: burial class, flags, CCD
             template, configuration check, exclusions and the selection
``dock``     a shard of selected copies through one group of arms:
             ``primary``   Vina seeds 1-3 and the crystal-pose control
             ``secondary`` Vinardo, AD4, the deprotonated ligand, kept metals
             ``pockets``   fpocket top-3 site finding and the decoy pocket
             ``alphafold`` docking into the superposed AlphaFold model
``report``   see scripts/redocking_report.py

Every copy is docked in a child process with a time limit; a shard stops
starting copies when its budget is spent and records the rest as not reached.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import multiprocessing as mp
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.database.ip_ligands import ip_series_label  # noqa: E402

LOGGER = logging.getLogger("redocking")

CCD_URL = "https://files.rcsb.org/ligands/download/{comp_id}.cif"
SYMMETRY_CONTACT = 4.5
SITE_RESIDUE_DISTANCE = 6.0
AF_MIN_IDENTITY = 0.90
AF_MIN_COVERAGE = 0.50
AF_MAX_LENGTH = 2700
TRUE_SITE_DISTANCE = 4.0
ARMS = ("primary", "secondary", "pockets", "alphafold")
#: The plan's exhaustiveness; only the synthetic tests lower it (``--exhaustiveness``).
EXHAUSTIVENESS = [32]


# ---------------------------------------------------------------- helpers
def shard_entries(entries: Sequence[str], index: int, count: int) -> List[str]:
    """Entries of one shard: sorted, round robin (an entry's copies stay together)."""
    ordered = sorted(set(entries))
    return [e for i, e in enumerate(ordered) if i % count == index]


def structure_path(directory: Path, pdb_id: str) -> Optional[Path]:
    for suffix in (".pdb", ".cif", ".pdb.gz", ".cif.gz", ".ent"):
        for name in (pdb_id.upper(), pdb_id.lower()):
            p = directory / f"{name}{suffix}"
            if p.exists():
                return p
    return None


def copy_key(pdb_id: str, chain: str, resseq: int, icode: str) -> str:
    return f"{pdb_id.upper()}:{chain}:{int(resseq)}{icode or ''}"


def find_copy(arrays, chain: str, resseq: int, icode: str):
    from cryptic_ip.validation.burial_metrics import find_ligand_instances

    for key, comp_id, atoms in find_ligand_instances(arrays):
        if str(key[1]) == str(chain) and int(key[2]) == int(resseq) and str(key[3] or "") == str(icode or ""):
            return key, comp_id, atoms
    raise LookupError(f"copy {chain}:{resseq}{icode} not found")


def ligand_extent(mol) -> float:
    """Largest heavy-atom distance of a conformer (Å)."""
    from rdkit import Chem

    xyz = Chem.RemoveHs(mol).GetConformer().GetPositions()
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    return float(d.max())


def box_side(ligand, seed: int = 1) -> float:
    """Plan section 4: max(22, d_max + 16) Å from the seed-1 ETKDG conformer."""
    from cryptic_ip.docking.ligand import start_pose

    pose = start_pose(ligand, [0.0, 0.0, 0.0], seed, crystal=None)
    return float(max(22.0, ligand_extent(pose.mol) + 16.0))


# ------------------------------------------------------------------ CCD
def fetch_ccd(comp_ids: Sequence[str], out_dir: Path, fetch=None) -> Dict[str, Path]:
    """CCD mmCIF files through the verified, retrying fetcher (cached on disk)."""
    from cryptic_ip.database.async_fetch import FetchJob, validate_nonempty

    out_dir.mkdir(parents=True, exist_ok=True)
    if fetch is None:
        from cryptic_ip.database.async_fetch import fetch_all

        def fetch(jobs):
            return fetch_all(jobs, concurrency=4, per_host=4)

    def validate(payload: bytes) -> None:
        validate_nonempty(payload)
        if b"_chem_comp" not in payload:
            from cryptic_ip.database.async_fetch import ValidationError

            raise ValidationError("not a CCD mmCIF")

    wanted = sorted({c.upper() for c in comp_ids})
    jobs = [FetchJob(key=c, url=CCD_URL.format(comp_id=c), dest=out_dir / f"{c}.cif", validator=validate)
            for c in wanted]
    results = fetch(jobs)
    out = {}
    for r in results:
        if r.ok:
            out[r.key] = out_dir / f"{r.key}.cif"
        else:
            LOGGER.warning("CCD %s: %s", r.key, r.error)
    return out


def ccd_template(path: Optional[Path]):
    """``(parent_mol, smiles, source)`` or raises LigandError."""
    from cryptic_ip.docking.ligand import LigandError, ccd_smiles, parent_molecule, parse_ccd_cif

    if path is None or not path.exists():
        raise LigandError("CCD entry could not be fetched")
    info = parse_ccd_cif(path.read_text())
    smiles, source = ccd_smiles(info["descriptors"])
    return parent_molecule(smiles), smiles, source


# --------------------------------------------------------------- census
def gemmi_flags(path: Path, chain: str, resseq: int, icode: str, ligand_xyz: np.ndarray) -> Dict[str, object]:
    """Altlocs of the copy, and polymer atoms of symmetry mates within 4.5 Å (X-ray only)."""
    import gemmi

    out: Dict[str, object] = {"altloc": False, "symmetry_contact": None}
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    model = st[0]
    for ch in model:
        if ch.name != chain:
            continue
        for res in ch:
            if res.seqid.num == int(resseq) and (res.seqid.icode.strip() or "") == (icode or ""):
                if any(a.altloc != "\0" for a in res):
                    out["altloc"] = True
    if st.cell.is_crystal() and st.cell.a > 1.5:
        ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
        contact = False
        for xyz in ligand_xyz:
            for mark in ns.find_atoms(gemmi.Position(*map(float, xyz)), "\0", radius=SYMMETRY_CONTACT):
                if mark.image_idx == 0:
                    continue
                cra = mark.to_cra(model)
                if cra.residue.entity_type == gemmi.EntityType.Polymer:
                    contact = True
                    break
            if contact:
                break
        out["symmetry_contact"] = contact
    return out


def census_entry(pdb_id: str, path: Path, meta: Dict[str, object], ccd_dir: Path) -> List[Dict[str, object]]:
    from cryptic_ip.docking.ligand import LigandError, configuration_matches, crystal_molecule
    from cryptic_ip.docking.receptor import nearby_metals
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.validation.burial_metrics import compute_ligand_burial, find_ligand_instances

    arrays = load_structure_arrays(path)
    found = find_ligand_instances(arrays)
    burial = {b.instance_id: b for b in compute_ligand_burial(path, n_points=256)}
    rows = []
    templates: Dict[str, object] = {}
    for key, comp_id, atoms in found:
        suffix = f"{key[2]}{key[3]}".strip()
        instance_id = f"{comp_id}_{key[1]}_{suffix}"
        b = burial.get(instance_id)
        heavy = atoms[arrays.elements[atoms] != "H"]
        xyz = arrays.coords[heavy]
        row: Dict[str, object] = {
            "pdb_id": pdb_id, "copy_key": copy_key(pdb_id, key[1], key[2], key[3]), "instance_id": instance_id,
            "comp_id": comp_id, "chain": key[1], "resseq": int(key[2]), "icode": key[3] or "",
            "n_heavy": int(len(heavy)), "n_phosphorus": int(np.sum(arrays.elements[heavy] == "P")),
            "species": ip_series_label(int(np.sum(arrays.elements[heavy] == "P"))),
            "centroid_x": float(xyz[:, 0].mean()), "centroid_y": float(xyz[:, 1].mean()),
            "centroid_z": float(xyz[:, 2].mean()),
            **{k: meta.get(k) for k in ("resolution", "experimental_method", "release_date",
                                        "homology_group", "homology_group_strict", "uniprot_ids")},
        }
        if b is not None:
            row.update({"burial_class": b.burial_class, "relative_sasa": b.relative_sasa,
                        "n_protein_contacts": b.n_protein_contacts, "n_contact_chains": b.n_contact_chains,
                        "interface": b.n_contact_chains >= 2, "mean_occupancy": b.mean_occupancy,
                        "partial_occupancy": b.mean_occupancy < 1.0})
        else:
            row.update({"burial_class": "unknown", "interface": None})
        metals = nearby_metals(arrays, xyz)
        row["metal"] = bool(metals)
        row["metal_ions"] = ";".join(sorted({m.element for m in metals}))
        try:
            row.update(gemmi_flags(path, key[1], key[2], key[3] or "", xyz))
        except Exception as exc:  # noqa: BLE001 - a flag must not lose the copy
            row["flag_error"] = f"{type(exc).__name__}: {exc}"[:200]
        status = "eligible"
        if row["burial_class"] == "crystal_artifact":
            status = "excluded: crystal artefact"
        elif row["burial_class"] == "unknown":
            status = "excluded: burial could not be measured"
        else:
            try:
                if comp_id not in templates:
                    templates[comp_id] = ccd_template(ccd_dir / f"{comp_id}.cif")
                parent, smiles, source = templates[comp_id]
                row["ccd_smiles"], row["smiles_source"] = smiles, source
                from rdkit import Chem

                row["template_heavy"] = Chem.RemoveHs(parent).GetNumAtoms()
                crystal, complete = crystal_molecule(arrays.elements[heavy], xyz, arrays.atom_names[heavy], parent)
                row["complete"] = bool(complete)
                row["configuration_match"] = bool(configuration_matches(crystal, parent)) if complete else None
                if complete and not row["configuration_match"]:
                    status = "excluded: configuration differs from the CCD"
            except LigandError as exc:
                status = f"excluded: {exc}"
            except Exception as exc:  # noqa: BLE001
                status = f"excluded: ligand check failed ({type(exc).__name__}: {exc})"[:200]
        row["status"] = status
        rows.append(row)
    return rows


def select_copies(rows: pd.DataFrame) -> pd.Series:
    """Plan section 1: at most one copy per (entry, burial class), complete copies first, then by key."""
    selected = pd.Series(False, index=rows.index)
    eligible = rows[rows["status"] == "eligible"].copy()
    if eligible.empty:
        return selected
    eligible["_incomplete"] = ~eligible["complete"].fillna(False).astype(bool)
    eligible = eligible.sort_values(["pdb_id", "burial_class", "_incomplete", "chain", "resseq", "icode"])
    first = eligible.groupby(["pdb_id", "burial_class"], sort=False).head(1)
    selected.loc[first.index] = True
    return selected


def cmd_census(args: argparse.Namespace) -> int:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.validation.burial_metrics import find_ligand_instances

    table_ids = set(pd.read_csv(args.table, usecols=["structure_id"])["structure_id"].str.upper())
    entries = pd.read_csv(args.entries, dtype=str)
    entries["pdb_id"] = entries["pdb_id"].str.upper()
    entries = entries[entries["pdb_id"].isin(table_ids)].drop_duplicates("pdb_id").set_index("pdb_id")
    mine = shard_entries(list(entries.index), args.shard_index, args.shard_count)
    LOGGER.info("%d entries in the table, %d in this shard", len(entries), len(mine))
    comp_ids = set()
    paths: Dict[str, Optional[Path]] = {}
    for pdb_id in mine:
        paths[pdb_id] = structure_path(args.structures_dir, pdb_id)
        if paths[pdb_id] is not None:
            try:
                comp_ids |= {c for _, c, _ in find_ligand_instances(load_structure_arrays(paths[pdb_id]))}
            except Exception:  # noqa: BLE001 - reported below
                pass
    fetch_ccd(sorted(comp_ids), args.ccd_dir)
    rows: List[Dict[str, object]] = []
    failures = []
    for pdb_id in mine:
        if paths[pdb_id] is None:
            failures.append({"pdb_id": pdb_id, "error": "structure not fetched"})
            continue
        try:
            got = census_entry(pdb_id, paths[pdb_id], entries.loc[pdb_id].to_dict(), args.ccd_dir)
            if not got:
                failures.append({"pdb_id": pdb_id, "error": "no IP copy found"})
            rows += got
        except Exception as exc:  # noqa: BLE001
            failures.append({"pdb_id": pdb_id, "error": f"{type(exc).__name__}: {exc}"[:300]})
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["selected"] = select_copies(frame)
        frame["primary_set"] = frame["selected"] & frame["complete"].fillna(False).astype(bool)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    pd.DataFrame(failures, columns=["pdb_id", "error"]).to_csv(args.output.with_name(
        args.output.stem.replace(".csv", "") + "_failures.csv"), index=False)
    print(json.dumps({"entries": len(mine), "copies": len(frame), "failures": len(failures),
                      "selected": int(frame["selected"].sum()) if not frame.empty else 0}))
    return 0


# ------------------------------------------------------------------ dock
class Context:
    """Everything one copy's arms share."""

    def __init__(self, row: Dict[str, object], structures_dir: Path, ccd_dir: Path, work: Path) -> None:
        from rdkit import Chem

        from cryptic_ip.analysis.structure_arrays import load_structure_arrays
        from cryptic_ip.docking.ligand import crystal_molecule, protonate

        self.row = row
        self.pdb_id = str(row["pdb_id"])
        self.path = structure_path(structures_dir, self.pdb_id)
        if self.path is None:
            raise FileNotFoundError(f"{self.pdb_id} not fetched")
        self.arrays = load_structure_arrays(self.path)
        self.key, self.comp_id, atoms = find_copy(self.arrays, str(row["chain"]), int(row["resseq"]),
                                                  str(row.get("icode") or "") if not _isnan(row.get("icode")) else "")
        self.heavy = atoms[self.arrays.elements[atoms] != "H"]
        self.xyz = self.arrays.coords[self.heavy]
        self.centre = self.xyz.mean(axis=0)
        self.parent, _, _ = ccd_template(ccd_dir / f"{self.comp_id}.cif")
        self.crystal, self.complete = crystal_molecule(self.arrays.elements[self.heavy], self.xyz,
                                                       self.arrays.atom_names[self.heavy], self.parent)
        self.crystal = Chem.RemoveHs(self.crystal, sanitize=False)
        self.ligands = {state: protonate(self.parent, state) for state in ("primary", "deprotonated")}
        self.side = box_side(self.ligands["primary"])
        self.work = work / self.pdb_id
        self.work.mkdir(parents=True, exist_ok=True)

    def receptor(self, metals: Sequence = (), name: str = "receptor"):
        from cryptic_ip.docking.engine import Receptor
        from cryptic_ip.docking.receptor import prepare_receptor

        cached = (self.work / "receptor_h.pdb", self.work / "receptor.pqr")
        protonated = cached if all(p.exists() for p in cached) else None
        prepared = prepare_receptor(self.arrays, self.work, keep_metals=metals, protonated=protonated, name=name)
        return Receptor(prepared.pdbqt), prepared


def _isnan(value) -> bool:
    return isinstance(value, float) and math.isnan(value)


def run_arm(ctx: Context, receptor, *, scoring: str = "vina", seed: int = 1, state: str = "primary",
            centre: Optional[Sequence[float]] = None, crystal=None, arm: str = "") -> Dict[str, object]:
    from cryptic_ip.docking import engine
    from cryptic_ip.docking.ligand import start_pose, to_pdbqt

    centre = ctx.centre if centre is None else np.asarray(centre, dtype=float)
    crystal = ctx.crystal if crystal is None else crystal
    pose = start_pose(ctx.ligands[state], centre, seed, crystal=ctx.crystal)
    size = [ctx.side] * 3
    result = engine.dock(receptor, to_pdbqt(pose.mol), centre, size, scoring=scoring, seed=seed,
                         exhaustiveness=EXHAUSTIVENESS[0])
    table = engine.pose_table(result, crystal, site_centroid=ctx.centre)
    ok = [r for r in table.rmsd if np.isfinite(r)]
    return {
        "arm": arm or f"{scoring}_{state}_s{seed}", "scoring": scoring, "seed": seed, "state": state,
        "box_centre": [float(c) for c in centre], "box_side": ctx.side, "start_rmsd": pose.start_rmsd,
        "start_seed": pose.seed, "n_poses": len(table.scores),
        "top_score": table.scores[0] if table.scores else None,
        "top_rmsd": table.rmsd[0] if table.rmsd else None,
        "top_rmsd_p": table.rmsd_p[0] if table.rmsd_p else None,
        "best_rmsd": float(min(ok)) if ok else None,
        "top_centroid_distance": table.centroid_distance[0] if table.centroid_distance else None,
        "top_pose_centroid": table.centroids[0] if table.centroids else None,
        "spearman": engine.spearman(table.scores, table.rmsd),
        "scores": table.scores, "rmsd": table.rmsd, "rmsd_p": table.rmsd_p,
        "centroid_distance": table.centroid_distance, "seconds": result.seconds,
    }


def arms_primary(ctx: Context) -> List[Dict[str, object]]:
    from cryptic_ip.docking import engine
    from cryptic_ip.docking.ligand import pose_on_crystal, poses_from_pdbqt, to_pdbqt
    from cryptic_ip.docking.rmsd import symmetric_rmsd

    receptor, prepared = ctx.receptor()
    out = []
    for seed in engine.SEEDS:
        out.append(run_arm(ctx, receptor, seed=seed, arm=f"vina_s{seed}"))
    out[0]["receptor"] = {"types": prepared.types, "his": prepared.his_states,
                          "removed_modified": prepared.strip.removed_modified_residues}
    if ctx.complete:
        placed = pose_on_crystal(ctx.ligands["primary"], ctx.crystal)
        before, after, pose_text = engine.score_and_minimise(receptor, to_pdbqt(placed), ctx.centre,
                                                             [ctx.side] * 3)
        minimised = poses_from_pdbqt(pose_text)[0]
        out.append({"arm": "crystal_control", "crystal_score": before, "crystal_minimised_score": after,
                    "minimised_rmsd": symmetric_rmsd(ctx.crystal, minimised)})
    return out


def arms_secondary(ctx: Context) -> List[Dict[str, object]]:
    from cryptic_ip.docking import engine
    from cryptic_ip.docking.ligand import start_pose, to_pdbqt
    from cryptic_ip.docking.receptor import nearby_metals

    receptor, _ = ctx.receptor()
    out = [run_arm(ctx, receptor, scoring="vinardo", arm="vinardo_s1"),
           run_arm(ctx, receptor, state="deprotonated", arm="vina_deprotonated_s1")]
    if engine.ad4_available():
        try:
            lig_text = to_pdbqt(start_pose(ctx.ligands["primary"], ctx.centre, 1, crystal=ctx.crystal).mol)
            rec_text = Path(receptor.pdbqt).read_text()
            gpf = ctx.work / "ad4.gpf"
            prefix = engine.write_gpf(gpf, Path(receptor.pdbqt), engine.pdbqt_types(rec_text),
                                      engine.pdbqt_types(lig_text), ctx.centre, [ctx.side] * 3)
            engine.run_autogrid(gpf)
            ad4 = engine.Receptor(receptor.pdbqt, ad4_maps=prefix)
            out.append(run_arm(ctx, ad4, scoring="ad4", arm="ad4_s1"))
        except Exception as exc:  # noqa: BLE001 - pre-declared: AD4 not run, with the reason
            out.append({"arm": "ad4_s1", "error": f"AD4 not run: {type(exc).__name__}: {exc}"[:300]})
    else:
        out.append({"arm": "ad4_s1", "error": "AD4 not run: autogrid4 not installed"})
    metals = nearby_metals(ctx.arrays, ctx.xyz)
    if metals:
        with_metals, _ = ctx.receptor(metals=metals, name="receptor_metals")
        for seed in engine.SEEDS:
            out.append(run_arm(ctx, with_metals, seed=seed, arm=f"vina_metals_s{seed}"))
    return out


def fpocket_pockets(ctx: Context) -> pd.DataFrame:
    """fpocket on the ligand-free structure (the benchmark's code path), labelled by the benchmark's rule."""
    import shutil
    import tempfile

    from cryptic_ip.analysis.analyzer import ProteinAnalyzer
    from cryptic_ip.analysis.labeling import LigandSite, assign_pocket_labels
    from cryptic_ip.analysis.structure_arrays import write_apo_structure
    from cryptic_ip.validation.burial_metrics import find_ligand_instances

    cache = ctx.work / "pockets.csv"
    if cache.exists():
        return pd.read_csv(cache)
    tmp = Path(tempfile.mkdtemp(prefix=f"fp_{ctx.pdb_id}_"))
    try:
        apo = write_apo_structure(ctx.path, tmp / f"{ctx.pdb_id}.pdb")
        analyzer = ProteinAnalyzer(str(apo), skip_electrostatics=True)
        analyzer.fpocket_timeout_s = 2700
        pockets = analyzer.detect_pockets(min_alpha_sphere=3)
        geometry = []
        for pid in pockets["pocket_id"]:
            row = pockets[pockets["pocket_id"] == pid].iloc[0]
            centre = np.array([row["center_x"], row["center_y"], row["center_z"]], dtype=float)
            spheres = analyzer._pocket_alpha_spheres(int(pid))
            geometry.append((int(pid), centre, spheres if spheres is not None and len(spheres) else centre[None]))
        sites = [LigandSite(f"{c}_{k[1]}_{k[2]}", c, ctx.arrays.coords[a[ctx.arrays.elements[a] != "H"]])
                 for k, c, a in find_ligand_instances(ctx.arrays)]
        labels = {a.pocket_id: a.label.value for a in assign_pocket_labels(geometry, sites)}
        out = pd.DataFrame({"pocket_id": [g[0] for g in geometry],
                            "center_x": [g[1][0] for g in geometry], "center_y": [g[1][1] for g in geometry],
                            "center_z": [g[1][2] for g in geometry],
                            "label": [labels.get(g[0], -1) for g in geometry]}).sort_values("pocket_id")
        out.to_csv(cache, index=False)
        return out
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def arms_pockets(ctx: Context, site_centroids: Sequence[np.ndarray], tertiary: bool) -> List[Dict[str, object]]:
    receptor, _ = ctx.receptor()
    pockets = fpocket_pockets(ctx)
    out = []
    decoy = pockets[pockets["label"] == 0].head(1)
    if len(decoy):
        c = decoy.iloc[0][["center_x", "center_y", "center_z"]].to_numpy(dtype=float)
        r = run_arm(ctx, receptor, centre=c, arm="decoy_s1")
        r["pocket_id"] = int(decoy.iloc[0]["pocket_id"])
        out.append(r)
    else:
        out.append({"arm": "decoy_s1", "error": "no negative fpocket pocket"})
    if tertiary:
        best = None
        top = pockets.head(3)
        for _, p in top.iterrows():
            c = p[["center_x", "center_y", "center_z"]].to_numpy(dtype=float)
            r = run_arm(ctx, receptor, centre=c, arm=f"pocket{int(p['pocket_id'])}_s1")
            r["pocket_id"] = int(p["pocket_id"])
            r["pocket_label"] = int(p["label"])
            out.append(r)
            if r["top_score"] is not None and (best is None or r["top_score"] < best["top_score"]):
                best = r
        summary: Dict[str, object] = {"arm": "site_finding", "n_pockets": int(len(top)),
                                      "any_positive_pocket": bool((top["label"] == 1).any())}
        if best is not None:
            summary["best_pocket_id"] = best["pocket_id"]
            summary["best_score"] = best["top_score"]
        out.append(summary)
    return out


def cmd_dock_one(args: argparse.Namespace) -> int:
    """Child process: one copy, one arm group; writes a JSON record."""
    row = json.loads(args.row_json.read_text())
    EXHAUSTIVENESS[0] = args.exhaustiveness
    record: Dict[str, object] = {"copy_key": row["copy_key"], "pdb_id": row["pdb_id"], "arm_group": args.arm}
    start = time.time()
    try:
        ctx = Context(row, args.structures_dir, args.ccd_dir, args.work_dir)
        record["complete"] = bool(ctx.complete)
        if args.arm == "primary":
            record["runs"] = arms_primary(ctx)
        elif args.arm == "secondary":
            record["runs"] = arms_secondary(ctx)
        elif args.arm == "pockets":
            sites = [np.array(s) for s in json.loads(row.get("_site_centroids", "[]"))]
            record["runs"] = arms_pockets_with_sites(ctx, sites, bool(row.get("_tertiary")))
        elif args.arm == "alphafold":
            record["runs"] = arms_alphafold(ctx, args.af_dir)
    except Exception as exc:  # noqa: BLE001 - recorded as the copy's failure reason
        record["error"] = f"{type(exc).__name__}: {exc}"[:500]
        record["traceback"] = traceback.format_exc()[-1500:]
    record["seconds"] = time.time() - start
    args.out_json.write_text(json.dumps(record, default=_json_default))
    return 0


def arms_pockets_with_sites(ctx: Context, sites: Sequence[np.ndarray], tertiary: bool) -> List[Dict[str, object]]:
    """Pocket arms, then for the site-finding runs the distance from each best pose to every true site."""
    out = arms_pockets(ctx, sites, tertiary)
    if not tertiary or not sites:
        return out
    runs = [r for r in out if str(r.get("arm", "")).startswith("pocket")]
    scored = [r for r in runs if r.get("top_score") is not None]
    summary = next(r for r in out if r.get("arm") == "site_finding")
    if scored:
        best = min(scored, key=lambda r: r["top_score"])
        centroid = np.asarray(best["top_pose_centroid"], dtype=float)
        dists = [float(np.linalg.norm(centroid - s)) for s in sites]
        summary["best_pose_distance_to_nearest_site"] = min(dists)
        summary["lands_in_true_site"] = bool(min(dists) <= TRUE_SITE_DISTANCE)
    return out


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)


# ------------------------------------------------------------- alphafold
def chain_sequence(arrays, chain: str) -> Tuple[str, List[int]]:
    from Bio.PDB.Polypeptide import three_to_index, index_to_one

    seq, numbers = [], []
    mask = (arrays.chain_ids == chain) & (arrays.atom_names == "CA") & arrays.is_polymer
    for i in np.flatnonzero(mask):
        try:
            seq.append(index_to_one(three_to_index(str(arrays.resnames[i]))))
        except (KeyError, ValueError):
            seq.append("X")
        numbers.append(int(arrays.resseqs[i]))
    return "".join(seq), numbers


def align_sequences(a: str, b: str) -> Tuple[List[Tuple[int, int]], float, float]:
    """Global alignment; returns aligned index pairs, identity over aligned pairs, coverage of ``a``."""
    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 2, -1
    aligner.open_gap_score, aligner.extend_gap_score = -5, -0.5
    aln = aligner.align(a, b)[0]
    pairs = []
    for (sa, ea), (sb, eb) in zip(*aln.aligned):
        pairs += list(zip(range(sa, ea), range(sb, eb)))
    same = sum(1 for i, j in pairs if a[i] == b[j])
    identity = same / len(pairs) if pairs else 0.0
    return pairs, identity, len(pairs) / max(1, len(a))


def kabsch(mobile: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t minimising |mobile @ R.T + t - target|."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, tc - mc @ rot.T


def arms_alphafold(ctx: Context, af_dir: Path) -> List[Dict[str, object]]:
    import dataclasses

    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.docking.engine import Receptor
    from cryptic_ip.docking.receptor import prepare_receptor
    from scipy.spatial import cKDTree

    arrays = ctx.arrays
    polymer = np.flatnonzero(arrays.is_polymer)
    tree = cKDTree(arrays.coords[polymer])
    near = sorted({int(polymer[i]) for idx in tree.query_ball_point(ctx.xyz, 4.5) for i in idx})
    if not near:
        return [{"arm": "alphafold_s1", "error": "no contacting chain"}]
    chains = pd.Series(arrays.chain_ids[near]).value_counts()
    chain = str(chains.index[0])
    seq, numbers = chain_sequence(arrays, chain)
    accessions = [a for a in str(ctx.row.get("uniprot_ids") or "").replace(",", ";").split(";") if a.strip()]
    best = None
    for acc in accessions:
        path = next(iter(sorted(af_dir.glob(f"AF-{acc.strip()}-F1-*.pdb"))), None)
        if path is None:
            continue
        model = load_structure_arrays(path)
        mseq, mnum = chain_sequence(model, str(model.chain_ids[0]))
        if len(mseq) > AF_MAX_LENGTH:
            continue
        pairs, identity, coverage = align_sequences(seq, mseq)
        if best is None or identity > best[2]:
            best = (acc, model, identity, coverage, pairs, mnum, path.name)
    if best is None:
        return [{"arm": "alphafold_s1", "error": "no AlphaFold model for the entry's accessions"}]
    acc, model, identity, coverage, pairs, mnum, name = best
    info = {"arm": "alphafold_s1", "accession": acc, "model": name, "chain": chain,
            "identity": identity, "coverage": coverage}
    if identity < AF_MIN_IDENTITY or coverage < AF_MIN_COVERAGE:
        info["error"] = f"model does not match the chain (identity {identity:.2f}, coverage {coverage:.2f})"
        return [info]
    ca = (arrays.chain_ids == chain) & (arrays.atom_names == "CA") & arrays.is_polymer
    ca_idx = np.flatnonzero(ca)
    site = set()
    site_tree = cKDTree(ctx.xyz)
    for i in np.flatnonzero((arrays.chain_ids == chain) & arrays.is_polymer):
        if site_tree.query_ball_point(arrays.coords[i], SITE_RESIDUE_DISTANCE):
            site.add(int(arrays.resseqs[i]))
    model_ca = np.flatnonzero((model.atom_names == "CA") & model.is_polymer)
    xs, ys = [], []
    for i, j in pairs:
        if numbers[i] in site:
            xs.append(model.coords[model_ca[j]])
            ys.append(arrays.coords[ca_idx[i]])
    info["site_ca_pairs"] = len(xs)
    if len(xs) < 3:
        info["error"] = "fewer than 3 site C-alpha pairs"
        return [info]
    rot, t = kabsch(np.array(xs), np.array(ys))
    moved = np.array(xs) @ rot.T + t
    info["site_ca_rmsd"] = float(np.sqrt(np.mean(np.sum((moved - np.array(ys)) ** 2, axis=1))))
    placed = dataclasses.replace(model, coords=model.coords @ rot.T + t)
    prepared = prepare_receptor(placed, ctx.work / f"af_{acc}", name="af_receptor")
    run = run_arm(ctx, Receptor(prepared.pdbqt), arm="alphafold_s1")
    run.update({k: v for k, v in info.items() if k != "arm"})
    return [run]


# ---------------------------------------------------------------- driver
def _child(argv: List[str]) -> None:
    main(argv)


def cmd_dock(args: argparse.Namespace) -> int:
    census = pd.read_csv(args.census, dtype={"icode": str}, keep_default_na=False, na_values=[""])
    census["icode"] = census["icode"].fillna("")
    selected = census[census["selected"].astype(str).str.lower() == "true"].copy()
    if args.arm in ("pockets", "alphafold"):
        selected = selected[selected["primary_set"].astype(str).str.lower() == "true"]
    mine = shard_entries(selected["pdb_id"].unique().tolist(), args.shard_index, args.shard_count)
    rows = selected[selected["pdb_id"].isin(mine)].sort_values(["pdb_id", "chain", "resseq"])
    if args.limit:
        rows = rows.head(args.limit)
    if args.skip_done is not None:
        done_keys = set()
        for path in args.skip_done.rglob(f"{args.arm}_*.jsonl"):
            for line in path.read_text().splitlines():
                if line.strip():
                    rec = json.loads(line)
                    if "runs" in rec:
                        done_keys.add(rec["copy_key"])
        rows = rows[~rows["copy_key"].isin(done_keys)]
        LOGGER.info("skipping %d copies docked in an earlier dispatch", len(done_keys))
    # Site-finding: the entry's first primary copy; true sites: every non-artefact copy of the entry.
    first = set(rows.groupby("pdb_id").head(1)["copy_key"])
    sites = {pdb: json.dumps(g[g["burial_class"] != "crystal_artifact"][["centroid_x", "centroid_y", "centroid_z"]]
                             .to_numpy(dtype=float).tolist()) for pdb, g in census.groupby("pdb_id")}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / f"{args.arm}_{args.shard_index}.jsonl"
    deadline = time.time() + args.budget_seconds
    done = not_reached = 0
    with open(out_path, "w") as fh:
        for _, row in rows.iterrows():
            data = {k: (None if _isnan(v) else v) for k, v in row.to_dict().items()}
            if time.time() > deadline:
                fh.write(json.dumps({"copy_key": data["copy_key"], "pdb_id": data["pdb_id"], "arm_group": args.arm,
                                     "error": "not reached: shard time budget"}) + "\n")
                not_reached += 1
                continue
            data["_tertiary"] = data["copy_key"] in first
            data["_site_centroids"] = sites.get(data["pdb_id"], "[]")
            key = data["copy_key"].replace(":", "_")
            row_json = args.work_dir / f"{key}_{args.arm}.row.json"
            out_json = args.work_dir / f"{key}_{args.arm}.out.json"
            args.work_dir.mkdir(parents=True, exist_ok=True)
            row_json.write_text(json.dumps(data, default=_json_default))
            child_args = ["dock-one", "--row-json", str(row_json), "--out-json", str(out_json), "--arm", args.arm,
                          "--structures-dir", str(args.structures_dir), "--ccd-dir", str(args.ccd_dir),
                          "--work-dir", str(args.work_dir), "--af-dir", str(args.af_dir),
                          "--exhaustiveness", str(args.exhaustiveness)]
            proc = mp.get_context("fork").Process(target=_child, args=(child_args,))
            proc.start()
            limit = min(args.copy_timeout, max(60.0, deadline + args.grace_seconds - time.time()))
            proc.join(limit)
            if proc.is_alive():
                proc.terminate()
                proc.join(10)
                record = {"copy_key": data["copy_key"], "pdb_id": data["pdb_id"], "arm_group": args.arm,
                          "error": f"timeout after {limit:.0f} s"}
            elif out_json.exists():
                record = json.loads(out_json.read_text())
            else:
                record = {"copy_key": data["copy_key"], "pdb_id": data["pdb_id"], "arm_group": args.arm,
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
    c = sub.add_parser("census")
    c.add_argument("--table", type=Path, required=True)
    c.add_argument("--entries", type=Path, required=True)
    c.add_argument("--structures-dir", type=Path, required=True)
    c.add_argument("--ccd-dir", type=Path, required=True)
    c.add_argument("--shard-index", type=int, default=0)
    c.add_argument("--shard-count", type=int, default=1)
    c.add_argument("--output", type=Path, required=True)
    for name in ("dock", "dock-one"):
        d = sub.add_parser(name)
        d.add_argument("--arm", choices=ARMS, required=True)
        d.add_argument("--structures-dir", type=Path, required=True)
        d.add_argument("--ccd-dir", type=Path, required=True)
        d.add_argument("--work-dir", type=Path, default=Path("/tmp/redock_work"))
        d.add_argument("--af-dir", type=Path, default=Path("af_models"))
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
            d.add_argument("--skip-done", type=Path, default=None,
                           help="arms directory of an earlier dispatch: copies docked there are skipped")
        else:
            d.add_argument("--row-json", type=Path, required=True)
            d.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "census":
        return cmd_census(args)
    if args.command == "dock":
        return cmd_dock(args)
    return cmd_dock_one(args)


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
