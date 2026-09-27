#!/usr/bin/env python3
"""Study K (docs/KINASE_PLAN.md): does restoring the nucleotide cosubstrate rescue docking?

Every inositol-phosphate kinase copy in study A's census fails, all as *scoring*
failures, often with a near-native pose sitting in the list. They are catalytic sites:
in the crystal the inositol phosphate is held against an ATP-analogue and its metals,
and the redocking protocol deletes both. This asks whether putting them back changes
the answer.

    python scripts/kinase.py cofactors --census census.csv --structures-dir S --out cofactors.json
    python scripts/kinase.py dock --census census.csv --cofactors cofactors.json \
        --structures-dir S --ccd-dir C --work-dir w --out-dir arms --shard 0 --shards 20
    python scripts/kinase.py report --census census.csv --cofactors cofactors.json \
        --copies copies.csv --apo-dir e32 --holo-dir arms --out-dir results/kinase

The apo arm is **not re-docked**: it is study F's pose lists at the same exhaustiveness,
seeds and starting poses. ``cryptic_ip/docking/receptor.py`` is not edited — appending
the cofactor to the receptor PDBQT happens here, with a test pinning the untouched path.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

LOGGER = logging.getLogger("kinase")

SEED = 20261004
N_BOOTSTRAP = 2000
EXHAUSTIVENESS = 32
N_POSES = 40
ENERGY_RANGE = 10.0
FROZEN_W = 0.1                    # study G's frozen re-ranking weight, reused unchanged

#: Nucleotide cosubstrates and their non-hydrolysable analogues, fixed in the plan.
COFACTORS = ("ATP", "ADP", "AMP", "ANP", "ACP", "AGS", "GTP", "GDP", "GNP", "ADX")
COFACTOR_DISTANCE = 6.0

#: The IPK superfamily, the plan's pre-specified stratum.
IPK_ACCESSIONS = {"N9UNA8": "EhIP6KA", "O43314": "PPIP5K2", "Q6PFW1": "PPIP5K1", "Q8NFU5": "IPMK",
                  "P23677": "ITPKA", "P27987": "ITPKB", "Q96DU7": "ITPKC", "Q13572": "ITPK1",
                  "Q92551": "IP6K1", "Q9UHH9": "IP6K2", "Q96PC2": "IP6K3", "Q9H8X2": "IPPK"}

#: Element -> AutoDock type for the appended cofactor. Nitrogen and oxygen are typed as
#: acceptors, which is what they overwhelmingly are in a nucleotide; aromatic carbons are
#: typed C rather than A. Stated as a limitation: Vina's scoring is typed, so this
#: slightly under-rewards stacking against the base.
ELEMENT_TYPE = {"C": "C", "N": "NA", "O": "OA", "P": "P", "S": "SA", "F": "F", "CL": "Cl",
                "BR": "Br", "I": "I", "MG": "Mg", "MN": "Mn", "ZN": "Zn", "CA": "Ca", "FE": "Fe"}


# ------------------------------------------------------------- cofactor census
def cofactor_atoms(arrays, ligand_xyz: np.ndarray, cutoff: float = COFACTOR_DISTANCE) -> List[dict]:
    """Heavy atoms of every nucleotide cofactor residue within ``cutoff`` of the ligand."""
    names = np.char.upper(arrays.resnames.astype(str))
    hits: Dict[Tuple[str, str, int, str], List[int]] = {}
    for i in np.flatnonzero(np.isin(names, list(COFACTORS))):
        if arrays.is_polymer[i] or str(arrays.elements[i]).upper() == "H":
            continue
        key = (str(arrays.resnames[i]), str(arrays.chain_ids[i]), int(arrays.resseqs[i]),
               str(arrays.icodes[i]))
        hits.setdefault(key, []).append(int(i))
    out = []
    for (resname, chain, resseq, icode), idx in sorted(hits.items()):
        coords = arrays.coords[idx]
        distance = float(np.min(np.linalg.norm(ligand_xyz[:, None, :] - coords[None, :, :], axis=2)))
        if distance > cutoff:
            continue
        out.append({"resname": resname, "chain": chain, "resseq": resseq, "icode": icode,
                    "min_distance": distance, "n_atoms": len(idx),
                    "atoms": [{"name": str(arrays.atom_names[i]), "element": str(arrays.elements[i]).upper(),
                               "xyz": [float(c) for c in arrays.coords[i]]} for i in idx]})
    return out


def cmd_cofactors(args: argparse.Namespace) -> int:
    from redocking import _isnan

    census = pd.read_csv(args.census, dtype={"icode": str}, keep_default_na=False, na_values=[""])
    rows = census[census["primary_set"].astype(str) == "True"] if "primary_set" in census else census
    out: Dict[str, object] = {}
    for _, row in rows.iterrows():
        key = str(row["copy_key"])
        entry: Dict[str, object] = {"pdb_id": str(row["pdb_id"]), "species": str(row.get("species")),
                                    "uniprot_ids": str(row.get("uniprot_ids") or ""),
                                    "burial_class": str(row.get("burial_class"))}
        try:
            from cryptic_ip.analysis.structure_arrays import load_structure_arrays
            from redocking import find_copy, structure_path

            path = structure_path(args.structures_dir, str(row["pdb_id"]))
            if path is None:
                raise FileNotFoundError(f"{row['pdb_id']} not fetched")
            arrays = load_structure_arrays(path)
            icode = str(row.get("icode") or "") if not _isnan(row.get("icode")) else ""
            _, _, atoms = find_copy(arrays, str(row["chain"]), int(row["resseq"]), icode)
            heavy = atoms[arrays.elements[atoms] != "H"]
            found = cofactor_atoms(arrays, arrays.coords[heavy])
            entry["cofactors"] = [{k: v for k, v in c.items() if k != "atoms"} for c in found]
            entry["qualifies"] = bool(found)
        except Exception as exc:  # noqa: BLE001 - recorded as this copy's reason
            entry["error"] = f"{type(exc).__name__}: {exc}"[:300]
            entry["qualifies"] = False
        entry["ipk"] = ipk_label(str(row.get("uniprot_ids") or ""))
        out[key] = entry
        if entry.get("qualifies"):
            LOGGER.info("%s (%s): %s", key, entry["ipk"] or "-",
                        [c["resname"] for c in entry.get("cofactors", [])])
    args.out.write_text(json.dumps(out, indent=2, default=str))
    qualifying = sum(1 for v in out.values() if v.get("qualifies"))
    print(json.dumps({"copies": len(out), "qualifying": qualifying,
                      "ipk": sum(1 for v in out.values() if v.get("ipk"))}))
    return 0


def ipk_label(uniprot_ids: str) -> Optional[str]:
    for acc in str(uniprot_ids).replace(";", ",").split(","):
        if acc.strip() in IPK_ACCESSIONS:
            return IPK_ACCESSIONS[acc.strip()]
    return None


# ------------------------------------------------------------- holo receptor
def cofactor_charges(resname: str, atoms: Sequence[dict], ccd_dir: Optional[Path]) -> Dict[str, float]:
    """Gasteiger charges from the CCD template, by atom name. Empty when unavailable.

    Vina ignores charges, so the primary arms are unaffected by an empty result; only
    K4, the electrostatic arm, needs them and it reports how many copies carried them.
    """
    if ccd_dir is None:
        return {}
    path = ccd_dir / f"{resname}.cif"
    if not path.exists():
        return {}
    try:
        from rdkit.Chem import AllChem

        from redocking import ccd_template

        mol, _, _ = ccd_template(path)
        AllChem.ComputeGasteigerCharges(mol)
        out = {}
        for atom in mol.GetAtoms():
            name = atom.GetPropsAsDict().get("molFileAlias") or atom.GetPropsAsDict().get("atomName")
            charge = float(atom.GetDoubleProp("_GasteigerCharge"))
            if name and np.isfinite(charge):
                out[str(name)] = charge
        return out
    except Exception as exc:  # noqa: BLE001 - charges are best-effort
        LOGGER.warning("no charges for %s: %s", resname, exc)
        return {}


def holo_receptor(apo_pdbqt: Path, cofactors: Sequence[dict], out_path: Path,
                  ccd_dir: Optional[Path] = None) -> Dict[str, object]:
    """The apo receptor PDBQT with the cofactor heavy atoms appended as rigid HETATM."""
    from cryptic_ip.docking.receptor import Atom, pdbqt_line

    lines = [ln for ln in apo_pdbqt.read_text().splitlines() if ln.strip()]
    serial = len(lines)
    appended, untyped, charged = 0, [], 0
    for residue in cofactors:
        charges = cofactor_charges(str(residue["resname"]), residue["atoms"], ccd_dir)
        for atom in residue["atoms"]:
            element = str(atom["element"]).upper()
            atype = ELEMENT_TYPE.get(element)
            if atype is None:
                untyped.append(element)
                continue
            serial += 1
            q = float(charges.get(str(atom["name"]), 0.0))
            charged += int(bool(charges))
            a = Atom("HETATM", str(atom["name"]), str(residue["resname"]), str(residue["chain"]),
                     int(residue["resseq"]), str(residue["icode"]), np.asarray(atom["xyz"], dtype=float),
                     element, q)
            lines.append(pdbqt_line(serial, a, atype))
            appended += 1
    out_path.write_text("\n".join(lines) + "\n")
    return {"appended": appended, "untyped": sorted(set(untyped)), "charged": charged,
            "residues": [f"{c['resname']}:{c['chain']}:{c['resseq']}" for c in cofactors]}


# -------------------------------------------------------------------- docking
def dock_copy(ctx, holo_pdbqt: Optional[Path]) -> List[Dict[str, object]]:
    """Study F's docking path exactly; only the receptor PDBQT differs."""
    from cryptic_ip.docking import engine
    from cryptic_ip.docking.ligand import start_pose, to_pdbqt
    from cryptic_ip.rescoring.electrostatics import ReceptorField, read_pdbqt, read_pdbqt_models

    if holo_pdbqt is None:
        receptor, _ = ctx.receptor()
    else:
        receptor = engine.Receptor(holo_pdbqt)
    field = ReceptorField(read_pdbqt(Path(receptor.pdbqt).read_text()))
    size = [ctx.side] * 3
    runs: List[Dict[str, object]] = []
    for seed in engine.SEEDS:
        pose = start_pose(ctx.ligands["primary"], ctx.centre, seed, crystal=ctx.crystal)
        result = engine.dock(receptor, to_pdbqt(pose.mol), ctx.centre, size, seed=seed,
                             exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES, energy_range=ENERGY_RANGE)
        table = engine.pose_table(result, ctx.crystal, site_centroid=ctx.centre)
        models = read_pdbqt_models(result.poses_pdbqt)
        if len(models) != len(table.scores):
            raise RuntimeError(f"{len(models)} PDBQT models but {len(table.scores)} scored poses")
        runs.append({"seed": seed, "vina": table.scores, "rmsd": table.rmsd,
                     "eel": [field.energy(m) for m in models], "seconds": result.seconds})
    return runs


def qualifying_copies(cofactors: Dict[str, dict]) -> List[str]:
    return sorted(k for k, v in cofactors.items() if v.get("qualifies"))


def cmd_dock(args: argparse.Namespace) -> int:
    from redocking import Context, _isnan, _json_default
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.docking.receptor import nearby_metals
    from redocking import find_copy, structure_path

    cofactors = json.loads(args.cofactors.read_text())
    census = pd.read_csv(args.census, dtype={"icode": str}, keep_default_na=False, na_values=[""])
    census["copy_key"] = census["copy_key"].astype(str)
    keys = qualifying_copies(cofactors)
    if args.shards > 1:
        keys = keys[args.shard::args.shards]
    LOGGER.info("shard %d of %d: %d qualifying copies", args.shard, args.shards, len(keys))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / f"holo_{args.shard}.jsonl"
    with out_path.open("w") as handle:
        for key in keys:
            row = census[census["copy_key"] == key]
            record: Dict[str, object] = {"copy_key": key}
            start = time.time()
            try:
                if row.empty:
                    raise KeyError(f"{key} is not in the census")
                row = row.iloc[0].to_dict()
                ctx = Context(row, args.structures_dir, args.ccd_dir, args.work_dir)
                apo, _ = ctx.receptor()
                path = structure_path(args.structures_dir, str(row["pdb_id"]))
                arrays = load_structure_arrays(path)
                icode = str(row.get("icode") or "") if not _isnan(row.get("icode")) else ""
                _, _, atoms = find_copy(arrays, str(row["chain"]), int(row["resseq"]), icode)
                heavy = atoms[arrays.elements[atoms] != "H"]
                residues = cofactor_atoms(arrays, arrays.coords[heavy])
                metals = nearby_metals(arrays, arrays.coords[heavy])
                holo = ctx.work / "receptor_holo.pdbqt"
                if metals:
                    apo, _ = ctx.receptor(metals=metals, name="receptor_metals")
                record["build"] = holo_receptor(Path(apo.pdbqt), residues, holo, args.ccd_dir)
                record["build"]["metals"] = len(metals)
                record["runs"] = dock_copy(ctx, holo)
            except Exception as exc:  # noqa: BLE001 - recorded as this copy's failure reason
                record["error"] = f"{type(exc).__name__}: {exc}"[:500]
                record["traceback"] = traceback.format_exc()[-1200:]
            record["seconds"] = time.time() - start
            handle.write(json.dumps(record, default=_json_default) + "\n")
            handle.flush()
            LOGGER.info("%s: %s", key, record.get("error") or f"{len(record.get('runs', []))} runs")
    print(json.dumps({"copies": len(keys)}))
    return 0


# --------------------------------------------------------------------- report
def build(census: pd.DataFrame, cofactors: Dict[str, dict], apo: Sequence[dict],
          holo: Sequence[dict], copies: Optional[pd.DataFrame] = None,
          n_bootstrap: int = N_BOOTSTRAP) -> Dict[str, object]:
    from sampling_report import arm_table, decide, per_copy
    from cryptic_ip.benchmark import protocol
    from cryptic_ip.docking import stats

    census = census.copy()
    census["copy_key"] = census["copy_key"].astype(str)
    group_of = dict(zip(census["copy_key"], census["homology_group_strict"].astype(str)))

    runs_apo, runs_holo = per_copy(apo), per_copy(holo)
    shared = sorted({k for k, _ in runs_apo} & {k for k, _ in runs_holo})
    report: Dict[str, object] = {
        "plan": "docs/KINASE_PLAN.md", "cofactor_distance": COFACTOR_DISTANCE,
        "qualifying": len(qualifying_copies(cofactors)),
        "cofactor_residues": _residue_counts(cofactors),
        "copies_compared": len(shared),
        "accounting": {"holo_records": len(holo), "holo_errors": sum(1 for r in holo if "error" in r)},
    }
    report["K0_family"] = family_table(census, copies)
    report["K3_pyrophosphates"] = species_table(census, copies, n_bootstrap)
    if not shared:
        report["K1"] = {"decision": "not evaluable: no copy has both arms"}
        report["notes"] = _notes()
        return report

    tables = {"apo": arm_table(runs_apo, shared), "holo": arm_table(runs_holo, shared)}
    groups = np.array([group_of.get(k, "?") for k in shared])
    n_groups = int(pd.Series(groups).nunique())
    report["groups"] = n_groups
    report["levels"] = {name: {m: stats.mean_estimates(t[m].to_numpy(), groups, null=0.5,
                                                       n_bootstrap=n_bootstrap, seed=SEED)
                               for m in ("ceiling", "vina", "reranked")}
                        for name, t in tables.items()}

    estimates, p_values = {}, {}
    for name, metric, question in (("K1", "vina", "top-pose success, holo - apo"),
                                   ("K2", "ceiling", "best-of-list ceiling, holo - apo")):
        est = stats.paired_difference(tables["holo"][metric].to_numpy(), tables["apo"][metric].to_numpy(),
                                      groups, n_bootstrap=n_bootstrap, seed=SEED)
        estimates[name] = (est, question, metric)
        p_values[name] = est["per_group"]["p_value"]
    holm = protocol.holm(p_values)
    for name, (est, question, _) in estimates.items():
        report[name] = {"question": question, "estimate": est["per_group"], "per_copy": est["per_copy"],
                        "copies": len(shared), "groups": n_groups, "holm_p": holm.get(name),
                        "decision": decide(est["per_group"], n_groups, holm.get(name),
                                           good="better with the cofactor",
                                           null_label="no material change")}
    report["holm"] = holm

    # K4: the electrostatic arm, which is the only one that reads the cofactor's charges.
    charged = sum(1 for r in holo if (r.get("build") or {}).get("charged"))
    k4 = stats.paired_difference(tables["holo"]["reranked"].to_numpy(), tables["holo"]["vina"].to_numpy(),
                                 groups, n_bootstrap=n_bootstrap, seed=SEED)
    report["K4_electrostatics"] = {
        "question": "re-ranked (w = 0.1) minus Vina top-pose, both on the holo receptor",
        "estimate": k4["per_group"], "copies_with_cofactor_charges": charged,
        "note": ("the cofactor carries Gasteiger charges from its CCD template" if charged else
                 "no copy carried cofactor charges, so this repeats study G's re-ranking on a "
                 "receptor whose cofactor is uncharged and is not evidence about electrostatics")}

    # The metals-only decomposition, quoted from study A rather than re-docked.
    if copies is not None and "metals_success" in copies.columns:
        m = copies[copies["copy_key"].astype(str).isin(shared)]
        report["metals_only"] = {
            "copies": int(m["metals_success"].notna().sum()),
            "mean": float(pd.to_numeric(m["metals_success"], errors="coerce").mean()),
            "note": "study A's metals-only arm on the same copies, quoted not re-docked, so the "
                    "cofactor's contribution can be separated from the metals'"}
    report["notes"] = _notes()
    return report


def _notes() -> List[str]:
    return [
        "The IPK slice of study A was read before this plan was written, so K0 and K3 are a "
        "re-analysis of data in hand, not a test. K1, K2 and K4 are blind.",
        "Vina's scoring function does not read partial charges, so in K1 and K2 the cofactor is a "
        "shaped, typed occluder: it restores the site's shape and hydrogen bonding, not its "
        "electrostatics.",
        "The holo arm adds the cofactor and the catalytic metals together. Study A's metals-only "
        "arm is quoted on the same copies so the two contributions are not confounded.",
        "A win here would narrow where the protocol may be trusted; it would not rescue the "
        "proteome ranking, whose candidate pockets are not catalytic sites.",
    ]


def _residue_counts(cofactors: Dict[str, dict]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for v in cofactors.values():
        for c in v.get("cofactors", []) or []:
            counts[str(c["resname"])] = counts.get(str(c["resname"]), 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))


def family_table(census: pd.DataFrame, copies: Optional[pd.DataFrame]) -> Dict[str, object]:
    """K0, descriptive: the IPK superfamily's outcome in study A, read before the plan."""
    if copies is None:
        return {"note": "copies.csv not supplied"}
    frame = copies.copy()
    frame["ipk"] = frame.get("uniprot_ids", pd.Series(dtype=str)).astype(str).map(ipk_label)
    fam = frame[frame["ipk"].notna()]
    success = pd.to_numeric(fam.get("success"), errors="coerce")
    return {"copies": int(len(fam)), "with_outcome": int(success.notna().sum()),
            "successes": int(np.nansum(success.to_numpy())),
            "proteins": sorted(fam["ipk"].dropna().unique().tolist()),
            "burial_classes": fam.get("burial_class", pd.Series(dtype=str)).value_counts().to_dict(),
            "failure_kinds": fam.get("failure_kind", pd.Series(dtype=str)).value_counts().to_dict(),
            "note": "descriptive: this slice was read before the plan was written"}


def species_table(census: pd.DataFrame, copies: Optional[pd.DataFrame], n_bootstrap: int) -> Dict[str, object]:
    """K3, descriptive: the diphosphoinositols pooled as one PP-IP stratum."""
    from cryptic_ip.docking import stats

    if copies is None or "species" not in copies.columns:
        return {"note": "copies.csv not supplied"}
    frame = copies.copy()
    frame["success"] = pd.to_numeric(frame.get("success"), errors="coerce")
    frame = frame[frame["success"].notna()]
    out: Dict[str, object] = {}
    for label, mask in (("PP-IP (InsP7 + InsP8)", frame["species"].isin(["InsP7", "InsP8"])),
                        ("InsP6", frame["species"] == "InsP6")):
        sub = frame[mask]
        if sub.empty:
            continue
        out[label] = {"copies": int(len(sub)),
                      "estimate": stats.mean_estimates(sub["success"].to_numpy(),
                                                       sub["homology_group_strict"].astype(str).to_numpy(),
                                                       null=0.0, n_bootstrap=n_bootstrap, seed=SEED)}
    out["note"] = "descriptive: read before the plan was written"
    return out


def markdown(r: Dict[str, object]) -> str:
    def fmt(e):
        return "-" if not e or "point" not in e else f"{e['point']:.3f} [{e['low']:.3f}, {e['high']:.3f}]"

    k0, k1, k2 = r.get("K0_family", {}), r.get("K1", {}), r.get("K2", {})
    lines = ["## The IP kinases and the missing cosubstrate (docs/KINASE_PLAN.md)", "",
             f"{r.get('qualifying')} of the primary-set copies carry a nucleotide cofactor within "
             f"{r.get('cofactor_distance')} Å of the ligand: {json.dumps(r.get('cofactor_residues'))}.", "",
             f"**K0 (descriptive):** the IPK superfamily is {k0.get('successes')} of "
             f"{k0.get('with_outcome')} on top-pose success — {', '.join(k0.get('proteins') or [])}. "
             f"Failure kinds: {json.dumps(k0.get('failure_kinds'))}.", "",
             f"**K1 (primary):** {k1.get('decision')} — {fmt(k1.get('estimate'))} over "
             f"{k1.get('copies')} copies in {k1.get('groups')} strict groups, Holm p "
             f"{k1.get('holm_p')}.", "",
             f"**K2 (ceiling):** {k2.get('decision')} — {fmt(k2.get('estimate'))}.", "",
             "| arm | ceiling | Vina top pose | re-ranked |", "|---|---|---|---|"]
    for name, v in (r.get("levels") or {}).items():
        lines.append(f"| {name} | {fmt((v.get('ceiling') or {}).get('per_group'))} | "
                     f"{fmt((v.get('vina') or {}).get('per_group'))} | "
                     f"{fmt((v.get('reranked') or {}).get('per_group'))} |")
    k3 = r.get("K3_pyrophosphates") or {}
    lines += ["", "| ligand stratum | copies | top-pose success |", "|---|---|---|"]
    for label, v in k3.items():
        if isinstance(v, dict) and "estimate" in v:
            lines.append(f"| {label} | {v.get('copies')} | "
                         f"{fmt((v['estimate'] or {}).get('per_group'))} |")
    k4 = r.get("K4_electrostatics") or {}
    lines += ["", f"**K4 (electrostatics on the holo receptor):** {fmt(k4.get('estimate'))}; "
              f"{k4.get('copies_with_cofactor_charges')} copies carried cofactor charges. "
              f"{k4.get('note')}", ""]
    if r.get("metals_only"):
        m = r["metals_only"]
        lines += [f"**Metals only (quoted from study A):** mean {m.get('mean')} over {m.get('copies')} "
                  f"copies. {m.get('note')}", ""]
    lines += [f"- {n}" for n in r.get("notes", [])]
    return "\n".join(lines) + "\n"


def cmd_report(args: argparse.Namespace) -> int:
    from sampling_report import read_records

    census = pd.read_csv(args.census, dtype={"icode": str}, keep_default_na=False, na_values=[""])
    cofactors = json.loads(args.cofactors.read_text())
    copies = pd.read_csv(args.copies, low_memory=False) if args.copies and args.copies.exists() else None
    apo = read_records(args.apo_dir, "*.jsonl")
    holo = read_records(args.holo_dir, "*.jsonl")
    report = build(census, cofactors, apo, holo, copies, args.n_bootstrap)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "kinase.json").write_text(json.dumps(report, indent=2, default=str))
    text = markdown(report)
    (args.out_dir / "KINASE.md").write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    c = sub.add_parser("cofactors")
    c.add_argument("--census", type=Path, required=True)
    c.add_argument("--structures-dir", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)

    d = sub.add_parser("dock")
    d.add_argument("--census", type=Path, required=True)
    d.add_argument("--cofactors", type=Path, required=True)
    d.add_argument("--structures-dir", type=Path, required=True)
    d.add_argument("--ccd-dir", type=Path, required=True)
    d.add_argument("--work-dir", type=Path, required=True)
    d.add_argument("--out-dir", type=Path, required=True)
    d.add_argument("--shard", type=int, default=0)
    d.add_argument("--shards", type=int, default=1)

    r = sub.add_parser("report")
    r.add_argument("--census", type=Path, required=True)
    r.add_argument("--cofactors", type=Path, required=True)
    r.add_argument("--copies", type=Path)
    r.add_argument("--apo-dir", type=Path)
    r.add_argument("--holo-dir", type=Path)
    r.add_argument("--out-dir", type=Path, required=True)
    r.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)

    args = parser.parse_args(argv)
    return {"cofactors": cmd_cofactors, "dock": cmd_dock, "report": cmd_report}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
