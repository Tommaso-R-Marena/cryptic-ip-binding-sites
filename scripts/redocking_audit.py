#!/usr/bin/env python3
"""Audit the redocking census's configuration exclusions (results/redocking/ledger.jsonl).

For every copy the census excluded because its configuration differs from the
CCD template, report which stereocentres differ and how well defined each is in
the deposited coordinates: the signed chiral volume at the centre (Å³) computed
from its three heavy-atom neighbours. A centre with a volume near zero is nearly
planar - its configuration is not determined by the model - whereas a large
volume of the wrong sign is a genuinely different configuration.

    python scripts/redocking_audit.py --census census.csv --structures-dir s --ccd-dir ccd --output audit.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


def chiral_volume(xyz: np.ndarray, centre: int, neighbours: Sequence[int]) -> float:
    a, b, c = (xyz[n] - xyz[centre] for n in neighbours[:3])
    return float(np.dot(a, np.cross(b, c)))


def centre_report(crystal, template) -> List[Dict[str, object]]:
    """Per template stereocentre (carbon): template and crystal CIP labels and the crystal chiral volume."""
    from rdkit import Chem

    from rdkit.Chem import rdCIPLabeler

    tmpl = Chem.RemoveHs(template)
    cry = Chem.Mol(crystal)
    matches = cry.GetSubstructMatches(tmpl, useChirality=False, uniquify=False, maxMatches=100000)
    if not matches:
        return [{"error": "no substructure match"}]
    rdCIPLabeler.AssignCIPLabels(tmpl)
    rdCIPLabeler.AssignCIPLabels(cry)
    xyz = cry.GetConformer().GetPositions()
    centres = [(a.GetIdx(), a.GetProp("_CIPCode")) for a in tmpl.GetAtoms()
               if a.GetAtomicNum() == 6 and a.HasProp("_CIPCode")]

    def rows_for(match):
        out = []
        for t_idx, label in centres:
            atom = cry.GetAtomWithIdx(match[t_idx])
            heavy = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
            cry_label = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "?"
            out.append({"template_atom": t_idx, "template_cip": label, "crystal_cip": cry_label,
                        "differs": label != cry_label,
                        "abs_chiral_volume": abs(chiral_volume(xyz, match[t_idx], heavy))
                        if len(heavy) >= 3 else float("nan")})
        return out

    # Symmetry-equivalent atoms of a meso compound carry opposite labels, so an
    # arbitrary mapping can report differences that are not there: take the
    # mapping with the fewest.
    return min((rows_for(m) for m in matches), key=lambda rows: sum(r["differs"] for r in rows))


def audit_copy(row: Dict[str, object], structures_dir: Path, ccd_dir: Path) -> Dict[str, object]:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.docking.ligand import crystal_molecule
    from redocking import ccd_template, find_copy, structure_path

    path = structure_path(structures_dir, str(row["pdb_id"]))
    arrays = load_structure_arrays(path)
    icode = "" if pd.isna(row.get("icode")) else str(row.get("icode") or "")
    _, comp, atoms = find_copy(arrays, str(row["chain"]), int(row["resseq"]), icode)
    heavy = atoms[arrays.elements[atoms] != "H"]
    parent, _, _ = ccd_template(ccd_dir / f"{comp}.cif")
    crystal, complete = crystal_molecule(arrays.elements[heavy], arrays.coords[heavy], arrays.atom_names[heavy],
                                         parent)
    centres = centre_report(crystal, parent)
    differing = [c for c in centres if c.get("differs")]
    return {"copy_key": row["copy_key"], "comp_id": comp, "resolution": row.get("resolution"),
            "method": row.get("experimental_method"), "n_centres": len(centres), "n_differing": len(differing),
            "differing_abs_volumes": json.dumps([round(c["abs_chiral_volume"], 3) for c in differing]),
            "min_abs_volume_all": float(np.nanmin([c["abs_chiral_volume"] for c in centres]))
            if centres and "abs_chiral_volume" in centres[0] else float("nan"),
            "labels": json.dumps([(c["template_cip"], c["crystal_cip"]) for c in centres])}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--structures-dir", type=Path, required=True)
    parser.add_argument("--ccd-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    census = pd.read_csv(args.census, dtype={"icode": str})
    mismatched = census[census["status"] == "excluded: configuration differs from the CCD"]
    rows = []
    for row in mismatched.to_dict(orient="records"):
        try:
            rows.append(audit_copy(row, args.structures_dir, args.ccd_dir))
        except Exception as exc:  # noqa: BLE001 - an audit row, not a result
            rows.append({"copy_key": row["copy_key"], "error": f"{type(exc).__name__}: {exc}"[:200]})
    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(out.drop(columns=[c for c in ("labels",) if c in out]).to_string()[:20000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
