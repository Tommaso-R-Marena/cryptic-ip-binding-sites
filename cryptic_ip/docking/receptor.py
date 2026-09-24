"""Receptor preparation for redocking (docs/REDOCKING_PLAN.md).

1. **Strip** every non-polymer atom - waters, ions, every ligand, every other
   inositol phosphate copy, cofactors - so nothing templates the pose. All
   polymer chains of the deposited asymmetric unit are kept, because sites sit
   at interfaces. Selenomethionine is written as methionine; other modified
   polymer residues are removed and counted.
2. **Protonate** with PDB2PQR (AMBER force field) using PROPKA at pH 7.4. PROPKA
   assigns titratable residues, including histidine (HID/HIE/HIP), from its
   pKa estimates; PDB2PQR then optimises the hydrogen-bond network.
3. **Type** atoms for AutoDock (:func:`write_pdbqt`): polar hydrogens kept as
   ``HD``, non-polar hydrogens merged (their charges added to the parent atom);
   aromatic ring carbons ``A``; acceptor nitrogens without a hydrogen in
   histidine and nucleobase rings ``NA``; all oxygens ``OA``; sulfur ``SA``;
   metals by element. Charges are PDB2PQR's AMBER charges (used only by AD4).

The typer is ours, not Meeko's receptor builder, so that every structure goes
through one deterministic path: Meeko's builder rejects residues with
non-ideal geometry, which would drop copies non-randomly.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

#: Metals flagged within this distance of a ligand heavy atom (Å).
METAL_DISTANCE = 3.0
METAL_TYPES: Dict[str, str] = {"MG": "Mg", "ZN": "Zn", "CA": "Ca", "MN": "Mn", "FE": "Fe"}
METAL_CHARGES: Dict[str, float] = {"MG": 2.0, "ZN": 2.0, "CA": 2.0, "MN": 2.0, "FE": 2.0}
PH = 7.4

AROMATIC_CARBONS: Dict[str, frozenset] = {
    "PHE": frozenset({"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "TYR": frozenset({"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "TRP": frozenset({"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"}),
}
HISTIDINES = frozenset({"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"})
HIS_RING_CARBONS = frozenset({"CG", "CD2", "CE1"})
HIS_RING_NITROGENS = frozenset({"ND1", "NE2"})
NUCLEOTIDES = frozenset({
    "A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DU", "DI", "RA", "RC", "RG", "RU",
    "A3", "A5", "C3", "C5", "G3", "G5", "U3", "U5", "DA3", "DA5", "DC3", "DC5", "DG3", "DG5",
    "DT3", "DT5", "RA3", "RA5", "RC3", "RC5", "RG3", "RG5", "RU3", "RU5",
})
BASE_RING_CARBONS = frozenset({"C2", "C4", "C5", "C6", "C8"})
BASE_RING_NITROGENS = frozenset({"N1", "N3", "N7", "N9"})
H_BOND_MAX = 1.3
CHAIN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


def _glycosidic_nitrogen(resname: str) -> str:
    """N9 for purines, N1 for pyrimidines: bonded to the sugar, neither donor nor acceptor."""
    base = resname.lstrip("DR").rstrip("35") or resname
    return "N9" if base[:1] in ("A", "G", "I") else "N1"


class ReceptorError(RuntimeError):
    """The receptor cannot be prepared; the copy is excluded with this reason."""


@dataclass
class Atom:
    record: str
    name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    xyz: np.ndarray
    element: str
    charge: float = 0.0


@dataclass
class StripReport:
    polymer_atoms: int = 0
    removed_nonpolymer_atoms: int = 0
    removed_modified_residues: List[str] = field(default_factory=list)
    selenomethionines: int = 0


def _pdb_atom_line(serial: int, a: Atom, record: str = "ATOM") -> str:
    name = a.name if len(a.name) == 4 or len(a.element) == 2 else f" {a.name}"
    return (f"{record:<6s}{serial % 100000:5d} {name:<4s} {a.resname:>3s} {a.chain[:1]:1s}{a.resseq:4d}"
            f"{a.icode[:1]:1s}   {a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}"
            f"          {a.element:>2s}")


def write_polymer_pdb(arrays, out_path: Path) -> StripReport:
    """Polymer atoms only (all chains), MSE as MET, other modified residues removed."""
    from ..analysis.structure_arrays import CANONICAL_AMINO_ACIDS

    report = StripReport()
    lines: List[str] = []
    removed = set()
    serial = 0
    last_chain = None
    # PDB format holds one-character chain ids and four-digit residue numbers.
    # Each polymer chain gets its own character, so chains are never merged; an
    # entry that cannot be written faithfully is refused, not approximated.
    chains = list(dict.fromkeys(str(c) for c in arrays.chain_ids[arrays.is_polymer]))
    if len(chains) > len(CHAIN_LETTERS):
        raise ReceptorError(f"{len(chains)} polymer chains: more than PDB format can hold")
    letter = {c: (c if len(c) == 1 and all(len(x) == 1 for x in chains) else CHAIN_LETTERS[i])
              for i, c in enumerate(chains)}
    for i in range(arrays.n_atoms):
        if not arrays.is_polymer[i]:
            report.removed_nonpolymer_atoms += 0 if arrays.is_solvent[i] else 1
            continue
        resname = str(arrays.resnames[i])
        name = str(arrays.atom_names[i]).strip()
        element = str(arrays.elements[i]).upper()
        if element == "H":
            continue
        if resname == "MSE":
            resname = "MET"
            if name == "SE":
                name, element = "SD", "S"
            report.selenomethionines += 1
        elif resname not in CANONICAL_AMINO_ACIDS and resname not in NUCLEOTIDES:
            removed.add(f"{arrays.chain_ids[i]}:{resname}{arrays.resseqs[i]}{arrays.icodes[i]}".strip())
            continue
        chain = letter[str(arrays.chain_ids[i])]
        if not -999 <= int(arrays.resseqs[i]) <= 9999:
            raise ReceptorError(f"residue number {int(arrays.resseqs[i])} does not fit PDB format")
        if last_chain is not None and chain != last_chain:
            lines.append("TER")
        last_chain = chain
        serial += 1
        atom = Atom("ATOM", name, resname, chain, int(arrays.resseqs[i]), str(arrays.icodes[i]),
                    arrays.coords[i], element)
        lines.append(_pdb_atom_line(serial, atom))
        report.polymer_atoms += 1
    lines += ["TER", "END"]
    report.removed_modified_residues = sorted(removed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    return report


def nearby_metals(arrays, ligand_xyz: np.ndarray, cutoff: float = METAL_DISTANCE) -> List[Atom]:
    """Mg, Zn, Ca, Mn and Fe ions within ``cutoff`` of any ligand heavy atom."""
    out = []
    for i in np.flatnonzero(np.isin(np.char.upper(arrays.elements.astype(str)), list(METAL_TYPES))):
        if arrays.is_polymer[i]:
            continue
        d = np.min(np.linalg.norm(ligand_xyz - arrays.coords[i], axis=1))
        if d <= cutoff:
            el = str(arrays.elements[i]).upper()
            out.append(Atom("HETATM", el, str(arrays.resnames[i]), str(arrays.chain_ids[i]),
                            int(arrays.resseqs[i]), str(arrays.icodes[i]), arrays.coords[i], el,
                            METAL_CHARGES[el]))
    return out


def protonate(pdb_in: Path, workdir: Path, ph: float = PH, timeout: int = 1800) -> Tuple[Path, Path]:
    """PDB2PQR with PROPKA at ``ph``; returns ``(pdb_out, pqr_out)``."""
    exe = shutil.which("pdb2pqr") or shutil.which("pdb2pqr30")
    if exe is None:
        raise ReceptorError("pdb2pqr is not installed")
    pqr = workdir / "receptor.pqr"
    pdb = workdir / "receptor_h.pdb"
    cmd = [exe, "--ff=AMBER", "--ffout=AMBER", "--titration-state-method=propka", f"--with-ph={ph}",
           "--keep-chain", "--pdb-output", str(pdb), str(pdb_in), str(pqr)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0 or not pdb.exists() or not pqr.exists():
        tail = (proc.stderr or proc.stdout).strip().splitlines()[-3:]
        raise ReceptorError(f"pdb2pqr failed: {' | '.join(tail)}"[:400])
    return pdb, pqr


def read_protonated(pdb: Path, pqr: Path) -> List[Atom]:
    """Atoms of PDB2PQR's PDB output with the charges of its PQR (same order)."""
    atoms: List[Atom] = []
    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        element = line[76:78].strip().upper() or name.lstrip("0123456789")[:1]
        atoms.append(Atom(line[:6].strip(), name, line[17:20].strip(), line[21:22].strip() or "A",
                          int(line[22:26]), line[26:27].strip(),
                          np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), element))
    charges = []
    for line in pqr.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            fields = line.split()
            charges.append(float(fields[-2]))
    if len(charges) != len(atoms):
        raise ReceptorError(f"PDB2PQR outputs disagree: {len(atoms)} atoms, {len(charges)} charges")
    for atom, q in zip(atoms, charges):
        atom.charge = q
    return atoms


def autodock_types(atoms: Sequence[Atom]) -> List[Tuple[Atom, str]]:
    """AutoDock atom types, polar hydrogens kept, non-polar hydrogens merged.

    Returns the atoms to write, each with its type; merged hydrogen charges
    are added to their parent atom's charge.
    """
    from scipy.spatial import cKDTree

    heavy_idx = [i for i, a in enumerate(atoms) if a.element != "H"]
    h_idx = [i for i, a in enumerate(atoms) if a.element == "H"]
    parent: Dict[int, int] = {}
    if heavy_idx and h_idx:
        tree = cKDTree(np.array([atoms[i].xyz for i in heavy_idx]))
        dist, pos = tree.query(np.array([atoms[i].xyz for i in h_idx]))
        for h, d, p in zip(h_idx, dist, pos):
            if d <= H_BOND_MAX:
                parent[h] = heavy_idx[int(p)]
    has_h = {p for p in parent.values()}
    charge = {i: atoms[i].charge for i in range(len(atoms))}
    keep_h = set()
    for h, p in parent.items():
        if atoms[p].element in ("N", "O"):
            keep_h.add(h)
        else:
            charge[p] += charge[h]
    out: List[Tuple[Atom, str]] = []
    for i, a in enumerate(atoms):
        if a.element == "H":
            if i not in keep_h:
                continue
            t = "HD"
        elif a.element == "C":
            aromatic = (a.name in AROMATIC_CARBONS.get(a.resname, ()) or
                        (a.resname in HISTIDINES and a.name in HIS_RING_CARBONS) or
                        (a.resname in NUCLEOTIDES and a.name in BASE_RING_CARBONS))
            t = "A" if aromatic else "C"
        elif a.element == "N":
            ring_n = ((a.resname in HISTIDINES and a.name in HIS_RING_NITROGENS) or
                      (a.resname in NUCLEOTIDES and a.name in BASE_RING_NITROGENS
                       and a.name != _glycosidic_nitrogen(a.resname)))
            t = "NA" if ring_n and i not in has_h else "N"
        elif a.element == "O":
            t = "OA"
        elif a.element in ("S", "SE"):
            t = "SA"
        elif a.element == "P":
            t = "P"
        elif a.element in METAL_TYPES:
            t = METAL_TYPES[a.element]
        else:
            raise ReceptorError(f"no AutoDock type for element {a.element!r} ({a.resname} {a.name})")
        typed = Atom(a.record, a.name, a.resname, a.chain, a.resseq, a.icode, a.xyz, a.element, charge[i])
        out.append((typed, t))
    return out


def pdbqt_line(serial: int, a: Atom, atype: str) -> str:
    record = "HETATM" if a.record == "HETATM" else "ATOM"
    name = a.name[:4]
    name = name if len(name) == 4 else f" {name:<3s}"
    return (f"{record:<6s}{serial % 100000:5d} {name:<4s} {a.resname[:3]:>3s} {a.chain[:1]:1s}{a.resseq % 10000:4d}"
            f"{a.icode[:1]:1s}   {a.xyz[0]:8.3f}{a.xyz[1]:8.3f}{a.xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}"
            f"    {a.charge:6.3f} {atype:<2s}")


def write_pdbqt(typed: Sequence[Tuple[Atom, str]], out_path: Path, extra: Sequence[Atom] = ()) -> Dict[str, int]:
    """Rigid receptor PDBQT; ``extra`` atoms (kept metals) are appended as HETATM."""
    lines = []
    counts: Dict[str, int] = {}
    for serial, (atom, t) in enumerate(typed, start=1):
        lines.append(pdbqt_line(serial, atom, t))
        counts[t] = counts.get(t, 0) + 1
    for k, atom in enumerate(extra, start=len(typed) + 1):
        t = METAL_TYPES[atom.element]
        lines.append(pdbqt_line(k, atom, t))
        counts[t] = counts.get(t, 0) + 1
    out_path.write_text("\n".join(lines) + "\n")
    return counts


@dataclass
class PreparedReceptor:
    pdbqt: Path
    types: Dict[str, int]
    strip: StripReport
    his_states: Dict[str, int]
    metals_kept: int = 0


def prepare_receptor(arrays, workdir: Path, *, keep_metals: Sequence[Atom] = (), ph: float = PH,
                     protonated: Optional[Tuple[Path, Path]] = None, name: str = "receptor") -> PreparedReceptor:
    """Strip, protonate and type one structure; ``protonated`` reuses an earlier PDB2PQR run."""
    workdir.mkdir(parents=True, exist_ok=True)
    strip = StripReport()
    if protonated is None:
        strip = write_polymer_pdb(arrays, workdir / "polymer.pdb")
        if strip.polymer_atoms == 0:
            raise ReceptorError("no polymer atoms")
        protonated = protonate(workdir / "polymer.pdb", workdir, ph=ph)
    atoms = read_protonated(*protonated)
    his: Dict[str, int] = {}
    for a in atoms:
        if a.resname in HISTIDINES and a.name == "CA":
            his[a.resname] = his.get(a.resname, 0) + 1
    typed = autodock_types(atoms)
    out = workdir / f"{name}.pdbqt"
    types = write_pdbqt(typed, out, extra=keep_metals)
    return PreparedReceptor(out, types, strip, his, len(keep_metals))
