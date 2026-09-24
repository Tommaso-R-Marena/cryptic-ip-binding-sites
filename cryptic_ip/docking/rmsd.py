"""Symmetry-corrected ligand RMSD without superposition (docs/REDOCKING_PLAN.md).

A docked pose is compared with the crystal pose **in the receptor frame**: no
alignment is done, because aligning the two ligands would remove exactly the
placement error docking is judged on.

Inositol phosphates are highly symmetric, and a naive atom-order RMSD is wrong
for them in two ways:

* the heavy-atom graph has automorphisms (the ring and its six phosphates can be
  relabelled around the ring), so the same pose written with atoms in another
  order would score > 0 Å;
* each phosphate's terminal oxygens are equivalent (P=O and P-O(-) are one
  resonance structure), so a pose that differs only by which terminal oxygen is
  called O1P scores > 0 Å.

:func:`symmetric_rmsd` handles both exactly. The **core** (every heavy atom
except the terminal oxygens on phosphorus) is matched by all graph
automorphisms, found by RDKit substructure matching on an element-only graph
(bond orders, charges and hydrogens removed, so protonation states cannot break
a match). For each core mapping, the terminal oxygens of each phosphorus are
assigned by exhaustive permutation. The squared deviations of different
phosphorus atoms' oxygens are independent, so assigning each set separately is
exact. RDKit's own ``symmetrizeConjugatedTerminalGroups`` is not relied on:
whether it covers P-O depends on the RDKit version, and a test checks this
implementation against brute-force enumeration.

The reference may be **incomplete** (a crystal copy with unmodelled atoms): it is
matched as a substructure of the probe, and the RMSD runs over the atoms it has.
"""

from __future__ import annotations

from itertools import permutations
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

#: Largest number of core automorphisms enumerated. An inositol phosphate core
#: has 12 (the dihedral group of the ring); this bound only guards against a
#: pathological input.
MAX_CORE_MATCHES = 100_000


def _element_graph(mol):
    """Heavy atoms only, all bonds single, no charges, no hydrogens, no stereo."""
    from rdkit import Chem

    heavy = Chem.RemoveHs(mol, sanitize=False)
    rw = Chem.RWMol()
    index = {}
    for atom in heavy.GetAtoms():
        new = Chem.Atom(atom.GetAtomicNum())
        new.SetNoImplicit(True)
        index[atom.GetIdx()] = rw.AddAtom(new)
    for bond in heavy.GetBonds():
        rw.AddBond(index[bond.GetBeginAtomIdx()], index[bond.GetEndAtomIdx()], Chem.BondType.SINGLE)
    graph = rw.GetMol()
    graph.UpdatePropertyCache(strict=False)
    return graph, heavy


def terminal_phosphate_oxygens(graph) -> Dict[int, List[int]]:
    """Phosphorus atom -> its terminal oxygens (oxygens bonded to nothing else).

    Works on any molecule graph; hydrogens are ignored, so a protonated P-OH
    oxygen is terminal too.
    """
    out: Dict[int, List[int]] = {}
    for atom in graph.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        terminal = [
            n.GetIdx() for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 8 and sum(1 for x in n.GetNeighbors() if x.GetAtomicNum() > 1) == 1
        ]
        if terminal:
            out[atom.GetIdx()] = sorted(terminal)
    return out


def _core(graph, terminal: Dict[int, List[int]]):
    """The graph without terminal phosphate oxygens, with maps between indices."""
    from rdkit import Chem

    drop = {o for oxygens in terminal.values() for o in oxygens}
    keep = [a.GetIdx() for a in graph.GetAtoms() if a.GetIdx() not in drop]
    rw = Chem.RWMol(graph)
    for idx in sorted(drop, reverse=True):
        rw.RemoveAtom(idx)
    core = rw.GetMol()
    core.UpdatePropertyCache(strict=False)
    return core, keep  # core index i corresponds to graph index keep[i]


def heavy_coordinates(mol, conf_id: int = -1) -> Tuple[object, np.ndarray]:
    """The element graph of ``mol`` and its heavy-atom coordinates, in graph order."""
    graph, heavy = _element_graph(mol)
    coords = np.asarray(heavy.GetConformer(conf_id).GetPositions(), dtype=float)
    return graph, coords


def atom_mappings(reference_graph, probe_graph) -> List[Tuple[Dict[int, int], Dict[int, List[int]], Dict[int, List[int]]]]:
    """Every core mapping of reference onto probe, with each side's terminal oxygens.

    Returns a list of ``(core_map, ref_terminal, probe_terminal)`` where
    ``core_map`` maps reference graph indices to probe graph indices.
    """
    ref_terminal = terminal_phosphate_oxygens(reference_graph)
    probe_terminal = terminal_phosphate_oxygens(probe_graph)
    ref_core, ref_keep = _core(reference_graph, ref_terminal)
    probe_core, probe_keep = _core(probe_graph, probe_terminal)
    matches = probe_core.GetSubstructMatches(
        ref_core, uniquify=False, useChirality=False, maxMatches=MAX_CORE_MATCHES
    )
    out = []
    for match in matches:
        core_map = {ref_keep[i]: probe_keep[j] for i, j in enumerate(match)}
        out.append((core_map, ref_terminal, probe_terminal))
    return out


def _best_terminal_assignment(ref_xyz: np.ndarray, probe_xyz: np.ndarray) -> float:
    """Least total squared deviation over injective assignments of ref oxygens to probe oxygens."""
    n = len(ref_xyz)
    if n == 0:
        return 0.0
    best = np.inf
    for chosen in permutations(range(len(probe_xyz)), n):
        d = float(np.sum((ref_xyz - probe_xyz[list(chosen)]) ** 2))
        if d < best:
            best = d
    return best


def symmetric_rmsd(
    reference,
    probe,
    *,
    reference_conf: int = -1,
    probe_conf: int = -1,
    phosphorus_only: bool = False,
) -> float:
    """Heavy-atom RMSD of ``probe`` against ``reference``, no superposition, symmetry-corrected.

    Args:
        reference: RDKit molecule (e.g. the crystal copy); may lack atoms.
        probe: RDKit molecule of the same compound (e.g. a docked pose).
        reference_conf, probe_conf: Conformer ids.
        phosphorus_only: RMSD over phosphorus atoms only (still under the
            graph automorphisms, so it is not fooled by relabelling).

    Raises:
        ValueError: if the reference graph is not a substructure of the probe's.
    """
    ref_graph, ref_xyz = heavy_coordinates(reference, reference_conf)
    probe_graph, probe_xyz = heavy_coordinates(probe, probe_conf)
    mappings = atom_mappings(ref_graph, probe_graph)
    if not mappings:
        raise ValueError("reference is not a substructure of the probe: different compounds")
    phosphorus = [a.GetIdx() for a in ref_graph.GetAtoms() if a.GetAtomicNum() == 15]
    best = np.inf
    n_atoms = 0
    for core_map, ref_terminal, probe_terminal in mappings:
        if phosphorus_only:
            pairs = [(i, core_map[i]) for i in phosphorus]
            total = float(sum(np.sum((ref_xyz[i] - probe_xyz[j]) ** 2) for i, j in pairs))
            n_atoms = len(pairs)
        else:
            total = float(sum(np.sum((ref_xyz[i] - probe_xyz[j]) ** 2) for i, j in core_map.items()))
            n_atoms = len(core_map)
            feasible = True
            for p_ref, oxygens in ref_terminal.items():
                candidates = probe_terminal.get(core_map[p_ref], [])
                if len(candidates) < len(oxygens):
                    feasible = False
                    break
                total += _best_terminal_assignment(ref_xyz[oxygens], probe_xyz[candidates])
                n_atoms += len(oxygens)
            if not feasible:
                continue
        if total < best:
            best = total
    if not np.isfinite(best) or n_atoms == 0:
        raise ValueError("no feasible atom mapping")
    return float(np.sqrt(best / n_atoms))


def naive_rmsd(reference, probe, *, reference_conf: int = -1, probe_conf: int = -1) -> float:
    """Atom-order heavy-atom RMSD, no symmetry and no superposition (for tests and contrast only)."""
    _, a = heavy_coordinates(reference, reference_conf)
    _, b = heavy_coordinates(probe, probe_conf)
    if a.shape != b.shape:
        raise ValueError("different atom counts")
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def rmsd_matrix(poses: Sequence, conf_ids: Optional[Sequence[int]] = None) -> np.ndarray:
    """Pairwise symmetric RMSD between poses (list of molecules, or one molecule's conformers)."""
    if conf_ids is not None:
        items = [(poses[0], c) for c in conf_ids]
    else:
        items = [(m, -1) for m in poses]
    n = len(items)
    out = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            out[i, j] = out[j, i] = symmetric_rmsd(items[i][0], items[j][0],
                                                   reference_conf=items[i][1], probe_conf=items[j][1])
    return out
