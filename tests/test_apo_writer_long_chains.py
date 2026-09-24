"""The apo copy of an mmCIF entry with multi-character chain identifiers."""

from __future__ import annotations

import copy

import numpy as np
from Bio.PDB import MMCIFIO, PDBParser

from cryptic_ip.analysis.structure_arrays import load_structure_arrays, write_apo_structure
from cryptic_ip.testing.synthetic import default_benchmark_specs, write_synthetic_structure


def test_long_chain_ids_are_renamed_without_losing_or_merging_chains(tmp_path):
    pdb = write_synthetic_structure(default_benchmark_specs(n_buried=1, n_surface=0, n_decoy=0, seed=0)[0], tmp_path)
    structure = PDBParser(QUIET=True).get_structure("x", str(pdb))
    model = structure[0]
    original = list(model)[0]
    twin = copy.deepcopy(original)
    for atom in twin.get_atoms():
        atom.coord = atom.coord + 60.0  # a second copy, far away
    original.id, twin.id = "AA", "AB"
    model.add(twin)
    cif = tmp_path / "long.cif"
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif))

    source = load_structure_arrays(cif)
    apo = load_structure_arrays(write_apo_structure(cif, tmp_path / "apo.pdb"))
    assert len(set(apo.chain_ids)) == 2
    assert apo.n_atoms == int(source.is_polymer.sum())
    # Coordinates survive unchanged, so coordinate-based labels still line up.
    np.testing.assert_allclose(
        np.sort(apo.coords, axis=0), np.sort(source.coords[source.is_polymer], axis=0), atol=1e-3
    )
