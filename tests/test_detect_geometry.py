"""Ideal-geometry fixtures, bulk-solvent burial, and side-chain rotamers.

These three belong together because the first exists for the third. The
project's existing synthetic structures place atoms on a lattice - their own
header says NOT A REAL STRUCTURE - and in them a CA-CB "bond" can be 14.7
angstrom long. Rotating a side chain about such an axis flings its atoms
through the protein, every conformer is rejected as a clash, and the result
looks exactly like a tightly packed pocket rather than like invalid input.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import cavity as CV
from cryptic_ip.detect import geometry_fixtures as GF
from cryptic_ip.detect import rotamer as RO


# ------------------------------------------------------- ideal geometry

def test_atom_placement_reproduces_the_requested_internal_coordinates():
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([1.5, 0.0, 0.0])
    c = np.array([2.0, 1.4, 0.0])
    d = GF.place_atom(a, b, c, bond=1.52, angle=114.0, torsion=-60.0)
    assert np.linalg.norm(d - c) == pytest.approx(1.52, abs=1e-9)
    # The B-C-D angle must be what was asked for.
    v1, v2 = b - c, d - c
    cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    assert np.rad2deg(np.arccos(cos)) == pytest.approx(114.0, abs=1e-6)


def test_atom_placement_reproduces_the_requested_torsion():
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([1.5, 0.0, 0.0])
    c = np.array([2.0, 1.4, 0.0])
    for wanted in (-180.0, -60.0, 0.0, 60.0, 120.0):
        d = GF.place_atom(a, b, c, 1.52, 114.0, wanted)
        assert _torsion(a, b, c, d) == pytest.approx(wanted, abs=1e-6)


def _torsion(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1, v), w))
    return float(np.degrees(np.arctan2(y, x)))


def test_three_collinear_atoms_give_no_frame_and_are_refused():
    with pytest.raises(GF.FixtureError, match="collinear"):
        GF.place_atom(np.zeros(3), np.array([1.0, 0, 0]), np.array([2.0, 0, 0]), 1.5, 110.0, 60.0)


def test_a_built_backbone_has_ideal_bond_lengths():
    residues = GF.build_backbone(20)
    GF.add_carbonyl(residues)
    for residue in residues:
        GF.add_side_chain(residue, "ALA")
    report = GF.bond_length_report(GF.assemble(residues))
    assert report["N_CA_mean"] == pytest.approx(GF.N_CA, abs=1e-3)
    assert report["CA_C_mean"] == pytest.approx(GF.CA_C, abs=1e-3)
    assert report["CA_CB_mean"] == pytest.approx(GF.CA_CB, abs=1e-3)
    # Every residue identical, so the deviation is numerical only.
    assert report["CA_CB_max_deviation"] < 1e-6


def test_a_helix_has_the_rise_a_helix_has():
    residues = GF.build_backbone(20, phi_psi=(GF.HELIX,))
    ca = np.array([r.atom("CA") for r in residues])
    # 1.5 angstrom per residue along the axis, so 19 residues span about 28.5.
    span = float(np.linalg.norm(ca[-1] - ca[0]))
    assert 24.0 < span < 32.0


def test_a_strand_is_more_extended_than_a_helix():
    helix = GF.build_backbone(20, phi_psi=(GF.HELIX,))
    strand = GF.build_backbone(20, phi_psi=(GF.STRAND,))
    extent = lambda rs: float(np.linalg.norm(rs[-1].atom("CA") - rs[0].atom("CA")))  # noqa: E731
    assert extent(strand) > 1.8 * extent(helix)


def test_too_short_a_backbone_is_refused():
    with pytest.raises(GF.FixtureError, match="at least 3"):
        GF.build_backbone(2)


def test_an_unknown_residue_is_refused_rather_than_approximated():
    residues = GF.build_backbone(5)
    with pytest.raises(GF.FixtureError, match="no build plan"):
        GF.add_side_chain(residues[0], "XYZ")


def test_the_barrel_has_valid_geometry_and_a_channel():
    fixture = GF.helix_barrel(n_helices=6, residues_per_helix=18, barrel_radius=9.0)
    report = GF.bond_length_report(fixture)
    assert report["CA_CB_mean"] == pytest.approx(GF.CA_CB, abs=1e-3)
    found = CV.detect_cavities(fixture.coords, [str(e) for e in fixture.elements], spacing=0.8)
    assert found.cavities, "a helix barrel must have a channel down its axis"


def test_the_barrel_places_the_named_lining_residues():
    fixture = GF.helix_barrel(n_helices=4, residues_per_helix=14, lining={0: "TRP", 2: "ARG"})
    names = set(str(n) for n in fixture.resnames)
    assert {"TRP", "ARG", "ALA"} <= names


def test_too_few_helices_is_refused():
    with pytest.raises(GF.FixtureError, match="at least 3"):
        GF.helix_barrel(n_helices=2)


# --------------------------------------------------- bulk-solvent burial

def _shell(n=400, radius=9.0, cut=None):
    i = np.arange(n) + 0.5
    polar = np.arccos(1.0 - 2.0 * i / n)
    azimuth = np.pi * (1.0 + 5.0 ** 0.5) * i
    points = radius * np.column_stack([
        np.sin(polar) * np.cos(azimuth), np.sin(polar) * np.sin(azimuth), np.cos(polar),
    ])
    return points if cut is None else points[points[:, 2] < cut]


def _slotted_block(width, depth=10.0, step=2.0, half=12.0):
    grid = np.arange(-half, half + 0.1, step)
    return np.array([
        [x, y, z] for x in grid for y in grid for z in grid
        if not (abs(x) < width / 2 and z > half - depth and abs(y) < 8.0)
    ], dtype=float)


def test_a_closed_shell_is_fully_sealed_with_no_mouth():
    atoms = _shell()
    cavity = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8).largest()
    assert cavity.sealed_fraction == pytest.approx(1.0)
    assert cavity.mouth_area == 0.0
    assert cavity.is_buried and not cavity.reaches_bulk


def test_opening_a_hole_in_a_shell_makes_it_reachable():
    atoms = _shell(cut=7.5)
    cavity = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8).largest()
    assert cavity.sealed_fraction < 0.2
    assert cavity.mouth_area > 0.0
    assert not cavity.is_buried and cavity.reaches_bulk


def test_a_wide_slot_is_an_open_groove_not_a_buried_cavity():
    # The distinction the enclosure count cannot make: a groove is enclosed in
    # five of seven directions along most of its length and reads as buried.
    atoms = _slotted_block(12.0)
    cavity = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8).largest()
    assert cavity.mean_enclosure >= 5.0, "the groove really is enclosed by the directional count"
    assert not cavity.is_buried, "but it is not buried, because solvent reaches it"


def test_a_throat_narrower_than_a_water_molecule_still_counts_as_sealed():
    # A 1.8 angstrom clearance is impassable to water, so the cavity behind it
    # is buried however much it looks like a slot on paper.
    atoms = _slotted_block(6.0)
    cavity = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8).largest()
    assert cavity.sealed_fraction == pytest.approx(1.0)


def test_widening_the_slot_monotonically_reduces_the_sealed_fraction():
    sealed = []
    for width in (6.0, 10.0, 14.0):
        atoms = _slotted_block(width)
        found = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8)
        sealed.append(found.largest().sealed_fraction if found.cavities else 0.0)
    assert sealed[0] > sealed[1] >= sealed[2]


def test_the_mouth_area_grows_as_the_opening_widens():
    areas = []
    for width in (10.0, 14.0):
        atoms = _slotted_block(width)
        areas.append(CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8).largest().mouth_area)
    assert areas[1] > areas[0]


def test_sealed_volume_is_the_sealed_part_of_the_volume():
    atoms = _slotted_block(12.0)
    cavity = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8).largest()
    assert cavity.sealed_volume == pytest.approx(cavity.volume * cavity.sealed_fraction)


def test_the_map_reports_buried_cavities_separately():
    atoms = _shell()
    found = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8)
    assert found.buried_volume > 0.0
    assert found.largest_buried() is found.largest()
    assert found.as_dict()["n_buried"] == 1


def test_raw_connectivity_without_a_probe_finds_almost_nothing_buried():
    # With probe_steps=0 a single-voxel gap is a channel, which is why the
    # default erodes by about a water's radius first.
    atoms = _shell(cut=8.2)
    strict = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8, probe_steps=0)
    lenient = CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=0.8, probe_steps=3)
    strict_sealed = strict.largest().sealed_fraction if strict.cavities else 0.0
    lenient_sealed = lenient.largest().sealed_fraction if lenient.cavities else 0.0
    assert lenient_sealed >= strict_sealed


# ------------------------------------------------------------- rotamers

def _barrel(radius=7.0, lining=None):
    return GF.helix_barrel(
        n_helices=6, residues_per_helix=18, barrel_radius=radius,
        sequence="ALA", lining=lining if lining is not None else {0: "TRP", 2: "TRP", 4: "TRP"},
    )


def test_rotatable_torsions_are_found_on_valid_geometry():
    fixture = _barrel()
    chains = RO.rotatable_side_chains(fixture)
    kinds = {(c.resname, c.torsion) for c in chains}
    assert ("TRP", "chi1") in kinds and ("TRP", "chi2") in kinds


def test_alanine_and_glycine_have_no_rotatable_chi():
    fixture = GF.helix_barrel(n_helices=4, residues_per_helix=12, sequence="ALA")
    assert RO.rotatable_side_chains(fixture) == []


def test_an_unphysical_torsion_axis_is_refused_with_its_length():
    # The lattice fixtures produce axes up to 14.7 angstrom long. Rotating
    # about them produced contacts of 0.17 angstrom - two atoms superimposed -
    # and every conformer was silently rejected as a clash.
    fixture = _barrel()
    broken = np.array(fixture.coords, copy=True)
    cb = (np.asarray(fixture.atom_names, dtype=object) == "CB")
    broken[cb] += np.array([12.0, 0.0, 0.0])

    class _Bad:
        coords = broken
        atom_names = fixture.atom_names
        resnames = fixture.resnames
        residue_index = fixture.residue_index
        is_polymer = fixture.is_polymer

    with pytest.raises(RO.RotamerError, match="no covalent geometry"):
        RO.rotatable_side_chains(_Bad())


def test_rotation_preserves_every_bond_length_in_the_side_chain():
    fixture = _barrel()
    chain = next(c for c in RO.rotatable_side_chains(fixture) if c.torsion == "chi1")
    origin = fixture.coords[chain.axis_from]
    axis = fixture.coords[chain.axis_to] - origin
    moved = RO.rotate_about_axis(fixture.coords, origin, axis, chain.moving, 73.0)
    # Distances within the moved set and from the axis atoms are rigid.
    for i in chain.moving:
        assert np.linalg.norm(moved[i] - moved[chain.axis_to]) == pytest.approx(
            float(np.linalg.norm(fixture.coords[i] - fixture.coords[chain.axis_to])), abs=1e-9,
        )


def test_rotation_leaves_the_backbone_untouched():
    fixture = _barrel()
    chain = next(c for c in RO.rotatable_side_chains(fixture) if c.torsion == "chi1")
    origin = fixture.coords[chain.axis_from]
    axis = fixture.coords[chain.axis_to] - origin
    moved = RO.rotate_about_axis(fixture.coords, origin, axis, chain.moving, 120.0)
    backbone = np.asarray([
        i for i, name in enumerate(fixture.atom_names) if str(name) in ("N", "CA", "C", "O")
    ])
    assert np.allclose(moved[backbone], fixture.coords[backbone])


def test_rotating_by_the_requested_angle_changes_the_torsion_by_that_much():
    fixture = _barrel()
    chain = next(c for c in RO.rotatable_side_chains(fixture) if c.torsion == "chi1")
    mask = fixture.residue_index == chain.residue
    names = {
        str(name): index
        for index, name in zip(np.nonzero(mask)[0], fixture.atom_names[mask])
    }
    n, ca, cb, cg = names["N"], names["CA"], names["CB"], names["CG"]
    before = _torsion(*(fixture.coords[i] for i in (n, ca, cb, cg)))
    origin = fixture.coords[chain.axis_from]
    axis = fixture.coords[chain.axis_to] - origin
    moved = RO.rotate_about_axis(fixture.coords, origin, axis, chain.moving, 60.0)
    after = _torsion(*(moved[i] for i in (n, ca, cb, cg)))
    delta = (after - before + 540.0) % 360.0 - 180.0
    assert abs(delta) == pytest.approx(60.0, abs=1e-6)


def test_a_zero_length_axis_is_refused():
    with pytest.raises(RO.RotamerError, match="zero length"):
        RO.rotate_about_axis(np.zeros((3, 3)), np.zeros(3), np.zeros(3), [0], 60.0)


def test_clash_screening_is_relative_to_the_input_not_absolute():
    # A predicted model routinely carries marginal contacts. Against a fixed
    # cutoff every rotation of such a side chain is a clash, including the one
    # that opens the pocket, and the ensemble comes back empty with no reason.
    coords = np.array([[0.0, 0, 0], [2.0, 0, 0]])
    residue_index = np.array([0, 5])
    # The input already has the two atoms 2.0 angstrom apart, inside the 2.9 cutoff.
    assert RO.clashes(coords, [0], residue_index, 0, reference=None) is True
    assert RO.clashes(coords, [0], residue_index, 0, reference=2.0) is False


def test_a_rotation_that_makes_contact_worse_is_still_rejected():
    residue_index = np.array([0, 5])
    # The input had the two atoms 2.0 A apart; this brings them to 1.0.
    closer = np.array([[1.0, 0, 0], [2.0, 0, 0]])
    assert RO.clashes(closer, [0], residue_index, 0, reference=2.0) is True


def test_sequence_neighbours_are_exempt_from_the_clash_test():
    coords = np.array([[0.0, 0, 0], [1.5, 0, 0]])
    assert RO.closest_contact(coords, [0], np.array([3, 4]), 3) == float("inf")


def test_the_ensemble_records_which_inputs_were_already_strained():
    fixture = _barrel(radius=6.5)
    ensemble = RO.rotamer_conformers(fixture, RO.rotatable_side_chains(fixture))
    payload = ensemble.as_dict()
    assert "strained_inputs" in payload
    assert payload["n_accepted"] + payload["n_rejected"] > 0


def test_rotamers_open_a_channel_a_side_chain_occludes():
    # The capability the backbone ensemble cannot have: elastic-network modes
    # move atoms rigidly with their residue, so no amount of low-frequency
    # motion rotates a tryptophan out of a channel.
    fixture = _barrel(radius=7.0)
    elements = [str(e) for e in fixture.elements]
    apo = CV.detect_cavities(fixture.coords, elements, spacing=0.8)
    apo_volume = apo.largest().volume if apo.cavities else 0.0
    ensemble = RO.rotamer_conformers(fixture, RO.rotatable_side_chains(fixture))
    assert ensemble.n_accepted > 0
    best = apo_volume
    for conformer in ensemble.conformers:
        found = CV.detect_cavities(conformer, elements, spacing=0.8)
        best = max(best, found.largest().volume if found.cavities else 0.0)
    assert best > apo_volume, "no rotamer opened the channel at all"


def test_a_tighter_channel_gives_a_larger_rotamer_effect():
    # The physical expectation: the tighter the channel, the more one side
    # chain decides. A method whose effect did not scale that way would be
    # measuring something else.
    def gain(radius):
        fixture = _barrel(radius=radius)
        elements = [str(e) for e in fixture.elements]
        apo = CV.detect_cavities(fixture.coords, elements, spacing=0.8)
        apo_volume = apo.largest().volume if apo.cavities else 0.0
        ensemble = RO.rotamer_conformers(fixture, RO.rotatable_side_chains(fixture))
        best = apo_volume
        for conformer in ensemble.conformers:
            found = CV.detect_cavities(conformer, elements, spacing=0.8)
            best = max(best, found.largest().volume if found.cavities else 0.0)
        return 0.0 if best == 0 else 1.0 - apo_volume / best

    assert gain(7.5) > gain(9.0)


def test_lining_residues_does_not_require_a_cavity_to_exist():
    # A side chain that occludes a cryptic site is not lining a cavity,
    # because there is no cavity yet.
    fixture = _barrel(radius=7.0)
    residues = RO.lining_residues(fixture, (0.0, 0.0, 0.0), radius=10.0)
    assert residues
