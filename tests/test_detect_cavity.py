"""The in-process cavity detector and the conformer-ensemble consensus.

The detector is tested against geometries whose answer is known by
construction - a closed shell has one cavity, a flat slab has none - because
there is no reference implementation here to compare against and fpocket's
volumes are not the same quantity.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import cavity as CV
from cryptic_ip.detect import ensemble as EN


def sphere_shell(n=300, radius=9.0):
    """A closed shell of atoms: exactly one cavity, at the origin."""
    i = np.arange(n) + 0.5
    polar = np.arccos(1.0 - 2.0 * i / n)
    azimuth = np.pi * (1.0 + 5.0 ** 0.5) * i
    return radius * np.column_stack([
        np.sin(polar) * np.cos(azimuth),
        np.sin(polar) * np.sin(azimuth),
        np.cos(polar),
    ])


def flat_slab(step=2.0, half=10.0):
    grid = np.arange(-half, half + step, step)
    return np.stack(np.meshgrid(grid, grid, [0.0]), -1).reshape(-1, 3)


def slot(gap, step=2.0, half=8.0):
    """Two plates `gap` apart, walled on three sides."""
    grid = np.arange(-half, half + step, step)
    heights = np.arange(-gap / 2, gap / 2 + 0.1, step)
    top = np.stack(np.meshgrid(grid, grid, [gap / 2]), -1).reshape(-1, 3)
    bottom = np.stack(np.meshgrid(grid, grid, [-gap / 2]), -1).reshape(-1, 3)
    walls = [[-half, y, z] for y in grid for z in heights]
    walls += [[x, s * half, z] for x in grid for s in (-1, 1) for z in heights]
    return np.vstack([top, bottom, np.asarray(walls, dtype=float)])


# ------------------------------------------------------------------- the grid

def test_the_occupancy_grid_marks_the_inside_of_an_atom():
    occupied, origin = CV.occupancy_grid(np.zeros((1, 3)), np.array([1.7]), spacing=0.5)
    assert occupied.any()
    # The marked region is about a sphere of radius 1.7 + 1.4.
    volume = float(occupied.sum()) * 0.5 ** 3
    expected = 4.0 / 3.0 * np.pi * (1.7 + CV.PROBE_RADIUS) ** 3
    assert 0.8 * expected < volume < 1.2 * expected


def test_an_empty_structure_is_refused():
    with pytest.raises(CV.CavityError, match="no atoms"):
        CV.occupancy_grid(np.zeros((0, 3)), np.zeros(0))


def test_a_nonpositive_spacing_is_refused():
    with pytest.raises(CV.CavityError, match="spacing"):
        CV.occupancy_grid(np.zeros((1, 3)), np.array([1.7]), spacing=0.0)


def test_an_absurd_grid_is_refused_rather_than_allocated():
    far = np.array([[0.0, 0, 0], [4000.0, 4000.0, 4000.0]])
    with pytest.raises(CV.CavityError, match="million points"):
        CV.occupancy_grid(far, np.array([1.7, 1.7]), spacing=0.5)


def test_radii_fall_back_to_carbon_for_an_unknown_element():
    radii = CV.atom_radii(["C", "N", "XX"])
    assert radii[0] == CV.VDW_RADII["C"]
    assert radii[2] == CV.DEFAULT_RADIUS


def test_a_shift_larger_than_the_grid_is_all_false_not_an_error():
    # A structure thinner than the scan reach produces exactly this, and the
    # slice pair became inconsistent before the guard existed - numpy raised
    # on a flat slab.
    grid = np.ones((3, 3, 3), dtype=bool)
    assert not CV._shifted(grid, (5, 0, 0)).any()


# --------------------------------------------------------------- the cavities

def test_a_closed_shell_has_exactly_one_fully_enclosed_cavity():
    shell = sphere_shell()
    found = CV.detect_cavities(shell, ["C"] * len(shell), spacing=0.8)
    assert len(found.cavities) == 1
    cavity = found.largest()
    assert np.allclose(cavity.centre, 0.0, atol=0.5)
    assert cavity.max_enclosure == len(CV.SCAN_DIRECTIONS)
    assert cavity.volume > CV.MIN_CAVITY_VOLUME


def test_a_flat_slab_has_no_cavity():
    slab = flat_slab()
    assert CV.detect_cavities(slab, ["C"] * len(slab), spacing=0.8).cavities == ()


def test_cutting_a_hole_in_a_shell_lowers_its_enclosure():
    shell = sphere_shell()
    closed = CV.detect_cavities(shell, ["C"] * len(shell), spacing=0.8).largest()
    holed_atoms = shell[shell[:, 2] < 6.0]
    holed = CV.detect_cavities(holed_atoms, ["C"] * len(holed_atoms), spacing=0.8).largest()
    assert holed is not None
    assert holed.mean_enclosure < closed.mean_enclosure


def test_a_narrow_slot_holds_no_cavity_and_a_wide_one_does():
    narrow = slot(5.0)
    wide = slot(9.0)
    assert CV.detect_cavities(narrow, ["C"] * len(narrow), spacing=1.0).cavities == ()
    assert CV.detect_cavities(wide, ["C"] * len(wide), spacing=1.0).cavities != ()


def test_requiring_more_enclosure_finds_fewer_cavities():
    shell = sphere_shell()
    loose = CV.detect_cavities(shell, ["C"] * len(shell), spacing=1.0, min_enclosure=3)
    tight = CV.detect_cavities(shell, ["C"] * len(shell), spacing=1.0, min_enclosure=7)
    assert loose.total_volume >= tight.total_volume


def test_cavities_come_back_largest_first():
    shell = sphere_shell()
    found = CV.detect_cavities(shell, ["C"] * len(shell), spacing=1.0, min_enclosure=3)
    volumes = [c.volume for c in found.cavities]
    assert volumes == sorted(volumes, reverse=True)


def test_a_map_with_no_cavities_reports_zero_rather_than_failing():
    slab = flat_slab()
    found = CV.detect_cavities(slab, ["C"] * len(slab), spacing=1.0)
    assert found.total_volume == 0.0 and found.largest() is None
    assert found.as_dict()["n_cavities"] == 0


def test_the_detector_is_deterministic():
    shell = sphere_shell()
    a = CV.detect_cavities(shell, ["C"] * len(shell), spacing=1.0)
    b = CV.detect_cavities(shell, ["C"] * len(shell), spacing=1.0)
    assert [c.as_dict() for c in a.cavities] == [c.as_dict() for c in b.cavities]


# ---------------------------------------------------------------- the tracks

def _maps(gaps):
    out = []
    for gap in gaps:
        atoms = slot(gap)
        out.append(CV.detect_cavities(atoms, ["C"] * len(atoms), spacing=1.0))
    return out


def test_a_site_absent_in_the_apo_and_open_later_is_fully_cryptic():
    # This is the quantity the project has never had: the apo structure has no
    # cavity at all, so measuring its burial would find nothing.
    maps = _maps([5.0, 5.0, 7.0, 9.0, 7.0])
    report = EN.analyse_ensemble(maps, tolerance=6.0)
    assert len(report.tracks) == 1
    track = report.tracks[0]
    assert track.apo_volume == 0.0
    assert track.crypticity == pytest.approx(1.0)
    assert track.classification == "cryptic"
    assert 0.0 < track.persistence < 1.0


def test_a_site_open_in_every_conformer_is_not_cryptic():
    maps = _maps([9.0, 9.0, 9.0, 9.0])
    track = EN.analyse_ensemble(maps, tolerance=6.0).tracks[0]
    assert track.persistence == pytest.approx(1.0)
    assert track.crypticity < 0.25
    assert track.classification == "open"


def test_crypticity_is_the_fraction_of_the_best_volume_the_apo_lacks():
    maps = _maps([7.0, 9.0])
    track = EN.analyse_ensemble(maps, tolerance=6.0).tracks[0]
    assert track.crypticity == pytest.approx(1.0 - track.apo_volume / track.max_volume)


def test_the_apo_conformer_is_the_first_one():
    # Every crypticity is measured against it, so the ordering is part of the
    # contract rather than a convention.
    open_first = EN.analyse_ensemble(_maps([9.0, 5.0, 5.0]), tolerance=6.0).tracks[0]
    closed_first = EN.analyse_ensemble(_maps([5.0, 9.0, 9.0]), tolerance=6.0).tracks[0]
    assert open_first.crypticity == pytest.approx(0.0)
    assert closed_first.crypticity == pytest.approx(1.0)


def test_two_distant_sites_stay_two_tracks():
    left = slot(9.0) + np.array([-40.0, 0.0, 0.0])
    right = slot(9.0) + np.array([40.0, 0.0, 0.0])
    both = np.vstack([left, right])
    found = CV.detect_cavities(both, ["C"] * len(both), spacing=1.0)
    assert len(found.cavities) == 2
    report = EN.analyse_ensemble([found, found], tolerance=4.0)
    assert len(report.tracks) == 2


def test_an_empty_ensemble_is_refused():
    with pytest.raises(EN.EnsembleError, match="no conformers"):
        EN.analyse_ensemble([])


def test_a_nonpositive_tolerance_is_refused():
    with pytest.raises(EN.EnsembleError, match="tolerance"):
        EN.match_cavities(_maps([9.0]), tolerance=0.0)


def test_the_report_round_trips_to_a_dict():
    payload = EN.analyse_ensemble(_maps([5.0, 9.0]), tolerance=6.0).as_dict()
    assert payload["n_conformers"] == 2
    assert payload["n_cryptic"] == 1
    assert payload["tracks"][0]["classification"] == "cryptic"


def test_scan_with_nodes_puts_the_apo_first_and_one_map_per_conformer():
    import cryptic_ip.detect.openability as OA

    atoms = slot(9.0)
    nodes = np.arange(len(atoms))
    modes = OA.normal_modes(atoms, cutoff=9.0, n_modes=6)
    report = EN.scan_with_nodes(
        atoms, ["C"] * len(atoms), atoms, nodes, modes,
        amplitudes=(-1.0, 1.0), which_modes=(0,), spacing=1.2, tolerance=6.0,
    )
    # One apo plus two displaced conformers.
    assert report.n_conformers == 3


def test_a_mismatched_node_map_is_refused():
    import cryptic_ip.detect.openability as OA

    atoms = slot(9.0)
    modes = OA.normal_modes(atoms, cutoff=9.0, n_modes=5)
    with pytest.raises(EN.EnsembleError, match="one node index per atom"):
        EN.scan_with_nodes(atoms, ["C"] * len(atoms), atoms, np.arange(3), modes)


def test_a_node_index_outside_the_network_is_refused():
    import cryptic_ip.detect.openability as OA

    atoms = slot(9.0)
    modes = OA.normal_modes(atoms, cutoff=9.0, n_modes=5)
    bad = np.full(len(atoms), len(atoms) + 5)
    with pytest.raises(EN.EnsembleError, match="outside the network"):
        EN.scan_with_nodes(atoms, ["C"] * len(atoms), atoms, bad, modes)
