"""The elastic-network opening cost.

The decisive correctness check on an anisotropic network model is its null
space: a free elastic network has exactly six zero modes, three translations
and three rotations, and any implementation error in the Hessian shows up as a
different number.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import openability as OA


def helix(n=60):
    """An alpha helix: 1.5 A rise, 100 degrees of turn, 2.3 A radius."""
    t = np.arange(n)
    return np.column_stack([
        2.3 * np.cos(np.deg2rad(100 * t)),
        2.3 * np.sin(np.deg2rad(100 * t)),
        1.5 * t,
    ])


class _Arrays:
    def __init__(self, coords, names=None):
        self.coords = np.asarray(coords, dtype=float)
        n = len(self.coords)
        self.is_polymer = np.ones(n, dtype=bool)
        self.atom_names = np.array(names if names is not None else ["CA"] * n, dtype=object)
        self.chain_ids = np.array(["A"] * n, dtype=object)
        self.resseqs = np.arange(n)


# ------------------------------------------------------------------- Hessian

def test_the_hessian_is_symmetric():
    dense = OA.anm_hessian(helix()).toarray()
    assert np.allclose(dense, dense.T)


def test_the_hessian_has_exactly_six_rigid_body_modes():
    # Three translations and three rotations cost no energy in a free network.
    # A different count means the construction is wrong, and this is the single
    # check that catches almost every way of getting it wrong.
    dense = OA.anm_hessian(helix()).toarray()
    values = np.linalg.eigvalsh(dense)
    scale = float(values.max())
    assert int(np.sum(values <= 1e-8 * scale)) == 6


def test_a_rigid_translation_costs_no_energy():
    coords = helix()
    hessian = OA.anm_hessian(coords).toarray()
    translation = np.tile(np.array([1.0, 0.0, 0.0]), len(coords))
    assert abs(float(translation @ hessian @ translation)) < 1e-8 * float(np.abs(hessian).max())


def test_a_rigid_rotation_costs_no_energy():
    coords = helix()
    centred = coords - coords.mean(axis=0)
    hessian = OA.anm_hessian(coords).toarray()
    # An infinitesimal rotation about z: v = omega x r.
    rotation = np.cross(np.array([0.0, 0.0, 1.0]), centred).reshape(-1)
    energy = float(rotation @ hessian @ rotation)
    assert abs(energy) < 1e-6 * float(np.abs(hessian).max()) * float(np.sum(rotation ** 2))


def test_stretching_a_contact_costs_energy():
    coords = helix()
    hessian = OA.anm_hessian(coords).toarray()
    displacement = np.zeros(3 * len(coords))
    displacement[0] = 1.0            # move one node along x
    assert float(displacement @ hessian @ displacement) > 0.0


def test_the_spring_constant_scales_every_eigenvalue():
    # Gamma cancels from any ranking, which is why it is allowed to be 1.
    a = np.linalg.eigvalsh(OA.anm_hessian(helix(), gamma=1.0).toarray())
    b = np.linalg.eigvalsh(OA.anm_hessian(helix(), gamma=3.0).toarray())
    assert np.allclose(b, 3.0 * a)


def test_a_wider_cutoff_stiffens_the_network():
    short = OA.normal_modes(helix(), cutoff=8.0, n_modes=5)
    long = OA.normal_modes(helix(), cutoff=15.0, n_modes=5)
    assert long.eigenvalues[0] > short.eigenvalues[0]


@pytest.mark.parametrize("bad", [{"cutoff": 0.0}, {"cutoff": -1.0}])
def test_a_nonpositive_cutoff_is_refused(bad):
    with pytest.raises(OA.OpenabilityError, match="cutoff"):
        OA.anm_hessian(helix(), **bad)


def test_too_few_nodes_is_refused():
    with pytest.raises(OA.OpenabilityError, match="at least 4"):
        OA.anm_hessian(np.zeros((3, 3)))


# --------------------------------------------------------------------- modes

def test_normal_modes_skip_the_rigid_body_subspace():
    modes = OA.normal_modes(helix(), n_modes=10)
    assert modes.n_rigid_found == 6
    assert modes.n_modes == 10
    assert np.all(modes.eigenvalues > 0)
    assert np.all(np.diff(modes.eigenvalues) >= -1e-12)


def test_a_disconnected_network_is_refused_rather_than_analysed():
    # Two lobes with no contact have twelve zero modes, not six, and no
    # collective motion is defined across them. The guard failed at first
    # because it scaled its threshold by the largest *returned* eigenvalue -
    # which are the smallest ones - so every zero mode looked large next to
    # its peers.
    coords = helix()
    split = np.vstack([coords[:30], coords[30:] + np.array([500.0, 0.0, 0.0])])
    with pytest.raises(OA.OpenabilityError, match="disconnected"):
        OA.normal_modes(split, n_modes=5)


def test_a_mode_reshapes_to_one_displacement_per_node():
    modes = OA.normal_modes(helix(), n_modes=4)
    assert modes.displacement(0).shape == (60, 3)


def test_asking_for_a_mode_outside_the_set_is_refused():
    modes = OA.normal_modes(helix(), n_modes=3)
    with pytest.raises(OA.OpenabilityError, match="outside"):
        modes.displacement(5)


# ------------------------------------------------------------------ overlap

def test_the_opening_overlap_is_a_cosine_and_cannot_exceed_one():
    # It could, at first: the normalisation divided by the displacement norm
    # alone, when Cauchy-Schwarz bounds the summed projection by that norm
    # times the square root of the residue count. The result was "overlaps"
    # above 1.5 going into a ranking as though they were cosines.
    coords = helix()
    modes = OA.normal_modes(coords, n_modes=12)
    centre = coords[30] + np.array([6.0, 0.0, 0.0])
    lining = list(range(25, 36))
    outward = OA.outward_directions(coords[lining], centre)
    overlap = OA.opening_overlap(modes, lining, outward)
    assert np.all(overlap >= 0.0)
    assert np.all(overlap <= 1.0 + 1e-12)


def test_a_mode_that_moves_the_lining_straight_out_scores_one():
    # Hand-built: a mode whose displacement is exactly the outward direction.
    coords = np.array([[1.0, 0, 0], [-1.0, 0, 0], [0, 1.0, 0], [0, -1.0, 0]])
    outward = OA.outward_directions(coords, (0, 0, 0))

    class _Modes:
        n_modes = 1
        n_nodes = 4

        def displacement(self, mode):
            return outward.copy()

    overlap = OA.opening_overlap(_Modes(), [0, 1, 2, 3], outward)
    assert overlap[0] == pytest.approx(1.0)


def test_a_mode_perpendicular_to_the_opening_scores_zero():
    coords = np.array([[1.0, 0, 0], [-1.0, 0, 0]])
    outward = OA.outward_directions(coords, (0, 0, 0))

    class _Modes:
        n_modes = 1
        n_nodes = 2

        def displacement(self, mode):
            return np.array([[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]])

    assert OA.opening_overlap(_Modes(), [0, 1], outward)[0] == pytest.approx(0.0)


def test_an_empty_lining_is_refused():
    modes = OA.normal_modes(helix(), n_modes=3)
    with pytest.raises(OA.OpenabilityError, match="no lining"):
        OA.opening_overlap(modes, [], np.zeros((0, 3)))


def test_a_residue_on_the_pocket_centre_is_refused():
    with pytest.raises(OA.OpenabilityError, match="exactly on the pocket centre"):
        OA.outward_directions(np.zeros((1, 3)), (0, 0, 0))


# ------------------------------------------------------------------- report

def test_describe_pocket_reports_an_openability_and_a_cost():
    coords = helix()
    report = OA.describe_pocket(_Arrays(coords), coords[30] + np.array([6.0, 0, 0]), n_modes=10)
    assert 0.0 <= report.openability <= 1.0
    assert report.opening_cost_gamma_a2 > 0.0
    assert report.n_lining > 0
    assert 0 <= report.best_mode < 10


def test_the_cost_is_not_labelled_kt_because_it_is_not_in_kt():
    # The field was called opening_cost_kt until it was noticed that the spring
    # constant has never been calibrated against a measured opening free
    # energy. This pipeline has already read fpocket's mean local hydrophobic
    # density as a burial depth in angstrom; a name that implies a unit the
    # number does not have is the same mistake.
    coords = helix()
    report = OA.describe_pocket(_Arrays(coords), coords[30] + np.array([6.0, 0, 0]), n_modes=6)
    assert hasattr(report, "opening_cost_gamma_a2")
    assert not hasattr(report, "opening_cost_kt")
    assert OA.COST_UNITS == "gamma*A^2"


def test_the_cost_grows_as_the_square_of_the_opening():
    coords = helix()
    arrays = _Arrays(coords)
    centre = coords[30] + np.array([6.0, 0, 0])
    one = OA.describe_pocket(arrays, centre, n_modes=6, opening=1.0).opening_cost_gamma_a2
    two = OA.describe_pocket(arrays, centre, n_modes=6, opening=2.0).opening_cost_gamma_a2
    assert two == pytest.approx(4.0 * one, rel=1e-9)


def test_explicit_lining_keys_are_used_when_given():
    coords = helix()
    report = OA.describe_pocket(
        _Arrays(coords), coords[30] + np.array([6.0, 0, 0]),
        lining_keys=["A:28", "A:29", "A:30"], n_modes=5,
    )
    assert report.n_lining == 3


def test_lining_keys_that_match_nothing_are_refused():
    coords = helix()
    with pytest.raises(OA.OpenabilityError, match="lining residue keys"):
        OA.describe_pocket(
            _Arrays(coords), coords[30] + np.array([6.0, 0, 0]),
            lining_keys=["Z:999"], n_modes=5,
        )


def test_a_pocket_centre_far_from_the_protein_is_refused():
    coords = helix()
    with pytest.raises(OA.OpenabilityError, match="within"):
        OA.describe_pocket(_Arrays(coords), (500.0, 0.0, 0.0), n_modes=5)


def test_a_structure_with_no_alpha_carbons_is_refused():
    coords = helix(10)
    with pytest.raises(OA.OpenabilityError, match="alpha carbons"):
        OA.describe_pocket(
            _Arrays(coords, names=["CB"] * 10), coords[5] + np.array([6.0, 0, 0]), n_modes=3,
        )


def test_the_report_round_trips_to_a_dict():
    coords = helix()
    payload = OA.describe_pocket(_Arrays(coords), coords[30] + np.array([6.0, 0, 0]), n_modes=5).as_dict()
    assert payload["opening_cost_gamma_a2"] > 0
    assert len(payload["per_mode_overlap"]) == 5


# ------------------------------------------------------------- the ensemble

def test_perturbed_conformers_move_by_the_requested_amplitude():
    coords = helix()
    modes = OA.normal_modes(coords, n_modes=5)
    conformers = OA.perturbed_conformers(coords, modes, amplitudes=(2.0,), which_modes=(0,))
    assert len(conformers) == 1
    moved = float(np.max(np.linalg.norm(conformers[0] - coords, axis=1)))
    assert moved == pytest.approx(2.0, rel=1e-9)


def test_both_signs_of_a_mode_are_generated():
    # A pocket that opens one way usually closes the other, so an ensemble that
    # only went one direction would miss half the openings.
    coords = helix()
    modes = OA.normal_modes(coords, n_modes=5)
    conformers = OA.perturbed_conformers(coords, modes, amplitudes=(-1.0, 1.0), which_modes=(0,))
    assert len(conformers) == 2
    assert not np.allclose(conformers[0], conformers[1])
    assert np.allclose(conformers[0] + conformers[1], 2.0 * coords)


def test_the_default_ensemble_is_twelve_conformers():
    coords = helix()
    modes = OA.normal_modes(coords, n_modes=5)
    assert len(OA.perturbed_conformers(coords, modes)) == 12


def test_perturbing_a_mismatched_structure_is_refused():
    modes = OA.normal_modes(helix(), n_modes=3)
    with pytest.raises(OA.OpenabilityError, match="different structures"):
        OA.perturbed_conformers(helix(30), modes)
