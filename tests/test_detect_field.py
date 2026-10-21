"""The rotation-invariant electrostatic descriptor.

The central claim is that the shell-resolved power spectrum does not change
when the structure is rotated about the pocket centre. That claim is exact for
the continuous integral and approximate for a sampled one, so the tests check
both the algebra it rests on (an orthonormal basis within each degree) and the
rate at which the sampled version approaches it.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import field as F


def _random_charges(seed=0, n=60, scale=6.0):
    rng = np.random.default_rng(seed)
    coords = rng.normal(size=(n, 3)) * scale
    charges = rng.choice([-1.0, -0.5, 0.5, 1.0], size=n)
    return F.ChargeSet(coords, charges)


# --------------------------------------------------------------- the sphere

def test_the_sampling_sphere_is_near_uniform():
    points = F.fibonacci_sphere(2000)
    assert np.allclose(np.linalg.norm(points, axis=1), 1.0)
    # A uniform set has a centroid at the origin and an isotropic second moment.
    assert np.linalg.norm(points.mean(axis=0)) < 0.02
    second = points.T @ points / len(points)
    assert np.allclose(second, np.eye(3) / 3.0, atol=0.02)


def test_the_sampling_sphere_is_deterministic():
    # Two spectra are only comparable if they were sampled in the same
    # directions, so this must not depend on any random state.
    assert np.array_equal(F.fibonacci_sphere(64), F.fibonacci_sphere(64))


def test_too_few_points_is_refused():
    with pytest.raises(F.FieldError, match="at least 4"):
        F.fibonacci_sphere(3)


# ------------------------------------------------------------- the harmonics

@pytest.mark.parametrize("max_degree", [0, 1, 2, 3, 4])
def test_the_basis_is_orthonormal_within_each_degree(max_degree):
    # This is the test that was missing. The first implementation used
    # hand-written Cartesian forms; l=0..3 were exactly orthonormal and l=4 had
    # diagonal entries as low as 0.25, so every obvious check passed and the
    # only symptom was a rotation invariance that would not converge.
    points = F.fibonacci_sphere(20000)
    weight = 4.0 * np.pi / len(points)
    for degree, basis in enumerate(F.real_spherical_harmonics(points, max_degree)):
        assert basis.shape == (len(points), 2 * degree + 1)
        gram = weight * (basis.T @ basis)
        assert np.max(np.abs(gram - np.eye(gram.shape[0]))) < 1e-4, f"degree {degree} is not orthonormal"


def test_the_monopole_is_the_constant_function():
    points = F.fibonacci_sphere(500)
    basis = F.real_spherical_harmonics(points, 0)[0]
    assert np.allclose(basis, 0.5 / np.sqrt(np.pi))


def test_higher_degrees_integrate_to_zero_against_a_constant():
    points = F.fibonacci_sphere(20000)
    weight = 4.0 * np.pi / len(points)
    for degree, basis in enumerate(F.real_spherical_harmonics(points, 4)):
        if degree == 0:
            continue
        assert np.max(np.abs(weight * basis.sum(axis=0))) < 1e-4


def test_a_zero_direction_is_refused():
    with pytest.raises(F.FieldError, match="cannot be zero"):
        F.real_spherical_harmonics(np.zeros((1, 3)), 2)


def test_a_negative_degree_is_refused():
    with pytest.raises(F.FieldError, match="cannot be negative"):
        F.real_spherical_harmonics(F.fibonacci_sphere(10), -1)


# -------------------------------------------------------------- the potential

def test_a_single_positive_charge_gives_a_positive_screened_field():
    charges = F.ChargeSet(np.zeros((1, 3)), np.array([1.0]))
    near = F.screened_potential(charges, np.array([[3.0, 0.0, 0.0]]))
    far = F.screened_potential(charges, np.array([[12.0, 0.0, 0.0]]))
    assert near[0] > far[0] > 0.0


def test_screening_makes_the_field_fall_faster_than_coulomb():
    charges = F.ChargeSet(np.zeros((1, 3)), np.array([1.0]))
    r1, r2 = 4.0, 16.0
    screened = F.screened_potential(charges, np.array([[r1, 0, 0], [r2, 0, 0]]))
    # Unscreened would fall as r1/r2 = 1/4; screening must be steeper.
    assert screened[1] / screened[0] < 0.25


def test_the_field_is_capped_near_a_charge_rather_than_diverging():
    charges = F.ChargeSet(np.zeros((1, 3)), np.array([1.0]))
    on_top = F.screened_potential(charges, np.zeros((1, 3)))
    assert np.isfinite(on_top[0])


def test_a_net_neutral_but_separated_pair_has_a_dipole_and_no_monopole():
    charges = F.ChargeSet(np.array([[3.0, 0, 0], [-3.0, 0, 0]]), np.array([1.0, -1.0]))
    spectrum = F.field_spectrum(charges, (0, 0, 0), points_per_shell=1536)
    # The monopole is not exactly zero because the quadrature is not exact; it
    # is zero to the precision the sampling supports, and the dipole dominates
    # it by ten orders of magnitude.
    assert abs(spectrum.monopole) < 1e-3
    assert spectrum.power[0, 1] > 1e6 * spectrum.power[0, 0]


def test_no_charges_gives_a_flat_field_rather_than_an_error():
    empty = F.ChargeSet(np.zeros((0, 3)), np.zeros(0))
    assert np.allclose(F.screened_potential(empty, F.fibonacci_sphere(10)), 0.0)


def test_a_nonpositive_debye_length_is_refused():
    with pytest.raises(F.FieldError, match="debye_length"):
        F.screened_potential(_random_charges(), F.fibonacci_sphere(10), debye_length=0.0)


# ------------------------------------------------------- rotation invariance

@pytest.mark.parametrize(
    "axis,angle",
    [((0, 0, 1), 0.7), ((1, 0, 0), 2.1), ((0.3, -0.7, 0.5), 1.1), ((-0.5, 0.1, 0.8), 3.0)],
)
def test_the_power_spectrum_survives_a_rotation(axis, angle):
    charges = _random_charges()
    before = F.field_spectrum(charges, (0, 0, 0), points_per_shell=1536)
    turned = F.ChargeSet(F.rotate(charges.coords, axis, angle), charges.charges)
    after = F.field_spectrum(turned, (0, 0, 0), points_per_shell=1536)
    assert F.spectrum_discrepancy(before.power, after.power) < 5e-3


def test_the_invariance_improves_as_the_sampling_is_refined():
    # An invariance that does not converge is a basis error, not a quadrature
    # error. That distinction is what found the l=4 bug, so it is pinned.
    charges = _random_charges()
    errors = []
    for points in (192, 768, 3072):
        before = F.field_spectrum(charges, (0, 0, 0), points_per_shell=points)
        turned = F.ChargeSet(F.rotate(charges.coords, (0.3, -0.7, 0.5), 1.1), charges.charges)
        after = F.field_spectrum(turned, (0, 0, 0), points_per_shell=points)
        errors.append(F.spectrum_discrepancy(before.power, after.power))
    assert errors[0] > errors[1] > errors[2]
    assert errors[-1] < 1e-3


def test_translating_the_whole_system_with_the_centre_changes_nothing():
    charges = _random_charges()
    shift = np.array([11.0, -4.0, 7.0])
    before = F.field_spectrum(charges, (0, 0, 0), points_per_shell=768)
    moved = F.ChargeSet(charges.coords + shift, charges.charges)
    after = F.field_spectrum(moved, shift, points_per_shell=768)
    assert np.allclose(before.power, after.power, rtol=1e-9, atol=1e-12)


def test_the_spectrum_is_not_invariant_to_moving_the_centre_alone():
    # A descriptor that ignored where the pocket is would be useless, so this
    # checks the invariance is the intended one and not something broader.
    charges = _random_charges()
    here = F.field_spectrum(charges, (0, 0, 0), points_per_shell=768)
    there = F.field_spectrum(charges, (6.0, 0, 0), points_per_shell=768)
    assert not np.allclose(here.power, there.power, rtol=0.05)


def test_doubling_every_charge_quadruples_the_power():
    charges = _random_charges()
    single = F.field_spectrum(charges, (0, 0, 0), points_per_shell=768)
    double = F.field_spectrum(
        F.ChargeSet(charges.coords, charges.charges * 2.0), (0, 0, 0), points_per_shell=768,
    )
    assert np.allclose(double.power, 4.0 * single.power, rtol=1e-9)


# ----------------------------------------------------------- the cutoff

def test_the_screening_cutoff_changes_the_field_by_less_than_it_claims():
    # The cutoff is an approximation and its error must stay below the
    # quadrature error it sits beside, or it is the dominant error.
    rng = np.random.default_rng(1)
    coords = rng.normal(size=(600, 3)) * 25
    charges = F.ChargeSet(coords, rng.choice([-1.0, 1.0], size=600))
    exact = F.field_spectrum(charges, (0, 0, 0), cutoff_debye=None, points_per_shell=768)
    cut = F.field_spectrum(charges, (0, 0, 0), cutoff_debye=5.0, points_per_shell=768)
    # Measured in spectrum_discrepancy, which is the metric the module's own
    # figures are quoted in: a component-wise relative error is meaningless on
    # a spectrum with near-cancelling components, and measuring it that way
    # was how this test first failed for the wrong reason.
    assert F.spectrum_discrepancy(exact.power, cut.power) < 2e-3


def test_a_tighter_cutoff_is_a_worse_approximation():
    rng = np.random.default_rng(1)
    coords = rng.normal(size=(600, 3)) * 25
    charges = F.ChargeSet(coords, rng.choice([-1.0, 1.0], size=600))
    exact = F.field_spectrum(charges, (0, 0, 0), cutoff_debye=None, points_per_shell=384)

    def error(cutoff):
        got = F.field_spectrum(charges, (0, 0, 0), cutoff_debye=cutoff, points_per_shell=384)
        return F.spectrum_discrepancy(exact.power, got.power)

    assert error(2.0) > error(5.0) > error(8.0)


# ------------------------------------------------------------ charge reading

class _Arrays:
    """The minimum of StructureArrays that protein_charges reads."""

    def __init__(self, rows):
        self.coords = np.array([r[0] for r in rows], dtype=float)
        self.resnames = np.array([r[1] for r in rows], dtype=object)
        self.atom_names = np.array([r[2] for r in rows], dtype=object)
        self.is_polymer = np.array([r[3] for r in rows], dtype=bool)


def test_arginine_and_lysine_are_read_as_positive():
    arrays = _Arrays([
        ([0, 0, 0], "ARG", "NH1", True),
        ([1, 0, 0], "ARG", "NH2", True),
        ([2, 0, 0], "LYS", "NZ", True),
    ])
    charges = F.protein_charges(arrays)
    assert charges.net_charge == pytest.approx(2.0)


def test_aspartate_and_glutamate_are_read_as_negative():
    arrays = _Arrays([
        ([0, 0, 0], "ASP", "OD1", True),
        ([1, 0, 0], "ASP", "OD2", True),
        ([2, 0, 0], "GLU", "OE1", True),
        ([3, 0, 0], "GLU", "OE2", True),
    ])
    assert F.protein_charges(arrays).net_charge == pytest.approx(-2.0)


def test_a_bound_ligands_own_charge_is_not_part_of_the_pockets_field():
    # A structure crystallised with a sulfate must not look like a site that
    # attracts one because the sulfate is in the field.
    arrays = _Arrays([
        ([0, 0, 0], "LYS", "NZ", True),
        ([1, 0, 0], "SO4", "O1", False),
    ])
    charges = F.protein_charges(arrays)
    assert charges.net_charge == pytest.approx(1.0)
    assert len(charges.charges) == 1


def test_a_backbone_carbon_carries_no_charge():
    arrays = _Arrays([([0, 0, 0], "ALA", "CB", True)])
    assert len(F.protein_charges(arrays).charges) == 0


def test_histidine_is_only_partly_protonated_at_physiological_ph():
    arrays = _Arrays([([0, 0, 0], "HIS", "ND1", True), ([1, 0, 0], "HIS", "NE2", True)])
    net = F.protein_charges(arrays).net_charge
    assert 0.0 < net < 0.5


# ------------------------------------------------------- ligand templates

def test_every_template_reaches_its_own_net_charge():
    for name, template in F.LIGAND_TEMPLATES.items():
        assert template.shell_charge[-1] == pytest.approx(template.net_charge), name


def test_template_charge_accumulates_monotonically_outward():
    for name, template in F.LIGAND_TEMPLATES.items():
        cumulative = np.array(template.shell_charge)
        assert np.all(np.diff(cumulative) <= 1e-9), f"{name} loses charge going outward"


def test_the_panel_holds_the_competitors_the_project_could_not_separate():
    # 30 of 75 candidates are flagged as explained by another phosphate- or
    # sulfate-rich ligand, so those ligands have to be in the panel for a
    # specificity claim to mean anything.
    assert {"IP6", "ATP", "PAPS", "UDP_GLUCOSE", "SULFATE", "PHOSPHATE"} <= set(F.LIGAND_TEMPLATES)


def test_ip6_and_a_heparin_fragment_have_similar_charge_but_different_shape():
    # The pair that net charge alone cannot separate: the descriptor has to.
    ip6 = F.LIGAND_TEMPLATES["IP6"]
    heparin = F.LIGAND_TEMPLATES["HEPARIN_TETRA"]
    assert abs(ip6.net_charge - heparin.net_charge) <= 1.0
    assert heparin.radius_of_gyration > ip6.radius_of_gyration * 1.5


def test_a_template_realises_as_a_charge_set_with_the_right_total():
    charges = F.LIGAND_TEMPLATES["IP6"].as_charge_set()
    assert charges.net_charge == pytest.approx(-9.0, abs=1e-6)


# ------------------------------------------------------- complementarity

def _basic_shell(radius=4.5, n=12, charge=1.0):
    """A ring of positive charge at one radius: an idealised sequestering site."""
    directions = F.fibonacci_sphere(n)
    return F.ChargeSet(radius * directions, np.full(n, charge))


def test_a_positive_shell_at_the_phosphate_radius_has_a_favourable_depth():
    spectrum = F.field_spectrum(_basic_shell(), (0, 0, 0), points_per_shell=768)
    result = F.complementarity(spectrum, F.LIGAND_TEMPLATES["IP6"])
    assert result["depth_kt"] > 0.0


def test_a_negative_pocket_is_unfavourable_for_a_polyanion():
    spectrum = F.field_spectrum(_basic_shell(charge=-1.0), (0, 0, 0), points_per_shell=768)
    result = F.complementarity(spectrum, F.LIGAND_TEMPLATES["IP6"])
    assert result["depth_kt"] < 0.0


def test_depth_and_shape_are_reported_separately():
    # Folding them into one number is what lets a shallow surface patch with
    # many arginines outscore a buried site.
    spectrum = F.field_spectrum(_basic_shell(), (0, 0, 0), points_per_shell=768)
    result = F.complementarity(spectrum, F.LIGAND_TEMPLATES["IP6"])
    assert set(result) == {"shape", "depth_kt", "net_charge", "template_rg"}


def test_shape_is_scale_free_in_the_charge_magnitude():
    weak = F.field_spectrum(_basic_shell(charge=0.2), (0, 0, 0), points_per_shell=768)
    strong = F.field_spectrum(_basic_shell(charge=2.0), (0, 0, 0), points_per_shell=768)
    a = F.complementarity(weak, F.LIGAND_TEMPLATES["IP6"])
    b = F.complementarity(strong, F.LIGAND_TEMPLATES["IP6"])
    assert a["shape"] == pytest.approx(b["shape"], abs=1e-6)
    assert b["depth_kt"] > a["depth_kt"]


def test_depth_scales_with_the_charge_magnitude():
    weak = F.field_spectrum(_basic_shell(charge=0.5), (0, 0, 0), points_per_shell=768)
    strong = F.field_spectrum(_basic_shell(charge=1.0), (0, 0, 0), points_per_shell=768)
    a = F.complementarity(weak, F.LIGAND_TEMPLATES["IP6"])["depth_kt"]
    b = F.complementarity(strong, F.LIGAND_TEMPLATES["IP6"])["depth_kt"]
    assert b == pytest.approx(2.0 * a, rel=1e-6)


def test_comparing_spectra_of_different_shape_is_refused():
    a = F.field_spectrum(_basic_shell(), (0, 0, 0), max_degree=2, points_per_shell=384)
    b = F.field_spectrum(_basic_shell(), (0, 0, 0), max_degree=4, points_per_shell=384)
    with pytest.raises(F.FieldError, match="same shape"):
        F._correlation(a.feature_vector(), b.feature_vector())


# ------------------------------------------------------------- the report

def test_describe_pocket_scores_the_whole_panel():
    arrays = _Arrays(
        [(list(4.5 * d), "ARG", "NH1", True) for d in F.fibonacci_sphere(12)]
    )
    report = F.describe_pocket(arrays, (0, 0, 0), points_per_shell=384)
    assert set(report.per_ligand) == set(F.LIGAND_TEMPLATES)
    name, shape = report.best_ligand()
    assert name in F.LIGAND_TEMPLATES
    assert -1.0 <= shape <= 1.0


def test_the_report_round_trips_to_a_dict():
    arrays = _Arrays([(list(4.5 * d), "LYS", "NZ", True) for d in F.fibonacci_sphere(8)])
    payload = F.describe_pocket(arrays, (0, 0, 0), points_per_shell=192).as_dict()
    assert "spectrum" in payload and "ligands" in payload
    assert len(payload["spectrum"]["power"]) == len(F.DEFAULT_SHELLS)


def test_best_ligand_on_an_empty_report_is_refused():
    spectrum = F.field_spectrum(_basic_shell(), (0, 0, 0), points_per_shell=192)
    with pytest.raises(F.FieldError, match="no ligands"):
        F.PocketFieldReport(spectrum=spectrum).best_ligand()


def test_the_feature_labels_match_the_feature_vector():
    spectrum = F.field_spectrum(_basic_shell(), (0, 0, 0), points_per_shell=192)
    labels = F.iter_shell_labels(spectrum.shells, spectrum.max_degree)
    assert len(labels) == spectrum.feature_vector().size
