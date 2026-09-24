"""Tests for the geometry primitives, checked against analytic ground truth.

These assertions are exact statements about geometry, not tolerances chosen to
fit the implementation: an isolated sphere has a known surface area, a point at
the centre of a closed shell is a known distance from it, and so on. An
implementation that fails them is wrong regardless of what it produces on real
structures.
"""

import numpy as np
import pytest

from cryptic_ip.analysis.geometry import (
    PROBE_RADIUS,
    burial_depth,
    convex_hull_volume,
    count_within,
    element_radius,
    enclosure_fraction,
    fibonacci_sphere,
    pairwise_min_distance,
    point_cloud_shape,
    shrake_rupley_sasa,
    total_sasa,
)


def test_fibonacci_sphere_points_are_unit_vectors():
    points = fibonacci_sphere(500)
    assert points.shape == (500, 3)
    assert np.allclose(np.linalg.norm(points, axis=1), 1.0)


def test_fibonacci_sphere_is_well_distributed():
    """The centroid of a uniform spherical covering sits at the origin."""
    points = fibonacci_sphere(2000)
    assert np.allclose(points.mean(axis=0), 0.0, atol=1e-2)


def test_fibonacci_sphere_rejects_non_positive_counts():
    with pytest.raises(ValueError):
        fibonacci_sphere(0)


def test_isolated_atom_sasa_matches_the_analytic_sphere_area():
    radius = 1.7
    sasa = total_sasa(np.array([[0.0, 0.0, 0.0]]), np.array([radius]), n_points=4096)
    assert sasa == pytest.approx(4 * np.pi * (radius + PROBE_RADIUS) ** 2, rel=1e-6)


def test_atom_enclosed_by_a_shell_has_zero_sasa():
    shell = fibonacci_sphere(200) * 4.0
    coords = np.vstack([[[0.0, 0.0, 0.0]], shell])
    radii = np.full(len(coords), 1.7)
    assert shrake_rupley_sasa(coords, radii, subset=[0], n_points=2048)[0] == pytest.approx(0.0)


def test_sasa_subset_is_measured_in_the_full_context():
    """The occluding set and the measured set are independent, by design."""
    coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    radii = np.array([1.7, 1.7])
    in_context = total_sasa(coords, radii, subset=[0], n_points=2048)
    isolated = total_sasa(coords[:1], radii[:1], n_points=2048)
    assert in_context < isolated


def test_sasa_rejects_mismatched_shapes():
    with pytest.raises(ValueError):
        shrake_rupley_sasa(np.zeros((3, 3)), np.zeros(2))
    with pytest.raises(ValueError):
        shrake_rupley_sasa(np.zeros((3, 2)), np.zeros(3))


def test_burial_depth_is_the_distance_to_the_nearest_exposed_atom():
    coords = np.array([[10.0, 0.0, 0.0], [-10.0, 0.0, 0.0]])
    sasa = np.array([50.0, 50.0])
    assert burial_depth([0.0, 0.0, 0.0], coords, sasa) == pytest.approx(10.0)


def test_burial_depth_ignores_atoms_below_the_exposure_threshold():
    coords = np.array([[2.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    sasa = np.array([0.5, 50.0])
    # The nearby atom is not solvent exposed, so the surface is the distant one.
    assert burial_depth([0.0, 0.0, 0.0], coords, sasa, min_atom_sasa=5.0) == pytest.approx(10.0)


def test_burial_depth_without_exposed_atoms_is_nan():
    coords = np.array([[1.0, 0.0, 0.0]])
    assert np.isnan(burial_depth([0.0, 0.0, 0.0], coords, np.array([0.0])))


def test_enclosure_is_one_inside_a_shell_and_zero_far_outside():
    shell = fibonacci_sphere(400) * 6.0
    radii = np.full(len(shell), 1.7)
    assert enclosure_fraction([0.0, 0.0, 0.0], shell, radii) == pytest.approx(1.0)
    assert enclosure_fraction([500.0, 0.0, 0.0], shell, radii) == pytest.approx(0.0)


def test_enclosure_decreases_monotonically_with_distance_from_a_shell():
    """Moving out of a cavity must reduce enclosure at every step.

    The absolute value just outside a hollow shell is not one half: rays leaving
    tangentially still strike the far side of the shell. Monotonicity is the
    property the descriptor actually relies on - buried sites score higher than
    surface sites - so that is what is asserted.
    """
    shell = fibonacci_sphere(800) * 10.0
    radii = np.full(len(shell), 1.7)
    values = [
        enclosure_fraction([offset, 0.0, 0.0], shell, radii)
        for offset in (0.0, 11.0, 20.0, 40.0)
    ]
    assert values[0] == pytest.approx(1.0)
    assert values == sorted(values, reverse=True)
    assert values[-1] < 0.15


def test_point_cloud_shape_of_a_sphere_is_isotropic():
    shape = point_cloud_shape(fibonacci_sphere(500) * 5.0)
    assert shape.radius_of_gyration == pytest.approx(5.0, rel=1e-3)
    assert shape.asphericity == pytest.approx(0.0, abs=1e-2)


def test_point_cloud_shape_of_a_line_is_maximally_aspherical():
    line = np.column_stack([np.linspace(-10, 10, 50), np.zeros(50), np.zeros(50)])
    assert point_cloud_shape(line).asphericity == pytest.approx(1.0, abs=1e-6)


def test_point_cloud_shape_of_degenerate_input_is_nan():
    assert np.isnan(point_cloud_shape(np.zeros((1, 3))).radius_of_gyration)


def test_convex_hull_volume_approximates_a_sphere():
    volume = convex_hull_volume(fibonacci_sphere(1000) * 4.0)
    assert volume == pytest.approx(4 / 3 * np.pi * 4.0**3, rel=0.02)


def test_convex_hull_volume_of_too_few_points_is_nan():
    assert np.isnan(convex_hull_volume(np.zeros((3, 3))))


def test_pairwise_min_distance_and_count_within():
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[3.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    assert pairwise_min_distance(a, b) == pytest.approx(3.0)
    count, indices = count_within(a, b, 5.0)
    assert count == 1
    assert indices.tolist() == [0]
    assert np.isnan(pairwise_min_distance(np.empty((0, 3)), b))


def test_element_radius_falls_back_for_unknown_elements():
    assert element_radius("P") == pytest.approx(1.80)
    assert element_radius("p") == pytest.approx(1.80)
    assert element_radius("Xx") == pytest.approx(1.70)
    assert element_radius(None) == pytest.approx(1.70)
