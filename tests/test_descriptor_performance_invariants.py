"""Guards for the descriptor optimisations: faster, never different.

The per-structure caches and the nearest-first SASA are pure optimisations.
These tests pin that they change no value: SASA against a plain reference
implementation, and every pocket's descriptors against a fresh extractor, in
any order.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.analysis.features import PocketFeatureExtractor
from cryptic_ip.analysis.geometry import PROBE_RADIUS, fibonacci_sphere, shrake_rupley_sasa
from cryptic_ip.analysis.structure_arrays import load_structure_arrays
from cryptic_ip.testing.synthetic import default_benchmark_specs, write_synthetic_structure


def _reference_sasa(coords, radii, n_points):
    """Shrake-Rupley written as plainly as possible."""
    expanded = radii + PROBE_RADIUS
    sphere = fibonacci_sphere(n_points)
    out = np.zeros(len(coords))
    for i in range(len(coords)):
        points = coords[i] + expanded[i] * sphere
        buried = np.zeros(len(points), dtype=bool)
        for j in range(len(coords)):
            if j != i:
                buried |= np.sum((points - coords[j]) ** 2, axis=1) < expanded[j] ** 2
        area = 4.0 * np.pi * expanded[i] * expanded[i]
        out[i] = area * float(np.count_nonzero(~buried)) / len(points)
    return out


@pytest.fixture(scope="module")
def arrays(tmp_path_factory):
    spec = default_benchmark_specs(n_buried=1, n_surface=0, n_decoy=0, seed=4)[0]
    return load_structure_arrays(write_synthetic_structure(spec, tmp_path_factory.mktemp("perf")))


def test_sasa_matches_the_reference_exactly(arrays):
    coords, radii = arrays.coords[:400], arrays.radii[:400]
    fast = shrake_rupley_sasa(coords, radii, n_points=96)
    assert np.array_equal(fast, _reference_sasa(coords, radii, 96))


def test_sasa_subset_matches_full(arrays):
    coords, radii = arrays.coords[:300], arrays.radii[:300]
    full = shrake_rupley_sasa(coords, radii, n_points=64)
    subset = np.arange(0, 300, 7)
    assert np.array_equal(shrake_rupley_sasa(coords, radii, subset=subset, n_points=64), full[subset])


def _features(extractor, centres):
    return [extractor.extract(k, c).features for k, c in enumerate(centres)]


def _same(a, b):
    for x, y in zip(a, b):
        assert x.keys() == y.keys()
        for key in x:
            if isinstance(x[key], float) and np.isnan(x[key]):
                assert np.isnan(y[key]), key
            else:
                assert x[key] == y[key], key


def test_cached_descriptors_do_not_depend_on_order(arrays):
    rng = np.random.default_rng(0)
    centres = arrays.coords[rng.choice(len(arrays.coords), 12, replace=False)]
    forward = _features(PocketFeatureExtractor(arrays, n_points=64), centres)
    shared = PocketFeatureExtractor(arrays, n_points=64)
    backward = [shared.extract(k, c).features for k, c in reversed(list(enumerate(centres)))][::-1]
    _same(forward, backward)
    # Each pocket alone, on a fresh extractor, gives the same values too.
    for k, centre in enumerate(centres):
        _same([forward[k]], [PocketFeatureExtractor(arrays, n_points=64).extract(k, centre).features])
