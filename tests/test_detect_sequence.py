"""Tests for the information-theoretic sequence descriptors.

The module they cover is built around three ways of getting this question
wrong, so the tests are organised the same way: the complexity proxies are
checked against textbook parse counts rather than against a notion of
Kolmogorov complexity, the entropy estimators are checked against a known
truth at the window sizes a pocket lining actually has, and the composition
null is checked by planting an ordering signal and then destroying it with a
shuffle that leaves the composition untouched.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.detect import sequence as sq


# --------------------------------------------------------------- fixtures

def helical_face(n: int = 36, period: float = 3.6, seed: int = 11) -> str:
    """A window with a basic residue every ``period`` positions.

    The filler alphabet holds no basic residue, so the only arginines are the
    planted ones and the composition is a function of ``n`` and ``period``
    alone. This is the positive control for a helical face lining a pocket.
    """
    filler = "ALVGST"
    marks = {int(round(i * period)) for i in range(int(n / period) + 1)}
    rng = np.random.default_rng(seed)
    return "".join(
        "R" if i in marks else filler[rng.integers(0, len(filler))]
        for i in range(n)
    )


def shuffled(sequence: str, seed: int) -> str:
    """The same letters in a different order: the composition-matched control."""
    letters = list(sequence)
    np.random.default_rng(seed).shuffle(letters)
    return "".join(letters)


# --------------------------------------------------------------- alphabets

def test_recode_identity_is_the_input_when_it_is_standard():
    assert sq.recode("ARNDCQ", "identity") == "ARNDCQ"


def test_recode_lowercases_in():
    assert sq.recode("arndcq", "identity") == "ARNDCQ"


def test_recode_charge_classes_histidine_on_its_own():
    assert sq.recode("RKHDEA", "charge") == "++h--0"


def test_recode_drops_rather_than_buckets_nonstandard_residues():
    # Assigning X a class would invent information; the length difference says
    # how many were dropped.
    assert sq.recode("RKXZ-U", "charge") == "++"


def test_recode_covers_every_amino_acid_in_every_alphabet():
    for name in sq.ALPHABETS:
        assert len(sq.recode(sq.AMINO_ACIDS, name)) == 20


def test_recode_rejects_an_unknown_alphabet():
    with pytest.raises(sq.SequenceError):
        sq.recode("ARND", "hydrophilicity")


def test_symbol_counts_are_descending():
    assert list(sq.symbol_counts("AABBC")) == [2, 2, 1]


def test_symbol_counts_of_nothing_is_empty():
    assert sq.symbol_counts("").size == 0


# ----------------------------------------------------------------- entropy

def test_a_single_symbol_has_no_entropy():
    assert sq.shannon_entropy("AAAA") == 0.0


def test_plugin_entropy_of_two_equal_symbols_is_one_bit():
    assert sq.shannon_entropy("ABAB", estimator="plugin") == pytest.approx(1.0)


def test_entropy_needs_symbols():
    with pytest.raises(sq.SequenceError):
        sq.shannon_entropy("")


def test_entropy_rejects_an_unknown_estimator():
    with pytest.raises(sq.SequenceError):
        sq.shannon_entropy("ARND", estimator="jackknife")


@pytest.mark.parametrize("n", [20, 100])
def test_grassberger_beats_miller_madow_beats_plugin_at_pocket_sizes(n):
    """The measurement that fixes the default estimator.

    A uniform source over the twenty amino acids has entropy log2(20) = 4.322
    bits. At the size of a pocket lining the plug-in estimate is badly low and
    - the part that ruins a comparison between two windows of different length
    - how low depends on N. The ordering of the three estimators' bias is the
    reason Grassberger is the default.
    """
    truth = float(np.log2(20))
    rng = np.random.default_rng(7)
    biases = {"plugin": [], "miller_madow": [], "grassberger": []}
    for _ in range(200):
        window = "".join(rng.choice(list(sq.AMINO_ACIDS), size=n))
        for name, values in biases.items():
            values.append(sq.shannon_entropy(window, estimator=name) - truth)
    mean = {name: float(np.mean(values)) for name, values in biases.items()}
    assert mean["plugin"] < 0.0
    assert abs(mean["grassberger"]) < abs(mean["plugin"])
    assert abs(mean["miller_madow"]) < abs(mean["plugin"])


def test_plugin_bias_depends_on_window_length():
    """Why plug-in entropy cannot compare a 12-residue lining with a 20-residue one."""
    truth = float(np.log2(20))
    rng = np.random.default_rng(13)
    bias = {}
    for n in (12, 20):
        values = [
            sq.shannon_entropy("".join(rng.choice(list(sq.AMINO_ACIDS), size=n)),
                               estimator="plugin") - truth
            for _ in range(200)
        ]
        bias[n] = float(np.mean(values))
    # Both low, and the shorter window lower: a difference of this size would
    # be reported as a signal by a plug-in comparison.
    assert bias[12] < bias[20] < 0.0
    assert bias[20] - bias[12] > 0.2


def test_block_entropy_needs_enough_symbols_for_the_block():
    with pytest.raises(sq.SequenceError):
        sq.block_entropy("ARN", 4)


def test_block_entropy_rejects_a_block_shorter_than_one():
    with pytest.raises(sq.SequenceError):
        sq.block_entropy("ARND", 0)


def test_block_entropy_saturates_when_every_kmer_is_unique():
    assert sq.block_entropy_is_saturated("ARNDCQEG", 3)
    assert not sq.block_entropy_is_saturated("ARARARAR", 2)


def test_saturation_is_reported_for_a_window_shorter_than_the_block():
    assert sq.block_entropy_is_saturated("AR", 3)


def test_entropy_rate_needs_k_of_at_least_two():
    with pytest.raises(sq.SequenceError):
        sq.entropy_rate("ARNDCQ", 1)


def test_a_repeat_has_a_lower_entropy_rate_than_its_own_shuffle():
    """The measure that separates ordering from composition, on its own."""
    repeat = "ARARARARARARARARARAR"
    scramble = shuffled(repeat, seed=4)
    # Same letters, so the same single-symbol entropy.
    assert sq.shannon_entropy(repeat) == pytest.approx(sq.shannon_entropy(scramble))
    assert sq.entropy_rate(repeat, 2) <= sq.entropy_rate(scramble, 2)


def test_excess_entropy_needs_k_max_of_at_least_two():
    with pytest.raises(sq.SequenceError):
        sq.excess_entropy("ARNDCQ", 1)


# -------------------------------------------------------------- complexity

@pytest.mark.parametrize("text,parses", [
    ("a", 1),
    ("aaaaaaaa", 2),
    ("abababab", 3),
    ("abcdefgh", 8),
])
def test_lempel_ziv_parse_counts_match_the_textbook(text, parses):
    assert sq.lempel_ziv_complexity(text) == parses


def test_lempel_ziv_needs_a_sequence():
    with pytest.raises(sq.SequenceError):
        sq.lempel_ziv_complexity("")


def test_lempel_ziv_is_lower_for_a_repeat_than_for_its_shuffle():
    repeat = "ARARARARARARARARARAR"
    assert sq.lempel_ziv_complexity(repeat) < sq.lempel_ziv_complexity(
        shuffled(repeat, seed=8))


def test_normalised_lz_is_zero_for_a_one_letter_window():
    assert sq.lz_normalised("AAAAAAAA") == 0.0


def test_normalised_lz_needs_two_symbols_of_length():
    with pytest.raises(sq.SequenceError):
        sq.lz_normalised("A")


def test_the_parse_obeys_its_own_definition():
    """Each phrase novel, each phrase's own prefix not: the definition itself.

    This is the test the textbook cases could not do. The first parser here
    passed ``aaaaaaaa`` to 2, ``abababab`` to 3 and ``abcdefgh`` to 8 while
    advancing one symbol per phrase on everything else, because those three
    strings are insensitive to the fault.
    """
    rng = np.random.default_rng(23)
    for alphabet in ("ab", "abc", sq.AMINO_ACIDS):
        window = "".join(rng.choice(list(alphabet), size=60))
        phrases = sq.lempel_ziv_phrases(window)
        assert "".join(phrases) == window
        start = 0
        for index, phrase in enumerate(phrases):
            history = window[:start + len(phrase) - 1]
            is_last = index == len(phrases) - 1
            if not is_last:
                # A completed phrase is novel in the history before its last symbol.
                assert phrase not in history
            if len(phrase) > 1:
                # And its own prefix was not, which is why it had to grow.
                assert phrase[:-1] in window[:start + len(phrase) - 2]
            start += len(phrase)


def test_a_random_window_parses_into_far_fewer_phrases_than_its_length():
    """The regression that caught the degenerate parser.

    A parser that advances one symbol per phrase reports close to ``n``. The
    correct parse of a random 500-residue window over twenty letters is about
    213 phrases against an asymptotic expectation of 241; the broken one
    reported 485.
    """
    rng = np.random.default_rng(19)
    window = "".join(rng.choice(list(sq.AMINO_ACIDS), size=500))
    parses = sq.lempel_ziv_complexity(window)
    assert 150 < parses < 280
    assert parses < 0.6 * len(window)


@pytest.mark.parametrize("n,expected", [(20, 0.88), (36, 0.84), (100, 0.83)])
def test_the_lz_normaliser_is_asymptotic_and_misses_at_pocket_sizes(n, expected):
    """Why `lz_normalised` must not be compared across window lengths.

    ``n / log_a(n)`` is the parse count of a random string in the limit. At the
    length of a pocket lining the ratio sits below 1 and is not even monotone
    in N, because the observed-alphabet correction shrinks as the finite-size
    deficit closes. That residual is the same kind of trap as the plug-in
    entropy bias, and the composition null is the answer to both.
    """
    rng = np.random.default_rng(29)
    values = [
        sq.lz_normalised("".join(rng.choice(list(sq.AMINO_ACIDS), size=n)))
        for _ in range(60)
    ]
    assert float(np.mean(values)) == pytest.approx(expected, abs=0.03)


def test_the_normaliser_is_not_monotone_in_window_length():
    """The measurement behind the warning, in one assertion."""
    rng = np.random.default_rng(29)
    mean = {}
    for n in (20, 100):
        mean[n] = float(np.mean([
            sq.lz_normalised("".join(rng.choice(list(sq.AMINO_ACIDS), size=n)))
            for _ in range(60)
        ]))
    assert mean[20] > mean[100]


def test_the_lz_normaliser_does_converge_on_a_binary_alphabet():
    rng = np.random.default_rng(31)
    window = "".join(rng.choice(list("01"), size=4000))
    assert sq.lz_normalised(window) == pytest.approx(1.0, abs=0.08)


def test_compressor_overhead_is_measured_not_assumed():
    for method, overhead in sq.COMPRESSOR_OVERHEAD.items():
        assert overhead > 0, method


def test_compressed_bits_subtracts_the_header():
    """On a pocket-sized window the header is most of the output."""
    window = sq.AMINO_ACIDS
    raw_bytes = len(__import__("zlib").compress(window.encode()))
    assert sq.compressed_bits(window, "zlib") == pytest.approx(
        (raw_bytes - sq.COMPRESSOR_OVERHEAD["zlib"]) * 8)
    # And the header really is the larger part at this length.
    assert sq.COMPRESSOR_OVERHEAD["zlib"] * 8 > 0.0


def test_compressed_bits_never_goes_negative():
    assert sq.compressed_bits("A", "zlib") >= 0.0


def test_compressed_bits_rejects_an_unknown_method():
    with pytest.raises(sq.SequenceError):
        sq.compressed_bits("ARND", "brotli")


def test_compressed_bits_needs_a_sequence():
    with pytest.raises(sq.SequenceError):
        sq.compressed_bits("", "zlib")


# ------------------------------------------------------------- periodicity

def test_charge_track_reads_the_charges():
    assert list(sq.as_track("RKDEA", "charge")) == [1.0, 1.0, -1.0, -1.0, 0.0]


def test_track_rejects_an_unknown_name():
    with pytest.raises(sq.SequenceError):
        sq.as_track("ARND", "polarity")


def test_track_needs_recognised_residues():
    with pytest.raises(sq.SequenceError):
        sq.as_track("XZB", "charge")


def test_periodicity_is_bounded_by_one_at_the_nyquist_period():
    """The bug this estimator was rewritten to fix.

    A pure alternating signal is entirely explained at period 2, but the sine
    quadrature vanishes there, so normalising by ``variance * n / 2`` returned
    2.0 - and an explained-variance fraction cannot exceed 1.
    """
    alternating = np.array([1.0, -1.0] * 10)
    assert sq.periodicity_power(alternating, 2.0) == pytest.approx(1.0)


def test_a_constant_series_explains_nothing():
    assert sq.periodicity_power(np.ones(10), 3.6) == 0.0


def test_periodicity_needs_a_period_above_one_residue():
    with pytest.raises(sq.SequenceError):
        sq.periodicity_power(np.arange(10.0), 1.0)


def test_periodicity_needs_a_few_points():
    with pytest.raises(sq.SequenceError):
        sq.periodicity_power(np.arange(3.0), 3.6)


def test_an_exact_probe_misses_a_period_four_face_that_the_band_finds():
    """Why the helical probe is a band and not the single value 3.6.

    A strictly period-4 charge train is *exactly* orthogonal to a 3.6-residue
    probe over 36 points, so the exact probe scores it at zero for a sequence
    that plainly has a periodic basic face.
    """
    track = sq.as_track("".join("R" if i % 4 == 0 else "A" for i in range(36)))
    assert sq.periodicity_power(track, 3.6) < 1e-6
    assert sq.band_periodicity(track, sq.HELIX_BAND) > 0.5


def test_the_helix_band_includes_a_three_ten_helix():
    """The lower edge is 3.0 because a 3-10 helix is exactly 3.0 per turn."""
    track = sq.as_track("".join("R" if i % 3 == 0 else "A" for i in range(36)))
    assert sq.HELIX_BAND[0] <= 3.0
    assert sq.band_periodicity(track, sq.HELIX_BAND) > 0.9


def test_band_periodicity_rejects_an_inverted_band():
    with pytest.raises(sq.SequenceError):
        sq.band_periodicity(np.arange(10.0), (4.2, 3.0))


def test_band_periodicity_needs_at_least_two_probes():
    with pytest.raises(sq.SequenceError):
        sq.band_periodicity(np.arange(10.0), sq.HELIX_BAND, probes=1)


def test_a_period_seven_pattern_has_helical_band_power_but_is_not_a_helix():
    """The harmonic confound, and why the dominant period is reported too.

    A sparse charge pattern repeating every 7 residues has content at 3.5 as
    well, so it scores in the helical band. The fundamental still wins the
    scan, which is what separates a real helical face from a harmonic of
    something else.
    """
    window = "".join("R" if i % 7 == 0 else "A" for i in range(42))
    assert sq.band_periodicity(sq.as_track(window), sq.HELIX_BAND) > 0.25
    dominant = sq.dominant_charge_period(window)
    assert dominant["period"] == pytest.approx(7.0, abs=0.3)
    assert dominant["in_helix_band"] == 0.0


def test_the_dominant_period_of_a_helical_face_is_in_the_helical_band():
    dominant = sq.dominant_charge_period(helical_face())
    assert dominant["period"] == pytest.approx(3.6, abs=0.3)
    assert dominant["in_helix_band"] == 1.0
    assert dominant["power"] > 0.3


def test_the_dominant_period_of_an_alternating_face_is_in_the_strand_band():
    window = "".join("R" if i % 2 == 0 else "A" for i in range(20))
    dominant = sq.dominant_charge_period(window)
    assert dominant["in_strand_band"] == 1.0
    assert dominant["in_helix_band"] == 0.0


def test_dominant_period_rejects_an_inverted_band():
    with pytest.raises(sq.SequenceError):
        sq.dominant_period(np.arange(10.0), (12.0, 2.0))


# ---------------------------------------------------------------- spacing

def test_spacing_of_helical_gaps_is_all_helical():
    out = sq.spacing_spectrum([10, 13, 17, 20, 24])
    assert out["fraction_helical"] == 1.0
    assert out["fraction_strand"] == 0.0
    assert out["median_gap"] == pytest.approx(3.5)
    assert out["span"] == 14.0
    assert out["n_positions"] == 5.0


def test_spacing_of_alternating_gaps_is_all_strand():
    out = sq.spacing_spectrum([4, 6, 8, 10])
    assert out["fraction_strand"] == 1.0
    assert out["fraction_helical"] == 0.0


def test_spacing_marks_distant_positions_as_not_local():
    out = sq.spacing_spectrum([1, 2, 200])
    assert out["fraction_local"] == pytest.approx(0.5)


def test_spacing_is_order_independent():
    assert sq.spacing_spectrum([24, 10, 17, 13, 20]) == sq.spacing_spectrum(
        [10, 13, 17, 20, 24])


def test_spacing_needs_two_positions():
    with pytest.raises(sq.SequenceError):
        sq.spacing_spectrum([7])


# ----------------------------------------------------- the composition null

def test_every_registered_statistic_runs_on_a_pocket_sized_window():
    window = helical_face(20)
    for name in sq.STATISTICS:
        result = sq.composition_z(window, name, n_shuffles=50, seed=1)
        assert result.statistic == name
        assert np.isfinite(result.z)


def test_the_composition_null_rejects_an_unknown_statistic():
    with pytest.raises(sq.SequenceError):
        sq.composition_z(helical_face(20), "fractal_dimension", n_shuffles=10)


def test_the_composition_null_needs_a_window():
    with pytest.raises(sq.SequenceError):
        sq.composition_z("ARN", "entropy", n_shuffles=10)


def test_the_composition_null_is_seeded():
    first = sq.composition_z(helical_face(20), "lz_normalised", n_shuffles=100, seed=3)
    second = sq.composition_z(helical_face(20), "lz_normalised", n_shuffles=100, seed=3)
    assert first.z == second.z
    assert first.p_two_sided == second.p_two_sided


def test_the_permutation_p_value_can_never_be_zero():
    result = sq.composition_z(helical_face(36), "helix_charge_periodicity",
                              n_shuffles=100, seed=3)
    assert result.p_two_sided >= 1.0 / 101.0


def test_a_planted_helical_face_is_not_explained_by_its_composition():
    """The decisive test: an ordering signal the shuffle null cannot reproduce."""
    result = sq.composition_z(helical_face(36), "helix_charge_periodicity",
                              n_shuffles=400, seed=5)
    assert result.z > 2.0
    assert not result.is_composition_explained


def test_the_same_letters_shuffled_show_no_ordering_signal():
    """The negative control, with composition held exactly fixed.

    These windows have the same length, the same letter counts and therefore
    the same entropy as the planted one. If the machinery reported a signature
    here it would be reporting composition.
    """
    planted = helical_face(36)
    scores = []
    for seed in range(12):
        control = shuffled(planted, seed=seed + 100)
        assert sq.shannon_entropy(control) == pytest.approx(
            sq.shannon_entropy(planted))
        scores.append(sq.composition_z(control, "helix_charge_periodicity",
                                       n_shuffles=200, seed=2).z)
    scores = np.asarray(scores)
    assert abs(float(scores.mean())) < 2.0
    assert float(np.mean(np.abs(scores) > 2.0)) < 0.25


def test_a_null_result_serialises_its_own_verdict():
    result = sq.composition_z(helical_face(20), "entropy", n_shuffles=50, seed=1)
    out = result.as_dict()
    assert out["statistic"] == "entropy"
    assert out["n_shuffles"] == 50
    assert out["composition_explained"] is result.is_composition_explained


# ------------------------------------------------------------- the bundle

def test_describe_sequence_calls_a_planted_helical_face_a_helix_face():
    features = sq.describe_sequence(helical_face(36), n_shuffles=100)
    assert features.raw["helix_face"] == 1.0
    assert features.raw["strand_face"] == 0.0
    assert features.raw["dominant_charge_period"] == pytest.approx(3.6, abs=0.3)


def test_describe_sequence_does_not_call_a_harmonic_a_helix_face():
    window = "".join("R" if i % 7 == 0 else "A" for i in range(42))
    features = sq.describe_sequence(window, n_shuffles=100)
    assert features.raw["helix_face"] == 0.0


def test_describe_sequence_shows_the_plugin_bias_rather_than_hiding_it():
    features = sq.describe_sequence(helical_face(20), n_shuffles=50)
    assert features.raw["entropy_bias_visible"] == pytest.approx(
        features.raw["entropy"] - features.raw["entropy_plugin"])
    assert features.raw["entropy_bias_visible"] > 0.0


def test_describe_sequence_reports_saturated_blocks():
    features = sq.describe_sequence("ARNDCQEGHILK", n_shuffles=20)
    # Twelve distinct residues: every 2-mer and 3-mer is unique.
    assert 2 in features.saturated_blocks
    assert 3 in features.saturated_blocks


def test_describe_sequence_fills_every_requested_alphabet():
    features = sq.describe_sequence(helical_face(36), n_shuffles=20,
                                    alphabets=tuple(sq.ALPHABETS))
    assert set(features.per_alphabet) == set(sq.ALPHABETS)
    for entry in features.per_alphabet.values():
        assert entry["alphabet_size"] >= 1.0


def test_describe_sequence_drops_nonstandard_residues_from_the_length():
    features = sq.describe_sequence("ARNDXXCQ", n_shuffles=20)
    assert features.length == 6
    assert "X" not in features.sequence


def test_describe_sequence_refuses_a_window_it_cannot_measure():
    with pytest.raises(sq.SequenceError):
        sq.describe_sequence("ARX", n_shuffles=10)


def test_ordering_signals_keeps_only_what_the_shuffle_could_not_reproduce():
    features = sq.describe_sequence(helical_face(36), n_shuffles=200)
    signals = features.ordering_signals()
    for name, z in signals.items():
        assert abs(z) >= 2.0, name
    for name, result in features.nulls.items():
        if result.is_composition_explained:
            assert name not in signals


def test_the_bundle_serialises():
    out = sq.describe_sequence(helical_face(20), n_shuffles=20).as_dict()
    assert out["length"] == 20
    assert set(out) >= {"raw", "nulls", "per_alphabet", "saturated_blocks",
                        "ordering_signals"}


# ------------------------------------------------------- association testing

def test_a_comparison_below_five_groups_is_not_evidence():
    """The project's standing floor, enforced before any difference is read."""
    comparison = sq.compare_feature(
        "x", {"a": 1.0, "b": 1.2, "c": 1.1}, {"d": 0.0, "e": 0.1, "f": 0.2},
        n_bootstrap=100)
    assert comparison.n_site_groups == 3
    assert not comparison.is_evidence


def test_a_comparison_at_five_groups_each_side_is_evidence():
    sites = {f"s{i}": v for i, v in enumerate([1.0, 1.2, 0.9, 1.1, 1.3])}
    controls = {f"c{i}": v for i, v in enumerate([0.1, 0.0, 0.2, -0.1, 0.1])}
    assert sq.compare_feature("x", sites, controls, n_bootstrap=100).is_evidence


def test_a_clear_separation_has_an_interval_that_excludes_zero():
    sites = {f"s{i}": v for i, v in enumerate([1.0, 1.2, 0.9, 1.1, 1.3, 1.0])}
    controls = {f"c{i}": v for i, v in enumerate([0.1, 0.0, 0.2, -0.1, 0.1, 0.0])}
    comparison = sq.compare_feature("x", sites, controls, n_bootstrap=400)
    assert comparison.difference > 0.0
    assert comparison.low > 0.0
    assert comparison.exceeds_power
    assert comparison.p_value < 0.05


def test_overlapping_arms_do_not_separate():
    values = [0.0, 0.2, -0.1, 0.1, 0.05, -0.05]
    sites = {f"s{i}": v for i, v in enumerate(values)}
    controls = {f"c{i}": v for i, v in enumerate(reversed(values))}
    comparison = sq.compare_feature("x", sites, controls, n_bootstrap=400)
    assert comparison.low <= 0.0 <= comparison.high
    assert comparison.p_value > 0.05
    assert not comparison.exceeds_power


def test_a_comparison_needs_both_arms():
    with pytest.raises(sq.SequenceError):
        sq.compare_feature("x", {}, {"c": 0.0}, n_bootstrap=10)


def test_a_comparison_serialises_its_floor_and_its_power():
    sites = {f"s{i}": float(i) for i in range(5)}
    controls = {f"c{i}": float(i) + 1.0 for i in range(5)}
    out = sq.compare_feature("x", sites, controls, n_bootstrap=100).as_dict()
    assert out["is_evidence"] is True
    assert len(out["interval"]) == 2
    assert out["minimum_detectable"] > 0.0


def test_the_minimum_detectable_difference_is_infinite_without_groups():
    assert sq.minimum_detectable_difference([1.0], [0.0]) == float("inf")


def test_the_minimum_detectable_difference_shrinks_as_groups_are_added():
    spread = [0.0, 1.0, -1.0, 0.5, -0.5]
    few = sq.minimum_detectable_difference(spread, spread)
    many = sq.minimum_detectable_difference(spread * 6, spread * 6)
    assert many < few


def test_only_the_implemented_power_is_accepted():
    with pytest.raises(sq.SequenceError):
        sq.minimum_detectable_difference([0.0, 1.0], [0.0, 1.0], power=0.9)


def test_holm_scales_the_smallest_p_by_the_family_size():
    assert sq.holm({"a": 0.001, "b": 0.02, "c": 0.3})["a"] == pytest.approx(0.003)


def test_holm_is_monotone_in_the_ordering():
    corrected = sq.holm({"a": 0.001, "b": 0.02, "c": 0.3, "d": 0.5})
    ordered = [corrected[k] for k in sorted(corrected, key=lambda k: corrected[k])]
    assert ordered == sorted(ordered)


def test_holm_leaves_a_family_of_one_alone():
    assert sq.holm({"only": 0.04})["only"] == pytest.approx(0.04)


def test_holm_never_exceeds_one():
    assert all(p <= 1.0 for p in sq.holm({"a": 0.4, "b": 0.5, "c": 0.9}).values())


def test_holm_covers_every_feature_tried():
    family = {f"f{i}": 0.01 * (i + 1) for i in range(6)}
    assert set(sq.holm(family)) == set(family)


# ------------------------------- the hydropathy track (study Q9)

def test_an_amphipathic_helix_shows_hydropathy_periodicity():
    """One face hydrophobic, the other not, at the helical repeat."""
    amphipathic = "LKKLLKKLLKKLLKKLLKKL"
    assert sq.helix_hydropathy_periodicity(amphipathic) > 0.8


def test_the_same_letters_in_blocks_do_not():
    """Composition is identical; only the arrangement differs."""
    amphipathic = "LKKLLKKLLKKLLKKLLKKL"
    blocky = "LLLLLLLLLLKKKKKKKKKK"
    assert sorted(amphipathic) == sorted(blocky)
    assert sq.helix_hydropathy_periodicity(blocky) < sq.helix_hydropathy_periodicity(amphipathic)


def test_the_hydropathy_z_separates_arrangement_from_letters():
    """The whole design of study Q9 in one assertion."""
    amphipathic = sq.composition_z("LKKLLKKLLKKLLKKLLKKL",
                                   "helix_hydropathy_periodicity", n_shuffles=400)
    blocky = sq.composition_z("LLLLLLLLLLKKKKKKKKKK",
                              "helix_hydropathy_periodicity", n_shuffles=400)
    assert amphipathic.z > 3.0
    assert blocky.z < 2.0


def test_the_strand_hydropathy_statistic_fires_on_alternation():
    alternating = "LKLKLKLKLKLKLKLKLKLK"
    assert sq.strand_hydropathy_periodicity(alternating) > 0.8
    assert sq.helix_hydropathy_periodicity(alternating) < \
        sq.strand_hydropathy_periodicity(alternating)


def test_the_hydropathy_statistics_are_registered_for_the_null():
    for name in ("helix_hydropathy_periodicity", "strand_hydropathy_periodicity",
                 "dominant_hydropathy_period_power"):
        assert name in sq.STATISTICS
        assert sq.composition_z("LKKLLKKLLKKLLKKLLKKL", name, n_shuffles=50).n_shuffles == 50


def test_hydropathy_and_charge_are_different_questions():
    """A hydrophobic face with no charge pattern must not fire the charge test."""
    face = "LAALLAALLAALLAALLAAL"
    assert sq.helix_hydropathy_periodicity(face) > 0.5
    assert sq.helix_charge_periodicity(face) < 0.5


# ------------------------------------- the vectorised periodogram is the same fit

def _lstsq_reference(series, period):
    """The least-squares fit `periodicity_power` was written as, verbatim.

    Kept here rather than in the module, so the vectorised form is pinned
    against the code it replaced instead of against itself.
    """
    series = np.asarray(series, dtype=float)
    n = series.size
    centred = series - series.mean()
    variance = float(np.sum(centred * centred))
    if variance == 0.0:
        return 0.0
    angles = 2.0 * np.pi * np.arange(n) / float(period)
    basis = np.column_stack([np.cos(angles), np.sin(angles)])
    keep = np.linalg.norm(basis, axis=0) > 1e-12
    basis = basis[:, keep]
    if basis.shape[1] == 0:
        return 0.0
    coefficients, *_ = np.linalg.lstsq(basis, centred, rcond=None)
    explained = float(np.sum((basis @ coefficients) ** 2))
    return float(min(1.0, explained / variance))


def test_the_vectorised_fit_is_the_least_squares_fit():
    """A performance change that moved a number would have moved every result."""
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(200):
        n = int(rng.integers(4, 60))
        series = rng.normal(size=n)
        for period in list(rng.uniform(1.01, 14.0, size=5)) + [2.0, 3.0, 3.6, 4.2]:
            worst = max(worst, abs(_lstsq_reference(series, period)
                                   - sq.periodicity_power(series, period)))
    assert worst < 1e-9


def test_the_nyquist_quadrature_is_still_dropped():
    """At period 2 the sine vanishes; fitting it would return more than one."""
    alternating = np.array([1.0, -1.0] * 12)
    assert sq.periodicity_power(alternating, 2.0) == pytest.approx(1.0)
    assert sq.periodicity_power(alternating, 2.0) <= 1.0


def test_a_grid_of_probes_agrees_with_the_probes_one_at_a_time():
    rng = np.random.default_rng(3)
    series = rng.normal(size=37)
    grid = [2.0, 2.5, 3.0, 3.6, 4.2, 7.0, 11.5]
    together = sq.probe_powers(series, grid)
    apart = [sq.periodicity_power(series, period) for period in grid]
    assert together == pytest.approx(apart, abs=1e-12)


def test_a_flat_series_has_no_periodicity_at_any_probe():
    flat = np.zeros(20)
    assert sq.probe_powers(flat, [2.0, 3.6, 9.0]) == pytest.approx([0.0, 0.0, 0.0])


def test_a_probe_grid_is_still_refused_where_a_single_probe_would_be():
    with pytest.raises(sq.SequenceError):
        sq.probe_powers(np.zeros(3), [3.6])
    with pytest.raises(sq.SequenceError):
        sq.probe_powers(np.zeros(20), [1.0])
    with pytest.raises(sq.SequenceError):
        sq.probe_powers(np.zeros(20), [])


def test_the_basis_cache_does_not_leak_between_window_lengths():
    """Same period grid, different lengths: each must get its own basis."""
    rng = np.random.default_rng(7)
    short, long = rng.normal(size=12), rng.normal(size=48)
    assert sq.probe_powers(short, [3.6])[0] == pytest.approx(
        _lstsq_reference(short, 3.6), abs=1e-9)
    assert sq.probe_powers(long, [3.6])[0] == pytest.approx(
        _lstsq_reference(long, 3.6), abs=1e-9)
