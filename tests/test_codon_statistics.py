import math

import numpy as np
import pytest

from grodon.codon import codon_counts, milc
from grodon.statistics import (
    SequenceFilterWarning,
    _trim_sequences,
    filter_sequences,
    get_codon_statistics,
    get_codon_statistics_i,
    weighted_median,
)


@pytest.fixture
def genes():
    return [
        "ATG" + "GCTGCCGCA" * 30 + "TAA",
        "ATG" + "GCCGCCGCT" * 30 + "TAA",
        "ATG" + "GCTGCTGCA" * 30 + "TAA",
    ]


def test_filtering_threshold_is_strict_and_metadata_stays_aligned():
    sequences = ["AAA" * 80, "AAA" * 81, "AAN" * 90, "CCC" * 82, "AAA" * 90 + "A"]
    with pytest.warns(SequenceFilterWarning, match="Filtered 3"):
        result = filter_sequences(
            sequences, [False, True, False, True, False], depth_of_coverage=[1, 2, 3, 4, 5]
        )
    assert result.filtered == 3
    np.testing.assert_array_equal(result.indices, [1, 3])
    np.testing.assert_array_equal(result.depth_of_coverage, [2, 4])
    assert result.genes == (sequences[1], sequences[3])


def test_missing_highly_expressed_genes_fail_early():
    with pytest.raises(ValueError, match="No highly expressed"):
        filter_sequences(["AAA" * 90], [False])
    with pytest.warns(SequenceFilterWarning), pytest.raises(ValueError, match="after filtering"):
        filter_sequences(["AAA"], [True])


@pytest.mark.parametrize("he", [[True], ["yes", "no"], [1, 0]])
def test_invalid_highly_expressed_annotations_rejected(he):
    with pytest.raises((ValueError, TypeError)):
        filter_sequences(["AAA" * 90, "CCC" * 90], he)


def test_weighted_median_is_interpolated():
    # Midpoints are 0.125 and 0.625, so their interpolation at 0.5 is 7.5.
    assert weighted_median([0, 10], [1, 3]) == pytest.approx(7.5)
    assert weighted_median([10, 0, 999], [3, 1, 0]) == pytest.approx(7.5)
    assert weighted_median([1, 2, 5, 10], [1, 1, 1, 1]) == pytest.approx(3.5)
    assert math.isnan(weighted_median([0, np.nan], [1, 1]))


@pytest.mark.parametrize("weights", [[0, 0], [1, -1], [1, float("inf")], [1]])
def test_weighted_median_invalid_weights(weights):
    with pytest.raises(ValueError):
        weighted_median([0, 1], weights)


def test_aggregates_and_weighted_gc(genes):
    he = [True, True, False]
    result = get_codon_statistics(genes, he)
    counts = codon_counts(genes)
    assert result["CUB"] == pytest.approx(np.median(milc(counts)))
    assert result["CUBHE"] == pytest.approx(np.median(milc(counts)[:2]))
    assert result["ConsistencyHE"] == pytest.approx(np.mean(milc(counts[:2])))
    assert result["FilteredSequences"] == 0
    assert result["nHE"] == 2
    weighted = get_codon_statistics(genes, he, depth_of_coverage=[1, 3, 0])
    gc = [(g.count("G") + g.count("C")) / len(g) for g in genes]
    assert weighted["GC"] == pytest.approx((gc[0] + 3 * gc[1]) / 4)
    assert weighted["CUBHE"] == pytest.approx(weighted_median(milc(counts)[:2], [1, 3]))
    assert math.isnan(weighted["CPB"])


def test_equal_coverage_reproduces_unweighted_statistics(genes):
    unweighted = get_codon_statistics(genes, [True, True, False])
    weighted = get_codon_statistics(genes, [True, True, False], depth_of_coverage=[1, 1, 1])
    for key in ("CUB", "CUBHE", "ConsistencyHE", "GC", "GCdiv"):
        assert weighted[key] == pytest.approx(unweighted[key])


def test_coverage_units_do_not_overflow_weighted_estimates(genes):
    ordinary = get_codon_statistics(genes, [True, True, False], depth_of_coverage=[1, 3, 1])
    large = get_codon_statistics(
        genes, [True, True, False], depth_of_coverage=[1e307, 3e307, 1e307]
    )
    for key in ("CUB", "CUBHE", "ConsistencyHE", "GC", "GCdiv"):
        assert large[key] == pytest.approx(ordinary[key])


def test_all_metrics_returns_r_names(genes):
    result = get_codon_statistics(genes, [True, True, False], all_metrics=True)
    assert set(result) == {"MILC", "ENCprime", "B", "SCUO", "MCB", "nHE"}
    assert all(np.isfinite(value) for value in result.values())


def test_fragments_and_disabled_threshold():
    gene = "AAA" * 60
    assert get_codon_statistics([gene], [True], fragments=True)["nHE"] == 1
    assert get_codon_statistics(["AAA"], [True], fragments=None)["nHE"] == 1
    assert get_codon_statistics([gene], [True], fragments=np.bool_(True))["nHE"] == 1


def test_end_trim_and_random_trim_preserve_frame_and_requested_length():
    gene = "AAA" * 30 + "CCC" * 30 + "GGG" * 30
    rng = np.random.default_rng(12)
    assert _trim_sequences([gene], 90, "end", rng) == ("GGG" * 30,)
    for length in (3, 90, 150, 267):
        trimmed = _trim_sequences([gene], length, "random", rng)[0]
        assert len(trimmed) == length
        assert gene.find(trimmed) % 3 == 0


def test_short_exclusion_only_applies_when_length_filter_disabled(genes):
    sequences = [genes[0], "AAA" * 30]
    with pytest.warns(SequenceFilterWarning):
        result = get_codon_statistics(
            sequences, [True, False], fragments=None, exclude_short=True, trimlen=150
        )
    assert result["FilteredSequences"] == 1


def test_individual_background_is_seeded_and_cpb_is_undefined(genes):
    first = get_codon_statistics_i(genes, [True, False, False], n_le=3, random_state=8)
    second = get_codon_statistics_i(genes, [True, False, False], n_le=3, random_state=8)
    for key in ("CUBHE", "CUB", "GC", "ConsistencyHE"):
        assert first[key] == second[key]
    assert math.isnan(first["CPB"])
    different = get_codon_statistics_i(genes, [True, False, False], n_le=3, random_state=9)
    assert first["CUBHE"] != different["CUBHE"]


def test_individual_all_metrics_and_weighting(genes):
    result = get_codon_statistics_i(
        genes,
        [True, False, False],
        n_le=2,
        random_state=4,
        all_metrics=True,
        depth_of_coverage=[1, 2, 3],
    )
    assert set(result) == {"MILC", "ENCprime", "B", "SCUO", "MCB", "nHE"}
    assert all(np.isfinite(value) for value in result.values())
    with pytest.raises(ValueError, match="non-highly-expressed"):
        get_codon_statistics_i(genes, [True, True, True])


def test_individual_metrics_preserve_upstream_standard_code_convention(genes):
    bacterial = get_codon_statistics_i(
        genes, [True, False, False], n_le=2, random_state=11, genetic_code=11
    )
    standard = get_codon_statistics_i(
        genes, [True, False, False], n_le=2, random_state=11, genetic_code=1
    )
    assert bacterial["CUBHE"] == standard["CUBHE"]
    assert bacterial["CUB"] == standard["CUB"]
    assert bacterial["ConsistencyHE"] != standard["ConsistencyHE"]


def test_individual_zero_coverage_genes_do_not_cause_random_sampling_failure(genes):
    # Only one low-expression gene has positive coverage, with n_le=1. Every
    # seed must select that eligible gene and yield finite aggregate statistics.
    for seed in range(8):
        result = get_codon_statistics_i(
            genes, [True, False, False], n_le=1, random_state=seed, depth_of_coverage=[2, 0, 3]
        )
        assert np.isfinite(result["CUB"])
        scaled = get_codon_statistics_i(
            genes, [True, False, False], n_le=1, random_state=seed, depth_of_coverage=[20, 0, 30]
        )
        assert scaled["CUB"] == pytest.approx(result["CUB"])
        assert scaled["CUBHE"] == pytest.approx(result["CUBHE"])
    with pytest.raises(ValueError, match="Non-highly-expressed genes need positive coverage"):
        get_codon_statistics_i(genes, [True, False, False], n_le=1, depth_of_coverage=[1, 0, 0])


@pytest.mark.parametrize(
    "options",
    [
        {"trimlen": 151},
        {"trimside": "left"},
        {"n_le": 0},
        {"background": "bad"},
        {"fragments": "yes"},
        {"depth_of_coverage": [1, 2]},
        {"depth_of_coverage": [0, 0, 1]},
    ],
)
def test_invalid_statistics_options(genes, options):
    with pytest.raises((ValueError, TypeError)):
        get_codon_statistics(genes, [True, True, False], **options)
