"""Analytical checks use artificial sequences, with no organism examples."""

import math

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from grodon.codon import (
    CODONS,
    codon_counts,
    codon_pair_bias,
    codon_usage,
    enc_prime,
    karlin_b,
    mcb,
    milc,
    scuo,
)


def test_counts_alphabetical_and_biopython_inputs():
    counts = codon_counts(["aaacccgggttt", Seq("AAATTT"), SeqRecord(Seq("GGG"))])
    assert counts.shape == (3, 64)
    assert CODONS == tuple(sorted(CODONS))
    np.testing.assert_array_equal(counts.sum(axis=1), [4, 2, 1])
    assert counts[0, CODONS.index("CCC")] == 1
    assert counts[1, CODONS.index("AAA")] == 1


def test_gene_mapping_uses_values_in_insertion_order():
    np.testing.assert_array_equal(
        codon_counts({"gene_a": "AAA", "gene_b": "CCC"}), codon_counts(["AAA", "CCC"])
    )


@pytest.mark.parametrize("genes", [["AA"], ["AAN"], ["AAé"]])
def test_low_level_counts_reject_invalid_sequences(genes):
    with pytest.raises(ValueError):
        codon_counts(genes)


def test_metrics_against_analytical_two_codon_distribution():
    counts = codon_counts(["AAA" * 80 + "AAG" * 20, "AAA" * 20 + "AAG" * 80])
    expected_milc = 2 * (0.8 * math.log(1.6) + 0.2 * math.log(0.4)) - 0.01 + 0.5
    expected_scuo = 1 + (0.8 * math.log(0.8) + 0.2 * math.log(0.2)) / math.log(2)
    np.testing.assert_allclose(milc(counts), expected_milc)
    np.testing.assert_allclose(karlin_b(counts), 0.6)
    np.testing.assert_allclose(scuo(counts), expected_scuo)
    np.testing.assert_allclose(mcb(counts), 0.72)


def test_uniform_synonymous_frequencies_have_zero_distance():
    counts = codon_counts(["".join(CODONS) * 10] * 2)
    np.testing.assert_allclose(karlin_b(counts), 0)
    np.testing.assert_allclose(mcb(counts), 0)
    np.testing.assert_allclose(scuo(counts), 0, atol=1e-15)
    np.testing.assert_allclose(enc_prime(counts), 61)
    # Code 11 alternative starts remove isoleucine as a separate family.
    np.testing.assert_allclose(milc(counts), 0.5 - (64 - 20) / 640)


def test_default_alternative_starts_are_global_not_only_initial():
    counts = codon_counts(["ATT" * 100, "ATA" * 100])
    # ATT and ATA both join the seven-codon methionine family in code 11.
    np.testing.assert_allclose(milc(counts), 2 * math.log(2) - 6 / 100 + 0.5)
    # Turning that convention off restores the three-codon isoleucine family.
    np.testing.assert_allclose(
        milc(counts, alternative_starts=False), 2 * math.log(2) - 2 / 100 + 0.5
    )


def test_separate_reference_and_zero_reference_counts():
    counts = codon_counts(["AAA" * 75 + "AAG" * 25])
    reference = codon_counts(["AAA" * 50 + "AAG" * 50])
    expected = 2 * (0.75 * math.log(1.5) + 0.25 * math.log(0.5)) + 0.49
    assert milc(counts, reference_counts=reference)[0] == pytest.approx(expected)
    assert np.isinf(milc(counts, reference_counts=codon_counts(["AAA" * 100]))[0])


def test_encprime_ignores_stops_but_milc_keeps_them():
    counts = codon_counts(["".join(CODONS) * 10])
    without_stops = counts.copy()
    for codon in ("TAA", "TAG", "TGA"):
        without_stops[:, CODONS.index(codon)] = 0
    np.testing.assert_allclose(enc_prime(counts), enc_prime(without_stops))
    assert milc(counts)[0] != milc(without_stops)[0]


def test_genetic_codes_change_codon_families():
    counts = codon_counts(["TGA" * 50 + "TGG" * 50, "TGA" * 10 + "TGG" * 90])
    assert not np.allclose(milc(counts, genetic_code=4), milc(counts, genetic_code=11))
    with pytest.raises(ValueError, match="genetic code"):
        milc(counts, genetic_code=999)


def test_pair_bias_includes_junctions_and_unique_pair_normalization():
    # Concatenation produces AAA-AAG twice and the junction AAG-AAA once.
    assert codon_pair_bias(["AAAAAG", "AAAAAG"]) == pytest.approx(math.log(32 / 9))


def test_pair_bias_translation_matches_r_initial_codon_asymmetry():
    # Isolated codons translate as M, but only first TTG is M in the genome.
    assert codon_pair_bias(["TTGATGCTGGTG"]) == pytest.approx(1.5 * math.log(4 / 3))


@pytest.mark.parametrize("genes", [[], ["AAA"], ["AAA" * 50]])
def test_pair_bias_degenerate_inputs_return_nan(genes):
    assert math.isnan(codon_pair_bias(genes))


@pytest.mark.parametrize("counts", [np.zeros((2, 64)), np.ones((2, 63)), np.full((1, 64), -1)])
def test_invalid_count_matrices_raise(counts):
    with pytest.raises(ValueError):
        codon_usage(counts)


def test_invalid_metric_raises():
    with pytest.raises(ValueError, match="method"):
        codon_usage(codon_counts(["AAA"]), method="invalid")


@pytest.mark.parametrize("code", [True, False, 11.9, 11.0, None, float("nan")])
def test_genetic_code_rejects_booleans_and_float_truncation(code):
    with pytest.raises(ValueError, match="genetic_code"):
        milc(codon_counts(["AAA"]), genetic_code=code)
