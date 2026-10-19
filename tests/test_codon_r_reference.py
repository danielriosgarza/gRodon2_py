"""Golden values produced by the original R package, not a Python reimplementation.

Regenerate with scripts/generate_codon_reference.R. R is not required to run the
tests. Fixtures contain synthetic sequences and record exact upstream versions.
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pytest

from grodon.codon import codon_counts, codon_pair_bias, codon_usage
from grodon.prediction import predict_growth
from grodon.statistics import get_codon_statistics, weighted_median

REFERENCE = json.loads((Path(__file__).parent / "fixtures" / "codon_reference.json").read_text())


def assert_reference(actual, expected):
    assert actual.keys() == expected.keys()
    for key, value in expected.items():
        if value is None:
            assert np.isnan(actual[key]), key
        else:
            assert actual[key] == pytest.approx(value, rel=2e-12, abs=5e-14), key


@pytest.mark.parametrize("case", REFERENCE["cases"], ids=lambda case: case["name"])
def test_genome_statistics_match_original_r(case):
    options = case["options"] or {}
    with (
        pytest.warns(UserWarning)
        if case["statistics"]["FilteredSequences"]
        else _no_warning_context()
    ):
        actual = get_codon_statistics(case["genes"], case["highly_expressed"], **options)
    assert_reference(actual, case["statistics"])


@pytest.mark.parametrize("case", REFERENCE["cases"], ids=lambda case: case["name"])
def test_all_metrics_match_original_r(case):
    options = case["options"] or {}
    with (
        pytest.warns(UserWarning)
        if case["statistics"]["FilteredSequences"]
        else _no_warning_context()
    ):
        actual = get_codon_statistics(
            case["genes"], case["highly_expressed"], all_metrics=True, **options
        )
    assert_reference(actual, case["all_metrics"])


def _no_warning_context():
    from contextlib import nullcontext

    return nullcontext()


@pytest.mark.parametrize("case", REFERENCE["low_level"])
def test_per_gene_metrics_match_cordon_r(case):
    counts = codon_counts(case["genes"])
    for method, expected in case["expected"].items():
        actual = codon_usage(
            counts,
            method=method,
            genetic_code=case["genetic_code"],
            alternative_starts=case["alternative_starts"],
        )
        np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=5e-14, err_msg=method)


@pytest.mark.parametrize("case", REFERENCE["fixed_backgrounds"])
def test_individual_metrics_match_r_with_identical_shuffled_background(case):
    counts = codon_counts([case["gene"]])
    for method, expected in case["expected"].items():
        actual = codon_usage(
            counts, method=method, reference_counts=case["reference_counts"], genetic_code=1
        )
        assert actual[0] == pytest.approx(expected, rel=2e-12, abs=5e-14), method


@pytest.mark.parametrize("case", REFERENCE["pair_cases"])
def test_pair_bias_matches_original_r(case):
    genes = case["genes"]
    genes = [genes] if isinstance(genes, str) else genes
    actual = codon_pair_bias(genes, genetic_code=case["genetic_code"])
    assert actual == pytest.approx(case["expected"], rel=2e-12, abs=5e-14)


@pytest.mark.parametrize("case", REFERENCE["weighted_medians"])
def test_weighted_median_matches_matrixstats_r(case):
    assert weighted_median(case["values"], case["weights"]) == pytest.approx(
        case["expected"], rel=2e-12
    )


@pytest.mark.parametrize("case", REFERENCE["predictions"], ids=lambda case: case["name"])
def test_complete_prediction_pipeline_matches_original_r(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        actual = predict_growth(case["genes"], case["highly_expressed"], **case["options"])
    assert_reference(dict(actual), case["expected"])


def test_r_fixture_records_provenance():
    provenance = REFERENCE["provenance"]
    assert provenance["R"].startswith("R version ")
    assert provenance["upstream_version"] == "2.7.3"
    assert "coRdon" in provenance["packages"]
    assert set(provenance["source_md5"]) == {
        "utilities.R",
        "codon_pair_bias.R",
        "codon_usage_bias.R",
        "predict_growth.R",
    }
    assert len(provenance["model_data_md5"]) == 32
