"""Prediction dispatch, uncertainty, model provenance, and regression checks."""

import dataclasses
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import t

from grodon.models import (
    LinearModel,
    available_models,
    boxcox_transform,
    load_model,
    model_provenance,
)
from grodon.prediction import (
    PREDICTION_MODES,
    GrowthPredictionWarning,
    predict_from_statistics,
    predict_growth,
    predictGrowth,
)


@pytest.fixture
def statistics():
    return {
        "CUBHE": 0.8,
        "ConsistencyHE": 0.5,
        "CPB": 0.02,
        "CUB": 0.4,
        "GCdiv": 0.1,
        "MILC": 0.8,
        "ENCprime": 40.0,
        "B": 0.3,
        "SCUO": 0.4,
        "MCB": 0.2,
        "FilteredSequences": 0,
    }


def test_all_original_models_are_bundled():
    names = available_models()
    assert len(names) == 36
    assert len(names) == len(set(names))
    for name in names:
        model = load_model(name)
        assert model.name == name
        assert model.predictors[0] == "(Intercept)"
        assert model.df_residual > 0
    assert model_provenance()["git_revision"] == "4b13053e7f49ca99bceeee36940f5e820c5b0c32"
    assert (
        model_provenance()["sha256"]
        == "befeda83bd205e3a9b29d96ea7a260031e549abd6b5aefa34e6da735664138a9"
    )


@pytest.mark.parametrize("name", available_models())
def test_every_model_matches_independent_r_predictions(name):
    """Five fixed feature rows per model, calculated by unmodified predict.lm."""
    references = pd.read_csv(Path(__file__).parent / "data" / "reference_growth_models.csv")
    rows = references[references["model"] == name]
    assert len(rows) == 5
    model = load_model(name)
    for _, row in rows.iterrows():
        np.testing.assert_allclose(
            model.predict_transformed(row),
            row[["fit", "lower", "upper"]].to_numpy(dtype=float),
            rtol=2e-11,
            atol=2e-12,
        )
        np.testing.assert_allclose(
            model.predict(row),
            row[["hours", "hours_lower", "hours_upper"]].to_numpy(dtype=float),
            rtol=5e-11,
            atol=2e-12,
            equal_nan=True,
        )
        assert model.boxcox_lambda == pytest.approx(row["lambda"], abs=1e-15)
        assert model.df_residual == row["df_residual"]


@pytest.mark.parametrize(
    "training_set,suffix",
    [("vs", ""), ("original", ""), ("madin", "_madin"), ("AOA_NOB", "_AOANOB")],
)
@pytest.mark.parametrize(
    "mode,stem", [("full", "base"), ("partial", "partial"), ("metagenome_v1", "meta")]
)
@pytest.mark.parametrize("temperature", [None, 30.0])
def test_training_set_dispatch(statistics, training_set, suffix, mode, stem, temperature):
    result = predict_from_statistics(
        statistics, mode=mode, training_set=training_set, temperature=temperature
    )
    if temperature is not None:
        stem = "temp" if stem == "base" else f"{stem}_temp"
    assert result.model_name == f"gRodon_model_{stem}{suffix}"
    assert result.lower_ci < result.doubling_time < result.upper_ci


@pytest.mark.parametrize("temperature", [None, 30.0])
@pytest.mark.parametrize("mode", ["eukaryote", "metagenome_euk"])
@pytest.mark.parametrize("background", ["all", "individual"])
def test_eukaryotes_use_their_own_transform(statistics, temperature, mode, background):
    result = predict_from_statistics(
        statistics, mode=mode, temperature=temperature, background=background
    )
    stem = "base" if temperature is None else "temp"
    name = f"gRodon_model_{stem}_euk" + ("_i" if background == "individual" else "")
    assert result.model_name == name
    assert result.training_set == "eukaryote"
    model = load_model(name)
    row = dict(statistics, OGT=30.0, dCUB=-1.0)
    assert result.doubling_time == pytest.approx(model.predict(row)[0])
    wrong_lambda = load_model("gRodon_model_base_madin").boxcox_lambda
    assert result.doubling_time != pytest.approx(
        boxcox_transform(model.predict_transformed(row)[0], wrong_lambda, back_transform=True)
    )


@pytest.mark.parametrize(
    "mode,stem", [("meta_testing", "newmeta"), ("meta_nogc_testing", "newmeta_nogc")]
)
@pytest.mark.parametrize("background", ["all", "individual"])
@pytest.mark.parametrize("temperature", [None, 25.0])
def test_experimental_model_dispatch(statistics, mode, stem, background, temperature):
    result = predict_from_statistics(
        statistics, mode=mode, temperature=temperature, background=background
    )
    suffix = "_temp" if temperature is not None else ""
    assert result.model_name == f"gRodon_model_{stem}{suffix}" + (
        "_i" if background == "individual" else ""
    )
    assert result["dCUB"] == -1


@pytest.mark.parametrize("mode,trimlen", [("metagenome_150bp", 150), ("metagenome_240bp", 240)])
@pytest.mark.parametrize("temperature", [None, 30.0])
def test_short_read_dispatch(statistics, mode, trimlen, temperature):
    result = predict_from_statistics(statistics, mode=mode, temperature=temperature)
    assert result.background == "individual"
    assert (
        result.model_name == f"gRodon_model_{'base' if temperature is None else 'temp'}_t{trimlen}"
    )


def test_metagenome_v2_requires_matching_background(statistics):
    pooled = predict_from_statistics(statistics, mode="metagenome_v2")
    assert pooled.model_name == "gRodon_model_meta_madin"
    statistics["ConsistencyHE"] = 0.6
    with pytest.raises(ValueError, match="recomputing"):
        predict_from_statistics(statistics, mode="metagenome_v2")
    individual = predict_from_statistics(statistics, mode="metagenome_v2", background="individual")
    assert individual.model_name == "gRodon_model_newmeta_i"


def test_original_training_sets_are_distinct(statistics):
    predictions = [
        predict_from_statistics(statistics, training_set=name).d
        for name in ("vs", "madin", "AOA_NOB")
    ]
    assert len(set(predictions)) == 3


def test_confidence_interval_is_for_mean_not_individual_observation():
    model = LinearModel(
        "synthetic", ("(Intercept)", "x"), (2.0, 3.0), ((0.25, 0.0), (0.0, 0.25)), 20, 10.0, 0.0
    )
    fit, lower, upper = model.predict_transformed({"x": 0.0})
    assert fit == 2.0
    assert upper - fit == pytest.approx(t.ppf(0.975, 20) * 0.5)
    assert fit - lower == pytest.approx(upper - fit)
    assert model.predict({"x": 0.0}) == pytest.approx(np.exp([fit, lower, upper]))


def test_wider_confidence_level_widens_interval(statistics):
    narrow = predict_from_statistics(statistics, confidence=0.8)
    wide = predict_from_statistics(statistics, confidence=0.99)
    assert narrow.d == wide.d
    assert wide.lower_ci < narrow.lower_ci < narrow.upper_ci < wide.upper_ci


def test_result_is_mapping_compatible_and_immutable(statistics):
    original = dict(statistics)
    result = predict_from_statistics(statistics)
    assert statistics == original
    assert result["d"] == result.doubling_time == result.d
    assert dict(result)["LowerCI"] == result.lower_ci
    assert result.to_dict() == dict(result)
    assert result == dict(result)
    with pytest.raises(TypeError):
        result.statistics["d"] = 10.0
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.model_name = "different"


@pytest.mark.parametrize("lam", [-0.1414141414, 0.0, 0.3])
def test_boxcox_roundtrip_including_logarithm(lam):
    values = np.array([0.01, 0.5, 1.0, 10.0, 1000.0])
    transformed = boxcox_transform(values, lam)
    np.testing.assert_allclose(
        boxcox_transform(transformed, lam, back_transform=True), values, rtol=1e-12
    )
    if lam == 0:
        np.testing.assert_allclose(transformed, np.log(values))


def test_boxcox_invalid_domain_is_not_fabricated():
    assert np.isnan(boxcox_transform(3.0, -0.5, back_transform=True))
    with pytest.raises(ValueError, match="positive"):
        boxcox_transform([0.0, 1.0], -0.5)


@pytest.mark.parametrize("degrees", [0, -1, 1.5, float("inf"), float("nan"), True])
def test_model_rejects_invalid_residual_degrees_of_freedom(degrees):
    with pytest.raises(ValueError, match="degrees of freedom"):
        LinearModel("invalid", ("(Intercept)",), (1.0,), ((1.0,),), degrees, 1.0, 0.0)


@pytest.mark.parametrize("genes", ["ATG" * 100, b"ATG" * 100, "genes.fna"])
def test_prediction_rejects_bare_strings_before_treating_bases_as_genes(genes):
    with pytest.raises(TypeError, match="collection of sequences"):
        predict_growth(genes, [True])


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"mode": "unknown"}, "Invalid mode"),
        ({"training_set": "unknown"}, "training_set"),
        ({"background": "unknown"}, "background"),
        ({"background": "individual"}, "no model"),
        ({"confidence": 1.0}, "confidence"),
        ({"temperature": np.inf}, "finite"),
        ({"temperature": True}, "temperature"),
    ],
)
def test_invalid_prediction_options(statistics, kwargs, match):
    with pytest.raises(ValueError, match=match):
        predict_from_statistics(statistics, **kwargs)


def test_missing_and_nonfinite_features_fail_clearly(statistics):
    del statistics["CPB"]
    with pytest.raises(ValueError, match="CPB"):
        predict_from_statistics(statistics)
    assert np.isfinite(predict_from_statistics(statistics, mode="partial").d)
    statistics["CUBHE"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        predict_from_statistics(statistics, mode="partial")


def test_invalid_model_name():
    with pytest.raises(ValueError, match="Unknown model"):
        load_model("unknown")


def test_metagenome_pipeline_recalculates_background(monkeypatch, statistics):
    import grodon.statistics as statistics_module

    calls = []

    def fake_statistics(genes, highly_expressed, **kwargs):
        calls.append(kwargs)
        return dict(statistics, ConsistencyHE=0.7)

    monkeypatch.setattr(statistics_module, "get_codon_statistics", fake_statistics)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GrowthPredictionWarning)
        result = predict_growth(
            ["ATG" * 100] * 12,
            [True] * 10 + [False] * 2,
            mode="metagenome_v2",
            random_state=42,
            n_le=7,
        )
    assert [call["background"] for call in calls] == ["all", "individual"]
    assert calls[1]["random_state"] == 42
    assert calls[1]["n_le"] == 7
    assert result.model_name == "gRodon_model_newmeta_i"


@pytest.mark.parametrize("mode,trimlen", [("metagenome_150bp", 150), ("metagenome_240bp", 240)])
def test_short_read_pipeline_preprocessing(monkeypatch, statistics, mode, trimlen):
    import grodon.statistics as statistics_module

    calls = []

    def fake_statistics(genes, highly_expressed, **kwargs):
        calls.append(kwargs)
        return statistics

    monkeypatch.setattr(statistics_module, "get_codon_statistics", fake_statistics)
    result = predict_growth(
        ["ATG" * 100] * 12, [True] * 10 + [False] * 2, mode=mode, exclude_short=True, random_state=3
    )
    assert result.background == "individual"
    assert calls[0]["fragments"] is None
    assert calls[0]["trimlen"] == trimlen
    assert calls[0]["exclude_short"] is True
    assert calls[0]["all_metrics"] is True


def test_individual_metagenome_v1_is_reachable(statistics):
    result = predict_from_statistics(statistics, mode="metagenome_v1", background="individual")
    assert result.model_name == "gRodon_model_meta_madin_i"


@pytest.fixture
def synthetic_sequences():
    """Random codon strings with no natural organism or sequence as a source."""
    from grodon.codon import CODONS

    rng = np.random.default_rng(8531)
    sense = [codon for codon in CODONS if codon not in {"TAA", "TAG", "TGA"}]
    weights = np.arange(1, len(sense) + 1, dtype=float)
    weights /= weights.sum()
    genes = [
        "ATG" + "".join(rng.choice(sense, size=160, p=weights if index < 12 else None)) + "TAA"
        for index in range(40)
    ]
    return genes, [index < 12 for index in range(40)]


@pytest.mark.parametrize("mode", PREDICTION_MODES)
def test_complete_sequence_to_prediction_pipeline(synthetic_sequences, mode):
    genes, he = synthetic_sequences
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GrowthPredictionWarning)
        result = predict_growth(genes, he, mode=mode, temperature=30, random_state=11, n_le=8)
    assert np.isfinite([result.d, result.lower_ci, result.upper_ci]).all()
    assert 0 < result.lower_ci < result.d < result.upper_ci
    assert result["nHE"] == 12
    assert result.mode == mode


def test_individual_sequence_predictions_are_reproducible(synthetic_sequences):
    genes, he = synthetic_sequences
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GrowthPredictionWarning)
        first = predict_growth(
            iter(genes),
            iter(he),
            mode="metagenome_v2",
            background="individual",
            random_state=29,
            n_le=8,
        )
        second = predict_growth(
            dict(enumerate(genes)),
            he,
            mode="metagenome_v2",
            background="individual",
            random_state=29,
            n_le=8,
        )
    assert first.d == second.d
    assert first.lower_ci == second.lower_ci
    assert first["CUBHE"] == second["CUBHE"]


def test_r_interface_positional_arguments(synthetic_sequences):
    genes, he = synthetic_sequences
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GrowthPredictionWarning)
        snake = predict_growth(genes, he, mode="partial", temperature=30, training_set="vs")
        camel = predictGrowth(genes, he, "partial", 30, "vs")
    assert camel.d == snake.d
    assert camel.model_name == "gRodon_model_partial_temp"
