from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from grodon.fitting import (
    estimate_boxcox_lambda,
    fit_gc_models,
    fit_linear_model,
    fit_models,
    fit_trim_models,
    fitGCModels,
    fitModels,
    fitTrimModels,
    get_statistics,
    get_statistics_batch,
)
from grodon.io import GeneSequences


@pytest.fixture
def training():
    rng = np.random.default_rng(846)
    data = pd.DataFrame(
        rng.normal(size=(100, 9)),
        columns=["CUBHE", "ConsistencyHE", "CPB", "OGT", "MILC", "ENCprime", "B", "SCUO", "MCB"],
    )
    data["CUB"] = rng.uniform(0.8, 2, size=len(data))
    data["GCdiv"] = rng.uniform(0, 0.3, size=len(data))
    data["d"] = np.exp(
        1
        + 0.3 * data.CUBHE
        + 0.2 * data.ConsistencyHE
        + 0.1 * data.CPB
        + rng.normal(scale=0.2, size=len(data))
    )
    return data


def test_fitted_coefficients_and_confidence_interval_against_closed_form():
    x = np.arange(1, 11, dtype=float)
    y = np.array([3.0, 3.2, 4.4, 4.8, 6.1, 7.4, 7.6, 8.5, 9.6, 10.4])
    model = fit_linear_model({"d": y, "x": x}, ("x",), boxcox_lambda=1)
    reference = stats.linregress(x, y - 1)
    np.testing.assert_allclose(model.coefficients, [reference.intercept, reference.slope])
    x0 = 5.5
    fit, lower, upper = model.predict_transformed({"x": x0})
    residuals = y - 1 - (reference.intercept + reference.slope * x)
    expected_se = np.sqrt(
        np.sum(residuals**2) / 8 * (1 / 10 + (x0 - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
    )
    expected_margin = stats.t.ppf(0.975, 8) * expected_se
    assert fit == pytest.approx(reference.intercept + reference.slope * x0)
    assert upper - fit == pytest.approx(expected_margin)
    assert fit - lower == pytest.approx(expected_margin)
    assert model.predict({"x": x0})[0] == pytest.approx(fit + 1)


def test_zero_lambda_uses_logarithmic_transform():
    x = np.arange(10, dtype=float)
    response = np.exp(1 + 0.1 * x)
    model = fit_linear_model({"d": response, "x": x}, ("x",), boxcox_lambda=0)
    np.testing.assert_allclose(model.coefficients, [1, 0.1], atol=1e-14)
    assert model.predict({"x": 4})[0] == pytest.approx(np.exp(1.4))


def test_all_model_families_and_predictor_choices(training):
    full = fit_models(training, training, lambda_value=0, eukaryotic_lambda=0)
    assert len(full) == 10
    assert full["gRodon_model_base"].predictors == ("(Intercept)", "CUBHE", "ConsistencyHE", "CPB")
    assert full["gRodon_model_partial_temp"].predictors == (
        "(Intercept)",
        "CUBHE",
        "ConsistencyHE",
        "OGT",
    )
    assert full["gRodon_model_meta"].predictors == ("(Intercept)", "CUBHE")
    trim = fit_trim_models(training, training, lambda_value=0)
    assert len(trim) == 3
    assert trim["gRodon_model_trim_temp"].predictors[-1] == "OGT"
    assert trim["gRodon_model_trim"].df_residual == len(training) - 6


def test_gc_models_fix_upstream_legacy_temperature_slot(training):
    original = training.copy(deep=True)
    result = fit_gc_models(training, training, lambda_value=0, legacy_lambda=1)
    assert result["meta_model_temp"].predictors == ("(Intercept)", "dCUB", "GCdiv", "OGT")
    assert result["meta_old_model_temp"].predictors == ("(Intercept)", "CUBHE", "OGT")
    assert result["meta_old_model_temp"].boxcox_lambda == 1
    pd.testing.assert_frame_equal(training, original)


def test_boxcox_search_and_r_list_wrappers(training):
    exponent = estimate_boxcox_lambda(training, ("CUBHE", "ConsistencyHE", "CPB"))
    assert exponent in np.linspace(-2, 2, 100)
    assert abs(exponent) < 0.5  # Synthetic log-linear response.
    assert fitModels(training, training)[6] == exponent
    assert len(fitTrimModels(training, training)) == 3
    legacy = fitGCModels(training, training)
    assert len(legacy) == 8
    assert legacy[6].predictors == ("(Intercept)", "CUBHE", "OGT")
    coarse = estimate_boxcox_lambda(training, ("CUBHE",), lambda_grid=[-1, 0, 1], interpolate=False)
    assert coarse in (-1, 0, 1)


def test_fitting_validation(training):
    with pytest.raises(ValueError, match="missing training columns"):
        fit_linear_model(training, ("absent",), boxcox_lambda=1)
    with pytest.raises(ValueError, match="rank deficient"):
        fit_linear_model(training.assign(x=1), ("x",), boxcox_lambda=1)
    with pytest.raises(ValueError, match="positive"):
        fit_linear_model(training.assign(d=0), ("CUBHE",), boxcox_lambda=1)
    with pytest.raises(ValueError, match="finite"):
        fit_linear_model(training.assign(d=np.nan), ("CUBHE",), boxcox_lambda=1)
    with pytest.raises(ValueError, match="more observations"):
        fit_linear_model(training.iloc[:2], ("CUBHE",), boxcox_lambda=1)
    with pytest.raises(ValueError, match="strictly increasing"):
        estimate_boxcox_lambda(training, ("CUBHE",), lambda_grid=[1, 0])
    with pytest.raises(ValueError, match="nonzero"):
        fit_gc_models(training.assign(CUB=0), training)


def test_boxcox_extreme_range_raises_explicit_error_without_runtime_warnings():
    import warnings

    with warnings.catch_warnings(record=True) as caught:
        with pytest.raises(ValueError, match="Box-Cox.*overflowed"):
            estimate_boxcox_lambda({"d": [1e-300, 1e-200, 1e200, 1e300], "x": [1, 2, 3, 4]}, ("x",))
    assert not caught


def test_annotated_statistics_uses_original_r_pattern(monkeypatch):
    observed = {}

    def mock_statistics(genes, he, **kwargs):
        observed.update(genes=genes, he=he, kwargs=kwargs)
        return {"CUBHE": 0.3}

    monkeypatch.setattr("grodon.statistics.get_codon_statistics", mock_statistics)
    genes = GeneSequences(
        ("ATG",) * 5,
        (
            "g1 30S ribosomal protein S1",
            "g2 50S ribosomal protein methyltransferase",
            "g3 hydroxylase 50S ribosomal protein",
            "g4 hypothetical protein",
            "g5 50s RIBOSOMAL PROTEIN L1",
        ),
    )
    result = get_statistics(genes=genes, bg="individual")
    assert observed["he"] == [True, False, False, False, True]
    assert observed["kwargs"]["background"] == "individual"
    assert result["File"] == "genome provided directly"
    with pytest.raises(ValueError, match="no highly expressed"):
        get_statistics(genes={"hypothetical protein": "ATG"})


def test_statistics_file_path_supports_gzip(tmp_path, monkeypatch):
    import gzip

    path = tmp_path / "genes.fna.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">gene1 30S ribosomal protein S1\nATGAAA\n>gene2 other protein\nATGCCC\n")
    monkeypatch.setattr(
        "grodon.statistics.get_codon_statistics",
        lambda genes, he, **kwargs: {"nHE": sum(he), "nGenes": len(genes)},
    )
    result = get_statistics(path)
    assert result == {"nHE": 1, "nGenes": 2, "File": "genes.fna.gz"}


def test_batch_sorting_parallel_seed_reproducibility_and_failures(tmp_path, monkeypatch):
    for name in ["b.fna", "a.fna", "notes.txt"]:
        (tmp_path / name).touch()

    def mocked(path, *, random_state=None, **kwargs):
        if Path(path).suffix == ".txt":
            raise ValueError("invalid FASTA")
        return {"File": Path(path).name, "random": np.random.default_rng(random_state).random()}

    monkeypatch.setattr("grodon.fitting.get_statistics", mocked)
    serial = get_statistics_batch(tmp_path, pattern="*.fna", random_state=7)
    parallel = get_statistics_batch(tmp_path, pattern="*.fna", random_state=7, workers=2)
    pd.testing.assert_frame_equal(serial, parallel)
    assert serial.File.tolist() == ["a.fna", "b.fna"]
    with pytest.raises(ValueError, match="invalid FASTA"):
        get_statistics_batch(tmp_path)
    with pytest.warns(RuntimeWarning, match="skipping notes.txt"):
        skipped = get_statistics_batch(tmp_path, errors="skip")
    assert len(skipped) == 2
    with pytest.raises(ValueError, match="workers"):
        get_statistics_batch(tmp_path, workers=0)
