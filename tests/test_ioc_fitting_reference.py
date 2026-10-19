"""Numerical fixtures executed by R, using synthetic traits and training data.

These tests need no R installation; regeneration uses the checked-in R script.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from grodon.copiotrophy import fit_ioc, predict_ioc
from grodon.fitting import fit_gc_models, fit_models, fit_trim_models

REFERENCE = Path(__file__).parent / "reference" / "ioc_fitting"


def test_ioc_fit_predict_and_completeness_match_r_execution():
    traits = pd.read_csv(REFERENCE / "ioc_traits.csv")
    heldout = pd.read_csv(REFERENCE / "ioc_heldout.csv")
    fitted_reference = pd.read_csv(REFERENCE / "ioc_fit.csv")
    predictions = pd.read_csv(REFERENCE / "ioc_predictions.csv")
    pca = pd.read_csv(REFERENCE / "ioc_pca.csv")
    fitted = fit_ioc(traits)
    np.testing.assert_allclose(fitted.values, fitted_reference.IoC, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(fitted.model.center, pca.center, rtol=1e-13)
    np.testing.assert_allclose(fitted.model.scale, pca.scale, rtol=1e-13)
    np.testing.assert_allclose(fitted.model.sdev**2, pca.variance, rtol=1e-12)
    complete = fit_ioc(traits, completeness=[90, 100, 50, 75, 100, 80])
    np.testing.assert_allclose(complete.values, fitted_reference.complete_IoC, rtol=1e-13)
    np.testing.assert_allclose(predict_ioc(heldout, model=fitted), predictions.custom, rtol=1e-13)
    for name in predictions.columns.drop("custom"):
        np.testing.assert_allclose(predict_ioc(heldout, model=name), predictions[name], rtol=1e-13)


@pytest.fixture(scope="module")
def models_and_data():
    training = pd.read_csv(REFERENCE / "training.csv")
    extremo = pd.read_csv(REFERENCE / "training_extremo.csv")
    models = {
        **fit_models(training, extremo),
        **fit_gc_models(training, extremo),
        **fit_trim_models(training, extremo),
    }
    training["dCUB"] = (training.CUB - training.CUBHE) / training.CUB
    return models, training


def test_all_fitted_models_coefficients_covariance_and_lambda_match_r(models_and_data):
    models, _ = models_and_data
    reference = pd.read_csv(REFERENCE / "fitted_parameters.csv", keep_default_na=False)
    for name, entries in reference.groupby("model"):
        fitted = models[name]
        params = entries.loc[entries.kind == "parameter"].set_index("row").value
        assert fitted.boxcox_lambda == pytest.approx(params["lambda"], abs=1e-14), name
        assert fitted.df_residual == params["df"], name
        assert fitted.sigma == pytest.approx(params["sigma"], rel=1e-11), name
        np.testing.assert_allclose(
            fitted.coefficients,
            params[list(fitted.predictors)],
            rtol=1e-10,
            atol=1e-10,
            err_msg=name,
        )
        covariance = (
            entries.loc[entries.kind == "covariance"]
            .pivot(index="row", columns="column", values="value")
            .loc[list(fitted.predictors), list(fitted.predictors)]
        )
        np.testing.assert_allclose(
            fitted.covariance, covariance, rtol=1e-9, atol=1e-10, err_msg=name
        )


def test_fitted_models_confidence_intervals_match_r(models_and_data):
    models, training = models_and_data
    reference = pd.read_csv(REFERENCE / "fitted_predictions.csv")
    for row in reference.itertuples():
        actual = models[row.model].predict_transformed(training.iloc[row.observation - 1])
        # scipy 1.11's t.ppf differs from R qt by ~1.8e-9 at 53 df;
        # this changes these confidence endpoints by up to ~2.7e-10.
        # Coefficients and covariance are checked more tightly above.
        np.testing.assert_allclose(
            actual, [row.fit, row.lwr, row.upr], rtol=1e-9, atol=1e-9, err_msg=row.model
        )
