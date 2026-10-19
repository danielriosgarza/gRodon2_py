import json

import numpy as np
import pandas as pd
import pytest

from grodon.copiotrophy import (
    IoCModel,
    available_ioc_models,
    fit_ioc,
    fitIoC,
    load_ioc_model,
    predict_ioc,
    predictIoC,
)


@pytest.fixture
def traits():
    return pd.DataFrame(
        {
            "dCUB": [-0.1, -0.2, -0.5, -0.3, -0.04, -0.7],
            "nCAZy": [10, 50, 80, 45, 3, 130],
            "nGenes": [1500, 2200, 4000, 3000, 800, 5100],
        }
    )


def test_fitted_training_indices_reproduce_and_have_fixed_origin(traits):
    fitted = fit_ioc(traits)
    np.testing.assert_allclose(fitted.values, predict_ioc(traits, model=fitted), atol=1e-14)
    assert fitted.values.min() == 1
    assert fitted.model.rotation[1, 0] * fitted.model.direction < 0
    assert fitted.model.variance_explained.sum() == pytest.approx(1)
    assert dict(fitted).keys() == {"IoC", "model"}


def test_pca_variance_matches_correlation_eigenvalues(traits):
    fitted = fit_ioc(traits)
    raw = np.column_stack((traits.nCAZy / traits.nGenes, traits.dCUB, traits.nGenes))
    expected = np.linalg.eigvalsh(np.corrcoef(raw.T))[::-1]
    np.testing.assert_allclose(fitted.model.sdev**2, expected, atol=1e-14)
    np.testing.assert_allclose(np.cov(fitted.model.x.T), np.diag(expected), atol=1e-14)


def test_completeness_preserves_cazy_fraction(traits):
    complete = fit_ioc(traits)
    incomplete = fit_ioc(traits, completeness=50)
    np.testing.assert_allclose(complete.values, incomplete.values, atol=1e-13)
    np.testing.assert_allclose(incomplete.model.center, complete.model.center * [1, 1, 2])
    np.testing.assert_allclose(incomplete.model.scale, complete.model.scale * [1, 1, 2])


def test_completeness_vector_is_applied_per_row(traits):
    completeness = np.array([90, 100, 50, 75, 100, 80])
    fitted = fit_ioc(traits, completeness=completeness)
    features = np.column_stack(
        (traits.nCAZy / traits.nGenes, traits.dCUB, traits.nGenes / (completeness / 100))
    )
    np.testing.assert_allclose(fitted.model.center, features.mean(axis=0))
    np.testing.assert_allclose(
        predict_ioc(traits, model=fitted, completeness=completeness), fitted.values
    )


def test_model_json_roundtrip_and_r_wrappers(traits):
    fitted = fitIoC(dCUB=traits.dCUB, nCAZy=traits.nCAZy, nGenes=traits.nGenes)
    portable = json.loads(json.dumps(fitted.model.to_dict(include_scores=True)))
    restored = IoCModel.from_dict(portable)
    np.testing.assert_allclose(predictIoC(traits, model=restored), fitted["IoC"])
    np.testing.assert_allclose(predict_ioc(traits, model=portable), fitted.values)
    assert not restored.rotation.flags.writeable


@pytest.mark.parametrize(
    ("model", "expected"),
    [
        ("soil", [4.08567438047225, 6.597473850375349, 1.7759270403046217]),
        ("permafrost", [3.0350574798655394, 5.799359200797587, 0.818234831923514]),
        ("pacific", [2.9796760762054975, 4.833526936939097, 1.8661161681315923]),
        ("human", [2.8149265160222106, 5.445159002286746, 0.6326342816496546]),
        ("all_habitat", [7.360090001389548, 8.136600812611634, 6.456018831565352]),
    ],
)
def test_bundled_model_reference_predictions(model, expected):
    # Independently evaluated from the unmodified R prcomp objects in
    # inst/extdata/IoC_Models.rda (gRodon 2.7.3), including their training scores.
    # These inputs are synthetic trait rows, not selected organisms.
    actual = predict_ioc(
        model=model, dcub=[-0.2, -0.5, -0.01], n_cazy=[30, 90, 0], n_genes=[2000, 4000, 500]
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)
    # A held-out batch must not move the fitted index's origin.
    np.testing.assert_allclose(
        predict_ioc(model=model, dcub=-0.2, n_cazy=30, n_genes=2000), actual[:1]
    )


@pytest.mark.parametrize("completeness", [0, -1, 101, np.nan, [50, 60], [[90]]])
def test_invalid_completeness(traits, completeness):
    with pytest.raises(ValueError, match="completeness"):
        fit_ioc(traits, completeness=completeness)


def test_invalid_inputs(traits):
    with pytest.raises(ValueError, match="constant"):
        fit_ioc(traits.assign(dCUB=1))
    with pytest.raises(ValueError, match="two observations"):
        fit_ioc(traits.iloc[:1])
    with pytest.raises(ValueError, match="finite"):
        fit_ioc(traits.assign(nCAZy=np.nan))
    with pytest.raises(ValueError, match="positive"):
        fit_ioc(traits.assign(nGenes=0))
    with pytest.raises(ValueError, match="same length"):
        fit_ioc(dcub=[1, 2], n_cazy=[1], n_genes=[10, 20])
    with pytest.raises(ValueError, match="both"):
        fit_ioc(traits, dcub=[1, 2])
    with pytest.raises(ValueError, match="unknown IoC model"):
        load_ioc_model("unknown")
    with pytest.raises(TypeError, match="model must"):
        predict_ioc(traits, model=42)


def test_bundled_models_are_complete():
    assert set(available_ioc_models()) == {"permafrost", "pacific", "soil", "human", "all_habitat"}
    assert load_ioc_model("permafrost").n_observations == 533
    assert load_ioc_model("pacific").n_observations == 6865
