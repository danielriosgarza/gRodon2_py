# Copiotrophy API

[Back to the API index](API.md)

The Index of Copiotrophy (IoC) combines three supplied traits using a standardized
principal component analysis (PCA): CAZyme fraction, differential codon usage bias
(`dCUB`), and completeness-corrected gene count. This API fits an index to a table
or places new observations on an existing index. It does not annotate CAZymes,
read genomes, calculate `dCUB` from sequences, or return a growth rate.

All names documented here are available from both `grodon` and
`grodon.copiotrophy`. The signatures use `numpy.typing.ArrayLike` (`npt.ArrayLike`),
`numpy.typing.NDArray` (`npt.NDArray`), `pandas.DataFrame` (`pd.DataFrame`), and
`collections.abc.Mapping`.

## Fit and predict with synthetic traits

The following examples use artificial numeric traits and require no sequence
files, R installation, network access, or optional genome dependencies. The
examples in this document can run in order in the same Python session.

```python
import json

import numpy as np
import pandas as pd

from grodon import (
    IoCModel,
    available_ioc_models,
    fit_ioc,
    fitIoC,
    load_ioc_model,
    predict_ioc,
    predictIoC,
)

traits = pd.DataFrame(
    {
        "dCUB": [-0.1, -0.2, -0.5, -0.3, -0.04, -0.7],
        "nCAZy": [10, 50, 80, 45, 3, 130],
        "nGenes": [1500, 2200, 4000, 3000, 800, 5100],
    }
)
completeness = np.array([90, 100, 85, 95, 100, 80])

fitted = fit_ioc(traits, completeness=completeness)
scores = predict_ioc(traits, model=fitted, completeness=completeness)

assert scores.shape == (6,)
np.testing.assert_allclose(scores, fitted.values)
assert fitted.values.min() == 1.0

new_traits = {
    "dCUB": [-0.25, -0.55],
    "nCAZy": [40, 100],
    "nGenes": [2500, 4200],
}
new_scores = predict_ioc(new_traits, model=fitted.model, completeness=[95, 90])
assert new_scores.shape == (2,)
```

The model retains its training center, scales, PCA loadings, and index origin.
Predicting a new batch does not refit those parameters or shift the batch minimum
to one. New scores can be less than one, including negative values; they are not
clipped. Scores from separately fitted models have different reference frames.

## Input contract and completeness

`fit_ioc` and `predict_ioc` accept either a table/mapping in `x`, or all three
keyword vectors. Mixing the two forms raises `ValueError`.

| Input | Table/mapping key | Vector keyword | Requirements |
| --- | --- | --- | --- |
| Differential codon usage bias | `dCUB` | `dcub` | Finite numeric values; no range restriction |
| Number of CAZyme genes | `nCAZy` | `n_cazy` | Finite, nonnegative numeric values |
| Observed gene count | `nGenes` | `n_genes` | Finite, strictly positive numeric values |
| Assembly completeness | Separate argument | `completeness` | Finite percentages greater than zero and at most 100 |

Column names are case-sensitive. Extra columns in `x` are ignored, including a
column named `completeness`; pass that column explicitly as the `completeness`
argument. A plain two-dimensional NumPy array is not a supported `x` table because
it has no named columns. Supply its columns through the three vector keywords.

Feature vectors must be one-dimensional, nonempty, and equally long. Lists,
tuples, NumPy arrays, and pandas Series are accepted when their values can be
converted to floating point. Scalars become length-one vectors. The functions
do not align independently supplied Series by index: they use their values in
positional order. Output arrays also preserve row order, not pandas index labels.
Missing values (`NaN`, `None` converted to `NaN`) and infinities are rejected.
The implementation accepts numeric counts without requiring integer values and
does not enforce `nCAZy <= nGenes`; callers are responsible for those data semantics.

`completeness` accepts a scalar, a length-one vector, or one value per observation.
It is a percentage: use `95` for 95% completeness. A value of `0.95` means 0.95%.
No clipping or missing-value imputation is applied.

For each row, feature preparation is exactly:

```text
rCAZy = nCAZy / nGenes
corrected_nGenes = nGenes / (completeness / 100)
features = [rCAZy, dCUB, corrected_nGenes]
```

The CAZyme fraction is calculated from the observed counts before correcting
gene count. Completeness does not change `rCAZy` or `dCUB`. When predicting, supply
the completeness appropriate to the observations being predicted. To reproduce
fitted training values, use the same training completeness values again.

## `fit_ioc`

```text
fit_ioc(
    x: pd.DataFrame | Mapping[str, Any] | None = None,
    *,
    dcub: npt.ArrayLike | None = None,
    n_cazy: npt.ArrayLike | None = None,
    n_genes: npt.ArrayLike | None = None,
    completeness: npt.ArrayLike = 100,
) -> IoCFit
```

`x` may be positional; the remaining parameters are keyword-only. See the shared
input contract above for meanings and validation.

Fitting requires at least two observations and nonzero variance in each prepared
feature. In particular, varying `nCAZy` and `nGenes` is insufficient if their ratio
is constant. Gene-count variance is checked after completeness correction.

The function centers the three prepared columns and scales them by their sample
standard deviations (`ddof=1`). It computes PCA using NumPy's singular value
decomposition. Features are rows of the resulting rotation matrix; principal
components are columns. The number of returned components is `min(n, 3)`, where
`n` is the number of training observations.

For the index, PC1 scores are oriented so that the effective `dCUB` loading
is negative. The stored rotation and PCA coordinates keep their original
orientation. The fitted index is:

```text
direction = -sign(rotation[1, 0])
oriented_PC1 = PCA_scores[:, 0] * direction
min_score = min(oriented_PC1 on the training observations)
IoC = oriented_PC1 - min_score + 1
```

The return value is an [`IoCFit`](#iocfit), containing a length-`n` index array and
the fitted [`IoCModel`](#iocmodel). Custom fitted models have `name="custom"`.

In addition to invalid inputs, `ValueError` is raised for too few observations,
constant prepared features, or a zero PC1 `dCUB` loading, whose orientation is
undefined. A numerical SVD failure can propagate as `numpy.linalg.LinAlgError`.

## `predict_ioc`

```text
predict_ioc(
    x: pd.DataFrame | Mapping[str, Any] | None = None,
    *,
    model: str | IoCModel | IoCFit | Mapping[str, Any],
    dcub: npt.ArrayLike | None = None,
    n_cazy: npt.ArrayLike | None = None,
    n_genes: npt.ArrayLike | None = None,
    completeness: npt.ArrayLike = 100,
) -> npt.NDArray[np.float64]
```

`model` is required and keyword-only. Input feature preparation follows the same
contract as `fit_ioc`, but one observation is sufficient and prediction columns
may be constant. Predictions use the supplied model's stored training parameters:

```text
PCA_scores = ((features - model.center) / model.scale) @ model.rotation
IoC = PCA_scores[:, 0] * model.direction - model.min_score + 1
```

| Accepted `model` value | Behavior |
| --- | --- |
| Built-in name such as `"soil"` | Calls `load_ioc_model(name)` |
| `IoCModel` | Uses that model directly |
| `IoCFit` | Uses its `.model` field |
| Raw model dictionary from `IoCModel.to_dict()` | Calls `IoCModel.from_dict(...)` |
| Mapping with a `"model"` entry containing an `IoCModel` | Uses that entry; for example, `dict(fitted)` |

A mapping with a `"model"` key takes the last route, even if it has other keys.
Nested serialized dictionaries such as `{"model": fitted.model.to_dict()}` are
not recursively reconstructed and raise `TypeError`. Pass the inner dictionary
directly or reconstruct an `IoCModel` first. A model filename is not a supported
argument: strings are interpreted as built-in model names.

The result is always a one-dimensional `float64` NumPy array of length `n`, even
for scalar input. It contains index values only, without confidence intervals.

```python
one_score = predict_ioc(
    model="soil", dcub=-0.2, n_cazy=30, n_genes=2000, completeness=100
)
assert one_score.shape == (1,)
assert np.isfinite(one_score).all()

from_vectors = predict_ioc(
    model=fitted,
    dcub=traits["dCUB"],
    n_cazy=traits["nCAZy"],
    n_genes=traits["nGenes"],
    completeness=completeness,
)
np.testing.assert_allclose(from_vectors, fitted.values)
np.testing.assert_allclose(
    predict_ioc(traits, model=dict(fitted), completeness=completeness),
    fitted.values,
)
```

Invalid feature inputs and unknown built-in names raise `ValueError`.
Unsupported model objects raise `TypeError`; omitting the required `model`
argument also raises `TypeError`. Dictionary reconstruction can raise the
exceptions described under serialization below.

## Bundled models

```text
available_ioc_models() -> tuple[str, ...]
load_ioc_model(name: str) -> IoCModel
```

`available_ioc_models()` returns the following tuple, in this order:

```text
("permafrost", "pacific", "soil", "human", "all_habitat")
```

`load_ioc_model` loads a PCA model exported from the original gRodon 2.7.3 model
data. Names are case-sensitive; other names raise `ValueError`. Loading uses
packaged JSON data, is cached, and does not fit a model, invoke R, or download
anything. These models retain their original training centers, scales, and
minimum oriented PC1 scores. Training scores are not bundled, so `.scores` and
`.x` are `None` for the distributed models.

```python
names = available_ioc_models()
assert names == ("permafrost", "pacific", "soil", "human", "all_habitat")

soil_model = load_ioc_model("soil")
assert soil_model.name == "soil"
assert soil_model.scores is None
np.testing.assert_allclose(
    predict_ioc(new_traits, model=soil_model),
    predict_ioc(new_traits, model="soil"),
)
```

Model provenance is recorded in `src/grodon/data/ioc_models.json`. See
[compatibility and validation](COMPATIBILITY.md) for the R comparison scope.

## `IoCModel`

```text
IoCModel(
    center: npt.NDArray[np.float64],
    scale: npt.NDArray[np.float64],
    rotation: npt.NDArray[np.float64],
    sdev: npt.NDArray[np.float64],
    min_score: float,
    n_observations: int,
    name: str = "custom",
    scores: npt.NDArray[np.float64] | None = None,
)
```

This is a frozen dataclass. Array-valued constructor inputs are copied into
`float64` arrays and marked read-only, so ordinary array assignment cannot change
the fitted parameters. Lists such as those read from JSON are accepted through
this conversion. Field reassignment is prevented by the dataclass.

Let `k` be the number of PCA components, from one to three.

| Field | Shape/type | Meaning and constructor validation |
| --- | --- | --- |
| `center` | `(3,)` float64 array | Finite training means of `[rCAZy, dCUB, corrected_nGenes]` |
| `scale` | `(3,)` float64 array | Finite, strictly positive training scales |
| `rotation` | `(3, k)` float64 array | Finite PCA loadings, features in rows and components in columns; `rotation[1, 0]` must be nonzero |
| `sdev` | `(k,)` float64 array | Finite, nonnegative principal-component standard deviations |
| `min_score` | `float` | Finite minimum of the oriented training PC1 scores; used to preserve the index origin |
| `n_observations` | `int` | Training observation count, at least two |
| `name` | `str`, default `"custom"` | Model label; it does not select a bundled model when constructing an object |
| `scores` | `(n_observations, k)` float64 array or `None` | Optional finite, **unoriented** training PCA coordinates; available after `fit_ioc` |

The constructor validates these shapes and numeric constraints, but does not
refit the PCA, prove that loadings are orthonormal, or recompute `min_score` from
`scores`. Use `fit_ioc` or an authentic exported model when those relationships
must hold. Shape or numeric-constraint violations raise `ValueError`; numeric
conversion errors may also propagate from NumPy.

### Properties

| Property | Return type | Meaning |
| --- | --- | --- |
| `direction` | `float` | `-sign(rotation[1, 0])`, either `-1.0` or `1.0` |
| `variance_explained` | `(k,)` float64 array | `sdev**2 / sum(sdev**2)` |
| `x` | float64 array or `None` | Alias of `scores`, matching R's `prcomp$x` |

For valid fitted models, `variance_explained` sums to one. A manually constructed
model with all-zero `sdev` is accepted by the constructor but produces undefined
(`NaN`) variance fractions because the denominator is zero.

### `transform`

```text
IoCModel.transform(features: npt.ArrayLike) -> npt.NDArray[np.float64]
```

Projects an already prepared two-dimensional array with exactly three columns,
in the order `[rCAZy, dCUB, corrected_nGenes]`. Returns an `(n, k)` array of raw
PCA coordinates. It performs centering and scaling, but neither calculates
CAZyme fractions/completeness correction nor applies the IoC direction and
offset. It does not look up DataFrame column names. Nonfinite values or an
incorrect array shape raise `ValueError`; a `(0, 3)` array is accepted and returns
an empty `(0, k)` projection.

```python
prepared = np.column_stack(
    (
        traits["nCAZy"] / traits["nGenes"],
        traits["dCUB"],
        traits["nGenes"] / (completeness / 100),
    )
)
coordinates = fitted.model.transform(prepared)
np.testing.assert_allclose(coordinates, fitted.model.x)
reconstructed = (
    coordinates[:, 0] * fitted.model.direction - fitted.model.min_score + 1
)
np.testing.assert_allclose(reconstructed, fitted.values)
```

### Serialization: `to_dict` and `from_dict`

```text
IoCModel.to_dict(self, *, include_scores: bool = False) -> dict[str, Any]

@classmethod
IoCModel.from_dict(cls, data: Mapping[str, Any]) -> IoCModel
```

`to_dict()` returns fresh Python lists for the arrays, suitable for JSON storage.
It includes the following keys; this individual-model dictionary has no
`schema_version` or outer `models` key.

| Key | Serialized value | Required by `from_dict`? |
| --- | --- | --- |
| `features` | `["rCAZy", "dCUB", "nGenes"]` | Optional; this exact feature order is assumed if absent |
| `center` | List of three numbers | Yes |
| `scale` | List of three numbers | Yes |
| `rotation` | Three lists, each containing `k` numbers | Yes |
| `sdev` | List of `k` numbers | Yes |
| `min_score` | Number | Yes |
| `n_observations` | Integer | Yes |
| `name` | String | Optional; defaults to `"custom"` |
| `scores` | List of `n_observations` rows, each containing `k` numbers, or `null` | Optional; defaults to `None` |

The serialized feature label `nGenes` denotes the completeness-corrected
gene-count feature. `to_dict()` omits `scores` by default. With
`include_scores=True`, it includes them only if the model actually has training
scores; bundled models do not acquire scores through serialization.

`from_dict()` checks the `features` order, removes that key, and passes every
remaining key to the `IoCModel` constructor. A wrong order raises `ValueError`.
Unknown extra keys and missing required constructor keys raise `TypeError`;
constructor shape and numeric validation still applies. It does not accept the
outer package model bundle, a JSON string, or an `IoCFit` result dictionary in
place of an individual model dictionary.

This example saves a model as JSON text in memory and loads it again without
writing any files:

```python
portable_model = fitted.model.to_dict()
assert "scores" not in portable_model
saved_json = json.dumps(portable_model, allow_nan=False)
restored = IoCModel.from_dict(json.loads(saved_json))

assert restored.scores is None
assert not restored.rotation.flags.writeable
np.testing.assert_allclose(
    predict_ioc(traits, model=restored, completeness=completeness),
    fitted.values,
)
np.testing.assert_allclose(
    predict_ioc(traits, model=portable_model, completeness=completeness),
    fitted.values,
)

saved_with_scores = json.dumps(fitted.model.to_dict(include_scores=True), allow_nan=False)
restored_with_scores = IoCModel.from_dict(json.loads(saved_with_scores))
np.testing.assert_allclose(restored_with_scores.x, fitted.model.x)
```

## `IoCFit`

```text
IoCFit(values: npt.NDArray[np.float64], model: IoCModel)
```

`fit_ioc` returns this frozen dataclass. Its fields are `values`, a length-`n`
float64 array of training IoC scores, and `model`, the fitted `IoCModel`.
Unlike `IoCModel`'s parameter arrays, the `values` array returned by `fit_ioc` is
writable. Freezing the dataclass prevents field reassignment, not mutation of
that array. The direct `IoCFit` constructor performs no additional validation or
array conversion; normally obtain instances through `fit_ioc`.

| Access or method | Result |
| --- | --- |
| `fitted.values` | Training IoC scores |
| `fitted.IoC` | Alias of `values` |
| `fitted.model` | Fitted `IoCModel` |
| `fitted["IoC"]` | The same array as `values` |
| `fitted["model"]` | The same object as `model` |
| `fitted.keys()` | `("IoC", "model")` |
| `iter(fitted)` | Iterator over the two keys, in the order above |
| `len(fitted)` | `2`, the number of result fields, not the number of observations |
| `dict(fitted)` | `{"IoC": fitted.values, "model": fitted.model}` |

Other keys raise `KeyError`. `IoCFit` implements this dictionary-style protocol,
but is not a `Mapping` subclass and does not provide `items()`, `get()`, or a
`to_dict()` method. `dict(fitted)` is not directly JSON-serializable because it
contains a NumPy array and an `IoCModel`; serialize the model through
`fitted.model.to_dict()` and values through `fitted.values.tolist()` if needed.

## R-style compatibility wrappers

These wrappers are exported alongside the Python-style API. Their signatures
retain the original mixed-case parameter names and allow positional arguments:

```text
fitIoC(x=None, dCUB=None, nCAZy=None, nGenes=None, completeness=100)

predictIoC(
    x=None,
    dCUB=None,
    nCAZy=None,
    nGenes=None,
    completeness=100,
    model=None,
)
```

`fitIoC` returns `IoCFit`; `predictIoC` returns the same one-dimensional float64
array as `predict_ioc`. They forward `dCUB` to `dcub`, `nCAZy` to `n_cazy`, and
`nGenes` to `n_genes`. All input, completeness, and model rules above still apply.
Although `predictIoC` has a default `model=None` for signature compatibility,
prediction requires a supported model; leaving it as `None` raises `TypeError`.
These wrappers do not execute R.

```python
compatible = fitIoC(
    dCUB=traits["dCUB"],
    nCAZy=traits["nCAZy"],
    nGenes=traits["nGenes"],
    completeness=completeness,
)
np.testing.assert_allclose(compatible["IoC"], fitted.values)
np.testing.assert_allclose(
    predictIoC(traits, model=compatible, completeness=completeness),
    fitted.values,
)
```

For sequence-derived `dCUB`, growth prediction, and the remaining package APIs,
return to the [API index](API.md).
