# Fitting, training statistics, and growth-model API

[Back to the API reference](API.md). These functions run in Python and do not require R.
All names below can be imported from `grodon`, except the explicitly identified
module-only compatibility alias. See [compatibility notes](COMPATIBILITY.md) for
the preserved R conventions and corrections to upstream fitting helpers.

## Training-data contract

Fitting functions accept a pandas `DataFrame` or a mapping of column names to
equal-length arrays. Each row is one training observation. `d` contains measured
doubling time in **hours**; `OGT`, when used, contains optimal growth temperature
in **degrees Celsius**. A `File` column may identify observations but is not fitted.

Only columns needed by a particular regression are inspected. They must be numeric
and finite, and every `d` must be positive. Missing rows are not silently dropped.
Each fit needs more observations than coefficients, counting the intercept, and a
full-rank design: constant or collinear predictor columns are rejected. Additional
unused columns are permitted. Rows in the ordinary and temperature tables are
fitted independently; the functions do not join them by index or `File`.

## `estimate_boxcox_lambda`

```text
estimate_boxcox_lambda(data, predictors, *, lambda_grid=None, interpolate=True) -> float
```

| Parameter | Meaning |
| --- | --- |
| `data` | Training table satisfying the contract above. |
| `predictors` | Sequence of numeric column names, excluding `d` and the intercept. |
| `lambda_grid` | Finite, strictly increasing, one-dimensional candidate exponents. Default: 41 equally spaced values from −2 to 2. |
| `interpolate` | Default `True`. When the grid has fewer than 100 points, interpolate its likelihood profile at 100 equally spaced positions. |

Returns the exponent at the first maximum of the regression profile likelihood.
The response is divided by its geometric mean for the likelihood search. The cubic
spline uses the endpoint conditions corresponding to the R `MASS::boxcox` workflow.
With `interpolate=False`, the returned exponent is an element of `lambda_grid`.
A grid of 100 or more points is used directly even when interpolation is enabled.

Raises `ValueError` for invalid training data, an empty/nonfinite/unordered grid,
fewer than four grid points when interpolation is needed, an exact fit with
undefined likelihood, or numerical overflow in scaling or residual calculation.
Use a scientifically justified exponent or search range; this helper does not
assess the scientific suitability of a training set.

## `fit_linear_model`

```text
fit_linear_model(data, predictors, *, boxcox_lambda, name="custom") -> LinearModel
```

Fits ordinary least squares to Box-Cox-transformed `d`. `predictors` supplies the
column order; the function inserts an intercept automatically. Supply a finite
`boxcox_lambda` explicitly, including `0.0` for a logarithmic transformation.
`name` labels the returned model and defaults to `"custom"`.

The result contains fitted coefficients, their covariance scaled by the residual
variance, `df_residual = n_observations - n_coefficients`, residual standard
deviation `sigma`, and `training_set="user supplied"`. It is not registered as a
bundled model. Call its `.predict(...)` method with already computed statistics;
`predict_growth(...)` continues to select the shipped models.

Raises `ValueError` for violations of the training contract, a nonfinite exponent,
or nonfinite transformed responses. Unlike exponent estimation, a fixed-exponent
fit can have zero residual variance. An unused extra column does not invalidate
otherwise complete training data.

This synthetic example is executable and illustrates the API, not a biological
training dataset:

```python
import numpy as np
from grodon import estimate_boxcox_lambda, fit_linear_model

data = {
    "x": np.arange(1.0, 11.0),
    "d": np.array([3.0, 3.2, 4.4, 4.8, 6.1, 7.4, 7.6, 8.5, 9.6, 10.4]),
}
exponent = estimate_boxcox_lambda(data, ("x",))
model = fit_linear_model(data, ("x",), boxcox_lambda=exponent, name="synthetic")
estimate, lower, upper = model.predict({"x": 5.5}, confidence=0.95)
assert model.predictors == ("(Intercept)", "x")
assert model.df_residual == 8
assert 0 < lower < estimate < upper
```

## `fit_models`

```text
fit_models(
    stat_data, stat_data_extremo, *, lambda_value=None, eukaryotic_lambda=None
) -> dict[str, LinearModel | float]
```

`stat_data` requires `d`, `CUBHE`, `ConsistencyHE`, and `CPB`.
`stat_data_extremo` requires those columns plus `OGT`. This function returns eight
models and two exponents:

| Dictionary key | Predictor columns, excluding the intercept |
| --- | --- |
| `gRodon_model_base` | `CUBHE`, `ConsistencyHE`, `CPB` |
| `gRodon_model_temp` | `CUBHE`, `ConsistencyHE`, `CPB`, `OGT` |
| `gRodon_model_partial` | `CUBHE`, `ConsistencyHE` |
| `gRodon_model_partial_temp` | `CUBHE`, `ConsistencyHE`, `OGT` |
| `gRodon_model_meta` | `CUBHE` |
| `gRodon_model_meta_temp` | `CUBHE`, `OGT` |
| `gRodon_model_euk` | `CUBHE` |
| `gRodon_model_euk_temp` | `CUBHE`, `OGT` |
| `lambda_milc` | Float exponent shared by the base, partial, and meta models. |
| `lambda_milc_euk` | Float exponent used by the euk models. |

`lambda_value=None` estimates `lambda_milc` on `stat_data` using the three full-model
predictors. `eukaryotic_lambda=None` separately estimates the euk exponent using
`CUBHE`. Explicit finite values skip those searches. Each temperature model uses
the corresponding exponent estimated on the ordinary table. No separate
temperature-table exponent is estimated.

Training and exponent-search validation errors propagate. The model names specify
the regression structures; passing a table does not establish that its organisms
are appropriate for every structure being fitted.

## `fit_gc_models`

```text
fit_gc_models(
    stat_data, stat_data_extremo, *, lambda_value=None, legacy_lambda=None
) -> dict[str, LinearModel | float]
```

Both tables require `d`, `CUB`, `CUBHE`, and `GCdiv`; the second also requires `OGT`.
The function copies the tables and recomputes `dCUB = (CUB - CUBHE) / CUB`, replacing
any supplied `dCUB` column. `CUB` must be nonzero. The input tables are not changed.

| Dictionary key | Predictor columns, excluding the intercept |
| --- | --- |
| `meta_model_base` | `dCUB`, `GCdiv` |
| `meta_model_temp` | `dCUB`, `GCdiv`, `OGT` |
| `meta_nogc_model_base` | `dCUB` |
| `meta_nogc_model_temp` | `dCUB`, `OGT` |
| `meta_old_model_base` | `CUBHE` |
| `meta_old_model_temp` | `CUBHE`, `OGT` |
| `lambda_meta` | Float exponent shared by the corrected and no-GC models. |
| `lambda_oldmeta` | Float exponent for the legacy models. |

If omitted, `lambda_value` is estimated from `stat_data` using `dCUB` and `GCdiv`,
and `legacy_lambda` using `CUBHE`. Their temperature variants reuse those exponents.
Missing required columns, zero `CUB`, and invalid regression inputs raise errors.
Supply numeric `CUB` and `CUBHE` columns because their arithmetic precedes fitting.
When predicting with a returned `LinearModel`, supply `dCUB` yourself; the model's
`.predict()` method does not derive it from `CUB` and `CUBHE`.

The returned legacy temperature model is the fitted legacy model. This corrects
an upstream R return-slot typo; it does not alter the already bundled R models.

## `fit_trim_models`

```text
fit_trim_models(stat_data, stat_data_extremo, *, lambda_value=None)
# Returns dict[str, LinearModel | float]
```

Both tables require `d`, `MILC`, `ENCprime`, `B`, `SCUO`, and `MCB`.
`stat_data_extremo` additionally requires `OGT`.

| Dictionary key | Contents |
| --- | --- |
| `gRodon_model_trim` | Model using `MILC`, `ENCprime`, `B`, `SCUO`, `MCB`. |
| `gRodon_model_trim_temp` | Model using those five predictors plus `OGT`. |
| `lambda_cub` | Float exponent shared by both models. |

`lambda_value=None` estimates the exponent on `stat_data` using all five metrics.
Otherwise the supplied finite exponent is used. This helper fits numeric tables;
it does not trim sequences. Compute statistics with the intended trimming settings
before fitting, and use the same definitions for subsequent predictions.

## `get_statistics`

```text
get_statistics(
    gene_file=None, *, genes=None, genetic_code=11, bg="all",
    fragments=False, exclude_short=False, trimlen=None, trimside="start",
    all_metrics=False, random_state=None
) -> dict[str, object]
```

Calculates statistics from annotated **nucleotide CDS**, selecting highly expressed
genes from annotations. It neither predicts genes from assemblies nor reads
amino-acid FASTA as codon data.

| Parameter | Accepted values and behavior |
| --- | --- |
| `gene_file` | Path string or `Path` to plain or gzip-compressed FASTA. Mutually exclusive with `genes`. |
| `genes` | `GeneSequences`, an ordered mapping of descriptions to sequences, or a sequence of annotated records such as Biopython `SeqRecord`. Bare strings without annotations cannot identify the HE set. |
| `genetic_code` | NCBI genetic-code identifier, integer or string; default `11`. |
| `bg` | `"all"` or `"individual"`; default `"all"`. |
| `fragments` | `False`: retain lengths >240 bp. `True`: >120 bp. `None`: no minimum-length cutoff. All retained sequences must still be in-frame A/C/G/T DNA. |
| `exclude_short` | Default `False`. With `fragments=None` and a `trimlen`, removes genes shorter than `trimlen`. |
| `trimlen` | `None`, or a positive integer divisible by three. Longer sequences are trimmed to this many bases. |
| `trimside` | `"start"`, `"end"`, or `"random"`; default `"start"`. Random windows remain in-frame. |
| `all_metrics` | Default `False`; chooses the return fields listed below. |
| `random_state` | Integer seed, NumPy `Generator`, or `None`. Controls individual backgrounds and random trimming. |

The HE annotation regex is case-insensitive
`^(?!.*(?:methyl|hydroxy)).*0S ribosomal protein`. It matches, for example,
`30S ribosomal protein S1` and `50S ribosomal protein L1`, while excluding
descriptions containing `methyl` or `hydroxy`. This is narrower than the default
`highly_expressed_mask` pattern. For explicit HE flags, coverage weighting, or a
custom `n_le`, use `get_codon_statistics` instead; this wrapper does not expose
those parameters. Its individual-background sampling uses the default `n_le=100`.

Returns `CUBHE`, `GC`, `GCdiv`, `ConsistencyHE`, `CUB`, `CPB`, `FilteredSequences`,
`nHE`, and `File`. With `all_metrics=True`, returns `MILC`, `ENCprime`, `B`, `SCUO`,
`MCB`, `nHE`, and `File` instead. `File` is the input basename or
`"genome provided directly"`. No measured `d` or `OGT` is added: join independently
measured training values to these statistics before fitting.

Invalid input combinations, no matching or surviving HE genes, invalid settings,
and unsuitable individual backgrounds raise errors. Removed sequences produce
`SequenceFilterWarning`. Individual backgrounds return NaN for `CPB`, so those
outputs cannot train a regression requiring finite `CPB`. Integer seeds reproduce
Python calculations with unchanged inputs/settings, not R's random stream.

## `get_statistics_batch`

```text
get_statistics_batch(
    directory, *, workers=1, pattern="*", errors="raise", **kwargs
) -> pandas.DataFrame
```

| Parameter | Meaning |
| --- | --- |
| `directory` | Existing directory path; otherwise `NotADirectoryError`. |
| `workers` | Positive integer, default `1`. More workers use threads. |
| `pattern` | `Path.glob` pattern, default `"*"`; use a FASTA-specific pattern such as `"*.fna.gz"` to avoid unrelated files. |
| `errors` | `"raise"` propagates a failed file's exception. `"skip"` emits `RuntimeWarning` and omits that file. |
| `**kwargs` | Options forwarded to `get_statistics`, such as `bg`, `trimlen`, and `all_metrics`. Batch `random_state` must be an integer or `None`, not a `Generator`. |

Files are sorted by path and non-files are ignored. Returns one row per successful
file, with `get_statistics` fields as columns, including `File`. Result order is
independent of worker count. No matched files, or all skipped files, produce an
empty DataFrame. A fixed integer seed creates separate deterministic per-file
streams; changing the matched file list can change those streams. This is a
statistics-table helper, separate from the assembly-based `run_genomes` API.

For example, `get_statistics_batch("annotated_cds", pattern="*.fna.gz", workers=4,
random_state=7)` requires an existing folder of annotated CDS files. Failures are
not omitted unless `errors="skip"` is explicitly requested.

## Inspecting bundled models

```text
available_models() -> tuple[str, ...]
load_model(name: str) -> LinearModel
model_provenance() -> dict[str, object]
```

- `available_models()` returns the sorted names of all 36 shipped growth models,
  including temperature and individual-background variants.
- `load_model(name)` retrieves a cached immutable model using its original R object
  name. An unknown name raises `ValueError` with a suggestion to list the names.
- `model_provenance()` returns a copy of source metadata: `package`, `version`,
  `file`, `sha256`, `git_revision`, and `export_method`. The source is gRodon 2.7.3;
  the bundled coefficients are exported fitted values, not refitted on installation.

```python
from grodon import available_models, load_model, model_provenance

assert "gRodon_model_base_madin" in available_models()
bundled = load_model("gRodon_model_base_madin")
assert bundled.predictors == ("(Intercept)", "CUBHE", "ConsistencyHE", "CPB")
assert model_provenance()["version"] == "2.7.3"
```

## `LinearModel`

```text
LinearModel(
    name, predictors, coefficients, covariance, df_residual, sigma,
    boxcox_lambda, training_set="", description=""
)
```

A frozen dataclass representing one regression. Prefer obtaining instances from
`load_model` or `fit_linear_model` unless importing an independently fitted model.

| Field | Type and interpretation |
| --- | --- |
| `name` | String identifier. |
| `predictors` | Tuple of unique column names in coefficient order. `"(Intercept)"`, when present, means a constant 1. |
| `coefficients` | Tuple of finite floats, one per predictor. |
| `covariance` | Square tuple-of-tuples coefficient covariance, already multiplied by residual variance; equivalent to R `vcov(lm)`, not an unscaled matrix. |
| `df_residual` | Positive integer residual degrees of freedom. |
| `sigma` | Finite, nonnegative residual standard deviation on the transformed scale. |
| `boxcox_lambda` | Finite transformation exponent. |
| `training_set` | Descriptive string, default `""`. |
| `description` | Descriptive string, default `""`. |

Construction converts predictor/parameter sequences to tuples. Inconsistent
dimensions, duplicate names, nonfinite coefficients/covariance, nonsymmetric or
non-positive-semidefinite covariance, invalid degrees of freedom, and invalid
`sigma` raise `ValueError`. The intercept is explicit in this constructor; only
the fitting helper inserts it automatically.

### Methods

```text
model.design_row(statistics) -> numpy.ndarray
model.predict_transformed(statistics, confidence=0.95) -> tuple[float, float, float]
model.predict(statistics, confidence=0.95) -> tuple[float, float, float]
```

`statistics` is a mapping from predictor names to numeric scalars, representing
one observation. Every non-intercept predictor must exist and be finite. Extra
keys are ignored. Missing or nonfinite required values raise `ValueError`.
The methods do not derive statistics, choose a different model, or accept a
whole DataFrame for vectorized prediction.

`design_row` returns a float array ordered as `predictors`, inserting `1.0` for
the intercept. `predict_transformed` returns `(fit, lower, upper)` on the Box-Cox
scale. Its margin is `t_quantile * sqrt(row @ covariance @ row)` with
`df_residual` degrees of freedom. `confidence` must satisfy `0 < confidence < 1`.

`predict` inverse-transforms each of those three values to hours. These are
confidence limits for the fitted regression mean on the transformed scale,
back-transformed without a bias correction; they are not prediction intervals
for a new observed culture. Extrapolation can yield an undefined/nonfinite
inverse transformation. Inspect all returned values before using them.

There is no dedicated `LinearModel.save()`/`load()` method. The public constructor
accepts a JSON round-trip of `dataclasses.asdict(model)`:

```python
import json
from dataclasses import asdict
from grodon import LinearModel, load_model

original = load_model("gRodon_model_base_madin")
encoded = json.dumps(asdict(original), allow_nan=False)
restored = LinearModel(**json.loads(encoded))
assert restored == original
```

## `boxcox_transform`

```text
boxcox_transform(x, lam, *, back_transform=False) -> float | numpy.ndarray
```

`x` is a scalar or numeric array-like input; `lam` is a finite scalar exponent.
Forward transformation computes `(x**lam - 1) / lam`, or `log(x)` when `lam == 0`,
and requires finite positive inputs. Violations raise `ValueError`.

With `back_transform=True`, computes `(1 + lam*x)**(1/lam)`, or `exp(x)` at zero.
The inverse uses SciPy's domain behavior: values with `1 + lam*x < 0` yield NaN;
the boundary may yield zero or infinity, and overflow can also yield infinity.
Inverse nonfinite results are returned rather than raised. Scalar input gives
a Python float; array input preserves its shape in a NumPy array.

```python
import numpy as np
from grodon import boxcox_transform

values = np.array([1.0, 2.0, 4.0])
transformed = boxcox_transform(values, 0.0)
np.testing.assert_allclose(transformed, np.log(values))
np.testing.assert_allclose(
    boxcox_transform(transformed, 0.0, back_transform=True), values
)
```

## R-style fitting and statistics aliases

```text
fitModels(stat_data, stat_data_extremo)
fitGCModels(stat_data, stat_data_extremo)
fitTrimModels(stat_data, stat_data_extremo)

getStatistics(
    gene_file=None, genes=None, genetic_code="11", bg="all", fragments=False,
    exclude_short=False, trimlen=None, trimside="start", all_metrics=False, **kwargs
)
getStatisticsBatch(
    directory, genetic_code="11", mc_cores=1, bg="all", fragments=False,
    exclude_short=False, trimlen=None, trimside="start", all_metrics=False, **kwargs
)
```

The three fitting aliases accept only the two tables; use snake-case functions
for explicit exponent overrides. They return lists in the following order, with
positions shown starting at 1:

| Position | `fitModels` | `fitGCModels` | `fitTrimModels` |
| ---: | --- | --- | --- |
| 1 | `gRodon_model_base` | `meta_model_base` | `gRodon_model_trim` |
| 2 | `gRodon_model_temp` | `meta_model_temp` | `gRodon_model_trim_temp` |
| 3 | `gRodon_model_partial` | `meta_nogc_model_base` | `lambda_cub` |
| 4 | `gRodon_model_partial_temp` | `meta_nogc_model_temp` | — |
| 5 | `gRodon_model_meta` | `lambda_meta` | — |
| 6 | `gRodon_model_meta_temp` | `meta_old_model_base` | — |
| 7 | `lambda_milc` | `meta_old_model_temp` | — |
| 8 | `gRodon_model_euk` | `lambda_oldmeta` | — |
| 9 | `gRodon_model_euk_temp` | — | — |
| 10 | `lambda_milc_euk` | — | — |

`getStatistics` retains positional R-style arguments and returns the same
dictionary as `get_statistics`; `random_state` may be passed through `**kwargs`.
`getStatisticsBatch` returns the same DataFrame as `get_statistics_batch`.
Its `mc_cores` controls workers; `workers=` overrides it. A dictionary expansion
such as `**{"mc.cores": 4}` is also accepted for R-name compatibility, with
precedence `workers`, then `mc.cores`, then `mc_cores`.

`grodon.fitting.getStatisticsBatchWindows` is an alias of `getStatisticsBatch`.
It is available in that module, not re-exported from the top-level `grodon`
namespace. Both batch names use the same portable thread implementation.
