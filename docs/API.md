# API guide

## Growth prediction

`predict_growth(genes, highly_expressed, *, mode="full", temperature=None,
training_set="madin", depth_of_coverage=None, fragments=False,
exclude_short=False, genetic_code=None, n_le=100, background="all",
random_state=None, confidence=0.95)` returns `PredictionResult`.

The result is an immutable mapping using R keys. Convenience properties are
`doubling_time`, `d`, `lower_ci`, and `upper_ci`. Metadata are `model_name`,
`mode`, `training_set`, `background`, and `confidence`. `to_dict()` copies
the numerical mapping. `bg=` is an alias for `background=`.

`predict_from_statistics(statistics, ...)` applies a bundled model to already
computed statistics. If automatic mode selection requires an individual
background, the supplied statistics must describe that background; this
function cannot reconstruct genes from aggregate statistics.

`available_models()`, `load_model(name)`, and `model_provenance()` inspect
the fitted regression assets. `LinearModel.predict_transformed(...)` returns
the fit and confidence bounds on the Box-Cox scale. `LinearModel.predict(...)`
returns hours. `boxcox_transform` supports forward/inverse transformation
and the logarithmic limiting case at lambda zero.

## Codon statistics

`get_codon_statistics(genes, highly_expressed, *, fragments=False,
exclude_short=False, depth_of_coverage=None, genetic_code=11, n_le=100,
trimlen=None, trimside="start", all_metrics=False, background="all",
random_state=None)` returns a dictionary.

- Default keys: `CUBHE`, `GC`, `GCdiv`, `ConsistencyHE`, `CUB`, `CPB`,
  `FilteredSequences`, `nHE`.
- `all_metrics=True`: `MILC`, `ENCprime`, `B`, `SCUO`, `MCB`, `nHE`.
- Weighted/individual backgrounds do not calculate CPB; its value is NaN.
- `trimside` supports `start`, `end`, and `random`, always in-frame.
- `get_codon_statistics_i` is the individual-background convenience wrapper.
- `filter_sequences` returns sequences, aligned flags/coverage, counts, and
  original row indices. No surviving highly expressed genes is an error.

In `grodon.codon`, `codon_counts` produces rows of 64 counts in `CODONS`
order. `codon_usage` computes a metric relative to supplied reference counts;
`milc`, `enc_prime`, `karlin_b`, `scuo`, and `mcb` are convenience functions.
`codon_pair_bias` implements the original R package's pair-count convention.

## Copiotrophy

`fit_ioc(x=None, *, dcub=None, n_cazy=None, n_genes=None,
completeness=100)` accepts a table with `dCUB`, `nCAZy`, `nGenes` columns or
three explicit vectors. The returned `IoCFit` has `values` and `model`
attributes and R-compatible `['IoC']` and `['model']` entries.

`predict_ioc(..., model=...)` accepts either a fitted `IoCModel` or a built-in
model name. Use `available_ioc_models()` and `load_ioc_model(name)` to inspect
the available models. Invalid counts, incomplete columns, nonfinite input,
and invalid completeness produce explicit exceptions.

## Training and batch work

`get_statistics(gene_file=None, *, genes=None, genetic_code=11, bg="all", ...)`
uses R's narrower `0S ribosomal protein` annotation pattern. It adds a `File`
label to the statistic dictionary. `get_statistics_batch(directory,
workers=1, pattern="*", errors="raise", ...)` returns a pandas DataFrame.
Use a pattern such as `*.fna.gz` to avoid unrelated directory contents.
`errors="skip"` emits warnings for failed inputs. Integer random seeds give
repeatable results independent of worker count.

`fit_models(stat_data, stat_data_extremo)` fits full, partial, metagenome-v1,
and eukaryotic regressions. Tables contain `d` (measured hours), `CUBHE`,
`ConsistencyHE`, `CPB`, and, for the temperature table, `OGT`.

`fit_gc_models` additionally uses `CUB` and `GCdiv`; it derives `dCUB`.
`fit_trim_models` uses `MILC`, `ENCprime`, `B`, `SCUO`, `MCB`, and requires
`OGT` in its second, temperature-training table. These functions return named dictionaries of `LinearModel` objects
and Box-Cox exponents. Explicit lambda overrides are available.

`fit_linear_model(data, predictors, boxcox_lambda=..., name=...)` fits one
model. `estimate_boxcox_lambda(data, predictors)` reproduces the MASS grid
and spline likelihood search. Training requires finite complete observations,
positive measured doubling times, full-rank predictors, and positive residual
degrees of freedom; missing rows are not silently dropped.

The aliases `fitModels`, `fitGCModels`, and `fitTrimModels` retain R's ordered
list results. `getStatistics` and `getStatisticsBatch` retain R-style names.
