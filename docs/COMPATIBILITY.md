# Compatibility with R gRodon 2.7.3

Source: `jlw-ecoevo/gRodon2`, revision
`4b13053e7f49ca99bceeee36940f5e820c5b0c32`.

## Covered functionality

All exported R entry points (`predictGrowth`, `fitIoC`, `predictIoC`) are
implemented, along with the documented statistics, batch processing, and
model fitting helpers. The port includes all ten source prediction modes,
36 fitted growth regressions, and five fitted IoC PCA models. The extensive
source training datasets and example genomes are not runtime dependencies;
the port includes the fitted parameters and synthetic validation fixtures.

## Preserved numerical conventions

- Model coefficients, Box-Cox exponents, residual degrees of freedom, and
  covariance come directly from the stored R objects. No model retraining
  occurs on installation or prediction.
- coRdon's alternative-initiation convention changes the amino acid grouping
  of all alternative-start codons for usage statistics. This is preserved,
  including its effect on noninitial codons.
- CPB counts the concatenated coding sequences, including pairs spanning gene
  boundaries. It sums the scores of distinct pairs without weighting by their
  counts, and divides by the number of distinct pairs minus one.
- CPB preserves R's asymmetric initiation-codon translation: isolated codons
  in pair annotations are treated as starts, whereas the concatenated genome
  has only one initial codon. This can produce undefined values for artificial
  inputs lacking the relevant amino acid counts.
- Weighted medians interpolate at cumulative-weight midpoints, following
  `matrixStats::weightedMedian`.
- Individual-background codon metrics use genetic code 1 because that is the
  convention used by R's `singleMILC` and the trained individual models.
  The requested genetic code is still used for consistency statistics.
- The default length cutoff is strict: **greater than** 240 bp, rather than
  greater than or equal to 240 bp.
- New IoC predictions retain the training PCA origin and offset, and can be
  below one. PCA component signs are oriented by the dCUB coefficient.
- The shipped individual-background metagenome-v1 temperature model is
  retained exactly, including its upstream duplication of a GC-corrected
  model. Refitting helpers correct that return-slot error for new models.

## Deliberate corrections and extensions

| R source behavior | Python behavior |
| --- | --- |
| Reversed training-set condition silently forces Madin for most non-short-read calls | Honors `vs` and `AOA_NOB` for their supported modes |
| Short-read setting can retain an irrelevant dataset label | Consistently selects Madin |
| Eukaryotic metagenome output uses bacterial inverse-transform exponent | Uses the corresponding eukaryotic exponent |
| Explicit individual-background metagenome paths rejected by a mode-name typo | Supports available individual metagenome models |
| `n_le` ignored by some individual-background paths | Honors the supplied sample count |
| End trimming includes an extra base; one random-trim path hardcodes 150 | Exact, in-frame requested lengths and valid sampling bounds |
| `fitGCModels` returns the wrong legacy temperature model | Returns the fitted legacy model in its proper slot |
| Missing temperature sometimes triggers ambiguous R condition errors | Handles `None`/NaN explicitly |
| Eukaryotic metagenome warning uses the bacterial saturation threshold | Uses the eukaryotic threshold |
| Box-Cox helper undefined at exponent zero | Implements the logarithmic limit |
| Invalid lengths/coverage/metadata can recycle or fail indirectly | Validates dimensions, finiteness, IDs, and nonnegative coverage |
| Random low-expression sampling may select only zero-coverage genes | Samples the positive-coverage background and rejects an absent background explicitly |
| Batch execution relies on Unix-style process forking | Portable threads and explicit per-file errors |

Stochastic methods use NumPy's random generators. An integer `random_state`
makes repeated Python calls reproducible, but matching R's integer seed does
not select identical samples or permutations. Python shares shuffled
background calculations across short-read metrics; R regenerates backgrounds
for each metric. Each metric uses the same statistical procedure, with a
different joint Monte Carlo realization. Bitwise cross-language equality is
therefore not claimed for stochastic end-to-end predictions.

## Independent reference validation

- `scripts/export_growth_reference.R` runs the original R `predict.lm` on five
  synthetic feature rows for each of the 36 stored growth models. The 180
  committed cases check transformed predictions, 95% confidence bounds,
  inverse-transformed hours, degrees of freedom, and exponents.
- `scripts/generate_ioc_fitting_reference.R` runs original IoC functions and
  the published regression formulas with R's `MASS::boxcox` and `stats::lm`
  on synthetic numeric data. Fixtures verify all
  five built-in IoC models, custom PCA fitting, completeness correction,
  Box-Cox selection, and all 16 fitted regressions' coefficients and covariance.
- Analytical tests cover codon metrics, filtering, alignment, weight scaling,
  genetic codes, interval behavior, and sampling reproducibility.
- `scripts/generate_codon_reference.R` executes original gRodon/coRdon functions
  on artificial sequences. It checks aggregate and coverage-weighted metrics,
  filtering, trimming, genetic codes, alternative initiation conventions, CPB,
  weighted medians, and individual metrics with fixed shuffled backgrounds.
  These numeric comparisons use a relative tolerance of `2e-12`.
- Fourteen complete deterministic `predictGrowth` cases compare the original
  R sequence-to-prediction pipeline with Python, including temperature and
  weighted metagenome cases. Known upstream dispatch/transform bugs are covered
  by explicit correction tests instead of equality assertions.
- The CLI, Python API, build, and installed distribution are tested separately.

The reference fixtures record R and dependency versions. Default unit tests
consume these committed numeric fixtures and do not require R. Directly
validated numerical equivalence is distinct from new biological validation:
this port does not expand the scientific domain validated by the original
publications.

Local release verification for version 0.2.0 on Windows, including the folder
validation tools and 40 new R model references, passed **685 tests** in each of two
environments: Python 3.13.5 with NumPy 2.5.4, SciPy 1.18.1, pandas 3.0.6,
and Biopython 1.88; and Python 3.11.12 with NumPy 1.26.4, SciPy 1.11.4,
pandas 2.1.4, and Biopython 1.83. Both environments ran against a clean copy
of the public Git index, including command-line subprocess tests and the optional
live BLAST/R integration. Earlier installed-wheel verification on Python 3.11
also passed before adding the public benchmark-download utility.
One filesystem symlink test was skipped because this Windows environment does
not permit creating symlinks. Ruff lint/format, mypy, source-distribution build,
and wheel build also passed.
The configured Linux/macOS/Windows CI matrix has not been run remotely.

The assembly pipeline was exercised with actual Pyrodigal 3.7.1 and BLASTP
2.17.0 on the CIRM-BIA2698 Lactococcus assembly. Its `lactococcus-v1` profile
reproduced the earlier 2,832 CDS, 54 ribosomal matches (42 after filtering),
and all 12 numerical statistics. The general `bacteria-swissprot-v1` profile
selected 53 ribosomal genes (41 after filtering), giving a slightly different
estimate as expected from the changed annotation reference. Both profiles were
also checked against the original R gRodon 2.7.3 on their identical prepared
inputs; maximum absolute differences were `3.33e-16` and `5.27e-16`, respectively.
See the [pipeline guide](GENOME_PIPELINE.md) and the saved
[reproduction audit](../validation/pipeline_reproduction.json).

The older SciPy release's Student-t quantile differs from R by about
`1.8e-9` for one reference degrees-of-freedom value, causing confidence limits
to differ by up to `2.7e-10`. Fitting-interval tests therefore allow `1e-9`
relative/absolute tolerance across supported versions. This does not change
the fitted coefficients or the model-selection logic.

The subsequent [40-assembly folder test](TEST_GENOMES_VALIDATION.md) passed
**480/480** finite comparisons against original R gRodon 2.7.3, with no omitted
assemblies or failed preparations. This test covers the default full/Madin
model, genetic code 11, no supplied temperature, and the broad bacterial
ribosomal reference. It shares predicted CDS and ribosomal IDs between languages
and independently calculates the codon statistics and predictions. The maximum
absolute difference across all reported statistics was `1.42e-14` (an upper
confidence bound in hours); the maximum doubling-time difference was `5.77e-15`
hours. Acceptance used `rtol=1e-10, atol=1e-12`; exact binary64 equality is
reported separately. Gene counts agreed exactly.

`tests/data/folder_growth_reference.json` retains all 40 original R outputs
with input, model, and artifact provenance. The normal test suite uses their
R-computed codon statistics to check Python model predictions and confidence
limits offline. Repeating annotation and sequence-statistic comparisons requires
the full workflow documented in the test case. These checks establish software
equivalence for the recorded inputs and settings, not new biological validation.

## Updating model assets

`scripts/export_growth_models.py` extracts growth parameters from the source
`R/sysdata.rda` with the development-only `rdata` dependency.
`scripts/export_ioc_models.py` exports the five PCA models from the same file.
The bundled JSON records source checksum, version, and revision; the IoC JSON
records its source checksum as well. Preserve provenance when updating either set of
assets, regenerate reference fixtures, and run the full test suite.

Example commands, run from the Python repository with R and its reference
dependencies available:

```console
python scripts/export_growth_models.py ../gRodon2/R/sysdata.rda
python scripts/export_ioc_models.py ../gRodon2/R/sysdata.rda
Rscript scripts/export_growth_reference.R ../gRodon2/R/sysdata.rda tests/data/reference_growth_models.csv
Rscript scripts/generate_ioc_fitting_reference.R ../gRodon2 tests/reference/ioc_fitting
Rscript scripts/generate_codon_reference.R ../gRodon2 tests/fixtures/codon_reference.json
pytest
```
