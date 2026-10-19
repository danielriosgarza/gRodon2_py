# Python API reference

The public functions and result types described here are imported from `grodon`,
unless another module is shown. `grodon.__version__` reports the installed package
version. Python predictions use the bundled models and do not require R.

Related references:

- [Codon counts and individual metrics](API_CODON.md): lower-level codon functions
  and numerical conventions.
- [Model inspection, fitting, and statistic batches](API_FITTING.md):
  `LinearModel`, `available_models`, `load_model`, `model_provenance`,
  `boxcox_transform`, `fit_linear_model`, `estimate_boxcox_lambda`, `fit_models`,
  `fit_gc_models`, `fit_trim_models`, `get_statistics`, and `get_statistics_batch`.
- [Copiotrophy](API_COPIOTROPHY.md): `IoCModel`, `IoCFit`, `fit_ioc`, `predict_ioc`,
  `available_ioc_models`, and `load_ioc_model`.
- [Compatibility notes](COMPATIBILITY.md): intentional corrections and retained
  conventions from the original R implementation.

Signatures omit some type annotations for readability, but preserve argument
order, keyword-only boundaries, and defaults. A `*` marks keyword-only parameters.
Paths accept strings or `os.PathLike` objects, including `pathlib.Path`.

## Predict from coding sequences

```text
predict_growth(
    genes, highly_expressed, *,
    mode="full", temperature=None, training_set="madin",
    depth_of_coverage=None, fragments=False, exclude_short=False,
    genetic_code=None, n_le=100, background="all", bg=None,
    random_state=None, confidence=0.95,
) -> PredictionResult
```

`genes` contains individual coding DNA sequences already oriented and in frame.
It may be a `GeneSequences` object, an iterable of strings or Biopython
`SeqRecord` objects, or a mapping whose **values** are sequences. Mapping insertion
order is preserved. A filename or a single DNA string is not a sequence collection;
use `read_fasta` for a CDS file or `run_genome` for an unannotated assembly.

`highly_expressed` is a one-dimensional boolean sequence with one entry per input
gene, in the same order. Integer zero/one flags are not accepted as booleans.
Both the input and the retained set must contain at least one selected gene.
Ribosomal protein genes are the usual selected set. Preserve gene order between
calculations: the original codon-pair statistic includes inter-gene junctions.

| Keyword | Default | Meaning and constraints |
|---|---|---|
| `mode` | `"full"` | Model/preprocessing choice; see the mode table below. |
| `temperature` | `None` | Optional temperature in °C. A finite numeric value selects a temperature model. `None` or `"none"` omits temperature; NaN also omits it and emits `GrowthPredictionWarning`. Infinity and booleans are rejected. |
| `training_set` | `"madin"` | `"madin"`, `"vs"`/`"original"`, or `"AOA_NOB"`; names are normalized case-insensitively. Eukaryotic modes use their own training set. Short-read modes and `metagenome_v2` select Madin, with a warning if another set was requested. |
| `depth_of_coverage` | `None` | One finite, nonnegative value per input gene, aligned before filtering. Used in metagenomic modes; otherwise ignored with a warning. Retained genes and retained HE genes must have positive total coverage. |
| `fragments` | `False` | Keep genes longer than 240 bp. `True` changes this to longer than 120 bp. `None` accepts any positive in-frame length. The short-read modes set `fragments=None` internally. |
| `exclude_short` | `False` | In the 150/240 bp modes, remove genes shorter than the target length instead of retaining their shorter sequences. |
| `genetic_code` | `None` | Defaults to NCBI code 11 for prokaryotic modes and 1 for eukaryotic modes. Accepts an integer code, a numeric string, or a genetic-code name recognized by Biopython. Floats and booleans are rejected. An explicit choice emits a warning about the training-code convention. |
| `n_le` | `100` | Positive integer number of non-HE genes sampled with replacement for individual-background calculations. It does not change the fixed 100 nucleotide shuffles used per gene. |
| `background` | `"all"` | Pooled codon background or `"individual"`. Full and partial models do not support an individual background. Individual `metagenome_v1` requires Madin. |
| `bg` | `None` | Compatibility alias for `background`. Overrides `background="all"`; conflicting values when `background` is not `"all"` raise `ValueError`. |
| `random_state` | `None` | Integer seed or NumPy `Generator` for stochastic calculations. A fixed seed reproduces Python sampling; it does not reproduce R's random-number stream. Pooled full/partial predictions are deterministic. |
| `confidence` | `0.95` | Confidence level strictly between zero and one. Bounds describe the fitted regression mean, not an observation prediction interval. |

`PREDICTION_MODES` is the exported tuple of supported mode names:

| Mode | Behavior |
|---|---|
| `full` | Uses HE codon usage, HE consistency, and codon-pair bias. |
| `partial` | Uses the partial regression without the codon-pair predictor. |
| `metagenome_v1` | Original metagenomic regression; coverage can weight aggregated statistics. |
| `metagenome_v2` | Uses pooled statistics when `ConsistencyHE < 0.6`; otherwise recomputes an individual background. |
| `metagenome_150bp`, `metagenome_240bp` | Trim from the start to the target length, retain shorter sequences unless excluded, calculate all five codon metrics, and use the individual background. |
| `eukaryote` | Uses the eukaryotic regression and default genetic code 1. |
| `metagenome_euk` | Eukaryotic metagenomic model; switches to an individual background when pooled `ConsistencyHE >= 0.6`. |
| `meta_testing`, `meta_nogc_testing` | Exposed experimental regression variants retained for compatibility; these names are not aliases for the standard metagenomic modes. |

Missing coverage in `metagenome_v1`, `metagenome_v2`, or `metagenome_euk` emits a
warning. Missing temperature in a eukaryotic mode also emits a warning.

### Self-contained example: read, select, calculate, predict

This executable example uses artificial CDS records only. Its model output has
no biological interpretation. It records warnings explicitly, including the
expected warning about the estimate's range for these synthetic inputs.

```python
from io import StringIO
from itertools import product
from random import Random
import warnings

from grodon import (
    GrowthPredictionWarning,
    get_codon_statistics,
    highly_expressed_mask,
    predict_growth,
    read_fasta,
)

rng = Random(7)
codons = [
    "".join(bases)
    for bases in product("ACGT", repeat=3)
    if "".join(bases) not in {"TAA", "TAG", "TGA"}
]
records = []
for index in range(24):
    label = "ribosomal protein" if index < 12 else "synthetic protein"
    sequence = "ATG" + "".join(rng.choices(codons, k=100)) + "TAA"
    records.append(f">gene_{index} {label}\n{sequence}\n")

genes = read_fasta(StringIO("".join(records)))
selected = highly_expressed_mask(genes)
statistics = get_codon_statistics(genes, selected)
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always", GrowthPredictionWarning)
    result = predict_growth(genes, selected)

assert len(genes) == 24 and statistics["nHE"] == result["nHE"] == 12
print(result.model_name)
print([str(item.message) for item in caught])
```

## Predict from existing statistics

```text
predict_from_statistics(
    statistics, *, mode="full", temperature=None, training_set="madin",
    background="all", confidence=0.95,
) -> PredictionResult
```

`statistics` is a mapping of names to numeric scalars. The selected model's
predictors must be present and finite. For example, the full model needs
`CUBHE`, `ConsistencyHE`, and `CPB`; a temperature model also uses `OGT`, supplied
through `temperature`. The input mapping is copied. If `CUB` and `CUBHE` are
present, `dCUB` is recalculated as `(CUB - CUBHE) / CUB`; zero `CUB` yields NaN.
Providing `temperature` sets or replaces `OGT`.

Precomputed statistics must match the requested background and preprocessing.
For pooled `metagenome_v2` or `metagenome_euk` with `ConsistencyHE >= 0.6`, this
function raises `ValueError`: aggregate values cannot supply the required
gene-level recalculation. Use `predict_growth`, or provide statistics already
calculated with `background="individual"` and select that background explicitly.

## Prediction results and serialization

```text
PredictionResult(
    statistics, model_name, mode, training_set, background, confidence=0.95
)
```

Normally obtain this object from a prediction function. Direct construction does
not fit a model or validate the presence of prediction fields.

`PredictionResult` is an immutable `Mapping[str, float | int]`. Its numerical
mapping is copied and made read-only when the result is created. Model metadata
is stored as attributes outside that mapping.

| Member | Meaning |
|---|---|
| `result["d"]`, `result.d`, `result.doubling_time` | Estimated doubling time in hours. |
| `result["LowerCI"]`, `result.lower_ci` | Lower confidence bound in hours. |
| `result["UpperCI"]`, `result.upper_ci` | Upper confidence bound in hours. |
| `result.statistics` | Read-only numerical mapping. |
| `result.model_name` | Name of the bundled fitted model actually used. |
| `result.mode`, `result.training_set`, `result.background` | Effective prediction choices; background may have changed automatically. |
| `result.confidence` | Confidence level used for the bounds. |
| `result.to_dict()` or `dict(result)` | Detached mutable copy of the numerical mapping; does not include model metadata. |

A default full prediction returns the eight aggregate-statistic keys documented
below, plus `dCUB`, `d`, `LowerCI`, and `UpperCI`. A temperature prediction adds
`OGT`. Short-read results contain `MILC`, `ENCprime`, `B`, `SCUO`, `MCB`, `nHE`,
the three prediction fields, and `OGT` when supplied. Precomputed-statistic
predictions retain additional keys supplied by the caller.

Nonfinite values are preserved by `PredictionResult.to_dict()`. An undefined
inverse Box–Cox transformation can produce NaN, including in confidence bounds.
An undefined point estimate emits `GrowthPredictionWarning`; inspect its value
before using it. Estimates above the warning thresholds of 5 hours for
prokaryotic modes or 40 hours for eukaryotic modes also emit that warning.

This second self-contained example illustrates strict JSON serialization using
artificial summary statistics. Convert nonfinite values to `None` explicitly:

```python
import json
import math
from grodon import predict_from_statistics

prediction = predict_from_statistics(
    {"CUBHE": 0.8, "ConsistencyHE": 0.3, "CPB": 0.1}
)
numeric = {
    key: value if math.isfinite(value) else None
    for key, value in prediction.to_dict().items()
}
payload = {"statistics": numeric, "model_name": prediction.model_name}
print(json.dumps(payload, allow_nan=False))
```

The compatibility function `predictGrowth` returns the same result type, but
retains R-style positional arguments and the `bg` argument:

```text
predictGrowth(
    genes, highly_expressed, mode="full", temperature="none",
    training_set="madin", depth_of_coverage=None, fragments=False,
    exclude_short=False, genetic_code=None, n_le=100, bg="all", **kwargs,
) -> PredictionResult
```

## FASTA, identifiers, and coverage

### `read_fasta` and `GeneSequences`

```text
read_fasta(source) -> GeneSequences
GeneSequences(sequences: tuple[str, ...], names: tuple[str, ...])
```

`source` is a path to plain FASTA or a `.gz` FASTA, or an already opened text
stream such as `StringIO`. Gzip handling is selected by the path's suffix.
Sequences are joined across lines, stripped of whitespace, and uppercased.
Complete descriptions after `>` are retained. The reader does not enforce a
DNA alphabet, reading frame, or minimum length; statistics filtering and assembly
validation perform those checks later.

`GeneSequences` is a frozen dataclass implementing `Sequence[str]`:

| Member | Type/behavior |
|---|---|
| `sequences` | Tuple of sequence strings in FASTA order. |
| `names` | Tuple of complete FASTA descriptions in the same order. |
| `ids` | Tuple containing each description's first whitespace-delimited token. |
| `len(genes)`, iteration | Number of records; iterate over sequence strings. |
| `genes[index]` | One sequence string. |
| `genes[start:stop]` | A tuple of strings, not a new `GeneSequences` object. |

Prefer `read_fasta` when loading files. It rejects an empty file, empty headers
or records, sequence text before the first header, and duplicate first-token
IDs with `ValueError`. The dataclass constructor itself only checks that the
two supplied tuples have equal lengths; it does not repeat parser validation.

### `highly_expressed_mask`

```text
highly_expressed_mask(
    genes, *, pattern="^(?!.*(?:methyl|hydroxy)).*ribosomal protein", ids=None,
) -> list[bool]
```

`genes` must provide `GeneSequences` descriptions and IDs. With no `ids`, the
case-insensitive regular expression is searched against full descriptions.
The default selects descriptions containing `ribosomal protein` while excluding
those containing `methyl` or `hydroxy`.

With `ids`, selection uses exact, case-sensitive first-token ID membership;
`pattern` is ignored. Repeated selection IDs have set semantics. Unknown IDs
raise `ValueError`. The returned list always follows FASTA order and can be all
false; prediction/statistics functions subsequently reject an empty HE set.
Invalid regular expressions raise `re.error`.

`get_statistics` has a deliberately narrower R-compatible annotation pattern
containing `0S ribosomal protein`; see [the fitting API](API_FITTING.md).

### `read_coverage`

```text
read_coverage(
    path, genes, *, id_column="id", coverage_column="coverage", delimiter="\t",
) -> list[float]
```

Reads a text table with a header and aligns coverage to `genes.ids`, regardless
of table row order. Every FASTA ID must be present exactly once. Additional table
IDs are allowed, but all rows must have valid coverage values and unique IDs.
Coverage must be finite and nonnegative; zero values are permitted here.
Statistics calculations require positive total coverage in the retained HE set.

| Keyword | Meaning |
|---|---|
| `id_column` | Header naming the FASTA ID column. |
| `coverage_column` | Header naming the numeric coverage column. |
| `delimiter` | Field separator; the default is tab. Use `","` for CSV. |

Missing columns, missing FASTA IDs, duplicate table IDs, or invalid numeric values
raise `ValueError`. This reader takes a filesystem path, not an open stream.

**File-based example:** replace every `path/to/...` path with your own inputs.
The FASTA must contain coding sequences and usable ribosomal descriptions.

```python
from pathlib import Path
from grodon import highly_expressed_mask, predict_growth, read_coverage, read_fasta

genes = read_fasta(Path("path/to/predicted_cds.fna.gz"))
selected = highly_expressed_mask(genes)
coverage = read_coverage(Path("path/to/coverage.tsv"), genes)
result = predict_growth(
    genes, selected, mode="metagenome_v1", depth_of_coverage=coverage
)
print(result.d)
```

## Aggregate codon statistics and filtering

```text
get_codon_statistics(
    genes, highly_expressed, *, fragments=False, exclude_short=False,
    depth_of_coverage=None, genetic_code=11, n_le=100, trimlen=None,
    trimside="start", all_metrics=False, background="all", random_state=None,
) -> dict[str, float | int]

get_codon_statistics_i(genes, highly_expressed, **kwargs) -> dict[str, float | int]
```

Input sequence/flag/coverage conventions are the same as for `predict_growth`.
The `fragments`, `depth_of_coverage`, `n_le`, `background`, and `random_state`
parameters have the meanings described above, with two differences: this function
always applies supplied coverage, and it accepts an individual background without
checking whether a prediction mode supports it. Here `genetic_code` defaults
explicitly to 11, regardless of any future model choice.

| Additional keyword | Default | Meaning |
|---|---|---|
| `trimlen` | `None` | Optional target length in bases; must be a positive integer divisible by three. Longer sequences are trimmed; shorter retained sequences are kept. |
| `trimside` | `"start"` | `"start"`, `"end"`, or `"random"`; trimming stays in frame. Random trimming uses `random_state`. |
| `exclude_short` | `False` | When `fragments=None` and `trimlen` is supplied, discard sequences shorter than `trimlen` before trimming. |
| `all_metrics` | `False` | Return the five HE codon metrics instead of the default aggregate-statistic schema. |

Trimming places genes eligible for trimming before shorter retained genes, with
aligned flags and coverage, matching the upstream ordering convention. The
individual-background wrapper sets `background="individual"` even if that key
appears in `kwargs`.

| Default return key | Meaning |
|---|---|
| `CUBHE` | Median HE codon usage bias, measured by MILC. |
| `CUB` | Median MILC across retained genes with the pooled background; the sampled non-HE aggregate for the individual background. |
| `ConsistencyHE` | Mean HE MILC relative to the HE codon background. |
| `GC` | Fraction of retained DNA bases that are G or C. |
| `GCdiv` | `abs(0.5 - GC)`. |
| `CPB` | Codon-pair bias for an unweighted pooled background; NaN for weighted or individual backgrounds. |
| `FilteredSequences` | Number of original genes removed by length/frame/ambiguity filtering. |
| `nHE` | Number of HE genes after filtering. |

With `all_metrics=True`, the return fields are exactly `MILC`, `ENCprime`, `B`,
`SCUO`, `MCB`, and `nHE`. The five metrics summarize HE genes; `CPB`, GC, and
`FilteredSequences` are absent from this schema.

With coverage, the medians use interpolated weighted medians, `ConsistencyHE`
uses a weighted mean, and `GC` weights both GC-base counts and sequence lengths.
Coverage weights do not change reference codon counts.
Individual-background sampling excludes zero-coverage non-HE candidates.
For model compatibility, individual-background metric calculations retain
upstream `singleMILC`'s code-1 convention; the requested code still controls HE
consistency. See [the codon API](API_CODON.md) for definitions and lower-level
functions, including `codon_counts`, `codon_usage`, `milc`, `enc_prime`,
`karlin_b`, `scuo`, `mcb`, and `codon_pair_bias`.

### Inspect filtering separately

```text
filter_sequences(
    genes, highly_expressed, *, length_threshold=240, depth_of_coverage=None,
) -> grodon.statistics.FilteredSequences
```

Retains only A/C/G/T sequences whose lengths are **strictly greater** than
`length_threshold` and divisible by three. A 240 bp CDS therefore fails the
default filter; a 243 bp CDS can pass. `length_threshold` must be a nonnegative
integer. Case normalization happens before filtering.

`FilteredSequences` is a frozen dataclass available from `grodon.statistics`:

| Attribute | Type/meaning |
|---|---|
| `genes` | Tuple of retained sequence strings. |
| `highly_expressed` | Aligned one-dimensional NumPy boolean array. |
| `depth_of_coverage` | Aligned NumPy floating array, or `None`. |
| `filtered` | Integer number of removed genes. |
| `indices` | NumPy integer array of zero-based original input positions. |

The array contents are not made read-only by the frozen dataclass. Filtering
emits `grodon.statistics.SequenceFilterWarning` when genes are removed and raises
`ValueError` if no HE genes survive. No selection, mismatched array lengths, or
invalid coverage also raise `ValueError`; nonboolean flags raise `TypeError`.

## Predict from an assembly

```text
run_genome(genome, output_dir, *, config=None, rscript=None) -> GenomePipelineResult
```

This offline workflow calls genes, selects ribosomal proteins by local BLASTP
homology, and runs the Python prediction. `genome` is a plain or gzipped DNA
assembly FASTA; multiple contigs in one file belong to one genome. Install the
`genome` extra and make BLASTP available, either on `PATH` or through `config`.
Input assembly DNA may contain IUPAC ambiguity symbols; gaps and protein-only
symbols are rejected. Ambiguous predicted CDS can later be removed by the
standard statistics filter.

`output_dir` must not already exist, even if it is empty. The single-genome
workflow checks the input file, dependencies, and reference before reserving the
directory. Once a run starts, its intermediates and failure manifest are retained
if a later step fails. It does not overwrite or resume an earlier run.

### `GenomePipelineConfig`

```text
GenomePipelineConfig(
    reference="bacteria-swissprot-v1", blastp="blastp",
    mode="full", training_set="madin", temperature=None, genetic_code=11,
    gene_mode="single", closed=True, mask=True,
    min_identity=70.0, min_query_coverage=0.8, min_reference_coverage=0.8,
    max_evalue=1e-8, min_ribosomal=10, confidence=0.95,
    threads=1, timeout=600.0,
)
```

This frozen dataclass validates configuration on construction. Passing `None`
to `run_genome(config=...)` creates the default configuration.

| Field | Meaning and accepted values |
|---|---|
| `reference` | Bundled profile name or path to a custom ribosomal protein FASTA, plain or gzip. The default is `bacteria-swissprot-v1`; `lactococcus-v1` is also bundled. |
| `blastp` | Executable name resolved through `PATH`, or explicit BLASTP path. |
| `mode` | `"full"` or `"partial"`; the assembly workflow does not expose every CDS-level prediction mode. |
| `training_set` | `"madin"`, `"vs"`, or `"AOA_NOB"`. This configuration does not accept the `"original"` alias. |
| `temperature` | Finite numeric °C or `None`; unlike the CDS API, this field does not accept `"none"` or NaN. |
| `genetic_code` | Integer 4 or 11, matching the gene-calling workflow. |
| `gene_mode` | `"single"` trains one gene-calling model on all contigs. `"meta"` uses pretrained bins for the requested code. This is a gene-calling choice, not a growth-prediction mode. |
| `closed` | Boolean: prevent genes from extending beyond contig ends. |
| `mask` | Boolean: enable gene-caller masking of unknown-base regions. |
| `min_identity` | BLAST amino-acid identity percentage, from 0 through 100. |
| `min_query_coverage` | Minimum query alignment-span fraction, from 0 through 1. |
| `min_reference_coverage` | Minimum reference alignment-span fraction, from 0 through 1. |
| `max_evalue` | Positive finite maximum BLAST E-value. |
| `min_ribosomal` | Positive integer minimum HE count **after** the default CDS filter. This is an operational gate, not a calibrated confidence guarantee. |
| `confidence` | Level strictly between zero and one; an R comparison requires exactly 0.95. |
| `threads` | Positive integer BLASTP thread count; it does not parallelize separate genomes. |
| `timeout` | Positive finite timeout in seconds for external BLAST/R steps. |

The pipeline always uses `fragments=False`, so the default `>240 bp` filter
applies. If no HE genes qualify, too few survive filtering, or the point estimate
is nonfinite/nonpositive, the run raises an error rather than returning a
completed `GenomePipelineResult`. A finite point estimate does not guarantee
that both confidence bounds are finite.

**File-based example:** this requires an assembly file, the `genome` extra, and
BLASTP. Replace the input and executable paths, and choose a new output path.

```python
from pathlib import Path
from grodon import GenomePipelineConfig, run_genome

config = GenomePipelineConfig(blastp=Path("path/to/blastp"), threads=2)
run = run_genome(
    Path("path/to/assembly.fna.gz"),
    Path("path/to/new_run_directory"),
    config=config,
)
print(run.prediction.d)
print(run.output_dir)
```

### `GenomePipelineResult` and run files

```text
GenomePipelineResult(prediction, output_dir, summary)
```

Normally returned by `run_genome`; direct construction does not run an analysis.

| Member | Meaning |
|---|---|
| `prediction` | The numerical `PredictionResult`. |
| `output_dir` | Absolute `Path` of the run directory. |
| `summary` | Mapping containing settings, statistics, counts, warnings, and provenance. |
| `to_dict()` | Detached JSON-compatible summary: paths become strings and nonfinite numbers become `None`. |

The result is a frozen dataclass, but it does not deeply freeze the underlying
`summary` mapping. `to_dict()` is the supported way to obtain a detached copy;
its conversion differs from `PredictionResult.to_dict()`.

Important summary keys include `status`, `doubling_time_hours`,
`doubling_time_minutes`, `growth_rate_per_hour` (`log(2) / doubling_time_hours`),
`lower_ci_hours`, `upper_ci_hours`, `confidence`, `predicted_genes`,
`ribosomal_genes`, `retained_genes`, `retained_ribosomal_genes`, `genome_bp`,
`contigs`, `statistics`, `settings`, `warnings`, `model`, `reference`,
`input_sha256`, `reference_sha256`, and `scientific_fingerprint`.

| Files | Contents |
|---|---|
| `predicted_cds.fna`, `predicted_proteins.faa`, `predicted_genes.tsv` | Called genes, proteins, and gene coordinates/IDs. |
| `gene_calling.json` | Gene-caller version, settings, counts, and provenance. |
| `ribosomal_hits.tsv`, `selected_ribosomal.tsv`, `highly_expressed_ids.txt` | Raw BLAST hits, selected best hits, and the exact selected ID list. |
| `ribosomal_annotation.json`, `blastp.log` | Annotation settings/provenance and external-process diagnostics. |
| `reference.faa.gz`, `reference.json` | The reference used for this run and its metadata. |
| `result.json`, `result.csv` | Prediction summary and numerical output. |
| `run.json` | Run status, software versions, settings, timestamps, and artifact hashes/sizes. |

The scientific fingerprint combines input/reference/CDS identities, selected IDs,
statistics, and model identity. Paths, timestamps, and logs are not part of that
fingerprint.

### Optional original-R comparison

Set `rscript` to an executable name or path for an R installation that already
contains the original gRodon package, its dependencies (including Biostrings),
and jsonlite for the comparison bridge. No installation or
network download is performed by `run_genome`. R receives exactly the generated
CDS and selected HE IDs. Additional files are `r_parameters.json`, `r_result.json`,
`r_comparison.json`, and `r_comparison.log`.

| Outcome | Behavior |
|---|---|
| No R requested | Successful result has `r_comparison_status="not_requested"`. |
| R completes and matches | `r_comparison_status="complete"`, `r_matches=True`. |
| R completes with differences | Run still completes; `r_matches=False`, details and a warning are retained. Check this flag separately from run status. |
| R execution/comparison fails | An exception is raised and `run.json` marks failure. Existing native `result.json`/CSV are retained with `status="prediction_complete_r_failed"` and an error description. |

The comparison uses relative tolerance `1e-10` and absolute tolerance `1e-12`.
Its JSON boundary represents nonfinite numbers as `null`; paired nulls count as
matching there. A matching flag is therefore not a guarantee that every reported
statistic is finite. The independent folder benchmark audit checks finite values
separately. See the compatibility notes before comparing corrected model choices
outside the default Madin path.

### Batch assemblies

```text
run_genomes(
    genomes, output_dir, *, config=None, rscript=None, continue_on_error=False,
) -> dict[str, Any]
```

`genomes` is a nonempty iterable of paths, one assembly per path; passing one
string/path instead of a collection raises `TypeError`. Processing is sequential.
Each genome gets a numbered child directory under `output_dir/genomes/`.
The batch root must be new, and `summary.json`/`summary.csv` are updated as work
finishes. Per-genome dependency/input checks happen after the batch root exists.

With `continue_on_error=False`, the first failure is recorded in the summaries
before `RuntimeError` is raised. With `True`, later inputs continue and the
returned JSON-compatible dictionary contains `requested`, `completed`, `failed`,
`results`, settings, timestamps, and output location. Final status is `"complete"`
or `"completed_with_errors"`. Numerical R mismatches do not raise, so inspect
each completed row's `r_matches` flag independently.

## Errors, warnings, and related APIs

| Situation | Behavior |
|---|---|
| Invalid arguments, unsupported mode/code/background, wrong dimensions, missing model predictors, invalid data | Usually `ValueError`; type-specific contracts can raise `TypeError` as described above. |
| Missing input/executable | `FileNotFoundError`; high-level assembly dependency checks can instead raise `RuntimeError` summarizing missing tools. |
| Existing assembly output directory or reserved annotation output | `FileExistsError`; previous files are preserved. |
| BLAST/R process failure or timeout | `RuntimeError` with diagnostic-log location when available. |
| Filtered CDS | `SequenceFilterWarning` from `grodon.statistics`; use retained counts and indices to inspect exclusions. |
| Model setting or estimate requiring attention | Exported `GrowthPredictionWarning`, a `UserWarning` subclass. |

Warnings do not by themselves invalidate a result. For example, fewer than ten
input HE genes warns in the CDS API, whereas the assembly pipeline enforces its
configured minimum after filtering. Use Python's `warnings.catch_warnings` to
record warnings, as in the self-contained example; the assembly pipeline records
them in its result summary automatically.

Training/statistic batch helpers (`get_statistics`, `get_statistics_batch`,
`fit_models`, `fit_gc_models`, `fit_trim_models`) and R-compatible aliases
`getStatistics`, `getStatisticsBatch`, `fitModels`, `fitGCModels`, and
`fitTrimModels` are documented in [API_FITTING.md](API_FITTING.md).
The aliases `fitIoC` and `predictIoC`, together with their modern names and
copiotrophy result/model types, are documented in
[API_COPIOTROPHY.md](API_COPIOTROPHY.md).
