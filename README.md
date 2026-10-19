# gRodon2 Python

Native Python implementation of **gRodon 2.7.3**, for estimating microbial
maximum growth potential from coding-sequence composition and calculating the
Index of Copiotrophy. Python performs the calculations directly; R is only
needed for optional direct R comparisons and reference-test fixtures.

This is an independent port, with attribution to the original gRodon authors.
The source reference is revision `4b13053e7f49ca99bceeee36940f5e820c5b0c32`.

## Install

Requires Python 3.11 or newer. From this checkout:

```console
python -m pip install .
```

For development:

```console
python -m pip install -e ".[dev,reference,genome]"
pytest
ruff check .
ruff format --check .
mypy
python -m build
```

An optional `uv.lock` records a reproducible development environment:
`uv sync --extra dev --extra reference --extra genome --frozen` installs it.

The distribution is named `grodon2-py`; the import and command are `grodon`.
It has not been published to PyPI. Runtime dependencies are NumPy, SciPy,
pandas, and Biopython. Model parameters are packaged as portable JSON.

## Start from a bacterial genome

The reproducible genome pipeline accepts a plain or gzipped assembly FASTA,
predicts coding sequences with Pyrodigal, identifies ribosomal proteins with
local BLASTP, and calculates growth potential. It runs without an LLM, API key,
or internet connection after installation. It records inputs, settings, software
versions, reference provenance, and file checksums alongside the results.

From this checkout in Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[genome]"
.\.venv\Scripts\python.exe -m grodon setup --tools-dir .tools\blast-2.17.0 --version 2.17.0
.\.venv\Scripts\python.exe -m grodon doctor --blastp .tools\blast-2.17.0\bin\blastp.exe
.\.venv\Scripts\python.exe -m grodon genome "assembly.fna.gz" --outdir results\sample --blastp .tools\blast-2.17.0\bin\blastp.exe
```

`setup` is the explicit download step and leaves system settings unchanged.
The analysis uses the bundled `bacteria-swissprot-v1` reference by default:
39,743 reviewed bacterial ribosomal proteins from UniProtKB/Swiss-Prot 2026_03.
For the exact reference set used in the archived Lactococcus comparison, add
`--reference lactococcus-v1`. Existing analysis outputs are never overwritten;
choose a new output directory for each run.

Several genomes can use the same settings:

```powershell
.\.venv\Scripts\python.exe -m grodon genome-batch "one.fna.gz" "two.fna.gz" --outdir results\batch --blastp .tools\blast-2.17.0\bin\blastp.exe --continue-on-error
```

See the [genome pipeline guide](docs/GENOME_PIPELINE.md) for the Python API,
exact Lactococcus command, optional R comparison, outputs, installation on
other platforms, and reference coverage limitations.

## Python quick start with coding sequences

The direct `predict_growth` function takes **in-frame nucleotide coding
sequences**. Use the genome pipeline above for unannotated whole-genome contigs.
Provide a boolean flag for every sequence
indicating membership in the highly expressed reference set. Ribosomal
protein annotations were used to train the original models.

```python
from grodon import read_fasta, highly_expressed_mask, predict_growth

genes = read_fasta("coding_sequences.fna.gz")
highly_expressed = highly_expressed_mask(genes)
result = predict_growth(genes, highly_expressed)

print(result.doubling_time)  # hours
print(result.lower_ci, result.upper_ci)
print(result.model_name)  # exact fitted model selected
print(dict(result))  # R-compatible numerical keys
```

The annotation helper matches `ribosomal protein` case-insensitively,
excluding annotations containing `methyl` or `hydroxy`. Use `pattern=...`
or `ids=[...]` to supply your own selection. A plain list of DNA strings,
Biopython `SeqRecord` objects, or an insertion-ordered mapping of names to
sequences can also be passed to `predict_growth`.

Run a complete, self-contained synthetic example:

```console
python examples/synthetic_workflow.py
```

### Prediction modes

| Mode | Application | Main settings |
| --- | --- | --- |
| `full` | Individual prokaryotic genome | MILC, consistency, codon-pair bias |
| `partial` | Genome with less pair information | MILC and consistency |
| `metagenome_v1` | Community estimate | Pooled or explicit individual background |
| `metagenome_v2` | Community estimate with automatic background selection | Switches at consistency ≥ 0.6 |
| `eukaryote` | Microbial eukaryotic genome | Temperature is strongly recommended |
| `metagenome_euk` | Microbial eukaryotic community | Automatic background selection |
| `metagenome_150bp` | Short coding fragments | 150-bp trimming; five codon metrics |
| `metagenome_240bp` | Short coding fragments | 240-bp trimming; five codon metrics |

The API also exposes the original experimental `meta_testing` and
`meta_nogc_testing` modes. They are retained for reproducibility.

```python
from grodon import read_coverage

depth = read_coverage("coverage.tsv", genes)  # header: id<TAB>coverage
result = predict_growth(
    genes,
    highly_expressed,
    mode="metagenome_v2",
    depth_of_coverage=depth,
    random_state=42,
)
```

Coverage is aligned by FASTA ID, regardless of table order. Pooled models
retain the original unweighted codon reference, while coverage weights the
summary statistics and GC fraction. Individual-background models sample and shuffle genes, so set
`random_state` for repeatable results. R and NumPy do not generate identical
random streams from an identical integer seed.

`training_set` supports `madin` (default), `vs` (original), and `AOA_NOB` for
the relevant prokaryotic modes. Metagenome-v2 and short-read modes use Madin;
eukaryotic modes use their dedicated training data. `temperature` is in degrees
Celsius. `genetic_code` defaults to 11 for prokaryotes and 1 for eukaryotes.

By default, sequences must be longer than 240 bp, divisible by three, and
contain only A/C/G/T. `fragments=True` lowers the strict threshold to 120 bp;
`fragments=None` removes the minimum-length restriction while retaining frame
and alphabet checks. Short-read modes automatically apply their own trimming
rules, with `exclude_short=True` optionally removing shorter input sequences.
Metadata are filtered and reordered together with their sequences.

Results include `d`, `LowerCI`, and `UpperCI`, plus the statistics used. These
are **confidence intervals for the fitted mean**, not prediction intervals
for individual observations. Estimates beyond the original model's useful
range produce warnings. The estimates describe maximum growth potential;
they are not measurements of current growth in a sample. See the original
papers for scientific assumptions and limitations.

### Index of Copiotrophy

```python
import pandas as pd
from grodon import fit_ioc, predict_ioc

traits = pd.DataFrame(
    {
        "dCUB": [-0.8, -0.4, -0.6, -0.2, -0.9],
        "nCAZy": [90, 40, 100, 30, 160],
        "nGenes": [3500, 2200, 4200, 1900, 5000],
    }
)
fitted = fit_ioc(traits, completeness=95)
print(fitted.values)
print(predict_ioc(traits, model=fitted.model, completeness=95))
print(predict_ioc(traits, model="permafrost", completeness=95))
```

Built-in models: `permafrost`, `pacific`, `soil`, `human`, `all_habitat`.
Completeness is a percentage, either scalar or one value per row. The PCA
model retains its training center, scale, orientation, and score offset.
New samples may legitimately receive scores below one. Models support
`to_dict()` / `IoCModel.from_dict(...)` for JSON storage.

### Statistics and model fitting

`get_codon_statistics` provides MILC, consistency, genomic bias, GC content,
codon-pair bias, and filtering counts. `all_metrics=True` returns MILC,
ENCprime, B, SCUO, and MCB for short-read work. Lower-level numerical routines
are available in `grodon.codon`.

`get_statistics` selects highly expressed genes from an annotated FASTA;
`get_statistics_batch` processes a directory with optional parallel workers.
`fit_models`, `fit_gc_models`, and `fit_trim_models` fit the original model
families from tables of measured doubling times and codon statistics.
`fit_linear_model` provides a reusable fitted model with coefficients,
covariance, Box-Cox exponent, and confidence intervals. See
[the API guide](docs/API.md) for details.

R-style entry points `predictGrowth`, `fitIoC`, and `predictIoC` are supported,
along with the fitting helper aliases. Prefer snake_case names in new code.

## Command line

```console
grodon predict coding_sequences.fna.gz --mode full -o prediction.json
grodon predict coding_sequences.fna.gz --mode metagenome_v2 --coverage coverage.tsv --seed 42
grodon statistics coding_sequences.fna.gz --all-metrics -o statistics.json
grodon ioc traits.csv --model permafrost -o ioc.json
grodon fit-ioc traits.csv --save-model custom_ioc.json -o training_scores.json
grodon ioc traits.csv --model custom_ioc.json
```

`--he-ids` accepts a text file of FASTA IDs. `--he-pattern` sets the annotation
regular expression. JSON goes to stdout or `--output`; warnings go to stderr.
Undefined numerical values become JSON `null`. `python -m grodon` provides
the same commands. Use `grodon <command> --help` for every option.

## Compatibility and validation

The port bundles 36 original growth regressions and five original IoC PCA
models. Coefficients and uncertainty calculations are extracted from fitted
R objects, rather than re-estimated. Committed synthetic reference fixtures
allow normal tests to run without R. Export scripts regenerate those fixtures
against the R implementation. See [compatibility notes](docs/COMPATIBILITY.md)
for calibration conventions, deliberate bug fixes, and validation scope.

The raw-genome pipeline with `lactococcus-v1` reproduces the archived
CIRM-BIA2698 analysis's gene counts and all 12 Python output statistics. See
the [pipeline guide](docs/GENOME_PIPELINE.md#reproduce-the-prior-lactococcus-preparation)
for the exact input, settings, and reproduction command.

The [40-assembly test case](docs/TEST_GENOMES_VALIDATION.md) compares every
assembly in the supplied `reconstruct_gsmm/test_genomes` folder with original
R gRodon 2.7.3. All **480 comparisons passed** (12 statistics per assembly),
using the default full/Madin model and the same prepared CDS and ribosomal IDs.
The maximum doubling-time difference was `5.77e-15` hours. This is numerical
agreement within `rtol=1e-10, atol=1e-12`, rather than bitwise equality.
The report includes all results, input checksums, software versions, and commands
to repeat the analysis. It accounts for the folder's 3 DNA files and 37
protein-only files by retrieving DNA for the exact matching assembly accessions.

CI is configured for Python 3.11–3.13 on Linux, macOS, and Windows, including
tests, formatting, lint, type checking, distribution builds, and an installed
command smoke test. Local verification does not imply the remote CI matrix
has already run.

## Attribution and citation

Please cite the original scientific methods, as applicable:

- [Prokaryotic gRodon, Weissman et al. (2021)](https://doi.org/10.1073/pnas.2016810118).
- [Eukaryotic gRodon2](https://doi.org/10.1101/2021.10.15.464604).
- [Metagenome mode v2, Weissman et al. (2022)](https://doi.org/10.1128/msystems.00745-22).
- [Short-read prediction](https://doi.org/10.1101/2025.09.18.677194).
- [AOA/NOB models](https://doi.org/10.1126/science.ado0742).
- [Index of Copiotrophy](https://doi.org/10.1101/2025.09.01.673550).
- [MILC, Supek and Vlahoviček (2005)](https://doi.org/10.1186/1471-2105-6-182),
  and [coRdon](https://bioconductor.org/packages/coRdon).

The Python implementation retains this repository's MIT license. See
[NOTICE](NOTICE) for upstream attribution and the unresolved license metadata
on the original model assets before redistributing them publicly. Bundled
UniProt reference data carry CC BY 4.0 attribution; the reference manifests
record source releases and modifications.
