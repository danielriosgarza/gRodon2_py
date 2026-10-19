# Reproducible bacterial genome pipeline

`grodon genome` takes a bacterial genome assembly, calls its coding sequences,
identifies ribosomal proteins, and estimates maximum growth potential. The
pipeline uses Pyrodigal, local BLASTP, a versioned protein reference snapshot,
and the native Python gRodon implementation. No LLM, API key, or online service
is involved in analysis. Installation and an explicitly requested tool download
are the only steps that need network access.

Supply one isolate genome per FASTA file, including its chromosome and plasmids
if desired. Both plain FASTA and `.gz` FASTA are accepted. Input files are read
without modification. The workflow is for bacterial assemblies; it does not
establish a validated general assembly-annotation pipeline for archaea,
eukaryotes, viruses, or mixed microbial communities.

## Install in this checkout

Python 3.11 or later is required. The package has not been published to PyPI;
install from this repository or a wheel built from it.

These commands run from the repository root in Windows PowerShell and do not
require shell activation:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[genome]"
.\.venv\Scripts\python.exe -m grodon setup --tools-dir .tools\blast-2.17.0 --version 2.17.0
.\.venv\Scripts\python.exe -m grodon doctor --blastp .tools\blast-2.17.0\bin\blastp.exe
```

The `genome` extra installs Pyrodigal. `setup` explicitly downloads the pinned
NCBI BLAST+ release into the chosen directory, checks the publisher's archive
checksum, records SHA256 hashes, and verifies the executable version. It does
not alter PATH or system settings. Repeating setup with an intact matching
installation verifies and reuses it. An unrelated or altered installation is
not overwritten. The download is an official BLAST archive, whose size can be
larger than the files retained for BLASTP.

`doctor` prints JSON showing dependency versions and a `ready` field; it does
not install or download anything. It exits with status 0 when ready and 1 when
a dependency is unavailable. If BLASTP is already installed, pass its
path to `--blastp` and skip setup. Use the explicit path consistently for all
subsequent commands.

On Linux or macOS, the corresponding commands are:

```console
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[genome]"
.venv/bin/python -m grodon setup --tools-dir .tools/blast-2.17.0 --version 2.17.0
.venv/bin/python -m grodon doctor --blastp .tools/blast-2.17.0/bin/blastp
```

The managed installer supports Windows x64, Linux x64/ARM64, and macOS
x64/ARM64 when the requested NCBI build is available. An externally installed
BLASTP can also be used. For exact reproduction of the archived Lactococcus
preparation, use Pyrodigal 3.7.1 and BLAST+ 2.17.0.

## Analyze one genome

```powershell
.\.venv\Scripts\python.exe -m grodon genome "assembly.fna.gz" --outdir results\sample --blastp .tools\blast-2.17.0\bin\blastp.exe
```

The command prints a JSON summary and saves the full analysis under `--outdir`.
Choose a fresh output directory for each run. Existing analysis files are not
overwritten. A failed run may retain diagnostics and intermediate products;
inspect those and use a new output directory for a retry.

Default processing is:

1. Train Pyrodigal on all contigs from this genome and predict coding DNA and
   proteins with translation table 11, `closed=True`, and `mask=True`.
2. Compare the predicted proteins with the bundled ribosomal reference using
   local BLASTP, low-complexity filtering disabled, composition-based score
   adjustment disabled, and at most one HSP per query/reference pair. The BLAST
   target limit is explicitly set to the full reference count so that a default
   top-hit cap does not remove candidates before the coverage and identity filters.
3. Retain matches with at least 70% amino-acid identity, 80% query coverage,
   80% reference coverage, and E-value at most `1e-8`. Coverage is the aligned
   coordinate span divided by the full protein length. Choose the greatest
   bit-score hit per query; ties use reference FASTA order.
4. Apply gRodon's coding-sequence filter: length strictly greater than 240 bp,
   length divisible by three, and only A/C/G/T. Require at least 10 ribosomal
   genes **after** this filter, then run the full Madin model with 95% confidence
   intervals. Ribosomal genes are a proxy for high expression, not an expression
   measurement.

Gene prediction occurs once. Every subsequent calculation uses those same
predicted coding sequences and ribosomal IDs. Reference proteins identify
homologs; they never substitute for the query genome's coding sequences.

## Reference selection

```powershell
.\.venv\Scripts\python.exe -m grodon references
```

| Profile | Intended use | Snapshot |
| --- | --- | --- |
| `bacteria-swissprot-v1` | Default broad bacterial reference | 39,743 complete canonical ribosomal proteins across 934 taxon IDs, UniProtKB/Swiss-Prot 2026_03 |
| `lactococcus-v1` | Reproduce the prior Lactococcus analysis | Exact 88 RefSeq proteins from MG1363, IL1403, and CIRM-BIA2553 |

Both snapshots ship with the package. They are resolved offline and checked
against compressed and uncompressed SHA256 hashes before alignment. The broad
snapshot derives from reviewed bacterial UniProt entries with keyword KW-0689;
the export excludes fragments, probable or associated proteins, unsupported
labels, and nonstandard amino-acid sequences. Its scope is wider than the
Lactococcus profile, but neither the reference collection nor the default
homology thresholds guarantee sufficient matches for every bacterial lineage.

Use a custom curated ribosomal protein FASTA when appropriate:

```powershell
.\.venv\Scripts\python.exe -m grodon genome "assembly.fna.gz" --outdir results\custom-reference --blastp .tools\blast-2.17.0\bin\blastp.exe --reference "curated_ribosomal.faa"
```

Custom FASTA IDs must be unique and sequences must be valid protein sequences.
The resolver normalizes amino-acid case and removes one terminal stop; internal
stops and gaps are rejected. It records both original-input and normalized
checksums. The pipeline checks syntax; the user supplies the biological
curation. See the [reference snapshot documentation](../src/grodon/data/references/README.md)
for export rules, attribution, and how to rebuild a snapshot explicitly.

## Reproduce the prior Lactococcus preparation

For assembly `GCA_964062575.1`, strain CIRM-BIA2698, supply its genomic FASTA:

```powershell
.\.venv\Scripts\python.exe -m grodon genome "path\to\GCA_964062575.1_genomic.fna.gz" --outdir results\lactococcus-reproduction --blastp .tools\blast-2.17.0\bin\blastp.exe --reference lactococcus-v1
```

The reference profile and default preparation settings match the
[archived reproduction audit](../validation/pipeline_reproduction.json).
The [40-assembly benchmark](TEST_GENOMES_VALIDATION.md#reproduce-this-test-case)
includes a download command for this exact assembly and the other 39 inputs.
The new pipeline reproduced all 12 archived Python output statistics using
the same genome and reference profile: 2,832 predicted CDS, 54 selected
ribosomal genes, and 2,600 CDS with 42 ribosomal genes retained by gRodon's
filter. The doubling-time estimate was `0.12051619413002182` hours. This checks
software and input reproduction; it does not establish a measured culture
doubling time. Use the same source genome bytes and tool versions when
assessing reproduction.

`lactococcus-v1` is deliberately specific to the earlier reference choice.
The broad default was also run on this genome: it selected 53 ribosomal genes,
retained 41 after filtering, and predicted `0.12129379522429953` hours
(7.27763 minutes). One S1 homolog selected with the Lactococcus reference did
not pass the broad profile's homology thresholds. This observed difference
shows why reference choice is recorded and why exact reproduction requires
`lactococcus-v1`.

## Compare the same inputs with R

Pass an Rscript executable whose library paths contain the original `gRodon`
package and its dependencies, plus `Biostrings` and `jsonlite`:

```powershell
.\.venv\Scripts\python.exe -m grodon genome "assembly.fna.gz" --outdir results\with-r --blastp .tools\blast-2.17.0\bin\blastp.exe --rscript "C:\path\to\R\bin\Rscript.exe"
```

R and its libraries are installed separately. R uses the exact same saved coding sequences,
ribosomal IDs, and matching model settings. The comparison records individual
statistics and numerical differences. The original R prediction interface
uses 95% confidence intervals, so R comparison requires `--confidence 0.95`.
Python-only predictions can use other supported confidence levels.

The Python prediction is saved before the optional R process starts. If R
cannot run, `result.json` and `result.csv` retain the completed Python result
with status `prediction_complete_r_failed`; `run.json` records the failed
comparison and the command exits with an error. This preserves the prediction
when an R dependency or executable needs repair. A completed comparison records
`r_matches` and per-statistic differences; numerical disagreement is reported
in those results rather than treated as an execution exception.

R is optional. Dependency installation and R package setup are separate from
analysis; the pipeline does not install R packages during a comparison.

## Batch processing

```powershell
.\.venv\Scripts\python.exe -m grodon genome-batch "one.fna.gz" "two.fna.gz" --outdir results\batch --blastp .tools\blast-2.17.0\bin\blastp.exe --continue-on-error
```

Each assembly receives a separate output directory with the same configuration.
Batch `summary.json` and `summary.csv` record outcomes. With
`--continue-on-error`, a failed assembly is recorded and processing continues;
without it, processing stops on an error. `--threads` sets BLASTP threads for
each assembly; it is not a request to infer one community-wide prediction from
multiple files.

A batch completed with failed inputs prints its summary and exits with status
1, including when `--continue-on-error` is used. Successful commands exit with
status 0; argument errors and execution exceptions exit with status 2. These
codes let external scripts detect failed work without parsing log text.

## Python API

```python
from pathlib import Path
from grodon.pipeline import GenomePipelineConfig, run_genome, run_genomes

config = GenomePipelineConfig(
    reference="bacteria-swissprot-v1",
    blastp=str(Path(".tools/blast-2.17.0/bin/blastp.exe").resolve()),
    genetic_code=11,
    gene_mode="single",
    closed=True,
    mask=True,
    min_identity=70,
    min_query_coverage=0.8,
    min_reference_coverage=0.8,
    max_evalue=1e-8,
    min_ribosomal=10,
    mode="full",
    training_set="madin",
    confidence=0.95,
    threads=1,
    timeout=600,
)

result = run_genome("assembly.fna.gz", "results/python-sample", config=config)
print(result.prediction.doubling_time)  # hours
print(result.output_dir)
print(result.to_dict())

summary = run_genomes(
    ["one.fna.gz", "two.fna.gz"],
    "results/python-batch",
    config=config,
    continue_on_error=True,
)
print(summary)
```

Pass `rscript="..."` to either function for the optional R comparison. The CLI
and Python API use the same implementation and produce the same output layout.

## Settings and limitations

- **Genetic code:** Table 11 is the default. Set `--genetic-code 4` when that is
  the known translation table of the input organism. The same table is used for
  gene calling, translation, and gRodon statistics; taxonomy is not inferred
  automatically to choose a code. This assembly pipeline supports tables 4 and 11.
- **Gene-calling mode:** `--gene-mode single` trains on the supplied assembly.
  For a short or fragmented assembly from one isolate, `--gene-mode meta`
  explicitly selects pretrained gene-calling models for code 4 or 11. This is
  a gene-caller setting, not a validated community growth prediction workflow.
- **Contig ends:** `--closed` is on by default and excludes partial genes that
  run beyond contig boundaries. `--no-closed` allows partial predictions. Circular
  contig-origin spanning genes are not reconstructed by this workflow.
- **Unknown sequence:** `--mask` is on by default; `--no-mask` changes the
  gene-caller's treatment of ambiguous sequence gaps. Invalid FASTA and sequence
  characters fail with a diagnostic.
- **Homology:** Identity is a percentage; both coverage thresholds are fractions
  from 0 to 1. The thresholds and `--min-ribosomal` are explicit and recorded.
  Insufficient retained ribosomal genes stops prediction instead of silently
  returning an estimate from an inadequate reference set. The 10-gene minimum
  is an operational quality gate, not a guarantee of biological accuracy.
- **Prediction:** `--mode full` and `--mode partial` are available for this
  bacterial genome workflow. `--training-set` selects `madin`, `vs`, or
  `AOA_NOB`; `--temperature` supplies optimal growth temperature in Celsius when
  supported by the selected model. Model suitability remains a scientific choice.
- **Reproducibility:** Pin the package, reference profile, Pyrodigal, BLASTP, and
  numerical-library versions for repeat work. Run manifests record the versions
  used. Changing reference content or tool versions can change gene predictions,
  homolog assignments, and resulting growth estimates. Manifests and logs contain
  run-specific paths and timing and need not be byte-identical between runs.
  The result's `scientific_fingerprint` summarizes genome/reference/CDS checksums,
  selected gene IDs, the model, and numerical statistics independently of output
  paths and timestamps.

Estimates describe maximum growth potential, not observed growth under a
particular culture condition. Confidence intervals describe uncertainty in the
fitted model mean; they do not capture all uncertainty introduced by an assembly,
gene calling, reference coverage, or annotation errors.

## Validate a complete genome folder against R

The [documented 40-assembly test case](TEST_GENOMES_VALIDATION.md) applies this
pipeline to every assembly represented in the supplied test folder. All 12
statistics matched original R gRodon 2.7.3 for every assembly, using the default
full/Madin model with `bacteria-swissprot-v1`. Both languages receive the same
ordered CDS and ribosomal IDs; each computes its own codon statistics and growth
prediction. The annotation stage is shared, so this tests R/Python calculation
equivalence given a reproducible preparation, not independent annotation methods.

The test case provides three reusable scripts: prepare and verify the exact
assembly inputs, run the paired predictions, and independently audit the saved
results. It includes all 40 outcomes and full reproduction commands. The 37
protein-only inputs require matching nucleotide assemblies because synonymous
codons cannot be recovered from amino-acid sequences. Raw runs and downloaded
genomes stay in the local validation cache; compact reports and numeric R
fixtures are included in the repository and source distribution.

## Output files

| File | Content |
| --- | --- |
| `run.json` | Run configuration, input identity, software versions, provenance, warnings, and outcome |
| `result.json` | Python prediction and its accompanying metadata |
| `result.csv` | Tabular numerical result |
| `reference.json` | Reference profile or custom-reference provenance and checksums |
| `reference.faa.gz` | Copy of the actual normalized ribosomal reference used |
| `predicted_cds.fna` | Query-derived nucleotide coding sequences |
| `predicted_proteins.faa` | Translated query proteins |
| `predicted_genes.tsv` | Gene IDs, contigs, coordinates, strand, length, and partial-boundary flags |
| `gene_calling.json` | Gene-calling settings, counts, input mapping, and hashes |
| `ribosomal_hits.tsv` | Raw BLASTP alignments |
| `selected_ribosomal.tsv` | Best qualifying hit and supporting evidence for each selected query |
| `highly_expressed_ids.txt` | Selected ribosomal query IDs before gRodon filtering |
| `ribosomal_annotation.json` | BLASTP settings, command, version, counts, and hashes |
| `blastp.log` | BLASTP execution diagnostics |
| `r_parameters.json` | Optional R comparison settings and shared input paths |
| `r_result.json` | Optional result from the original installed R package |
| `r_comparison.json` | Optional R/Python numerical differences and agreement checks |
| `r_comparison.log` | Optional R execution diagnostics |

The four `r_` files are produced when `--rscript` is supplied. Batch output
additionally has `summary.json`, `summary.csv`, and
the separate per-genome output directories. Keep the full directory with the
original genome when archiving an analysis.
