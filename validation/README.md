# Public benchmark evidence

The [40-assembly benchmark report](../docs/TEST_GENOMES_VALIDATION.md) records
480 passing comparisons with original R gRodon 2.7.3. It includes the complete
results table, model settings, software versions, tolerances, and commands to
download the exact genomic inputs and repeat the analysis locally.

| File | Purpose |
| --- | --- |
| `test_genomes_comparison.csv` | One row per assembly, with R/Python doubling times, confidence intervals, differences, and outcomes |
| `test_genomes_comparison.json` | Full numerical results and the recorded audit of input identity, settings, counts, and all 680 output artifacts |
| `test_genomes_inputs_manifest.json` | Exact assembly accessions, official download URLs, input checksums, and preparation provenance |
| `pipeline_reproduction.json` | Earlier Lactococcus reproduction audit for the exact and broad reference profiles |
| `../tests/data/folder_growth_reference.json` | Compact original-R numerical reference for the 40 offline model regression cases |

The JSON and CSV files preserve the original recorded evidence byte for byte.
Historical paths identify where the recorded analysis ran; they are provenance,
not installation requirements or portable paths to open. The public reproduction
commands create a new input manifest with paths for the current computer.
Git attributes preserve archived line endings so the manifest/report/fixture
SHA256 links remain valid on all platforms.

The original folder contained 3 genomic FASTAs and 37 protein-only FASTAs.
Preparation fetched the exact matching nucleotide assemblies for those 37
accessions. Public reproduction starts directly from all 40 genomic FASTAs,
verified against their recorded SHA256 values. Both approaches supply the same
nucleotide genome bytes to the gene-calling and growth-prediction pipeline.

Genome downloads, full proteomes, predicted coding sequences, alignments,
per-run outputs, installed tools, local environments, and one-off development
scripts are excluded from Git and the source distribution. The curated ribosomal
protein snapshots under `src/grodon/data/references` are required runtime data
and are included with their attribution and provenance.

The numeric model tests run without downloading genomes or installing R:

```console
python -m pytest tests/test_folder_growth_reference.py
```

These tests evaluate the Python model using recorded R codon statistics. To
repeat gene prediction, annotation, codon statistics, and the live R comparison,
follow the full benchmark recipe. The historical artifact hashes document the
completed audit; rechecking raw artifacts requires retaining or regenerating
the corresponding local run directories. Use fresh output names such as
`validation/reproduced_comparison.json` so the archived benchmark stays intact.
