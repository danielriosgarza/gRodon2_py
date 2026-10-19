# Offline ribosomal reference snapshots

These reference proteins identify homologs in query genomes. They do not replace
query sequences in growth predictions. `grodon.references.resolve_reference`
checks compressed and uncompressed SHA256 values and record counts at runtime.
No downloads occur during prediction.

| Profile | Source | Scope |
| --- | --- | --- |
| `bacteria-swissprot-v1` | UniProtKB/Swiss-Prot 2026_03 | 39,743 complete canonical bacterial ribosomal proteins; broad default |
| `lactococcus-v1` | NCBI RefSeq, archived 2026-10-19 | The exact 88 proteins used for the CIRM-BIA2698 R/Python comparison |

The default snapshot is an export of the official UniProt query
`(taxonomy_id:2) AND (keyword:KW-0689) AND (reviewed:true)`. The source had 39,880
entries. The reproducible name and sequence filters exclude fragments, probable
proteins, homolog labels, fusion proteins, ribosome-associated proteins, and
nonstandard amino-acid sequences. Full source descriptions and identifiers are
preserved. Each JSON manifest records the query, URL, release, retrieval date,
selection rule, counts, checksums, and attribution.

The broad profile does not guarantee coverage of every bacterial lineage.
Reviewed reference content and homology thresholds determine sensitivity.
The Lactococcus profile is intended for exact reproduction and related organisms;
it is not a general bacterial reference. Custom protein FASTA files are supported
when users need their own curated reference collection.

## Attribution

UniProt data: The UniProt Consortium, UniProtKB/Swiss-Prot. Copyrightable parts
are available under [CC BY 4.0](https://www.uniprot.org/help/license). This export
selects records, removes unsupported entries, sorts identifiers, normalizes FASTA
line wrapping, and compresses the file. Cite the source database and the snapshot
release in downstream work. See [UniProt citations](https://www.uniprot.org/citations).

NCBI data: NCBI RefSeq and the contributing genome projects listed in the
`lactococcus-v1.json` provenance. NCBI makes sequence data available for public
use; its policy notes that submitters may assert intellectual-property rights.
See the [NCBI data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/).

## Maintenance

The package source includes `scripts/build_ribosomal_references.py`. Explicit
`--fetch` downloads the current UniProt release; this is a maintenance action,
never a prediction dependency. Save the raw source and response headers using
`--save-source-directory` when making a new snapshot. An archived source can be
re-exported offline with `--raw-fasta`, `--raw-metadata`, and
`--expected-source-sha256`. The checksum rejects a changed upstream snapshot.
Existing profile files are not overwritten unless `--force` is given.

The Lactococcus snapshot and its complete provenance ship with the package;
the original raw validation archive is not distributed. Maintainers with that
archive can regenerate the profile by explicitly supplying
`--lactococcus-validation-directory path/to/archive`, containing
`ribosomal_reference.faa` and `ribosomal_reference_metadata.json`.
The exact FASTA-byte checksum remains
`075848596f8e2f356d8c6ad244cc24ab27ea376a5305b32a59f790e029bb8d80`.
New scientific reference content should receive a new named version rather than
silently changing an existing profile. `gzip` encoding can differ between Python
versions; the decompressed FASTA checksum identifies the scientific input.
