"""Explicitly export pinned ribosomal reference snapshots for offline use.

This maintenance script is never run by the prediction pipeline. Recreate a
Swiss-Prot snapshot from the saved raw FASTA and response metadata, or explicitly
fetch a current release. Updating a database can change annotations and requires
a new profile version. The archived Lactococcus validation inputs reproduce the
original 88-protein reference byte for byte.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from grodon.io import read_fasta

QUERY = "(taxonomy_id:2) AND (keyword:KW-0689) AND (reviewed:true)"
MAX_DOWNLOAD_BYTES = 20_000_000
CANONICAL_NAME = re.compile(
    r"(?:Small|Large) ribosomal subunit protein [ub][LS]\d+[A-Z]?(?:-\d+)?"
    r"(?:/[ub]?[LS]?\d+[A-Z]?(?:-\d+)?)*$"
    r"|(?:30S|50S) ribosomal protein [LS]\d+[A-Z]?(?:/[LS]?\d+[A-Z]?)*$"
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def export_uniprot(raw: bytes, source: dict[str, Any], name: str) -> tuple[bytes, dict]:
    """Select complete canonical ribosomal proteins in accession order."""
    proteins = read_fasta(io.StringIO(raw.decode("utf-8")))
    accepted: list[tuple[str, str]] = []
    excluded: Counter[str] = Counter()
    families: Counter[str] = Counter()
    taxa: set[str] = set()
    for description, sequence in zip(proteins.names, proteins.sequences, strict=True):
        label = description.split(" ", 1)[1].split(" OS=", 1)[0]
        if not CANONICAL_NAME.fullmatch(label):
            excluded[label] += 1
            continue
        if set(sequence).difference("ACDEFGHIKLMNPQRSTVWY"):
            excluded["noncanonical amino-acid sequence"] += 1
            continue
        accepted.append((description, sequence))
        families[label] += 1
        match = re.search(r" OX=(\d+)\b", description)
        if match:
            taxa.add(match[1])
    accepted.sort(key=lambda pair: pair[0].split()[0])
    fasta = "".join(f">{description}\n{sequence}\n" for description, sequence in accepted).encode()
    headers = {key.lower(): value for key, value in source.get("headers", {}).items()}
    if not headers.get("x-uniprot-release"):
        raise ValueError("Source metadata must record the X-UniProt-Release response header")
    metadata = {
        "schema_version": 1,
        "name": name,
        "version": "1",
        "kind": "bundled",
        "database": "UniProtKB/Swiss-Prot",
        "database_release": headers["x-uniprot-release"],
        "database_release_date": headers.get("x-uniprot-release-date"),
        "retrieved_utc": source.get("retrieved_utc", headers.get("date")),
        "source_url": source["url"],
        "source_query": QUERY,
        "source_fasta_sha256": sha256(raw),
        "source_reference_count": len(proteins),
        "reference_count": len(accepted),
        "taxon_count": len(taxa),
        "selection": (
            "Reviewed bacterial entries with keyword KW-0689; retain canonical small/large "
            "ribosomal L/S-numbered proteins matching the recorded regex and containing "
            "only the 20 standard amino acids. Exclude fragments, probable proteins, "
            "homolog labels, fusion proteins, and ribosome-associated proteins. "
            "Sort full FASTA identifiers ascending; preserve full source headers."
        ),
        "selection_regex": CANONICAL_NAME.pattern,
        "excluded_label_counts": dict(sorted(excluded.items())),
        "protein_label_counts": dict(sorted(families.items())),
        "scope": (
            "Broad bacterial structural ribosomal homology references. Coverage is limited "
            "by reviewed UniProt content; acceptance thresholds may miss divergent taxa. "
            "Reference proteins annotate query genes and never replace query sequences."
        ),
        "license": "CC-BY-4.0",
        "license_url": "https://www.uniprot.org/help/license",
        "attribution": "The UniProt Consortium, UniProtKB/Swiss-Prot",
        "citation_url": "https://www.uniprot.org/citations",
    }
    return fasta, metadata


def export_lactococcus(directory: Path, name: str) -> tuple[bytes, dict]:
    """Preserve the archived reference FASTA and its full source manifest."""
    raw = (directory / "ribosomal_reference.faa").read_bytes()
    source = json.loads((directory / "ribosomal_reference_metadata.json").read_text())
    if sha256(raw) != source["reference_fasta_sha256"]:
        raise ValueError("Archived Lactococcus FASTA does not match its manifest")
    metadata = {
        "schema_version": 1,
        "name": name,
        "version": "1",
        "kind": "bundled",
        "database": "NCBI RefSeq",
        "retrieved_utc": source["created_utc"],
        "reference_count": source["unique_reference_records"],
        "scope": (
            "Lactococcus only: original CIRM-BIA2698 R/Python validation snapshot from "
            "MG1363, IL1403, and CIRM-BIA2553. Do not use as a universal bacterial reference."
        ),
        "selection": source["selection"],
        "attribution": "NCBI RefSeq and the contributing genome projects",
        "license_url": "https://www.ncbi.nlm.nih.gov/home/about/policies/",
        "license_note": (
            "NCBI makes sequence data available for public use; its disclaimer notes "
            "that submitters may assert intellectual-property rights. See source policy."
        ),
        "original_provenance": source,
    }
    return raw, metadata


def write_bundle(output: Path, name: str, fasta: bytes, metadata: dict, force: bool) -> None:
    """Write deterministic compressed FASTA and its integrity manifest."""
    asset = output / f"{name}.faa.gz"
    manifest = output / f"{name}.json"
    if not force and (asset.exists() or manifest.exists()):
        raise FileExistsError(f"Refusing to overwrite existing version {name!r}; use --force")
    compressed = gzip.compress(fasta, compresslevel=9, mtime=0)
    metadata.update(
        asset_filename=asset.name,
        asset_sha256=sha256(compressed),
        fasta_sha256=sha256(fasta),
        fasta_size_bytes=len(fasta),
    )
    output.mkdir(parents=True, exist_ok=True)
    asset.write_bytes(compressed)
    manifest.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"name": name, "count": metadata["reference_count"], "sha256": sha256(fasta)}))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", choices=("bacteria-swissprot-v1", "lactococcus-v1"))
    parser.add_argument("--output-directory", type=Path, default=Path("src/grodon/data/references"))
    parser.add_argument("--raw-fasta", type=Path, help="Saved UniProt response FASTA")
    parser.add_argument("--raw-metadata", type=Path, help="Saved URL and response headers JSON")
    parser.add_argument("--fetch", action="store_true", help="Explicitly download current UniProt")
    parser.add_argument("--save-source-directory", type=Path)
    parser.add_argument("--expected-source-sha256", help="Reject a changed upstream snapshot")
    parser.add_argument(
        "--lactococcus-validation-directory",
        type=Path,
        help="Local original reference archive (not distributed with this repository)",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.profile == "lactococcus-v1":
        if args.lactococcus_validation_directory is None:
            parser.error("The Lactococcus exporter requires --lactococcus-validation-directory")
        fasta, metadata = export_lactococcus(args.lactococcus_validation_directory, args.profile)
    else:
        if args.fetch:
            if args.raw_fasta or args.raw_metadata:
                parser.error("--fetch cannot be combined with --raw-fasta or --raw-metadata")
            url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
                {"query": QUERY, "format": "fasta", "sort": "accession asc"}
            )
            with urllib.request.urlopen(url, timeout=120) as response:
                raw = response.read(MAX_DOWNLOAD_BYTES + 1)
                source = {
                    "url": url,
                    "headers": dict(response.headers),
                    "retrieved_utc": datetime.now(UTC).isoformat(),
                }
            if len(raw) > MAX_DOWNLOAD_BYTES:
                raise ValueError("Reference download exceeded the explicit 20 MB limit")
        elif args.raw_fasta and args.raw_metadata:
            raw = args.raw_fasta.read_bytes()
            source = json.loads(args.raw_metadata.read_text(encoding="utf-8"))
        else:
            parser.error("Provide --raw-fasta and --raw-metadata, or explicitly use --fetch")
        if args.expected_source_sha256 and sha256(raw) != args.expected_source_sha256:
            raise ValueError("Downloaded source SHA256 differs from the requested snapshot")
        if args.save_source_directory:
            args.save_source_directory.mkdir(parents=True, exist_ok=True)
            (args.save_source_directory / "uniprot-raw.faa").write_bytes(raw)
            (args.save_source_directory / "uniprot-raw.headers.json").write_text(
                json.dumps(source, indent=2) + "\n", encoding="utf-8"
            )
        fasta, metadata = export_uniprot(raw, source, args.profile)
    write_bundle(args.output_directory, args.profile, fasta, metadata, args.force)


if __name__ == "__main__":
    main()
