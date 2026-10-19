"""Versioned, offline protein reference sets for ribosomal annotation.

The bundled files are immutable snapshots. Resolving a profile never contacts a
remote service. Custom FASTA files must contain only trusted ribosomal proteins;
their biological annotation is the caller's responsibility.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import tempfile
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from importlib import resources
from os import PathLike
from pathlib import Path
from typing import Any

from .io import GeneSequences, read_fasta

DEFAULT_REFERENCE = "bacteria-swissprot-v1"
REFERENCE_NAMES = (DEFAULT_REFERENCE, "lactococcus-v1")
_PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def reference_catalog() -> list[dict[str, Any]]:
    """List bundled provenance summaries without unpacking reference proteins."""
    directory = resources.files("grodon").joinpath("data", "references")
    fields = (
        "name",
        "version",
        "database",
        "database_release",
        "reference_count",
        "taxon_count",
        "scope",
        "license",
        "license_note",
        "license_url",
        "attribution",
        "fasta_sha256",
    )
    return [
        {key: metadata[key] for key in fields if key in metadata}
        for name in REFERENCE_NAMES
        for metadata in [json.loads(directory.joinpath(f"{name}.json").read_text(encoding="utf-8"))]
    ]


@dataclass(frozen=True)
class ReferenceInfo:
    """A validated plain FASTA and its provenance, valid inside its context.

    ``sha256`` hashes the exact bytes in ``path``, which is the file supplied to
    the aligner. Descriptions map FASTA identifiers to complete header text.
    Custom compressed or noncanonical FASTA is normalized in a temporary file;
    the original input byte checksum is retained in ``source_metadata``.
    """

    path: Path
    name: str
    sha256: str
    descriptions: Mapping[str, str]
    version: str
    source_metadata: Mapping[str, Any]

    @property
    def count(self) -> int:
        """Number of reference proteins."""
        return len(self.descriptions)


def sha256_file(path: str | PathLike[str]) -> str:
    """Hash file bytes without reading the entire file into memory."""
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def validate_protein_fasta(path: str | PathLike[str]) -> GeneSequences:
    """Validate unique IDs, nonempty records, and IUPAC amino-acid letters.

    A single terminal stop is accepted and removed in the returned sequences;
    internal stops, gaps, digits, punctuation, and control characters in labels
    are errors. Ambiguous amino-acid symbols B, J, X, Z, U, and O are accepted.
    This validates file syntax, not whether proteins really are ribosomal.
    """
    proteins = read_fasta(path)
    sequences = []
    for identifier, description, sequence in zip(
        proteins.ids, proteins.names, proteins.sequences, strict=True
    ):
        if any(ord(character) < 32 or ord(character) == 127 for character in description):
            raise ValueError(f"Control character in protein FASTA description: {identifier!r}")
        sequence = sequence.removesuffix("*")
        if not sequence:
            raise ValueError(f"Empty protein sequence for {identifier!r}")
        invalid = sorted(set(sequence).difference(_PROTEIN_ALPHABET))
        if invalid:
            raise ValueError(f"Invalid protein sequence for {identifier!r}: {''.join(invalid)!r}")
        sequences.append(sequence)
    return GeneSequences(tuple(sequences), proteins.names)


def _canonical_fasta(proteins: GeneSequences) -> bytes:
    return "".join(
        f">{description}\n{sequence}\n"
        for description, sequence in zip(proteins.names, proteins.sequences, strict=True)
    ).encode("utf-8")


@contextmanager
def resolve_reference(
    name_or_path: str | PathLike[str] = DEFAULT_REFERENCE,
) -> Iterator[ReferenceInfo]:
    """Resolve a bundled profile or custom protein FASTA, without network access.

    The context owns any temporary file. Keep it open until alignment finishes.
    Bundled compressed bytes, decompressed FASTA bytes, and record counts are
    checked against the shipped provenance manifest before use.
    """
    name = str(name_or_path)
    if name in REFERENCE_NAMES:
        directory = resources.files("grodon").joinpath("data", "references")
        metadata = json.loads(directory.joinpath(f"{name}.json").read_text(encoding="utf-8"))
        if metadata.get("name") != name or metadata.get("schema_version") != 1:
            raise ValueError(f"Invalid reference metadata for {name!r}")
        asset_name = metadata.get("asset_filename")
        if asset_name != f"{name}.faa.gz":
            raise ValueError(f"Invalid reference asset filename for {name!r}")
        compressed = directory.joinpath(asset_name).read_bytes()
        if hashlib.sha256(compressed).hexdigest() != metadata["asset_sha256"]:
            raise ValueError(f"Reference asset SHA256 mismatch for {name!r}")
        fasta = gzip.decompress(compressed)
        checksum = hashlib.sha256(fasta).hexdigest()
        if checksum != metadata["fasta_sha256"]:
            raise ValueError(f"Reference FASTA SHA256 mismatch for {name!r}")
        with tempfile.TemporaryDirectory(prefix="grodon-reference-") as temporary:
            path = Path(temporary) / f"{name}.faa"
            path.write_bytes(fasta)
            proteins = validate_protein_fasta(path)
            if len(proteins) != metadata["reference_count"]:
                raise ValueError(f"Reference record count mismatch for {name!r}")
            yield ReferenceInfo(
                path=path,
                name=name,
                sha256=checksum,
                descriptions=dict(zip(proteins.ids, proteins.names, strict=True)),
                version=str(metadata["version"]),
                source_metadata=metadata,
            )
        return

    path = Path(name_or_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Reference {name!r} is neither a file nor a bundled profile "
            f"({', '.join(REFERENCE_NAMES)})"
        )
    proteins = validate_protein_fasta(path)
    fasta = _canonical_fasta(proteins)
    checksum = hashlib.sha256(fasta).hexdigest()
    metadata = {
        "schema_version": 1,
        "name": path.name,
        "kind": "custom",
        "source_path": str(path),
        "source_sha256": sha256_file(path),
        "fasta_sha256": checksum,
        "reference_count": len(proteins),
        "normalization": "Uppercase amino acids; remove one terminal stop; one sequence line; LF",
        "scope": "Caller-supplied ribosomal protein references; annotation not independently checked",
    }
    with tempfile.TemporaryDirectory(prefix="grodon-reference-") as temporary:
        normalized = Path(temporary) / "custom-reference.faa"
        normalized.write_bytes(fasta)
        yield ReferenceInfo(
            path=normalized,
            name=path.name,
            sha256=checksum,
            descriptions=dict(zip(proteins.ids, proteins.names, strict=True)),
            version="custom",
            source_metadata=metadata,
        )
