"""FASTA and aligned metadata readers with no dependency on an R installation."""

from __future__ import annotations

import csv
import gzip
import re
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import TextIO, overload

DEFAULT_HE_PATTERN = r"^(?!.*(?:methyl|hydroxy)).*ribosomal protein"


@dataclass(frozen=True)
class GeneSequences(Sequence[str]):
    """Ordered coding sequences and complete FASTA descriptions.

    ``ids`` are the first whitespace-delimited tokens of the descriptions.
    Passing this object to ``predict_growth`` preserves the FASTA order.
    """

    sequences: tuple[str, ...]
    names: tuple[str, ...]

    def __post_init__(self) -> None:
        if len(self.sequences) != len(self.names):
            raise ValueError("sequences and names must have the same length")

    def __len__(self) -> int:
        return len(self.sequences)

    @overload
    def __getitem__(self, index: int) -> str: ...

    @overload
    def __getitem__(self, index: slice) -> tuple[str, ...]: ...

    def __getitem__(self, index: int | slice) -> str | tuple[str, ...]:
        return self.sequences[index]

    def __iter__(self) -> Iterator[str]:
        return iter(self.sequences)

    @property
    def ids(self) -> tuple[str, ...]:
        return tuple(name.split()[0] for name in self.names)


def _parse_fasta(stream: TextIO) -> GeneSequences:
    names: list[str] = []
    sequences: list[str] = []
    chunks: list[str] = []
    for number, line in enumerate(stream, 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if names:
                if not chunks:
                    raise ValueError(f"Empty sequence for {names[-1]!r}")
                sequences.append("".join(chunks).upper())
            name = line[1:].strip()
            if not name:
                raise ValueError(f"Empty FASTA header on line {number}")
            names.append(name)
            chunks = []
        else:
            if not names:
                raise ValueError(f"Sequence before first FASTA header on line {number}")
            chunks.append("".join(line.split()))
    if not names:
        raise ValueError("FASTA contains no sequences")
    if not chunks:
        raise ValueError(f"Empty sequence for {names[-1]!r}")
    sequences.append("".join(chunks).upper())
    genes = GeneSequences(tuple(sequences), tuple(names))
    if len(set(genes.ids)) != len(genes):
        raise ValueError("FASTA identifiers must be unique")
    return genes


def read_fasta(source: str | PathLike[str] | TextIO) -> GeneSequences:
    """Read plain or gzip-compressed FASTA, normalizing bases to uppercase.

    Ambiguous bases and out-of-frame lengths are retained here; the statistics
    pipeline applies the documented sequence filter and reports removals.
    """
    if hasattr(source, "read"):
        return _parse_fasta(source)  # type: ignore[arg-type]
    path = Path(source)  # type: ignore[arg-type]
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as stream:
        return _parse_fasta(stream)


def highly_expressed_mask(
    genes: GeneSequences,
    *,
    pattern: str = DEFAULT_HE_PATTERN,
    ids: Sequence[str] | None = None,
) -> list[bool]:
    """Select annotations by case-insensitive regex, or an explicit set of IDs."""
    if ids is not None:
        selected = set(ids)
        unknown = selected.difference(genes.ids)
        if unknown:
            raise ValueError(f"Highly expressed IDs missing from FASTA: {sorted(unknown)[:5]}")
        return [identifier in selected for identifier in genes.ids]
    expression = re.compile(pattern, re.IGNORECASE)
    return [expression.search(name) is not None for name in genes.names]


def read_coverage(
    path: str | PathLike[str],
    genes: GeneSequences,
    *,
    id_column: str = "id",
    coverage_column: str = "coverage",
    delimiter: str = "\t",
) -> list[float]:
    """Read a table and align coverage by FASTA ID, never by table row order."""
    import math

    values: dict[str, float] = {}
    with open(path, encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        if not reader.fieldnames or not {id_column, coverage_column}.issubset(reader.fieldnames):
            raise ValueError(f"Coverage table must contain {id_column!r} and {coverage_column!r}")
        for row in reader:
            identifier = row[id_column]
            if identifier in values:
                raise ValueError(f"Duplicate coverage ID: {identifier}")
            try:
                coverage = float(row[coverage_column])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid coverage for {identifier}") from exc
            if not math.isfinite(coverage) or coverage < 0:
                raise ValueError(f"Coverage for {identifier} must be finite and nonnegative")
            values[identifier] = coverage
    missing = set(genes.ids).difference(values)
    if missing:
        raise ValueError(f"Coverage is missing FASTA IDs: {sorted(missing)[:5]}")
    return [values[identifier] for identifier in genes.ids]
