"""Local assembly gene calling and ribosomal-protein homology annotation.

Pyrodigal is imported only when gene calling is requested. BLAST runs as a local
process with explicit arguments; neither operation needs a network connection.
"""

from __future__ import annotations

import csv
import hashlib
import importlib
import json
import math
import re
import shutil
import subprocess
import tempfile
from numbers import Integral
from os import PathLike
from pathlib import Path
from typing import Any, TypedDict

from .io import GeneSequences, read_fasta

_DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
_PROTEIN_ALPHABET = frozenset("ABCDEFGHIKLMNPQRSTVWXYZJUO*")
_SAFE_ID = re.compile(r"[A-Za-z0-9_.:-]+\Z")
_OUTFMT = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
_HIT_COLUMNS = (
    "id",
    "reference",
    "description",
    "identity_percent",
    "query_coverage",
    "reference_coverage",
    "evalue",
    "bitscore",
)


class RibosomalHit(TypedDict):
    """Best qualifying reference hit for one query protein."""

    id: str
    reference: str
    description: str
    identity_percent: float
    query_coverage: float
    reference_coverage: float
    evalue: float
    bitscore: float


class GeneCallingResult(TypedDict):
    """Paths and provenance for a completed gene-calling run."""

    cds_path: Path
    proteins_path: Path
    gene_table_path: Path
    metadata_path: Path
    input: str
    input_sha256: str
    pyrodigal_version: str
    settings: dict[str, Any]
    counts_by_contig: dict[str, int]
    contig_id_map: dict[str, str]
    total_predicted_genes: int
    total_genome_bp: int
    cds_sha256: str
    proteins_sha256: str
    gene_table_sha256: str


class RibosomalAnnotationResult(TypedDict):
    """Selected query IDs, best hits and reproducible BLAST metadata."""

    selected_ids: list[str]
    selected_hits: list[RibosomalHit]
    hits_path: Path
    selected_hits_path: Path
    ids_path: Path
    metadata_path: Path
    log_path: Path
    blastp_version: str
    settings: dict[str, Any]
    query_count: int
    reference_count: int
    selected_count: int
    proteins_sha256: str
    reference_sha256: str
    hits_sha256: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, metadata: dict[str, Any]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(metadata, stream, indent=2, allow_nan=False)
        stream.write("\n")


def _prepare_outputs(output_dir: Path, names: tuple[str, ...]) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {name: output_dir / name for name in names}
    existing = [path.name for path in paths.values() if path.exists()]
    if existing:
        raise FileExistsError("Output files already exist: " + ", ".join(existing))
    return paths


def _load_pyrodigal() -> Any:
    try:
        return importlib.import_module("pyrodigal")
    except ImportError as exc:
        raise ImportError(
            "Assembly gene calling requires Pyrodigal. Install the optional dependency "
            "with: pip install 'grodon2-py[genome]'"
        ) from exc


def _validate_sequences(records: GeneSequences, alphabet: frozenset[str], kind: str) -> None:
    for identifier, sequence in zip(records.ids, records.sequences, strict=True):
        invalid = set(sequence).difference(alphabet)
        if invalid:
            symbols = ", ".join(repr(symbol) for symbol in sorted(invalid))
            raise ValueError(f"Invalid {kind} symbols in {identifier!r}: {symbols}")


def _contig_ids(identifiers: tuple[str, ...]) -> dict[str, str]:
    reserved = {identifier for identifier in identifiers if _SAFE_ID.fullmatch(identifier)}
    used = set(reserved)
    mapping = {}
    for index, identifier in enumerate(identifiers, 1):
        if identifier in reserved:
            mapped = identifier
        else:
            mapped = f"contig_{index:05d}"
            suffix = 1
            while mapped in used:
                mapped = f"contig_{index:05d}_{suffix}"
                suffix += 1
            used.add(mapped)
        mapping[identifier] = mapped
    return mapping


def call_genes(
    genome_path: str | PathLike[str],
    output_dir: str | PathLike[str],
    *,
    genetic_code: int = 11,
    meta: bool = False,
    closed: bool = True,
    mask: bool = True,
) -> GeneCallingResult:
    """Predict coding DNA and proteins from a plain or gzipped assembly FASTA.

    Defaults reproduce the validation workflow: Pyrodigal learns one model from
    all input contigs, with closed ends, unknown-base masking and code 11. Meta
    mode explicitly uses pretrained bins for the requested code (4 or 11).
    DNA IUPAC symbols are accepted; gaps and protein-only symbols are rejected.

    Safe contig IDs retain ``{contig}_CDS{number:05d}``; other IDs receive stable
    replacements recorded in the gene table and metadata. Output targets must
    not already exist. Coordinates are Pyrodigal's one-based inclusive bounds.
    """
    if isinstance(genetic_code, bool) or not isinstance(genetic_code, Integral):
        raise ValueError("genetic_code must be an integer NCBI translation-table number")
    if not all(isinstance(value, bool) for value in (meta, closed, mask)):
        raise TypeError("meta, closed and mask must be booleans")
    genetic_code = int(genetic_code)
    source = Path(genome_path).resolve(strict=True)
    output = Path(output_dir).resolve()
    records = read_fasta(source)
    _validate_sequences(records, _DNA_ALPHABET, "DNA IUPAC")
    pyrodigal = _load_pyrodigal()
    finder_options: dict[str, Any] = {"meta": meta, "closed": closed, "mask": mask}
    if meta:
        bins = [
            item
            for item in pyrodigal.METAGENOMIC_BINS
            if item.training_info.translation_table == genetic_code
        ]
        if not bins:
            raise ValueError(
                f"Pyrodigal has no metagenomic training bins for code {genetic_code}; use 4 or 11"
            )
        finder_options["metagenomic_bins"] = pyrodigal.MetagenomicBins(bins)
    finder = pyrodigal.GeneFinder(**finder_options)
    if not meta:
        try:
            finder.train(*records.sequences, translation_table=genetic_code)
        except (ValueError, RuntimeError) as exc:
            raise ValueError(
                f"Single-genome gene-caller training failed: {exc}. "
                "For a short assembly, explicitly choose meta=True."
            ) from exc
    names = (
        "predicted_cds.fna",
        "predicted_proteins.faa",
        "predicted_genes.tsv",
        "gene_calling.json",
    )
    paths = _prepare_outputs(output, names)
    mapping = _contig_ids(records.ids)
    counts: dict[str, int] = {}
    with (
        paths[names[0]].open("x", encoding="utf-8", newline="\n") as cds,
        paths[names[1]].open("x", encoding="utf-8", newline="\n") as proteins,
        paths[names[2]].open("x", encoding="utf-8", newline="") as table,
    ):
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["id", "contig", "start", "end", "strand", "length_bp", "partial_begin", "partial_end"]
        )
        for contig, sequence in zip(records.ids, records.sequences, strict=True):
            genes = finder.find_genes(sequence)
            counts[contig] = len(genes)
            for number, gene in enumerate(genes, 1):
                identifier = f"{mapping[contig]}_CDS{number:05d}"
                dna = gene.sequence()
                protein = gene.translate().rstrip("*")
                if not dna or len(dna) % 3:
                    raise RuntimeError(f"Gene caller returned an out-of-frame CDS for {identifier}")
                cds.write(f">{identifier}\n{dna}\n")
                proteins.write(f">{identifier}\n{protein}\n")
                writer.writerow(
                    [
                        identifier,
                        contig,
                        gene.begin,
                        gene.end,
                        gene.strand,
                        len(dna),
                        gene.partial_begin,
                        gene.partial_end,
                    ]
                )
    result: GeneCallingResult = {
        "cds_path": paths[names[0]],
        "proteins_path": paths[names[1]],
        "gene_table_path": paths[names[2]],
        "metadata_path": paths[names[3]],
        "input": str(source),
        "input_sha256": _sha256(source),
        "pyrodigal_version": str(pyrodigal.__version__),
        "settings": {
            "meta": meta,
            "closed": closed,
            "mask": mask,
            "translation_table": genetic_code,
            "training": "pretrained metagenomic bins for requested code"
            if meta
            else "all input contigs, single-genome mode",
        },
        "counts_by_contig": counts,
        "contig_id_map": mapping,
        "total_predicted_genes": sum(counts.values()),
        "total_genome_bp": sum(map(len, records.sequences)),
        "cds_sha256": _sha256(paths[names[0]]),
        "proteins_sha256": _sha256(paths[names[1]]),
        "gene_table_sha256": _sha256(paths[names[2]]),
    }
    persistent = {
        key: value.name if isinstance(value, Path) else value for key, value in result.items()
    }
    _write_json(paths[names[3]], persistent)
    return result


def _resolve_blastp(executable: str | PathLike[str]) -> Path:
    candidate = str(executable)
    located = shutil.which(candidate)
    if located is None:
        raise FileNotFoundError(
            f"BLAST protein search executable not found: {candidate!r}. "
            "Install NCBI BLAST+ and provide blastp or its explicit executable path."
        )
    return Path(located).resolve()


def _relative_input(path: Path, output: Path) -> str:
    try:
        return str(path.relative_to(output))
    except ValueError:
        return str(path)


def _run_blast(command: list[str], output: Path, log_path: Path, timeout: float) -> str:
    try:
        process = subprocess.run(
            command,
            cwd=output,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            check=False,
            shell=False,
        )
    except (subprocess.TimeoutExpired, OSError) as exc:
        with log_path.open("a", encoding="utf-8", newline="\n") as log:
            log.write(json.dumps(command) + "\n" + str(exc) + "\n")
            if isinstance(exc, subprocess.TimeoutExpired):
                for captured in (exc.stdout, exc.stderr):
                    if captured:
                        text = (
                            captured.decode("utf-8", errors="replace")
                            if isinstance(captured, bytes)
                            else captured
                        )
                        log.write(text + "\n")
        if isinstance(exc, subprocess.TimeoutExpired):
            raise RuntimeError(
                f"BLAST timed out after {timeout:g} seconds; see {log_path}"
            ) from exc
        raise RuntimeError(f"Could not execute BLAST; see {log_path}: {exc}") from exc
    with log_path.open("a", encoding="utf-8", newline="\n") as log:
        log.write(json.dumps(command) + f"\nexit_status={process.returncode}\n")
        log.write(process.stdout + "\n" + process.stderr + "\n")
    if process.returncode:
        raise RuntimeError(f"BLAST exited with status {process.returncode}; see {log_path}")
    return process.stdout.strip()


def _normalize_lf(path: Path) -> None:
    """Normalize an external program's text output with bounded memory."""
    temporary: Path | None = None
    try:
        with (
            path.open("r", encoding="utf-8", newline=None) as source,
            tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="\n",
                dir=path.parent,
                prefix=".blast-hits-",
                suffix=".tmp",
                delete=False,
            ) as target,
        ):
            temporary = Path(target.name)
            shutil.copyfileobj(source, target, length=1024 * 1024)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _select_ribosomal_hits(
    path: Path,
    proteins: GeneSequences,
    references: GeneSequences,
    *,
    min_identity: float,
    min_query_coverage: float,
    min_reference_coverage: float,
    max_evalue: float,
) -> list[RibosomalHit]:
    query_lengths = dict(zip(proteins.ids, map(len, proteins.sequences), strict=True))
    reference_lengths = dict(zip(references.ids, map(len, references.sequences), strict=True))
    descriptions = dict(zip(references.ids, references.names, strict=True))
    order = {identifier: index for index, identifier in enumerate(references.ids)}
    selected: dict[str, RibosomalHit] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for number, fields in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if not fields:
                continue
            if len(fields) != 12:
                raise ValueError(f"Malformed BLAST output at line {number}: expected 12 columns")
            query, subject = fields[:2]
            if query not in query_lengths or subject not in reference_lengths:
                raise ValueError(f"Unknown query or reference ID in BLAST output at line {number}")
            try:
                identity, length, qlen, slen, qstart, qend, sstart, send, evalue, bitscore = map(
                    float, fields[2:]
                )
            except ValueError as exc:
                raise ValueError(f"Nonnumeric BLAST output at line {number}") from exc
            numbers = (identity, length, qlen, slen, qstart, qend, sstart, send, evalue, bitscore)
            if not all(math.isfinite(value) for value in numbers):
                raise ValueError(f"Nonfinite BLAST output at line {number}")
            if (
                not 0 <= identity <= 100
                or evalue < 0
                or bitscore < 0
                or length <= 0
                or qlen != query_lengths[query]
                or slen != reference_lengths[subject]
                or not 1 <= qstart <= qend <= qlen
                or not 1 <= sstart <= send <= slen
                or any(
                    not value.is_integer()
                    for value in (length, qlen, slen, qstart, qend, sstart, send)
                )
            ):
                raise ValueError(f"Invalid BLAST alignment values at line {number}")
            qcover, scover = (qend - qstart + 1) / qlen, (send - sstart + 1) / slen
            if (
                identity < min_identity
                or qcover < min_query_coverage
                or scover < min_reference_coverage
                or evalue > max_evalue
            ):
                continue
            row: RibosomalHit = {
                "id": query,
                "reference": subject,
                "description": descriptions[subject],
                "identity_percent": identity,
                "query_coverage": qcover,
                "reference_coverage": scover,
                "evalue": evalue,
                "bitscore": bitscore,
            }
            previous = selected.get(query)
            if (
                previous is None
                or bitscore > previous["bitscore"]
                or (
                    bitscore == previous["bitscore"]
                    and order[subject] < order[previous["reference"]]
                )
            ):
                selected[query] = row
    return [selected[identifier] for identifier in proteins.ids if identifier in selected]


def annotate_ribosomal(
    proteins_path: str | PathLike[str],
    reference_path: str | PathLike[str],
    output_dir: str | PathLike[str],
    *,
    blastp: str | PathLike[str] = "blastp",
    min_identity: float = 70,
    min_query_coverage: float = 0.8,
    min_reference_coverage: float = 0.8,
    max_evalue: float = 1e-8,
    threads: int = 1,
    timeout: float = 600,
) -> RibosomalAnnotationResult:
    """Select ribosomal homologs using explicit local BLASTP and reference FASTA.

    Filters use alignment coordinate spans divided by full query/reference
    lengths, matching the validated selection workflow. The greatest bit score
    wins; exact score ties prefer reference FASTA order. Selected IDs retain
    query FASTA order. No minimum number of selected proteins is imposed here.
    The caller owns reference suitability and the downstream minimum-gene rule.
    """
    for name, value, maximum in (
        ("min_identity", min_identity, 100),
        ("min_query_coverage", min_query_coverage, 1),
        ("min_reference_coverage", min_reference_coverage, 1),
    ):
        if isinstance(value, bool) or not math.isfinite(value) or not 0 <= value <= maximum:
            raise ValueError(f"{name} must be finite and between 0 and {maximum}")
    if isinstance(max_evalue, bool) or not math.isfinite(max_evalue) or max_evalue <= 0:
        raise ValueError("max_evalue must be positive and finite")
    if isinstance(threads, bool) or not isinstance(threads, Integral) or threads <= 0:
        raise ValueError("threads must be a positive integer")
    if isinstance(timeout, bool) or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be positive and finite")
    source, reference = (
        Path(proteins_path).resolve(strict=True),
        Path(reference_path).resolve(strict=True),
    )
    output = Path(output_dir).resolve()
    proteins, references = read_fasta(source), read_fasta(reference)
    _validate_sequences(proteins, _PROTEIN_ALPHABET, "protein")
    _validate_sequences(references, _PROTEIN_ALPHABET, "protein")
    executable = _resolve_blastp(blastp)
    names = (
        "ribosomal_hits.tsv",
        "selected_ribosomal.tsv",
        "highly_expressed_ids.txt",
        "ribosomal_annotation.json",
        "blastp.log",
    )
    paths = _prepare_outputs(output, names)
    version = _run_blast([str(executable), "-version"], output, paths[names[4]], min(timeout, 10))
    # BLAST itself does not read gzip FASTA. Normalize input FASTAs into local
    # scratch files only when compressed; preserve original paths in provenance.
    query_argument, reference_argument = (
        _relative_input(source, output),
        _relative_input(reference, output),
    )
    scratch: list[Path] = []
    try:
        for path, records, label in (
            (source, proteins, "query"),
            (reference, references, "reference"),
        ):
            if path.suffix.lower() == ".gz":
                plain = output / f".blast_{label}.faa"
                with plain.open("x", encoding="utf-8", newline="\n") as stream:
                    for name, sequence in zip(records.names, records.sequences, strict=True):
                        stream.write(f">{name}\n{sequence}\n")
                scratch.append(plain)
                if label == "query":
                    query_argument = plain.name
                else:
                    reference_argument = plain.name
        command = [
            str(executable),
            "-query",
            query_argument,
            "-subject",
            reference_argument,
            "-out",
            paths[names[0]].name,
            "-outfmt",
            _OUTFMT,
            "-evalue",
            format(max_evalue, ".15g"),
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-max_hsps",
            "1",
            "-num_threads",
            str(int(threads)),
            "-max_target_seqs",
            str(len(references)),
        ]
        _run_blast(command, output, paths[names[4]], timeout)
    except Exception:
        paths[names[0]].unlink(missing_ok=True)
        raise
    finally:
        for path in scratch:
            path.unlink(missing_ok=True)
    if not paths[names[0]].is_file():
        raise RuntimeError(f"BLAST completed without creating its hit table; see {paths[names[4]]}")
    _normalize_lf(paths[names[0]])
    selected = _select_ribosomal_hits(
        paths[names[0]],
        proteins,
        references,
        min_identity=min_identity,
        min_query_coverage=min_query_coverage,
        min_reference_coverage=min_reference_coverage,
        max_evalue=max_evalue,
    )
    selected_ids = [hit["id"] for hit in selected]
    with paths[names[1]].open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=_HIT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(selected)
    with paths[names[2]].open("x", encoding="utf-8", newline="\n") as stream:
        stream.write("".join(identifier + "\n" for identifier in selected_ids))
    result: RibosomalAnnotationResult = {
        "selected_ids": selected_ids,
        "selected_hits": selected,
        "hits_path": paths[names[0]],
        "selected_hits_path": paths[names[1]],
        "ids_path": paths[names[2]],
        "metadata_path": paths[names[3]],
        "log_path": paths[names[4]],
        "blastp_version": version,
        "settings": {
            "min_identity": min_identity,
            "min_query_coverage": min_query_coverage,
            "min_reference_coverage": min_reference_coverage,
            "max_evalue": max_evalue,
            "threads": int(threads),
            "timeout": timeout,
            "seg": "no",
            "comp_based_stats": 0,
            "max_hsps": 1,
            "max_target_seqs": len(references),
            "outfmt": _OUTFMT,
            "tie_break": "reference FASTA order",
        },
        "query_count": len(proteins),
        "reference_count": len(references),
        "selected_count": len(selected),
        "proteins_sha256": _sha256(source),
        "reference_sha256": _sha256(reference),
        "hits_sha256": _sha256(paths[names[0]]),
    }
    persistent = {
        key: value.name if isinstance(value, Path) else value for key, value in result.items()
    }
    persistent.update(
        proteins_input=_relative_input(source, output),
        reference_input=_relative_input(reference, output),
        blastp_executable=str(executable),
        command=command,
    )
    _write_json(paths[names[3]], persistent)
    return result
