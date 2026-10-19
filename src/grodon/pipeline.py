"""Reproducible, offline bacterial assembly-to-growth workflows.

Each run owns a new output directory. Failed runs retain their intermediates
and a failure manifest; no output from an earlier run is overwritten.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import platform
import re
import warnings
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from importlib import metadata
from numbers import Integral, Real
from os import PathLike
from pathlib import Path
from typing import Any

from ._version import __version__
from .annotation import annotate_ribosomal, call_genes
from .io import highly_expressed_mask, read_fasta
from .models import model_provenance
from .prediction import PredictionResult, predict_growth
from .references import DEFAULT_REFERENCE, resolve_reference, sha256_file
from .statistics import filter_sequences
from .tools import doctor


@dataclass(frozen=True)
class GenomePipelineConfig:
    """Explicit settings shared by gene calling, homology, and prediction.

    ``min_ribosomal`` is an operational quality gate on genes remaining after
    the gRodon sequence filter, not a calibrated confidence guarantee.
    """

    reference: str | PathLike[str] = DEFAULT_REFERENCE
    blastp: str | PathLike[str] = "blastp"
    mode: str = "full"
    training_set: str = "madin"
    temperature: float | None = None
    genetic_code: int = 11
    gene_mode: str = "single"
    closed: bool = True
    mask: bool = True
    min_identity: float = 70.0
    min_query_coverage: float = 0.8
    min_reference_coverage: float = 0.8
    max_evalue: float = 1e-8
    min_ribosomal: int = 10
    confidence: float = 0.95
    threads: int = 1
    timeout: float = 600.0

    def __post_init__(self) -> None:
        if self.mode not in {"full", "partial"}:
            raise ValueError("Genome pipeline mode must be 'full' or 'partial'")
        if self.training_set not in {"madin", "vs", "AOA_NOB"}:
            raise ValueError("training_set must be madin, vs, or AOA_NOB")
        if self.gene_mode not in {"single", "meta"}:
            raise ValueError("gene_mode must be 'single' or 'meta'")
        # These are the translation tables supported by Prodigal gene calling.
        if (
            isinstance(self.genetic_code, bool)
            or not isinstance(self.genetic_code, Integral)
            or self.genetic_code not in {4, 11}
        ):
            raise ValueError("Genome gene calling supports genetic_code 4 or 11")
        for name in ("closed", "mask"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be a boolean")
        for name in ("min_ribosomal", "threads"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        for name, lower, upper, exclusive in (
            ("min_identity", 0, 100, False),
            ("min_query_coverage", 0, 1, False),
            ("min_reference_coverage", 0, 1, False),
            ("confidence", 0, 1, True),
            ("max_evalue", 0, math.inf, True),
            ("timeout", 0, math.inf, True),
        ):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not math.isfinite(value)
                or (
                    not lower < float(value) < upper
                    if exclusive
                    else not lower <= float(value) <= upper
                )
            ):
                raise ValueError(f"Invalid {name}: {value!r}")
        if self.temperature is not None and (
            isinstance(self.temperature, bool)
            or not isinstance(self.temperature, Real)
            or not math.isfinite(self.temperature)
        ):
            raise ValueError("temperature must be finite or None")


@dataclass(frozen=True)
class GenomePipelineResult:
    """A completed prediction and the location of its auditable run files."""

    prediction: PredictionResult
    output_dir: Path
    summary: Mapping[str, Any]

    def to_dict(self) -> dict[str, Any]:
        """Return a detached, JSON-compatible summary."""
        return _json_value(dict(self.summary))


def _json_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, PathLike):
        return str(value)
    if hasattr(value, "item"):
        return _json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _write_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(_json_value(value), stream, indent=2, allow_nan=False)
        stream.write("\n")
    temporary.replace(path)


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _software() -> dict[str, str]:
    versions = {"python": platform.python_version(), "grodon2-py": __version__}
    for package in ("numpy", "scipy", "pandas", "biopython", "pyrodigal"):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "unavailable"
    return versions


def _new_directory(path: str | PathLike[str]) -> Path:
    output = Path(path).expanduser().absolute()
    # mkdir(exist_ok=False) rejects existing directories, files, and symlinks.
    # In particular, an existing empty directory is not an invitation to reuse.
    try:
        output.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise FileExistsError(f"Output already exists; choose a new directory: {output}") from exc
    return output.resolve()


def _artifacts(output: Path) -> dict[str, Any]:
    return {
        path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
        for path in sorted(output.iterdir())
        if path.is_file() and path.name not in {"run.json", "run.json.tmp"}
    }


def _csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(_json_value(rows))


_SUMMARY_FIELDS = [
    "input",
    "output_dir",
    "status",
    "doubling_time_hours",
    "doubling_time_minutes",
    "growth_rate_per_hour",
    "lower_ci_hours",
    "upper_ci_hours",
    "predicted_genes",
    "ribosomal_genes",
    "retained_genes",
    "retained_ribosomal_genes",
    "reference",
    "model",
    "r_matches",
    "r_comparison_status",
    "error",
]


def _summary_row(result: Mapping[str, Any]) -> dict[str, Any]:
    return {name: result.get(name) for name in _SUMMARY_FIELDS}


def run_genome(
    genome: str | PathLike[str],
    output_dir: str | PathLike[str],
    *,
    config: GenomePipelineConfig | None = None,
    rscript: str | PathLike[str] | None = None,
) -> GenomePipelineResult:
    """Analyze one bacterial genome FASTA (plain/gzip) entirely locally.

    Optional ``rscript`` reruns the original R package on exactly the same CDS
    and selected IDs. The Python prediction itself never requires R.
    """
    settings = GenomePipelineConfig() if config is None else config
    if not isinstance(settings, GenomePipelineConfig):
        raise TypeError("config must be a GenomePipelineConfig")
    if rscript is not None and settings.confidence != 0.95:
        raise ValueError("The original R comparison supports only confidence=0.95")
    source = Path(genome).expanduser().resolve(strict=True)
    if not source.is_file():
        raise ValueError(f"Genome must be a FASTA file: {source}")
    health = doctor(settings.blastp)
    if not health["ready"]:
        problems = [
            str(health[name].get("error") or f"{name} unavailable")
            for name in ("pyrodigal", "blastp")
            if not health[name]["available"]
        ]
        raise RuntimeError("Genome pipeline dependencies unavailable: " + "; ".join(problems))
    # Resolve/check the reference before reserving the output directory.
    with resolve_reference(settings.reference) as reference:
        output = _new_directory(output_dir)
        manifest: dict[str, Any] = {
            "schema_version": 1,
            "status": "running",
            "started_at": _now(),
            "input": {"path": str(source), "sha256": sha256_file(source)},
            "settings": asdict(settings),
            "software": _software(),
            "platform": platform.platform(),
            "dependencies": health,
            "reference": {
                "name": reference.name,
                "version": reference.version,
                "sha256": reference.sha256,
                "proteins": reference.count,
            },
            "model_provenance": model_provenance(),
        }
        _write_json(output / "run.json", manifest)
        try:
            _write_json(output / "reference.json", reference.source_metadata)
            # Freeze the exact reference bytes used, including custom references.
            with (output / "reference.faa.gz").open("wb") as stream:
                with gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0) as zipped:
                    zipped.write(reference.path.read_bytes())
            with warnings.catch_warnings(record=True) as recorded:
                warnings.simplefilter("always")
                called = call_genes(
                    source,
                    output,
                    genetic_code=settings.genetic_code,
                    meta=settings.gene_mode == "meta",
                    closed=settings.closed,
                    mask=settings.mask,
                )
                if called["input_sha256"] != manifest["input"]["sha256"]:
                    raise RuntimeError("Genome changed during processing; rerun from a stable file")
                annotated = annotate_ribosomal(
                    called["proteins_path"],
                    output / "reference.faa.gz",
                    output,
                    blastp=health["blastp"]["path"],
                    min_identity=settings.min_identity,
                    min_query_coverage=settings.min_query_coverage,
                    min_reference_coverage=settings.min_reference_coverage,
                    max_evalue=settings.max_evalue,
                    threads=settings.threads,
                    timeout=settings.timeout,
                )
                genes = read_fasta(called["cds_path"])
                mask = highly_expressed_mask(genes, ids=annotated["selected_ids"])
                if not any(mask):
                    raise ValueError("No ribosomal proteins pass the homology thresholds")
                filtered = filter_sequences(genes, mask)
                n_he = int(filtered.highly_expressed.sum())
                if n_he < settings.min_ribosomal:
                    raise ValueError(
                        f"Only {n_he} ribosomal genes remain after filtering; "
                        f"at least {settings.min_ribosomal} are required. "
                        "Review the assembly and reference coverage."
                    )
                prediction = predict_growth(
                    genes,
                    mask,
                    mode=settings.mode,
                    training_set=settings.training_set,
                    temperature=settings.temperature,
                    genetic_code=None if settings.genetic_code == 11 else settings.genetic_code,
                    fragments=False,
                    confidence=settings.confidence,
                )
            messages = list(dict.fromkeys(str(item.message) for item in recorded))
            d = prediction.d
            if not math.isfinite(d) or d <= 0:
                raise ValueError("The fitted model did not produce a finite positive doubling time")
            summary: dict[str, Any] = {
                "schema_version": 1,
                "status": "complete",
                "input": str(source),
                "output_dir": str(output),
                "input_sha256": manifest["input"]["sha256"],
                "doubling_time_hours": d,
                "doubling_time_minutes": d * 60,
                "growth_rate_per_hour": math.log(2) / d,
                "lower_ci_hours": prediction.lower_ci,
                "upper_ci_hours": prediction.upper_ci,
                "confidence": settings.confidence,
                "predicted_genes": len(genes),
                "ribosomal_genes": len(annotated["selected_ids"]),
                "retained_genes": len(filtered.genes),
                "retained_ribosomal_genes": n_he,
                "genome_bp": called["total_genome_bp"],
                "contigs": len(called["counts_by_contig"]),
                "reference": reference.name,
                "reference_sha256": reference.sha256,
                "model": prediction.model_name,
                "statistics": dict(prediction),
                "settings": asdict(settings),
                "warnings": messages,
            }
            # A fingerprint describes scientific inputs/outputs, independent of output
            # paths, timestamps, logs, or the order of equally scoring BLAST hits.
            fingerprint = {
                "input_sha256": summary["input_sha256"],
                "reference_sha256": reference.sha256,
                "cds_sha256": sha256_file(called["cds_path"]),
                "selected_ids": annotated["selected_ids"],
                "statistics": dict(prediction),
                "model": prediction.model_name,
            }
            summary["scientific_fingerprint"] = hashlib.sha256(
                json.dumps(_json_value(fingerprint), sort_keys=True, allow_nan=False).encode()
            ).hexdigest()
            if rscript is not None:
                from .r_comparison import compare_with_r

                summary.update(
                    status="prediction_complete_r_pending", r_comparison_status="running"
                )
                _write_json(output / "result.json", summary)
                _csv(output / "result.csv", [_summary_row(summary)], _SUMMARY_FIELDS)
                try:
                    comparison = compare_with_r(
                        called["cds_path"],
                        annotated["ids_path"],
                        output,
                        rscript=rscript,
                        mode=settings.mode,
                        training_set=settings.training_set,
                        temperature=settings.temperature,
                        genetic_code=settings.genetic_code,
                        python_statistics=dict(prediction),
                        timeout=settings.timeout,
                    )
                except Exception as exc:
                    summary.update(
                        status="prediction_complete_r_failed",
                        r_comparison_status="failed",
                        r_comparison_error=str(exc),
                    )
                    _write_json(output / "result.json", summary)
                    _csv(output / "result.csv", [_summary_row(summary)], _SUMMARY_FIELDS)
                    raise
                summary.update(status="complete", r_comparison_status="complete")
                summary["r_matches"] = comparison["matches"]
                if not comparison["matches"]:
                    messages.append("Python and R results differ; inspect r_comparison.json")
            else:
                summary["r_comparison_status"] = "not_requested"
            _write_json(output / "result.json", summary)
            _csv(output / "result.csv", [_summary_row(summary)], _SUMMARY_FIELDS)
            manifest.update(
                status="complete",
                finished_at=_now(),
                scientific_fingerprint=summary["scientific_fingerprint"],
                warnings=messages,
                artifacts=_artifacts(output),
            )
            _write_json(output / "run.json", manifest)
            return GenomePipelineResult(prediction, output, _json_value(summary))
        except Exception as exc:
            manifest.update(
                status="failed",
                finished_at=_now(),
                error={"type": type(exc).__name__, "message": str(exc)},
                artifacts=_artifacts(output),
            )
            _write_json(output / "run.json", manifest)
            raise


def run_genomes(
    genomes: Iterable[str | PathLike[str]],
    output_dir: str | PathLike[str],
    *,
    config: GenomePipelineConfig | None = None,
    rscript: str | PathLike[str] | None = None,
    continue_on_error: bool = False,
) -> dict[str, Any]:
    """Run a sequential batch, preserving a summary even after a failed input.

    Each FASTA represents one genome (possibly with multiple contigs/plasmids).
    On default fail-fast behavior the summary is saved before raising. With
    ``continue_on_error=True`` failures remain explicit in the returned summary.
    """
    if isinstance(genomes, (str, PathLike)):
        raise TypeError("genomes must be an iterable of genome paths, not one path")
    sources = [Path(path).expanduser().absolute() for path in genomes]
    if not sources:
        raise ValueError("Provide at least one genome")
    settings = GenomePipelineConfig() if config is None else config
    if not isinstance(settings, GenomePipelineConfig):
        raise TypeError("config must be a GenomePipelineConfig")
    if not isinstance(continue_on_error, bool):
        raise TypeError("continue_on_error must be a boolean")
    if rscript is not None and settings.confidence != 0.95:
        raise ValueError("The original R comparison supports only confidence=0.95")
    output = _new_directory(output_dir)
    result: dict[str, Any] = {
        "schema_version": 1,
        "status": "running",
        "started_at": _now(),
        "output_dir": str(output),
        "settings": asdict(settings),
        "requested": len(sources),
        "completed": 0,
        "failed": 0,
        "results": [],
    }

    def save() -> None:
        _write_json(output / "summary.json", result)
        _csv(
            output / "summary.csv",
            [_summary_row(row) for row in result["results"]],
            _SUMMARY_FIELDS,
        )

    save()
    for index, source in enumerate(sources, 1):
        stem = source.name.removesuffix(".gz")
        stem = re.sub(r"[^A-Za-z0-9_.-]", "_", Path(stem).stem)[:80] or "genome"
        destination = output / "genomes" / f"{index:04d}_{stem}"
        try:
            row = run_genome(source, destination, config=settings, rscript=rscript).to_dict()
            result["completed"] += 1
        except Exception as exc:
            row = {
                "input": str(source),
                "output_dir": str(destination),
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
            }
            result["failed"] += 1
            result["results"].append(row)
            if not continue_on_error:
                result.update(status="failed", finished_at=_now())
                save()
                raise RuntimeError(
                    f"Batch stopped at {source.name}: {exc}; see {output / 'summary.json'}"
                ) from exc
            save()
            continue
        result["results"].append(row)
        save()
    result.update(
        status="complete" if result["failed"] == 0 else "completed_with_errors", finished_at=_now()
    )
    save()
    return _json_value(result)
