"""Audit existing paired Python/R genome runs and render a complete folder report.

This script never predicts genes or growth rates. Input paths in the inventory
are absolute, or relative to the inventory file. A passing row requires finite
values for all twelve full-mode statistics, matching input identities/settings/
counts, and (by default) verified SHA256 hashes for every recorded artifact.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import struct
from collections import Counter, defaultdict
from collections.abc import Mapping
from datetime import UTC, datetime
from numbers import Real
from pathlib import Path
from typing import Any

RTOL = 1e-10
ATOL = 1e-12
FULL_METRICS = (
    "CUBHE",
    "GC",
    "GCdiv",
    "ConsistencyHE",
    "CUB",
    "CPB",
    "FilteredSequences",
    "nHE",
    "dCUB",
    "d",
    "LowerCI",
    "UpperCI",
)
_STATUSES = ("pass", "mismatch", "error", "not_run")


def _sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _finite(value: Any) -> bool:
    if not isinstance(value, Real) or isinstance(value, bool):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


def _json_safe(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _read_object(path: Path, errors: list[str]) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
        if not isinstance(value, dict):
            raise ValueError("expected a JSON object")
        return value
    except (OSError, ValueError) as exc:
        errors.append(f"Cannot read {path.name}: {exc}")
        return None


def _resolve_path(value: Any, base: Path) -> Path | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        path = Path(value)
        return (path if path.is_absolute() else base / path).resolve()
    except (OSError, ValueError):
        return None


def _object(value: Any, label: str, errors: list[str]) -> dict[str, Any]:
    if not isinstance(value, dict):
        errors.append(f"{label} must be a JSON object")
        return {}
    return value


def _warnings(value: Any, label: str, errors: list[str]) -> list[str]:
    # jsonlite auto_unbox emits a string when R records exactly one warning.
    if isinstance(value, str):
        return [value]
    if isinstance(value, list) and all(isinstance(item, str) for item in value):
        return value
    errors.append(f"{label} warnings must be a string or a list of strings")
    return []


def _check(
    checks: list[dict[str, Any]], errors: list[str], name: str, expected: Any, actual: Any
) -> bool:
    matches = expected is not None and actual is not None and expected == actual
    checks.append({"check": name, "expected": expected, "observed": actual, "matches": matches})
    if not matches:
        errors.append(f"{name}: expected {expected!r}, observed {actual!r}")
    return matches


def _input_identity(entry: dict[str, Any], base: Path, row: dict[str, Any]) -> None:
    checks, errors = row["input_checks"], row["errors"]
    for prefix, path_key in (("source", "source_file"), ("genome", "genome_path")):
        path = _resolve_path(entry.get(path_key), base)
        row["inputs"][path_key] = str(path) if path is not None else None
        expected = entry.get(f"{prefix}_sha256")
        row["inputs"][f"{prefix}_sha256"] = expected
        observed = None
        if path is not None:
            try:
                observed = _sha256(path)
            except OSError as exc:
                errors.append(f"Cannot read {prefix} input {path}: {exc}")
        _check(checks, errors, f"{prefix} input SHA256", expected, observed)
    assembly = entry.get("assembly", {})
    if isinstance(assembly, dict) and "accession" in assembly:
        _check(
            checks, errors, "Assembly accession identity", entry["accession"], assembly["accession"]
        )
    download = entry.get("download", {})
    if isinstance(download, dict):
        for key in (
            "genome_md5_matches",
            "source_protein_content_matches",
            "source_protein_order_matches",
        ):
            if download.get(key) is False:
                errors.append(f"Input preparation verification failed: {key}")
        for prefix in ("dataset_report", "checksums", "official_protein"):
            if f"{prefix}_path" not in download and f"{prefix}_sha256" not in download:
                continue
            path = _resolve_path(download.get(f"{prefix}_path"), base)
            observed = None
            if path is not None:
                try:
                    observed = _sha256(path)
                except OSError as exc:
                    errors.append(f"Cannot read input provenance artifact {path}: {exc}")
            _check(checks, errors, f"{prefix} SHA256", download.get(f"{prefix}_sha256"), observed)


def _verify_artifacts(run_dir: Path, manifest: dict[str, Any], *, enabled: bool) -> dict[str, Any]:
    result: dict[str, Any] = {"enabled": enabled, "verified": False, "files": [], "errors": []}
    if not enabled:
        return result
    errors = result["errors"]
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        errors.append("run.json has no artifact SHA256 manifest")
        return result
    actual_names = {item.name for item in run_dir.iterdir() if item.name != "run.json"}
    expected_names = set(artifacts)
    for name in sorted(actual_names - expected_names):
        errors.append(f"Unrecorded run artifact: {name}")
    for name, metadata in artifacts.items():
        path = run_dir / name
        safe = (
            isinstance(name, str)
            and Path(name).name == name
            and name not in {".", "..", "run.json"}
            and "/" not in name
            and "\\" not in name
            and ":" not in name
            and path.resolve().parent == run_dir.resolve()
            and not path.is_symlink()
        )
        record = {"name": name, "sha256": None, "bytes": None, "matches": False}
        result["files"].append(record)
        if not safe or not isinstance(metadata, dict):
            errors.append(f"Unsafe or malformed artifact entry: {name}")
            continue
        try:
            record["sha256"] = _sha256(path)
            record["bytes"] = path.stat().st_size
        except OSError as exc:
            errors.append(f"Cannot read artifact {name}: {exc}")
            continue
        record["matches"] = record["sha256"] == metadata.get("sha256") and record[
            "bytes"
        ] == metadata.get("bytes")
        if not record["matches"]:
            errors.append(f"Artifact SHA256 or size mismatch: {name}")
    result["verified"] = not errors
    return result


def _metric_comparison(python: dict[str, Any], original: dict[str, Any]) -> list[dict[str, Any]]:
    names = list(FULL_METRICS) + sorted((set(python) | set(original)) - set(FULL_METRICS))
    rows = []
    for name in names:
        left, right = python.get(name), original.get(name)
        row = {
            "metric": name,
            "python": _json_safe(left),
            "r": _json_safe(right),
            "present_in_python": name in python,
            "present_in_r": name in original,
            "finite": False,
            "absolute_difference": None,
            "relative_difference": None,
            "allowed_absolute_difference": None,
            "matches": False,
            "bitwise_equal": False,
            "reason": None,
        }
        if name not in python or name not in original:
            row["reason"] = "Missing statistic"
        elif not _finite(left) or not _finite(right):
            row["reason"] = "Both values must be finite numeric scalars"
        else:
            assert left is not None and right is not None
            left, right = float(left), float(right)
            absolute = abs(left - right)
            scale = max(abs(left), abs(right))
            relative = absolute / scale if scale else 0.0
            allowed = max(ATOL, RTOL * scale)
            row.update(
                finite=True,
                absolute_difference=absolute if math.isfinite(absolute) else None,
                relative_difference=relative if math.isfinite(relative) else None,
                allowed_absolute_difference=allowed,
                matches=absolute <= allowed,
                bitwise_equal=struct.pack("!d", left) == struct.pack("!d", right),
            )
            if not row["matches"]:
                row["reason"] = "Outside numerical tolerance"
        rows.append(row)
    return rows


def _settings(settings: Any) -> dict[str, Any]:
    if not isinstance(settings, dict):
        return {}
    result = {
        key: settings.get(key)
        for key in ("mode", "training_set", "temperature", "genetic_code", "confidence")
    }
    if result["temperature"] == "none":
        result["temperature"] = None
    if result["training_set"] == "original":
        result["training_set"] = "vs"
    # Pipeline input consists of complete predicted CDS; fragments is not a
    # configurable pipeline setting, but is explicit at the R boundary.
    result["fragments"] = settings.get("fragments", False)
    return result


def _validate_counts_and_settings(
    run_dir: Path,
    row: dict[str, Any],
    run: dict[str, Any],
    python: dict[str, Any],
    original: dict[str, Any],
    parameters: dict[str, Any],
) -> None:
    errors, mismatches = row["errors"], row["mismatches"]
    checks: list[dict[str, Any]] = row["agreement_checks"]
    _check(
        checks,
        mismatches,
        "Python/R prediction settings",
        _settings(python.get("settings")),
        _settings(original.get("settings")),
    )
    _check(checks, errors, "Run/result settings", run.get("settings"), python.get("settings"))
    _check(
        checks,
        errors,
        "R requested/reported settings",
        _settings(parameters),
        _settings(original.get("settings")),
    )
    if _settings(python.get("settings")).get("mode") != "full":
        mismatches.append("This report requires full-mode results with all twelve statistics")
    for key, name in (("cds_path", "predicted_cds.fna"), ("ids_path", "highly_expressed_ids.txt")):
        observed = _resolve_path(parameters.get(key), run_dir)
        _check(
            checks,
            errors,
            f"R {key} identity",
            str((run_dir / name).resolve()),
            str(observed) if observed is not None else None,
        )
    cds_count: int | None
    id_count: int | None
    try:
        with (run_dir / "predicted_cds.fna").open(encoding="ascii") as stream:
            cds_count = sum(line.startswith(">") for line in stream)
        ids = [
            line.strip()
            for line in (run_dir / "highly_expressed_ids.txt")
            .read_text(encoding="utf-8")
            .splitlines()
            if line.strip()
        ]
        id_count = len(ids)
        if len(set(ids)) != id_count:
            errors.append("The prepared highly expressed ID file contains duplicates")
    except (OSError, UnicodeError) as exc:
        errors.append(f"Cannot count prepared inputs: {exc}")
        cds_count = id_count = None
    pstats, rstats = python["statistics"], original["statistics"]
    row["counts"] = {
        "prepared_cds": cds_count,
        "prepared_ribosomal_ids": id_count,
        "python_predicted_genes": python.get("predicted_genes"),
        "r_input_genes": original.get("input_genes"),
        "python_ribosomal_genes": python.get("ribosomal_genes"),
        "r_highly_expressed_before_filtering": original.get("highly_expressed_before_filtering"),
        "python_retained_genes": python.get("retained_genes"),
        "python_retained_ribosomal_genes": python.get("retained_ribosomal_genes"),
        "python_nHE": pstats.get("nHE"),
        "r_nHE": rstats.get("nHE"),
    }
    for label, expected, actual in (
        ("Prepared/Python CDS count", cds_count, python.get("predicted_genes")),
        ("Prepared/R CDS count", cds_count, original.get("input_genes")),
        ("Prepared/Python ribosomal count", id_count, python.get("ribosomal_genes")),
        ("Prepared/R ribosomal count", id_count, original.get("highly_expressed_before_filtering")),
        ("Python retained HE count", pstats.get("nHE"), python.get("retained_ribosomal_genes")),
        ("Python/R retained HE count", pstats.get("nHE"), rstats.get("nHE")),
    ):
        _check(checks, mismatches, label, expected, actual)
    filtered = pstats.get("FilteredSequences")
    if cds_count is not None and _finite(filtered):
        _check(
            checks,
            mismatches,
            "Python retained CDS count",
            cds_count - filtered,
            python.get("retained_genes"),
        )


def _summarize_entry(
    entry: dict[str, Any],
    inventory_dir: Path,
    runs_dir: Path,
    verify_artifacts: bool,
) -> dict[str, Any]:
    accession = entry["accession"]
    run_dir = runs_dir / accession
    organism = entry.get("organism", {})
    row: dict[str, Any] = {
        "accession": accession,
        "organism": organism.get("name", "") if isinstance(organism, dict) else organism,
        "inventory_entry": entry,
        "source_kind": entry.get("source_kind"),
        "status": "not_run",
        "run_directory": str(run_dir),
        "errors": [],
        "mismatches": [],
        "inputs": {},
        "input_checks": [],
        "agreement_checks": [],
        "counts": {},
        "statistics": {"python": {}, "r": {}},
        "metric_comparisons": [],
        "numeric_key_sets_equal": False,
        "full_metric_set": False,
        "numerically_matches": False,
        "bitwise_equal": False,
        "max_absolute_difference": None,
        "max_relative_difference": None,
        "warnings": {"python": [], "r": []},
        "software": {},
        "r_versions": {},
        "settings": {},
        "r_settings": {},
        "dependencies": {},
        "doubling_time": {},
        "artifact_verification": {
            "enabled": verify_artifacts,
            "verified": False,
            "files": [],
            "errors": [],
        },
    }
    _input_identity(entry, inventory_dir, row)
    errors, mismatches = row["errors"], row["mismatches"]
    if entry.get("status") != "ready":
        errors.append(f"Input preparation status {entry.get('status')!r}: {entry.get('error')}")
    if not run_dir.exists():
        row["status"] = "error" if errors else "not_run"
        return row
    if not run_dir.is_dir() or run_dir.is_symlink() or run_dir.resolve() != run_dir:
        errors.append("Run directory is not a directory inside the requested runs folder")
        row["status"] = "error"
        return row
    run = _read_object(run_dir / "run.json", errors)
    if run is None:
        row["status"] = "error"
        return row
    row.update(
        run_status=run.get("status"),
        software=_object(run.get("software", {}), "Python software versions", errors),
        dependencies=_object(run.get("dependencies", {}), "Runtime dependencies", errors),
        settings=run.get("settings", {}),
        reference=run.get("reference", {}),
        model_provenance=run.get("model_provenance", {}),
        started_at=run.get("started_at"),
        finished_at=run.get("finished_at"),
    )
    if run.get("status") != "complete":
        errors.append(f"Run status {run.get('status')!r}: {run.get('error', 'not completed')}")
    row["artifact_verification"] = _verify_artifacts(run_dir, run, enabled=verify_artifacts)
    errors.extend(row["artifact_verification"]["errors"])
    run_input = run.get("input", {})
    if not isinstance(run_input, dict):
        run_input = {}
    _check(
        row["input_checks"],
        errors,
        "Run genome SHA256",
        entry.get("genome_sha256"),
        run_input.get("sha256"),
    )
    actual_path = _resolve_path(run_input.get("path"), inventory_dir)
    _check(
        row["input_checks"],
        errors,
        "Run genome path",
        row["inputs"].get("genome_path"),
        str(actual_path) if actual_path is not None else None,
    )
    python = _read_object(run_dir / "result.json", errors)
    original = _read_object(run_dir / "r_result.json", errors)
    parameters = _read_object(run_dir / "r_parameters.json", errors)
    if python is not None:
        row.update(
            prediction_status=python.get("status"),
            r_comparison_status=python.get("r_comparison_status"),
            scientific_fingerprint=python.get("scientific_fingerprint"),
            model=python.get("model"),
        )
        row["warnings"]["python"] = _warnings(python.get("warnings", []), "Python", errors)
        if python.get("status") != "complete" or python.get("r_comparison_status") != "complete":
            errors.append(
                f"R/Python pipeline incomplete: {python.get('r_comparison_error', python.get('status'))}"
            )
        _check(
            row["input_checks"],
            errors,
            "Result genome SHA256",
            entry.get("genome_sha256"),
            python.get("input_sha256"),
        )
    if original is not None:
        row["r_versions"] = _object(original.get("versions", {}), "R software versions", errors)
        row["r_settings"] = original.get("settings", {})
        row["warnings"]["r"] = _warnings(original.get("warnings", []), "R", errors)
    if python is None or original is None or parameters is None:
        row["status"] = "error"
        return row
    for label, value in (("Python", python), ("R", original)):
        if not isinstance(value.get("statistics"), dict):
            errors.append(f"{label} statistics are missing or malformed")
    if not isinstance(python.get("statistics"), dict) or not isinstance(
        original.get("statistics"), dict
    ):
        row["status"] = "error"
        return row
    pstats, rstats = python["statistics"], original["statistics"]
    row["statistics"] = {"python": _json_safe(pstats), "r": _json_safe(rstats)}
    row["metric_comparisons"] = _metric_comparison(pstats, rstats)
    row["numeric_key_sets_equal"] = set(pstats) == set(rstats)
    row["full_metric_set"] = set(pstats) == set(rstats) == set(FULL_METRICS)
    if not row["full_metric_set"]:
        mismatches.append("Python and R must each contain exactly the twelve full-mode statistics")
    for metric in row["metric_comparisons"]:
        if not metric["matches"]:
            mismatches.append(f"{metric['metric']}: {metric['reason']}")
    row["numerically_matches"] = row["full_metric_set"] and all(
        metric["matches"] for metric in row["metric_comparisons"]
    )
    row["bitwise_equal"] = row["full_metric_set"] and all(
        metric["bitwise_equal"] for metric in row["metric_comparisons"]
    )
    for key in ("absolute_difference", "relative_difference"):
        values = [metric[key] for metric in row["metric_comparisons"] if metric[key] is not None]
        row[f"max_{key}"] = max(values, default=None)
    for label, stats in (("python", pstats), ("r", rstats)):
        d, lower, upper = (stats.get(name) for name in ("d", "LowerCI", "UpperCI"))
        row["doubling_time"].update(
            {
                f"{label}_hours": d,
                f"{label}_minutes": d * 60 if _finite(d) else None,
                f"{label}_lower_ci_hours": lower,
                f"{label}_upper_ci_hours": upper,
            }
        )
        if not all(_finite(value) for value in (d, lower, upper)) or not 0 < lower <= d <= upper:
            mismatches.append(
                f"{label}: doubling time and confidence bounds must be finite, positive, and ordered"
            )
    _validate_counts_and_settings(run_dir, row, run, python, original, parameters)
    if _finite(pstats.get("d")) and _finite(rstats.get("d")):
        difference = abs(float(pstats["d"]) - float(rstats["d"]))
        row["doubling_time"]["absolute_difference_hours"] = difference
        row["doubling_time"]["absolute_difference_minutes"] = difference * 60
    row["status"] = "error" if errors else "mismatch" if mismatches else "pass"
    return _json_safe(row)


def summarize(
    manifest_path: str | Path,
    runs_dir: str | Path,
    *,
    verify_artifacts: bool = True,
    expected_count: int = 40,
) -> dict[str, Any]:
    """Read every inventory entry, verify existing results, and classify each run.

    Incomplete inventories cannot produce ``all_pass=True``. ``expected_count``
    permits reuse for other explicit test suites. Floating point equality is
    symmetric: abs(x-y) <= max(1e-12, 1e-10 * max(abs(x), abs(y))). Exact equality
    separately compares binary64 encodings, including the sign of zero.
    """
    if (
        isinstance(expected_count, bool)
        or not isinstance(expected_count, int)
        or expected_count < 1
    ):
        raise ValueError("expected_count must be a positive integer")
    manifest_path, runs_dir = Path(manifest_path).resolve(), Path(runs_dir).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    entries = manifest.get("entries") if isinstance(manifest, dict) else manifest
    if not isinstance(entries, list):
        raise ValueError("Input manifest must be a list or an object containing an entries list")
    seen = set()
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("accession"), str):
            raise ValueError("Every input entry must have a string accession")
        accession = entry["accession"]
        if (
            not accession
            or any(
                character not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-"
                for character in accession
            )
            or accession in {".", ".."}
        ):
            raise ValueError(f"Unsafe input accession: {accession!r}")
        if accession in seen:
            raise ValueError(f"Duplicate input accession: {accession}")
        seen.add(accession)
    rows = [
        _summarize_entry(entry, manifest_path.parent, runs_dir, verify_artifacts)
        for entry in entries
    ]
    counts = {status: sum(row["status"] == status for row in rows) for status in _STATUSES}
    complete_inventory = len(entries) == expected_count
    metric_summary = []
    names = list(FULL_METRICS) + sorted(
        {m["metric"] for row in rows for m in row["metric_comparisons"]} - set(FULL_METRICS)
    )
    for name in names:
        metrics = [
            metric
            for row in rows
            for metric in row["metric_comparisons"]
            if metric["metric"] == name
        ]
        metric_summary.append(
            {
                "metric": name,
                "compared_genomes": len(metrics),
                "finite_pairs": sum(metric["finite"] for metric in metrics),
                "matching_genomes": sum(metric["matches"] for metric in metrics),
                "bitwise_equal_genomes": sum(metric["bitwise_equal"] for metric in metrics),
                "max_absolute_difference": max(
                    (
                        m["absolute_difference"]
                        for m in metrics
                        if m["absolute_difference"] is not None
                    ),
                    default=None,
                ),
                "max_relative_difference": max(
                    (
                        m["relative_difference"]
                        for m in metrics
                        if m["relative_difference"] is not None
                    ),
                    default=None,
                ),
            }
        )
    versions: dict[str, dict[str, set[str]]] = {"python": defaultdict(set), "r": defaultdict(set)}
    blast_versions = set()
    for row in rows:
        for label, key in (("python", "software"), ("r", "r_versions")):
            for name, value in row[key].items():
                versions[label][name].add(str(value))
        blast = row["dependencies"].get("blastp", {})
        if isinstance(blast, dict) and blast.get("version") is not None:
            blast_versions.add(str(blast["version"]))
    protein_entries = [entry for entry in entries if entry.get("source_kind") == "protein"]
    protein_audit = {
        "protein_inputs": len(protein_entries),
        "matching_compressed_archives": sum(
            entry.get("download", {}).get("source_protein_md5_matches") is True
            for entry in protein_entries
        ),
        "different_compressed_archives": sum(
            entry.get("download", {}).get("source_protein_md5_matches") is False
            for entry in protein_entries
        ),
        "matching_ids_sequences_and_order": sum(
            entry.get("download", {}).get("source_protein_content_matches") is True
            and entry.get("download", {}).get("source_protein_order_matches") is True
            for entry in protein_entries
        ),
        "annotation_header_only_differences": [
            entry["accession"]
            for entry in protein_entries
            if entry.get("download", {}).get("source_protein_difference_kind")
            == "annotation_headers_differ_protein_sequences_identical"
        ],
    }
    report = {
        "schema_version": 1,
        "generated_at": datetime.now(UTC).isoformat(),
        "manifest": {
            "path": str(manifest_path),
            "sha256": _sha256(manifest_path),
            "metadata": {key: value for key, value in manifest.items() if key != "entries"}
            if isinstance(manifest, dict)
            else {},
        },
        "runs_directory": str(runs_dir),
        "expected_assemblies": expected_count,
        "actual_assemblies": len(entries),
        "complete_inventory": complete_inventory,
        "source_kind_counts": dict(
            Counter(entry.get("source_kind", "unknown") for entry in entries)
        ),
        "protein_input_audit": protein_audit,
        "counts": counts,
        "all_pass": complete_inventory and counts["pass"] == expected_count,
        "all_bitwise_equal": complete_inventory
        and counts["pass"] == expected_count
        and all(row["bitwise_equal"] for row in rows),
        "tolerance": {
            "relative": RTOL,
            "absolute": ATOL,
            "criterion": "abs(python-r) <= max(atol, rtol * max(abs(python), abs(r)))",
        },
        "required_metrics": list(FULL_METRICS),
        "verify_artifacts": verify_artifacts,
        "metric_summary": metric_summary,
        "observed_versions": {
            label: {name: sorted(values) for name, values in names.items()}
            for label, names in versions.items()
        },
        "blastp_versions": sorted(blast_versions),
        "errors": []
        if complete_inventory
        else [f"Expected {expected_count} assemblies, found {len(entries)} inventory entries"],
        "genomes": rows,
    }
    return _json_safe(report)


def _cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def _number(value: Any) -> str:
    return f"{value:.12g}" if _finite(value) else "—"


def render_markdown(report: Mapping[str, Any]) -> str:
    """Render the audit without hiding failed, missing, or nonfinite results."""
    counts = report["counts"]
    outcome = "PASS" if report["all_pass"] else "INCOMPLETE OR FAILED"
    source = report["source_kind_counts"]
    protein_audit = report["protein_input_audit"]
    protein_note = (
        f"The input preparation audit recorded identical protein identifiers, sequences, and order "
        f"for {protein_audit['matching_ids_sequences_and_order']} of "
        f"{protein_audit['protein_inputs']} protein inputs when compared with the current files for "
        "those exact NCBI assembly accessions. "
        f"{protein_audit['matching_compressed_archives']} compressed protein archives matched "
        f"the official checksum; {protein_audit['different_compressed_archives']} differed."
    )
    header_only = protein_audit["annotation_header_only_differences"]
    if header_only:
        protein_note += (
            f" For {len(header_only)} of the differing files, only descriptive annotation headers "
            "changed: the protein IDs, sequences, and order were identical. These are annotation "
            "text differences, not evidence of different assembly genomes. Accessions: "
            + ", ".join(f"`{value}`" for value in header_only)
            + "."
        )
    lines = [
        "# Folder validation: Python against original R gRodon",
        "",
        f"**{outcome}: {counts['pass']} of {report['expected_assemblies']} assemblies passed.** "
        f"{counts['mismatch']} numerical/settings/count mismatches; {counts['error']} errors; "
        f"{counts['not_run']} not run. Inventory entries: {report['actual_assemblies']}.",
        "",
        "## Inputs and comparison method",
        "",
        f"The supplied inventory contains {source.get('genomic', 0)} genomic DNA files and "
        f"{source.get('protein', 0)} protein-only files. Protein sequences cannot recover synonymous "
        "codons. Protein-only inputs require genomic DNA fetched for the same assembly accessions; "
        "the proteins are not back-translated. Each accession and both source/genome SHA256 "
        "identities are retained in the machine-readable report.",
        "",
        protein_note,
        "",
        "Both implementations consume the same predicted coding sequences and the same selected "
        "ribosomal gene IDs. Each row compares all twelve full-mode statistics, prepared gene "
        "counts, retained highly expressed gene counts, and prediction settings. A pass requires "
        "finite values and positive, ordered doubling-time confidence bounds. This checks agreement "
        "with the original implementation, not biological accuracy against measured growth rates.",
        "",
        "Numerical acceptance: `abs(Python − R) <= max(1e-12, 1e-10 × max(abs(Python), abs(R)))`. "
        "The separate exact-equality count compares binary64 encodings after JSON parsing, including "
        "the sign of zero. Numerical agreement does not imply bitwise identity. R serializes its "
        "reported floating-point values to JSON, which can itself introduce final-digit rounding.",
        "",
        (
            "Every produced run artifact is checked against the SHA256 and byte-size manifest; "
            "missing, modified, or unrecorded artifacts fail the audit."
            if report["verify_artifacts"]
            else "**Artifact hash verification was explicitly disabled for this report.**"
        ),
        "",
        f"Inventory SHA256: `{report['manifest']['sha256']}`.",
        "",
        "## Software and settings",
        "",
    ]
    for label, versions in report["observed_versions"].items():
        lines.append(
            f"- {label.upper()}: "
            + "; ".join(
                f"{_cell(name)} {', '.join(map(_cell, values))}"
                for name, values in sorted(versions.items())
            )
        )
    lines.extend(["- BLASTP: " + (", ".join(report["blastp_versions"]) or "not recorded"), ""])
    configurations = {
        json.dumps(row["settings"], sort_keys=True) for row in report["genomes"] if row["settings"]
    }
    lines.append(
        "Observed pipeline settings (all distinct configurations). Recorded executable paths "
        "are historical provenance from the validation machine; use the portable reproduction "
        "commands below for a new checkout."
    )
    for configuration in sorted(configurations):
        lines.extend(["", "```json", json.dumps(json.loads(configuration), indent=2), "```"])
    lines.extend(
        [
            "",
            "## Every assembly",
            "",
            "Doubling times and their absolute difference are shown in minutes. `nHE` is the "
            "retained highly expressed gene count. The result column gives the run status, "
            "statistics within tolerance, and bitwise-equal statistics. Comparisons use unrounded "
            "values; only this display is rounded. Full confidence intervals and hour values are "
            "preserved in the JSON and CSV exports.",
            "",
            "| Assembly | Organism | Input | CDS | nHE Python/R | Python min | R min | Absolute Δ min | Result (tolerance; exact) |",
            "|---|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in report["genomes"]:
        values, number_counts = row["doubling_time"], row["counts"]
        exact = sum(metric["bitwise_equal"] for metric in row["metric_comparisons"])
        matching = sum(metric["matches"] for metric in row["metric_comparisons"])
        cells = [
            row["accession"],
            row["organism"],
            row["source_kind"],
            number_counts.get("python_predicted_genes", "—"),
            f"{_number(number_counts.get('python_nHE'))}/{_number(number_counts.get('r_nHE'))}",
            _number(values.get("python_minutes")),
            _number(values.get("r_minutes")),
            _number(values.get("absolute_difference_minutes")),
            f"{row['status']} ({matching}/12; {exact}/12)",
        ]
        lines.append("| " + " | ".join(map(_cell, cells)) + " |")
    lines.extend(
        [
            "",
            "## Agreement by statistic",
            "",
            "| Statistic | Compared | Finite pairs | Within tolerance | Bitwise equal | Max absolute Δ | Max relative Δ |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in report["metric_summary"]:
        lines.append(
            "| "
            + " | ".join(
                map(
                    _cell,
                    [
                        row["metric"],
                        row["compared_genomes"],
                        row["finite_pairs"],
                        row["matching_genomes"],
                        row["bitwise_equal_genomes"],
                        _number(row["max_absolute_difference"]),
                        _number(row["max_relative_difference"]),
                    ],
                )
            )
            + " |"
        )
    lines.extend(["", "## Problems and warnings", ""])
    for message in report["errors"]:
        lines.append(f"- Inventory: {_cell(message)}")
    for row in report["genomes"]:
        for message in row["errors"] + row["mismatches"]:
            lines.append(f"- `{row['accession']}` ({row['status']}): {_cell(message)}")
        if row["status"] == "not_run":
            lines.append(f"- `{row['accession']}`: no run directory was found.")
    if report["all_pass"]:
        lines.append("No input, artifact, execution, or numerical comparison failures were found.")
    warning_groups: dict[str, list[str]] = defaultdict(list)
    for row in report["genomes"]:
        for language, messages in row["warnings"].items():
            for message in messages:
                warning_groups[f"{language}: {message}"].append(row["accession"])
    if warning_groups:
        lines.extend(["", "Warnings are retained even for passing comparisons:", ""])
        for message, accessions in sorted(warning_groups.items()):
            lines.append(f"- {_cell(message)} — " + ", ".join(f"`{value}`" for value in accessions))
    else:
        lines.extend(["", "No warnings were recorded in the available results."])
    lines.extend(
        [
            "",
            "## Reproduce this test case",
            "",
            "The archived report, input inventory, and compact reference fixture preserve their "
            "original contents and checksum relationships. Historical absolute paths are provenance, "
            "not required locations on another machine. Raw input genomes, the original protein "
            "folder, local tools, and complete run directories are not distributed. See the "
            "[benchmark archive notes](../validation/README.md).",
            "",
            "The public reproduction below fetches the same forty genomic DNA files by their "
            "archived NCBI URLs and SHA256 checksums; it does not need the original source folder. "
            "The new inventory therefore contains forty DNA inputs, whereas the historical inventory "
            "records three DNA and thirty-seven protein inputs. The genomic sequence bytes used for "
            "prediction are the same. New reports are written separately from the archived evidence.",
            "",
            "Run these PowerShell commands from the repository root with Python 3.11 or later. "
            "The genome extra pins Pyrodigal 3.7.1, and setup installs BLASTP 2.17.0 into `.tools`. "
            "For the full comparison, separately install original R gRodon 2.7.3, its dependencies, "
            "and jsonlite, then replace the example Rscript path with your installation. The "
            "observed versions above describe the archived runs.",
            "",
            "```powershell",
            "python -m venv .venv",
            '.\\.venv\\Scripts\\python.exe -m pip install -e ".[dev,genome]"',
            ".\\.venv\\Scripts\\python.exe -m grodon setup --tools-dir .tools/blast-2.17.0 --version 2.17.0",
            ".\\.venv\\Scripts\\python.exe scripts/fetch_benchmark_genomes.py --manifest validation/test_genomes_inputs_manifest.json --output-dir validation/test_genomes_inputs/public-genomes --fetch",
            ".\\.venv\\Scripts\\python.exe scripts/prepare_folder_genomes.py --input-dir validation/test_genomes_inputs/public-genomes --cache-dir validation/test_genomes_inputs/ncbi --manifest validation/reproduced_inputs.json --fetch",
            '.\\.venv\\Scripts\\python.exe scripts/validate_genome_folder.py --manifest validation/reproduced_inputs.json --runs-dir validation/test_genomes_runs --blastp .tools/blast-2.17.0/bin/blastp.exe --rscript "C:\\path\\to\\R\\bin\\Rscript.exe" --workers 8',
            ".\\.venv\\Scripts\\python.exe scripts/report_folder_validation.py --manifest validation/reproduced_inputs.json --runs-dir validation/test_genomes_runs --output-json validation/reproduced_comparison.json --output-csv validation/reproduced_comparison.csv --output-md validation/reproduced_comparison.md",
            "```",
            "",
            "On Linux or macOS, use `python3 -m venv .venv`, `.venv/bin/python`, and "
            "`.tools/blast-2.17.0/bin/blastp`, with the appropriate Rscript path. "
            "For an offline repeat, omit `--fetch` from both genome fetching and preparation to "
            "verify the existing input cache. "
            "Add `--resume` to the runner to reuse completed runs only after their recorded inputs, "
            "settings, and artifact hashes have been checked. To regenerate every analysis, select a "
            "fresh `--runs-dir` for both the runner and reporter and omit `--resume`. After dependencies "
            "are installed, network access is confined to input preparation; prediction, R comparison, "
            "reporting, and cached-input verification run locally.",
            "",
            "After a new complete run passes, its numeric R results can be exported separately. "
            "The second command below runs the existing archived model regression tests offline; "
            "it needs neither R nor the downloaded genomes:",
            "",
            "```powershell",
            ".\\.venv\\Scripts\\python.exe scripts/export_folder_reference.py --report validation/reproduced_comparison.json --output validation/reproduced_growth_reference.json",
            ".\\.venv\\Scripts\\python.exe -m pytest tests/test_folder_growth_reference.py",
            "```",
            "",
            "That regression test compares fitted-model predictions with recorded original R results. "
            "It does not repeat gene prediction, BLAST, or codon-statistic computation; the full workflow "
            "above exercises those stages.",
            "",
            "Repository artifacts: [full comparison JSON](../validation/test_genomes_comparison.json), "
            "[per-assembly CSV](../validation/test_genomes_comparison.csv), "
            "[input inventory](../validation/test_genomes_inputs_manifest.json), and "
            "[compact model regression fixture](../tests/data/folder_growth_reference.json). "
            "Input provenance follows the official "
            "[NCBI genome download documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/) "
            "and [NCBI Datasets REST API documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/).",
        ]
    )
    lines.extend(
        [
            "",
            f"Generated at {report['generated_at']}. The JSON companion preserves every "
            "per-metric value, difference, exact-equality flag, setting, version, artifact hash, "
            "input identity check, and error.",
            "",
        ]
    )
    return "\n".join(lines)


def render_csv(report: Mapping[str, Any]) -> str:
    """Return one flat CSV row per inventory accession, including failed runs."""
    fields = [
        "accession",
        "organism",
        "source_kind",
        "status",
        "predicted_genes",
        "python_nHE",
        "r_nHE",
        "python_doubling_hours",
        "r_doubling_hours",
        "python_lower_ci_hours",
        "python_upper_ci_hours",
        "r_lower_ci_hours",
        "r_upper_ci_hours",
        "python_doubling_minutes",
        "r_doubling_minutes",
        "max_absolute_difference",
        "max_relative_difference",
        "matches",
        "numerically_matches",
        "bitwise_equal",
        "bitwise_equal_statistics",
        "source_sha256",
        "genome_sha256",
        "errors",
        "mismatches",
    ]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    for genome in report["genomes"]:
        times, counts = genome["doubling_time"], genome["counts"]
        row = {
            key: genome[key]
            for key in (
                "accession",
                "organism",
                "source_kind",
                "status",
                "max_absolute_difference",
                "max_relative_difference",
                "numerically_matches",
                "bitwise_equal",
            )
        }
        row.update(
            predicted_genes=counts.get("python_predicted_genes"),
            python_nHE=counts.get("python_nHE"),
            r_nHE=counts.get("r_nHE"),
            matches=genome["status"] == "pass",
            bitwise_equal_statistics=sum(
                metric["bitwise_equal"] for metric in genome["metric_comparisons"]
            ),
            source_sha256=genome["inputs"].get("source_sha256"),
            genome_sha256=genome["inputs"].get("genome_sha256"),
            errors="; ".join(genome["errors"]),
            mismatches="; ".join(genome["mismatches"]),
        )
        for language in ("python", "r"):
            row[f"{language}_doubling_hours"] = times.get(f"{language}_hours")
            row[f"{language}_doubling_minutes"] = times.get(f"{language}_minutes")
            for bound in ("lower", "upper"):
                row[f"{language}_{bound}_ci_hours"] = times.get(f"{language}_{bound}_ci_hours")
        writer.writerow(row)
    return output.getvalue()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--runs-dir", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path)
    parser.add_argument("--expected-count", type=int, default=40)
    parser.add_argument("--no-verify-artifacts", action="store_true")
    arguments = parser.parse_args(argv)
    outputs = [arguments.output_json.resolve(), arguments.output_md.resolve()]
    if arguments.output_csv is not None:
        outputs.append(arguments.output_csv.resolve())
    if len(set(outputs)) != len(outputs) or arguments.manifest.resolve() in outputs:
        parser.error("Report outputs must be different files and must not overwrite the inventory")
    if any(path.is_relative_to(arguments.runs_dir.resolve()) for path in outputs):
        parser.error("Report outputs must be outside the audited runs directory")
    report = summarize(
        arguments.manifest,
        arguments.runs_dir,
        expected_count=arguments.expected_count,
        verify_artifacts=not arguments.no_verify_artifacts,
    )
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    outputs[0].write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    outputs[1].write_text(render_markdown(report), encoding="utf-8")
    if arguments.output_csv is not None:
        outputs[2].write_text(render_csv(report), encoding="utf-8", newline="")
    print(
        json.dumps(
            {
                "counts": report["counts"],
                "all_pass": report["all_pass"],
                "output_json": str(outputs[0]),
                "output_md": str(outputs[1]),
                "output_csv": str(outputs[2]) if len(outputs) > 2 else None,
            }
        )
    )
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
