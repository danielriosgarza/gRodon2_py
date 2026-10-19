"""Export original-R model regression values from a complete audited folder report.

The compact fixture stores numeric statistics and provenance, never sequences.
It exercises model evaluation from archived R codon statistics; the separate
folder report documents the measured assembly-to-prediction comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

EXPECTED_COUNT = 40
METRICS = {
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
}
TOLERANCE = {"relative": 1e-10, "absolute": 1e-12}


def _sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[a-f0-9]{64}", value) is None:
        raise ValueError(f"Missing or malformed SHA256 for {label}")
    return value


def _settings(values: Any) -> dict[str, Any]:
    if not isinstance(values, dict):
        raise ValueError("Prediction settings are missing")
    settings = {
        key: values.get(key)
        for key in ("mode", "training_set", "temperature", "genetic_code", "confidence")
    }
    if settings["temperature"] == "none":
        settings["temperature"] = None
    settings["fragments"] = values.get("fragments", False)
    expected = {
        "mode": "full",
        "training_set": "madin",
        "temperature": None,
        "genetic_code": 11,
        "confidence": 0.95,
        "fragments": False,
    }
    if settings != expected:
        raise ValueError(
            "The folder fixture requires full/madin/code11/default-temperature settings"
        )
    return settings


def _record(row: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    accession = row.get("accession")
    if not isinstance(accession, str) or re.fullmatch(r"GC[AF]_\d{9}\.\d+", accession) is None:
        raise ValueError(f"Invalid accession: {accession!r}")
    if (
        row.get("status") != "pass"
        or row.get("numerically_matches") is not True
        or row.get("full_metric_set") is not True
        or row.get("errors")
        or row.get("mismatches")
    ):
        raise ValueError(f"The report contains an incomplete or failing row: {accession}")
    audit = row.get("artifact_verification", {})
    if audit.get("enabled") is not True or audit.get("verified") is not True or audit.get("errors"):
        raise ValueError(f"Artifacts were not independently verified: {accession}")
    artifacts = audit.get("files")
    if not isinstance(artifacts, list) or not artifacts:
        raise ValueError(f"Artifact hashes are missing: {accession}")
    hashes = {}
    for artifact in artifacts:
        name = artifact.get("name")
        if not isinstance(name, str) or name in hashes or artifact.get("matches") is not True:
            raise ValueError(f"Artifact audit is incomplete or duplicated: {accession}")
        hashes[name] = _sha(artifact.get("sha256"), f"{accession}/{name}")
    required = (
        "predicted_cds.fna",
        "highly_expressed_ids.txt",
        "r_result.json",
        "r_comparison.json",
        "reference.json",
        "reference.faa.gz",
    )
    if any(name not in hashes for name in required):
        raise ValueError(f"Required prepared/R/reference artifact hashes are missing: {accession}")
    statistics = row.get("statistics", {}).get("r")
    if not isinstance(statistics, dict) or set(statistics) != METRICS:
        raise ValueError(f"All twelve original R statistics are required: {accession}")
    if not all(
        isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
        for value in statistics.values()
    ):
        raise ValueError(f"Original R statistics must be finite numeric scalars: {accession}")
    if not 0 < statistics["LowerCI"] <= statistics["d"] <= statistics["UpperCI"]:
        raise ValueError(f"Original R confidence bounds are invalid: {accession}")
    for name in ("nHE", "FilteredSequences"):
        if statistics[name] < 0 or int(statistics[name]) != statistics[name]:
            raise ValueError(f"Original R gene counts must be nonnegative integers: {accession}")
    prediction_settings = _settings(row.get("settings"))
    if _settings(row.get("r_settings")) != prediction_settings:
        raise ValueError(f"Python/R settings differ: {accession}")
    reference = row.get("reference")
    if not isinstance(reference, dict) or not reference.get("name"):
        raise ValueError(f"Reference provenance is missing: {accession}")
    _sha(reference.get("sha256"), f"{accession} reference")
    model_provenance = row.get("model_provenance")
    versions = row.get("r_versions")
    if not isinstance(model_provenance, dict) or not model_provenance:
        raise ValueError(f"Model provenance is missing: {accession}")
    if not isinstance(versions, dict) or not {"R", "gRodon"} <= set(versions):
        raise ValueError(f"Original R software versions are missing: {accession}")
    profile = {
        "model": row["model"],
        "prediction_settings": prediction_settings,
        "pipeline_settings": row["settings"],
        "reference": reference,
        "model_provenance": model_provenance,
        "software": {
            "python": row.get("software", {}),
            "r": versions,
            "blastp": row.get("dependencies", {}).get("blastp", {}).get("version"),
        },
    }
    inputs = row.get("inputs", {})
    record = {
        "accession": accession,
        "source_sha256": _sha(inputs.get("source_sha256"), f"{accession} source"),
        "genome_sha256": _sha(inputs.get("genome_sha256"), f"{accession} genome"),
        "artifact_sha256": {name: hashes[name] for name in required},
        "r_statistics": statistics,
    }
    return record, profile


def export_reference(report_path: Path, output_path: Path) -> dict[str, Any]:
    """Validate an audited forty-assembly report, then save portable R values."""
    report_path, output_path = report_path.resolve(), output_path.resolve()
    if report_path == output_path:
        raise ValueError("The fixture must not overwrite its source report")
    raw = report_path.read_bytes()
    report = json.loads(raw)
    if not isinstance(report, dict) or report.get("schema_version") != 1:
        raise ValueError("Unsupported folder report schema")
    if report.get("all_pass") is not True or report.get("complete_inventory") is not True:
        raise ValueError("Only complete, independently audited, all-pass reports can be exported")
    if report.get("verify_artifacts") is not True:
        raise ValueError("Artifact verification must be enabled")
    if (
        report.get("expected_assemblies") != EXPECTED_COUNT
        or report.get("actual_assemblies") != EXPECTED_COUNT
        or report.get("counts") != {"pass": EXPECTED_COUNT, "mismatch": 0, "error": 0, "not_run": 0}
        or report.get("errors")
    ):
        raise ValueError("The fixture requires exactly forty completed, passing assemblies")
    tolerance = report.get("tolerance", {})
    if any(tolerance.get(key) != value for key, value in TOLERANCE.items()):
        raise ValueError("The audited comparison tolerance differs from the regression tolerance")
    rows = report.get("genomes")
    if not isinstance(rows, list) or len(rows) != EXPECTED_COUNT:
        raise ValueError("The report must contain forty per-assembly audit records")
    records: list[dict[str, Any]] = []
    profiles: dict[str, dict[str, Any]] = {}
    profile_ids: dict[str, str] = {}
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Malformed per-assembly audit record")
        record, profile = _record(row)
        encoded = json.dumps(profile, sort_keys=True, allow_nan=False)
        if encoded not in profile_ids:
            profile_id = f"provenance_{len(profiles) + 1}"
            profile_ids[encoded] = profile_id
            profiles[profile_id] = profile
        record["provenance_id"] = profile_ids[encoded]
        records.append(record)
    if len({row["accession"] for row in records}) != EXPECTED_COUNT:
        raise ValueError("Duplicate accessions in audited report")
    bundle = {
        "schema_version": 1,
        "scope": "Model evaluation from archived original R codon statistics; assembly-to-prediction equivalence is measured by the separate audited folder report.",
        "generator": "scripts/export_folder_reference.py",
        "source_report": {
            "name": report_path.name,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "generated_at": report.get("generated_at"),
            "manifest_sha256": _sha(report.get("manifest", {}).get("sha256"), "input manifest"),
        },
        "tolerance": TOLERANCE,
        "provenance": profiles,
        "records": sorted(records, key=lambda row: row["accession"]),
    }
    encoded = json.dumps(bundle, indent=2, allow_nan=False) + "\n"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(encoded)
    return bundle


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args(argv)
    try:
        bundle = export_reference(arguments.report, arguments.output)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"reference export error: {exc}\n")
    print(f"Exported {len(bundle['records'])} audited original-R model reference records")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
