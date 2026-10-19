"""Offline model evaluation against forty archived actual original-R results.

These tests supply R's already calculated codon statistics to the Python model.
They do not test sequence statistics, gene calling, or ribosomal annotation; full
assembly-to-prediction equivalence is measured by the separate folder audit.
"""

from __future__ import annotations

import importlib.util
import json
import re
import warnings
from pathlib import Path

import pytest

from grodon import GrowthPredictionWarning, predict_from_statistics

REFERENCE = Path(__file__).with_name("data") / "folder_growth_reference.json"
BUNDLE = json.loads(REFERENCE.read_text(encoding="utf-8"))
CASES = BUNDLE["records"]


@pytest.fixture(scope="module")
def exporter():
    path = Path(__file__).resolve().parents[1] / "scripts" / "export_folder_reference.py"
    spec = importlib.util.spec_from_file_location("grodon_folder_reference_exporter", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_archive_has_forty_unique_records_and_audited_provenance():
    assert BUNDLE["schema_version"] == 1
    assert len(BUNDLE["records"]) == 40
    assert len({record["accession"] for record in BUNDLE["records"]}) == 40
    assert BUNDLE["tolerance"] == {"relative": 1e-10, "absolute": 1e-12}
    assert "Model evaluation" in BUNDLE["scope"]
    assert re.fullmatch(r"[a-f0-9]{64}", BUNDLE["source_report"]["sha256"])
    assert re.fullmatch(r"[a-f0-9]{64}", BUNDLE["source_report"]["manifest_sha256"])
    for record in BUNDLE["records"]:
        assert re.fullmatch(r"[a-f0-9]{64}", record["genome_sha256"])
        assert re.fullmatch(r"[a-f0-9]{64}", record["source_sha256"])
        for digest in record["artifact_sha256"].values():
            assert re.fullmatch(r"[a-f0-9]{64}", digest)
        profile = BUNDLE["provenance"][record["provenance_id"]]
        assert profile["model"] == "gRodon_model_base_madin"
        assert profile["reference"]["name"]
        assert {"R", "gRodon"} <= set(profile["software"]["r"])
        assert profile["model_provenance"]


@pytest.mark.parametrize("record", CASES, ids=lambda record: record["accession"])
def test_python_model_reproduces_archived_original_r_predictions(record):
    profile = BUNDLE["provenance"][record["provenance_id"]]
    settings = profile["prediction_settings"]
    original = record["r_statistics"]
    statistics = {
        key: value for key, value in original.items() if key not in {"d", "LowerCI", "UpperCI"}
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GrowthPredictionWarning)
        prediction = predict_from_statistics(
            statistics,
            mode=settings["mode"],
            training_set=settings["training_set"],
            temperature=settings["temperature"],
            confidence=settings["confidence"],
        )
    assert prediction.model_name == profile["model"]
    for name in ("d", "LowerCI", "UpperCI"):
        assert prediction[name] == pytest.approx(original[name], rel=1e-10, abs=1e-12), name


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("all_pass", False, "all-pass"),
        ("complete_inventory", False, "all-pass"),
        ("verify_artifacts", False, "verification"),
        ("actual_assemblies", 39, "forty"),
        ("counts", {"pass": 39, "mismatch": 1, "error": 0, "not_run": 0}, "forty"),
    ],
)
def test_export_refuses_partial_or_unaudited_reports_without_touching_output(
    exporter, tmp_path, field, value, message
):
    # This is only a deliberately incomplete audit envelope, not fabricated R values.
    report = {
        "schema_version": 1,
        "all_pass": True,
        "complete_inventory": True,
        "verify_artifacts": True,
        "expected_assemblies": 40,
        "actual_assemblies": 40,
        "counts": {"pass": 40, "mismatch": 0, "error": 0, "not_run": 0},
    }
    report[field] = value
    source = tmp_path / "report.json"
    source.write_text(json.dumps(report), encoding="utf-8")
    output = tmp_path / "existing-reference.json"
    output.write_bytes(b"preserve existing reference")
    with pytest.raises(ValueError, match=message):
        exporter.export_reference(source, output)
    assert output.read_bytes() == b"preserve existing reference"


def test_export_cannot_overwrite_its_source_report(exporter, tmp_path):
    source = tmp_path / "report.json"
    source.write_bytes(b"original audit")
    with pytest.raises(ValueError, match="source report"):
        exporter.export_reference(source, source)
    assert source.read_bytes() == b"original audit"
