"""Audit-boundary tests; no genomes, BLAST, R, or worker processes are executed."""

from __future__ import annotations

import importlib.util
import json
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import pytest

from grodon.pipeline import GenomePipelineConfig
from grodon.references import sha256_file


@pytest.fixture(scope="module")
def runner():
    path = Path(__file__).resolve().parents[1] / "scripts" / "validate_genome_folder.py"
    spec = importlib.util.spec_from_file_location("grodon_folder_validation_runner", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def entry(tmp_path):
    source = tmp_path / "source metadata.txt"
    genome = tmp_path / "synthetic input.fna"
    source.write_text("Synthetic source metadata\n", encoding="utf-8")
    genome.write_text(">artificial_sequence\nACGTACGT\n", encoding="utf-8")
    return {
        "accession": "GCF_000000001.1",
        "status": "ready",
        "source_file": str(source),
        "source_sha256": sha256_file(source),
        "genome_path": str(genome),
        "genome_sha256": sha256_file(genome),
    }


@pytest.fixture
def config():
    return GenomePipelineConfig(blastp="fake-blastp", timeout=60)


@pytest.fixture
def result():
    return {
        "status": "complete",
        "r_comparison_status": "complete",
        "r_matches": True,
        "doubling_time_hours": 2.0,
        "scientific_fingerprint": "a" * 64,
        "statistics": {
            "CUBHE": 0.7,
            "GC": 0.5,
            "GCdiv": 0.0,
            "ConsistencyHE": 0.4,
            "CUB": 0.5,
            "CPB": 0.1,
            "FilteredSequences": 2,
            "nHE": 12,
            "dCUB": -0.4,
            "d": 2.0,
            "LowerCI": 1.5,
            "UpperCI": 2.5,
        },
    }


@pytest.fixture
def pipeline_mock(monkeypatch, runner, result):
    calls = []

    def run_genome(genome, output, **kwargs):
        calls.append((genome, output, kwargs))
        return SimpleNamespace(to_dict=lambda: deepcopy(result))

    monkeypatch.setattr(runner, "run_genome", run_genome)
    return calls


def _write_json(path, value):
    path.write_text(json.dumps(value), encoding="utf-8")


def _refresh_artifacts(folder):
    manifest = json.loads((folder / "run.json").read_text(encoding="utf-8"))
    manifest["artifacts"] = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in folder.iterdir()
        if path.is_file() and path.name != "run.json"
    }
    _write_json(folder / "run.json", manifest)


@pytest.fixture
def completed_run(tmp_path, entry, config, result):
    folder = tmp_path / "completed"
    folder.mkdir()
    _write_json(folder / "result.json", result)
    _write_json(folder / "r_result.json", {"statistics": result["statistics"]})
    _write_json(
        folder / "r_comparison.json",
        {"matches": True, "python_statistics": result["statistics"]},
    )
    _write_json(folder / "r_parameters.json", {"mode": "full", "genetic_code": 11})
    (folder / "predicted_cds.fna").write_text(">artificial_gene\nATGAAATAA\n", encoding="utf-8")
    (folder / "highly_expressed_ids.txt").write_text("artificial_gene\n", encoding="utf-8")
    (folder / "r_comparison.log").write_text("Mock R completed\n", encoding="utf-8")
    _write_json(
        folder / "run.json",
        {
            "status": "complete",
            "input": {"sha256": entry["genome_sha256"]},
            "settings": asdict(config),
        },
    )
    _refresh_artifacts(folder)
    return folder


def test_one_passes_verified_input_and_exact_settings_to_real_pipeline_boundary(
    runner, tmp_path, entry, config, result, pipeline_mock
):
    output = tmp_path / "new output"
    record = runner._one(entry, str(output), asdict(config), "explicit-rscript", False)
    assert record["status"] == "pass"
    assert record["all_python_statistics_finite"] is True
    assert record["full_metric_set"] is True
    assert record["pipeline_complete"] is True
    assert record["r_matches"] is True
    assert record["reused"] is False
    assert record["scientific_fingerprint"] == result["scientific_fingerprint"]
    assert record["elapsed_seconds"] >= 0
    assert pipeline_mock == [
        (entry["genome_path"], output, {"config": config, "rscript": "explicit-rscript"})
    ]


@pytest.mark.parametrize("input_key", ["source_file", "genome_path"])
def test_changed_input_is_an_error_before_pipeline_or_output(
    runner, tmp_path, entry, config, pipeline_mock, input_key
):
    Path(entry[input_key]).write_text("changed after preparation", encoding="utf-8")
    output = tmp_path / "output"
    record = runner._one(entry, str(output), asdict(config), "rscript", False)
    assert record["status"] == "error"
    assert "checksum mismatch" in record["error"]
    assert not output.exists()
    assert not pipeline_mock


def test_input_not_ready_does_not_launch_pipeline(runner, tmp_path, entry, config, pipeline_mock):
    entry.update(status="failed", error="preparation failed")
    record = runner._one(entry, str(tmp_path / "output"), asdict(config), "rscript", False)
    assert record["status"] == "error"
    assert "preparation failed" in record["error"]
    assert not pipeline_mock


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf"), None, True, "2"])
def test_nonfinite_or_nonnumeric_statistic_never_passes(
    runner, tmp_path, entry, config, result, pipeline_mock, value
):
    result["statistics"]["UpperCI"] = value
    record = runner._one(entry, str(tmp_path / "output"), asdict(config), "rscript", False)
    assert record["status"] == "mismatch"
    assert record["all_python_statistics_finite"] is False
    assert record["r_matches"] is True


def test_empty_statistics_never_pass(runner, tmp_path, entry, config, result, pipeline_mock):
    result["statistics"] = {}
    record = runner._one(entry, str(tmp_path / "output"), asdict(config), "rscript", False)
    assert record["status"] == "mismatch"
    assert record["all_python_statistics_finite"] is False


@pytest.mark.parametrize("change", ["only_d", "missing_cpb", "extra_metric", "renamed_metric"])
def test_exact_full_metric_set_is_required_even_with_finite_values_and_r_agreement(
    runner, tmp_path, entry, config, result, pipeline_mock, change
):
    if change == "only_d":
        result["statistics"] = {"d": 2.0}
    elif change == "missing_cpb":
        del result["statistics"]["CPB"]
    elif change == "extra_metric":
        result["statistics"]["unexpected"] = 1.0
    else:
        result["statistics"]["cpb"] = result["statistics"].pop("CPB")
    record = runner._one(entry, str(tmp_path / "output"), asdict(config), "rscript", False)
    assert record["status"] == "mismatch"
    assert record["all_python_statistics_finite"] is True
    assert record["full_metric_set"] is False
    assert record["r_matches"] is True


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("status", "prediction_complete_r_pending"),
        ("status", "failed"),
        ("status", None),
        ("r_comparison_status", "running"),
        ("r_comparison_status", "failed"),
        ("r_comparison_status", None),
    ],
)
def test_incomplete_pipeline_or_comparison_cannot_pass_despite_complete_statistics(
    runner, tmp_path, entry, config, result, pipeline_mock, field, value
):
    result[field] = value
    record = runner._one(entry, str(tmp_path / "output"), asdict(config), "rscript", False)
    assert record["status"] == "mismatch"
    assert record["all_python_statistics_finite"] is True
    assert record["full_metric_set"] is True
    assert record["pipeline_complete"] is False
    assert record["r_matches"] is True


@pytest.mark.parametrize("r_matches", [False, None, 1, "true"])
def test_only_explicit_r_agreement_can_pass(
    runner, tmp_path, entry, config, result, pipeline_mock, r_matches
):
    result["r_matches"] = r_matches
    record = runner._one(entry, str(tmp_path / "output"), asdict(config), "rscript", False)
    assert record["status"] == "mismatch"
    assert record["all_python_statistics_finite"] is True


def test_pipeline_exception_is_captured_with_accession_and_elapsed_time(
    runner, monkeypatch, tmp_path, entry, config
):
    def fail(*args, **kwargs):
        raise RuntimeError("external comparison failed")

    monkeypatch.setattr(runner, "run_genome", fail)
    record = runner._one(entry, str(tmp_path / "output"), asdict(config), "rscript", False)
    assert record["status"] == "error"
    assert record["error"] == "RuntimeError: external comparison failed"
    assert record["accession"] == entry["accession"]
    assert record["elapsed_seconds"] >= 0


def test_existing_output_is_preserved_when_pipeline_refuses_overwrite(
    runner, monkeypatch, tmp_path, entry, config
):
    output = tmp_path / "existing"
    output.mkdir()
    protected = output / "user work.txt"
    protected.write_bytes(b"preserve these exact bytes")

    def refuse_overwrite(genome, path, **kwargs):
        assert path == output and protected.read_bytes() == b"preserve these exact bytes"
        raise FileExistsError("Refusing to overwrite existing output")

    monkeypatch.setattr(runner, "run_genome", refuse_overwrite)
    record = runner._one(entry, str(output), asdict(config), "rscript", False)
    assert record["status"] == "error"
    assert "FileExistsError" in record["error"]
    assert protected.read_bytes() == b"preserve these exact bytes"
    assert list(output.iterdir()) == [protected]


def test_valid_resume_verifies_artifacts_and_skips_pipeline(
    runner, completed_run, entry, config, result, pipeline_mock
):
    assert runner._check_resume(completed_run, entry, config) == result
    record = runner._one(entry, str(completed_run), asdict(config), "rscript", True)
    assert record["status"] == "pass"
    assert record["reused"] is True
    assert not pipeline_mock


@pytest.mark.parametrize("changed", ["input", "settings"])
def test_resume_rejects_different_input_or_configuration(
    runner, completed_run, entry, config, changed
):
    if changed == "input":
        entry["genome_sha256"] = "0" * 64
    else:
        config = GenomePipelineConfig(**{**asdict(config), "min_identity": 75})
    with pytest.raises(ValueError, match="different inputs/settings"):
        runner._check_resume(completed_run, entry, config)


@pytest.mark.parametrize(
    ("filename", "field", "value"),
    [
        ("run.json", "status", "failed"),
        ("result.json", "status", "prediction_complete_r_failed"),
        ("result.json", "r_comparison_status", "running"),
        ("result.json", "r_comparison_status", "failed"),
    ],
)
def test_resume_rejects_incomplete_runs_even_when_r_matches_is_true(
    runner, completed_run, entry, config, filename, field, value
):
    path = completed_run / filename
    saved = json.loads(path.read_text(encoding="utf-8"))
    saved[field] = value
    _write_json(path, saved)
    _refresh_artifacts(completed_run)
    with pytest.raises(ValueError, match="incomplete run"):
        runner._check_resume(completed_run, entry, config)


def test_completed_numerical_mismatch_is_resumed_as_mismatch(
    runner, completed_run, entry, config, pipeline_mock
):
    path = completed_run / "result.json"
    saved = json.loads(path.read_text(encoding="utf-8"))
    saved["r_matches"] = False
    _write_json(path, saved)
    _refresh_artifacts(completed_run)
    record = runner._one(entry, str(completed_run), asdict(config), "rscript", True)
    assert record["status"] == "mismatch"
    assert record["reused"] is True
    assert not pipeline_mock


@pytest.mark.parametrize("change_size", [False, True])
def test_resume_rejects_changed_artifact_bytes(runner, completed_run, entry, config, change_size):
    path = completed_run / "predicted_cds.fna"
    original = path.read_bytes()
    replacement = original + b"\n" if change_size else original.replace(b"AAA", b"CCC")
    assert replacement != original
    path.write_bytes(replacement)
    with pytest.raises(ValueError, match="checksum mismatch"):
        runner._check_resume(completed_run, entry, config)


@pytest.mark.parametrize("change", ["extra", "missing"])
def test_resume_rejects_changed_artifact_inventory(runner, completed_run, entry, config, change):
    if change == "extra":
        (completed_run / "unexpected.txt").write_text("extra", encoding="utf-8")
    else:
        (completed_run / "r_comparison.log").unlink()
    with pytest.raises(ValueError, match="inventory changed"):
        runner._check_resume(completed_run, entry, config)


@pytest.mark.parametrize("filename", ["r_result.json", "r_comparison.json", "r_parameters.json"])
def test_resume_requires_original_r_evidence_even_if_manifest_omits_it(
    runner, completed_run, entry, config, filename
):
    (completed_run / filename).unlink()
    _refresh_artifacts(completed_run)
    with pytest.raises(ValueError, match="Missing comparison artifacts"):
        runner._check_resume(completed_run, entry, config)


@pytest.mark.parametrize("absolute", [False, True])
def test_resume_rejects_artifact_paths_outside_run_directory(
    runner, completed_run, entry, config, tmp_path, absolute
):
    outside = tmp_path / "outside.txt"
    outside.write_text("outside bytes", encoding="utf-8")
    manifest_path = completed_run / "run.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    unsafe = str(outside.resolve()) if absolute else "../outside.txt"
    manifest["artifacts"][unsafe] = {
        "bytes": outside.stat().st_size,
        "sha256": sha256_file(outside),
    }
    _write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="inventory changed|Unsafe artifact path"):
        runner._check_resume(completed_run, entry, config)
    assert outside.read_text(encoding="utf-8") == "outside bytes"


def test_resume_rejects_symlink_artifact_resolving_outside_run(
    runner, completed_run, entry, config, tmp_path
):
    outside = tmp_path / "outside.txt"
    outside.write_text("external evidence", encoding="utf-8")
    link = completed_run / "linked-evidence.txt"
    try:
        link.symlink_to(outside)
    except (OSError, NotImplementedError):
        pytest.skip("Creating file symlinks is unavailable on this host")
    _refresh_artifacts(completed_run)
    with pytest.raises(ValueError, match="Unsafe artifact path"):
        runner._check_resume(completed_run, entry, config)


@pytest.mark.parametrize("kind", ["empty", "duplicate", "traversal", "absolute"])
def test_invalid_manifest_is_rejected_before_dependencies_or_output(
    runner, monkeypatch, tmp_path, entry, kind
):
    entries = [deepcopy(entry)]
    if kind == "empty":
        entries = []
    elif kind == "duplicate":
        entries.append(deepcopy(entry))
    else:
        entries[0]["accession"] = (
            "../GCF_000000001.1" if kind == "traversal" else str(tmp_path / "escaped-output")
        )
    manifest = tmp_path / "manifest.json"
    _write_json(manifest, {"entries": entries})

    def unexpected(*args, **kwargs):
        pytest.fail("Invalid inventory must fail before external dependency checks")

    monkeypatch.setattr(runner, "_resolved_program", unexpected)
    output = tmp_path / "runs"
    with pytest.raises(ValueError, match="unique assembly accessions|Invalid accession"):
        runner.run_suite(manifest, output, blastp="blastp", rscript="Rscript")
    assert not output.exists()


@pytest.mark.parametrize("workers", [0, -1, True, 1.5, "2", None])
def test_invalid_worker_count_fails_before_any_output(runner, tmp_path, workers):
    output = tmp_path / "runs"
    with pytest.raises(ValueError, match="workers"):
        runner.run_suite(
            tmp_path / "not-read.json", output, blastp="blastp", rscript="Rscript", workers=workers
        )
    assert not output.exists()


@pytest.fixture
def suite_boundary(runner, monkeypatch, tmp_path, entry):
    manifest = tmp_path / "input-manifest.json"
    _write_json(manifest, {"entries": [entry]})
    output = tmp_path / "suite"
    output.mkdir()
    monkeypatch.setattr(runner, "_resolved_program", lambda value: str(Path(value).resolve()))
    monkeypatch.setattr(runner, "doctor", lambda executable: {"ready": True})

    def unexpected(*args, **kwargs):
        pytest.fail("Invalid resume provenance must fail before launching any workers")

    monkeypatch.setattr(runner, "ProcessPoolExecutor", unexpected)
    return manifest, output


def test_resume_with_different_r_executable_preserves_prior_suite_and_inventory(
    runner, tmp_path, suite_boundary
):
    manifest, output = suite_boundary
    previous = {
        "schema_version": 1,
        "status": "complete",
        "rscript": str(tmp_path / "original-Rscript.exe"),
        "finished": 40,
    }
    _write_json(output / "suite.json", previous)
    (output / "input_manifest.json").write_bytes(b"previous inventory bytes")
    before = {path.name: path.read_bytes() for path in output.iterdir()}
    with pytest.raises(ValueError, match="different R executable"):
        runner.run_suite(
            manifest,
            output,
            blastp=str(tmp_path / "blastp.exe"),
            rscript=str(tmp_path / "different-Rscript.exe"),
            resume=True,
        )
    assert {path.name: path.read_bytes() for path in output.iterdir()} == before


@pytest.mark.parametrize("has_existing_run", [False, True])
def test_resume_without_suite_provenance_fails_without_writing_any_files(
    runner, tmp_path, suite_boundary, has_existing_run
):
    manifest, output = suite_boundary
    if has_existing_run:
        existing = output / "GCF_000000001.1"
        existing.mkdir()
        (existing / "run.json").write_bytes(b"existing run evidence")
    before = {
        str(path.relative_to(output)): path.read_bytes()
        for path in output.rglob("*")
        if path.is_file()
    }
    with pytest.raises(ValueError, match="without suite.json"):
        runner.run_suite(
            manifest,
            output,
            blastp=str(tmp_path / "blastp.exe"),
            rscript=str(tmp_path / "Rscript.exe"),
            resume=True,
        )
    after = {
        str(path.relative_to(output)): path.read_bytes()
        for path in output.rglob("*")
        if path.is_file()
    }
    assert after == before
    assert not (output / "suite.json").exists()
    assert not (output / "input_manifest.json").exists()


@pytest.mark.parametrize("previous", [{}, {"rscript": "Rscript.exe"}, {"rscript": None}, []])
def test_resume_rejects_missing_or_ambiguous_recorded_r_executable_without_overwriting(
    runner, tmp_path, suite_boundary, previous
):
    manifest, output = suite_boundary
    _write_json(output / "suite.json", previous)
    before = (output / "suite.json").read_bytes()
    with pytest.raises(ValueError, match="absolute recorded R executable"):
        runner.run_suite(
            manifest,
            output,
            blastp=str(tmp_path / "blastp.exe"),
            rscript=str(tmp_path / "Rscript.exe"),
            resume=True,
        )
    assert (output / "suite.json").read_bytes() == before
    assert sorted(path.name for path in output.iterdir()) == ["suite.json"]


def test_resume_preserves_original_r_path_when_requested_path_resolves_to_same_executable(
    runner, tmp_path, suite_boundary
):
    _, output = suite_boundary
    original = str(tmp_path / "r-install" / ".." / "Rscript.exe")
    _write_json(output / "suite.json", {"rscript": original})
    before = (output / "suite.json").read_bytes()
    recorded = runner._check_suite_resume(output, str(tmp_path / "Rscript.exe"))
    assert recorded == original
    assert (output / "suite.json").read_bytes() == before
