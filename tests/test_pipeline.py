"""Exercise orchestration using artificial CDS and the actual prediction model."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pytest

import grodon.pipeline as pipeline
import grodon.r_comparison as r_comparison
from grodon.pipeline import GenomePipelineConfig, run_genome, run_genomes
from grodon.prediction import PredictionResult


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def fake_pipeline(tmp_path, monkeypatch):
    fixture = json.loads((Path(__file__).parent / "fixtures" / "codon_reference.json").read_text())
    sequences = fixture["cases"][0]["genes"]
    reference = tmp_path / "reference.faa"
    reference.write_text(">synthetic_reference ribosomal test protein\nMAKRTA\n")
    source = tmp_path / "assembly.fna"
    source.write_text(">synthetic_contig\nACGTACGTACGT\n")
    state: dict[str, Any] = {
        "genes": list(sequences),
        "selected_ids": [f"gene_{i:03d}" for i in range(12)],
        "fail_inputs": set(),
        "calls": [],
        "annotations": [],
        "source": source,
        "reference": reference,
    }

    def doctor(executable):
        return {
            "ready": True,
            "pyrodigal": {"available": True, "version": "3.7.1", "error": None},
            "blastp": {
                "available": True,
                "version": "2.17.0+",
                "path": str(tmp_path / "mock-blastp"),
                "error": None,
            },
        }

    def call_genes(genome, output, **kwargs):
        genome, output = Path(genome), Path(output)
        state["calls"].append((genome, kwargs))
        if genome.name in state["fail_inputs"]:
            (output / "partial_gene_call.log").write_text("synthetic gene-caller failure\n")
            raise RuntimeError("synthetic gene-caller failure")
        cds, proteins = output / "predicted_cds.fna", output / "predicted_proteins.faa"
        records = state["genes"]
        cds.write_text(
            "".join(f">gene_{i:03d}\n{sequence}\n" for i, sequence in enumerate(records)),
            newline="\n",
        )
        proteins.write_text(
            "".join(f">gene_{i:03d}\nMAKRTA\n" for i in range(len(records))), newline="\n"
        )
        return {
            "cds_path": cds,
            "proteins_path": proteins,
            "input_sha256": _sha(genome),
            "total_genome_bp": sum(map(len, records)),
            "counts_by_contig": {"synthetic_contig": len(records)},
        }

    def annotate(proteins, reference_path, output, **kwargs):
        output = Path(output)
        state["annotations"].append((Path(reference_path), kwargs))
        selected = state["selected_ids"].copy()
        ids = output / "highly_expressed_ids.txt"
        ids.write_text("".join(identifier + "\n" for identifier in selected), newline="\n")
        (output / "ribosomal_hits.tsv").write_text("synthetic test hit table\n", newline="\n")
        (output / "blastp.log").write_text("synthetic BLASTP stub\n", newline="\n")
        return {"selected_ids": selected, "ids_path": ids}

    monkeypatch.setattr(pipeline, "doctor", doctor)
    monkeypatch.setattr(pipeline, "call_genes", call_genes)
    monkeypatch.setattr(pipeline, "annotate_ribosomal", annotate)
    state["config"] = GenomePipelineConfig(reference=reference)
    return state


def test_complete_pipeline_predicts_and_hashes_artifacts(tmp_path, fake_pipeline):
    state = fake_pipeline
    result = run_genome(state["source"], tmp_path / "run", config=state["config"])
    summary = result.to_dict()
    assert summary["status"] == "complete"
    assert summary["predicted_genes"] == 24
    assert summary["ribosomal_genes"] == summary["retained_ribosomal_genes"] == 12
    assert summary["retained_genes"] == 24
    assert result.prediction.model_name == "gRodon_model_base_madin"
    assert math.isfinite(result.prediction.d) and result.prediction.d > 0
    assert summary["doubling_time_minutes"] == pytest.approx(result.prediction.d * 60)
    assert summary["growth_rate_per_hour"] == pytest.approx(math.log(2) / result.prediction.d)
    assert json.loads((result.output_dir / "result.json").read_text()) == summary
    manifest = json.loads((result.output_dir / "run.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["input"]["sha256"] == _sha(state["source"])
    for name, item in manifest["artifacts"].items():
        path = result.output_dir / name
        assert item["sha256"] == _sha(path)
        assert item["bytes"] == path.stat().st_size
    with (result.output_dir / "result.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1 and rows[0]["status"] == "complete"
    assert float(rows[0]["doubling_time_hours"]) == result.prediction.d
    assert state["annotations"][0][0] == result.output_dir / "reference.faa.gz"
    assert (
        hashlib.sha256(
            gzip.decompress((result.output_dir / "reference.faa.gz").read_bytes())
        ).hexdigest()
        == manifest["reference"]["sha256"]
    )


def test_repeat_runs_have_same_scientific_fingerprint(tmp_path, fake_pipeline):
    state = fake_pipeline
    first = run_genome(state["source"], tmp_path / "first", config=state["config"])
    second = run_genome(state["source"], tmp_path / "second", config=state["config"])
    assert first.output_dir != second.output_dir
    assert first.summary["scientific_fingerprint"] == second.summary["scientific_fingerprint"]
    assert first.summary["statistics"] == second.summary["statistics"]
    assert (first.output_dir / "reference.faa.gz").read_bytes() == (
        second.output_dir / "reference.faa.gz"
    ).read_bytes()
    # Changed reference bytes affect the scientific identity, even when mocked
    # homology happens to select the same genes and the point estimate is equal.
    state["reference"].write_text(">synthetic_reference ribosomal test protein\nMAKRTAA\n")
    changed = run_genome(state["source"], tmp_path / "changed", config=state["config"])
    assert changed.summary["scientific_fingerprint"] != first.summary["scientific_fingerprint"]


@pytest.mark.parametrize("existing", ["directory", "file"])
def test_existing_outputs_are_never_reused(tmp_path, fake_pipeline, existing):
    destination = tmp_path / "existing"
    if existing == "directory":
        destination.mkdir()
        sentinel = destination / "keep.txt"
    else:
        sentinel = destination
    sentinel.write_bytes(b"preserve exactly")
    with pytest.raises(FileExistsError, match="Output already exists"):
        run_genome(fake_pipeline["source"], destination, config=fake_pipeline["config"])
    assert sentinel.read_bytes() == b"preserve exactly"
    assert fake_pipeline["calls"] == []


def test_dependency_preflight_does_not_create_outputs(tmp_path, fake_pipeline, monkeypatch):
    monkeypatch.setattr(
        pipeline,
        "doctor",
        lambda executable: {
            "ready": False,
            "pyrodigal": {"available": False, "error": "missing Pyrodigal"},
            "blastp": {"available": True, "error": None},
        },
    )
    output = tmp_path / "new" / "output"
    with pytest.raises(RuntimeError, match="dependencies unavailable.*Pyrodigal"):
        run_genome(fake_pipeline["source"], output, config=fake_pipeline["config"])
    assert not output.parent.exists()
    assert fake_pipeline["calls"] == []


def test_reference_and_input_preflight_do_not_create_outputs(tmp_path, fake_pipeline):
    output = tmp_path / "output"
    with pytest.raises(FileNotFoundError):
        run_genome(
            fake_pipeline["source"],
            output,
            config=GenomePipelineConfig(reference=tmp_path / "absent.faa"),
        )
    assert not output.exists()
    with pytest.raises(FileNotFoundError):
        run_genome(tmp_path / "absent.fna", output, config=fake_pipeline["config"])
    assert not output.exists()


def test_failure_manifest_keeps_intermediates_and_exception(tmp_path, fake_pipeline):
    state = fake_pipeline
    state["fail_inputs"].add(state["source"].name)
    output = tmp_path / "failed"
    with pytest.raises(RuntimeError, match="synthetic gene-caller failure"):
        run_genome(state["source"], output, config=state["config"])
    manifest = json.loads((output / "run.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["error"] == {"type": "RuntimeError", "message": "synthetic gene-caller failure"}
    assert "partial_gene_call.log" in manifest["artifacts"]
    assert (output / "reference.faa.gz").is_file()
    assert not (output / "result.json").exists()


@pytest.mark.parametrize("short_count,passes", [(2, True), (3, False)])
def test_minimum_ribosomal_gate_is_applied_after_sequence_filtering(
    tmp_path, fake_pipeline, short_count, passes
):
    state = fake_pipeline
    state["genes"][:short_count] = ["AAA" * 40] * short_count
    output = tmp_path / "output"
    if passes:
        result = run_genome(state["source"], output, config=state["config"])
        assert result.summary["ribosomal_genes"] == 12
        assert result.summary["retained_ribosomal_genes"] == 10
        assert result.summary["retained_genes"] == 22
    else:
        with pytest.raises(ValueError, match="Only 9 ribosomal genes remain after filtering"):
            run_genome(state["source"], output, config=state["config"])
        assert json.loads((output / "run.json").read_text())["status"] == "failed"
        assert not (output / "result.json").exists()


def test_no_homology_hits_fail_with_retained_evidence(tmp_path, fake_pipeline):
    state = fake_pipeline
    state["selected_ids"] = []
    with pytest.raises(ValueError, match="No ribosomal proteins"):
        run_genome(state["source"], tmp_path / "output", config=state["config"])
    assert (tmp_path / "output" / "ribosomal_hits.tsv").is_file()
    assert json.loads((tmp_path / "output" / "run.json").read_text())["status"] == "failed"


@pytest.mark.parametrize("doubling", [float("nan"), float("inf"), 0, -1])
def test_nonfinite_or_nonpositive_prediction_fails(tmp_path, fake_pipeline, monkeypatch, doubling):
    prediction = PredictionResult(
        {"d": doubling, "LowerCI": 1, "UpperCI": 2}, "synthetic", "full", "madin", "all"
    )
    monkeypatch.setattr(pipeline, "predict_growth", lambda *args, **kwargs: prediction)
    with pytest.raises(ValueError, match="finite positive doubling time"):
        run_genome(fake_pipeline["source"], tmp_path / "output", config=fake_pipeline["config"])
    manifest = json.loads((tmp_path / "output" / "run.json").read_text())
    assert manifest["status"] == "failed"
    assert not (tmp_path / "output" / "result.json").exists()


def test_batch_continue_isolates_failures_and_writes_every_row(tmp_path, fake_pipeline):
    state = fake_pipeline
    sources = []
    for name in ("good.fna", "bad.fna", "good second; name.fna"):
        source = tmp_path / name
        source.write_text(">synthetic\nACGTACGT\n")
        sources.append(source)
    state["fail_inputs"].add("bad.fna")
    result = run_genomes(
        sources, tmp_path / "batch", config=state["config"], continue_on_error=True
    )
    assert result["status"] == "completed_with_errors"
    assert (result["requested"], result["completed"], result["failed"]) == (3, 2, 1)
    assert [row["status"] for row in result["results"]] == ["complete", "failed", "complete"]
    assert [Path(row["output_dir"]).name for row in result["results"]] == [
        "0001_good",
        "0002_bad",
        "0003_good_second__name",
    ]
    with (tmp_path / "batch" / "summary.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 3
    assert "synthetic gene-caller failure" in rows[1]["error"]
    assert float(rows[0]["doubling_time_hours"]) > 0
    assert (
        json.loads((Path(result["results"][1]["output_dir"]) / "run.json").read_text())["status"]
        == "failed"
    )


def test_batch_fail_fast_saves_summary_before_raising(tmp_path, fake_pipeline):
    state = fake_pipeline
    bad, later = tmp_path / "bad.fna", tmp_path / "later.fna"
    for source in (bad, later):
        source.write_text(">synthetic\nACGTACGT\n")
    state["fail_inputs"].add(bad.name)
    with pytest.raises(RuntimeError, match="Batch stopped at bad.fna"):
        run_genomes([bad, later], tmp_path / "batch", config=state["config"])
    result = json.loads((tmp_path / "batch" / "summary.json").read_text())
    assert result["status"] == "failed"
    assert result["requested"] == 2 and result["failed"] == 1 and result["completed"] == 0
    assert len(result["results"]) == 1
    assert [source.name for source, _ in state["calls"]] == ["bad.fna"]


def test_batch_success_and_existing_output_protection(tmp_path, fake_pipeline):
    state = fake_pipeline
    destination = tmp_path / "batch"
    result = run_genomes([state["source"]], destination, config=state["config"])
    assert result["status"] == "complete" and result["failed"] == 0
    before = (destination / "summary.json").read_bytes()
    with pytest.raises(FileExistsError):
        run_genomes([state["source"]], destination, config=state["config"])
    assert (destination / "summary.json").read_bytes() == before


@pytest.mark.parametrize(
    "options",
    [
        {"genetic_code": 11.0},
        {"genetic_code": True},
        {"min_ribosomal": 0},
        {"threads": 1.2},
        {"temperature": float("nan")},
        {"confidence": 1},
        {"mode": "eukaryote"},
    ],
)
def test_invalid_pipeline_configuration_rejected(options):
    with pytest.raises((ValueError, TypeError)):
        GenomePipelineConfig(**options)


@pytest.mark.parametrize("config", [{}, False, "bad"])
def test_wrong_config_type_rejected_before_output(tmp_path, fake_pipeline, config):
    with pytest.raises(TypeError, match="GenomePipelineConfig"):
        run_genome(fake_pipeline["source"], tmp_path / "one", config=config)
    with pytest.raises(TypeError, match="GenomePipelineConfig"):
        run_genomes([fake_pipeline["source"]], tmp_path / "batch", config=config)
    assert not (tmp_path / "one").exists() and not (tmp_path / "batch").exists()


def test_r_comparison_confidence_is_validated_before_output(tmp_path, fake_pipeline):
    settings = GenomePipelineConfig(reference=fake_pipeline["reference"], confidence=0.9)
    with pytest.raises(ValueError, match="confidence=0.95"):
        run_genome(fake_pipeline["source"], tmp_path / "output", config=settings, rscript="unused")
    assert not (tmp_path / "output").exists()


def test_failed_r_comparison_preserves_completed_native_prediction(
    tmp_path, fake_pipeline, monkeypatch
):
    state = fake_pipeline
    output = tmp_path / "r_failed"

    def fail_comparison(cds_path, ids_path, destination, **kwargs):
        assert Path(cds_path).is_file() and Path(ids_path).is_file()
        pending = json.loads((destination / "result.json").read_text())
        assert pending["status"] == "prediction_complete_r_pending"
        assert pending["r_comparison_status"] == "running"
        assert pending["statistics"] == kwargs["python_statistics"]
        assert (destination / "result.csv").is_file()
        (destination / "r_run.log").write_text("synthetic R installation error\n")
        raise RuntimeError("synthetic R installation error")

    monkeypatch.setattr(r_comparison, "compare_with_r", fail_comparison)
    with pytest.raises(RuntimeError, match="synthetic R installation error"):
        run_genome(state["source"], output, config=state["config"], rscript="mock-Rscript")
    native = json.loads((output / "result.json").read_text())
    assert native["status"] == "prediction_complete_r_failed"
    assert native["r_comparison_status"] == "failed"
    assert native["r_comparison_error"] == "synthetic R installation error"
    assert math.isfinite(native["doubling_time_hours"]) and native["doubling_time_hours"] > 0
    assert native["statistics"]["d"] == native["doubling_time_hours"]
    assert native["retained_ribosomal_genes"] == 12
    with (output / "result.csv").open(newline="") as stream:
        row = next(csv.DictReader(stream))
    assert row["status"] == "prediction_complete_r_failed"
    assert float(row["doubling_time_hours"]) == native["doubling_time_hours"]
    manifest = json.loads((output / "run.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["error"]["message"] == "synthetic R installation error"
    for name in ("result.json", "result.csv", "r_run.log"):
        assert manifest["artifacts"][name]["sha256"] == _sha(output / name)


def test_successful_r_comparison_completes_run_without_changing_prediction(
    tmp_path, fake_pipeline, monkeypatch
):
    state = fake_pipeline
    native = run_genome(state["source"], tmp_path / "native", config=state["config"])
    calls = []

    def matching_comparison(cds_path, ids_path, destination, **kwargs):
        calls.append(kwargs)
        assert Path(cds_path).is_file() and Path(ids_path).is_file()
        pending = json.loads((destination / "result.json").read_text())
        assert pending["status"] == "prediction_complete_r_pending"
        assert kwargs["python_statistics"] == native.summary["statistics"]
        value = {"matches": True}
        (destination / "r_comparison.json").write_text(json.dumps(value) + "\n")
        return value

    monkeypatch.setattr(r_comparison, "compare_with_r", matching_comparison)
    result = run_genome(
        state["source"], tmp_path / "compared", config=state["config"], rscript="mock-Rscript"
    )
    assert len(calls) == 1
    assert calls[0]["rscript"] == "mock-Rscript"
    assert calls[0]["mode"] == "full" and calls[0]["genetic_code"] == 11
    assert result.summary["status"] == "complete"
    assert result.summary["r_comparison_status"] == "complete"
    assert result.summary["r_matches"] is True
    assert result.summary["statistics"] == native.summary["statistics"]
    assert result.summary["scientific_fingerprint"] == native.summary["scientific_fingerprint"]
    manifest = json.loads((result.output_dir / "run.json").read_text())
    assert manifest["status"] == "complete"
    assert "r_comparison.json" in manifest["artifacts"]
