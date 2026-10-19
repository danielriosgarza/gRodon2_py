"""Synthetic audit fixtures; no genome computation or R installation required."""

from __future__ import annotations

import csv
import importlib.util
import io
import json
import re
from copy import deepcopy
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "report_folder_validation.py"
_SPEC = importlib.util.spec_from_file_location("folder_validation_report", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
reporter = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(reporter)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def refresh_artifacts(case):
    run = json.loads(case["run"].read_text())
    run["artifacts"] = {
        path.name: {"sha256": reporter._sha256(path), "bytes": path.stat().st_size}
        for path in case["run_dir"].iterdir()
        if path.name != "run.json"
    }
    write_json(case["run"], run)


@pytest.fixture
def case(tmp_path):
    inputs = tmp_path / "input files"
    inputs.mkdir()
    genome = inputs / "GCF_000000001.1_genomic.fna"
    genome.write_text(">synthetic\nATGAAATAA\n", encoding="ascii")
    manifest = tmp_path / "inventory.json"
    entry = {
        "accession": "GCF_000000001.1",
        "source_file": str(genome),
        "source_kind": "genomic",
        "source_sha256": reporter._sha256(genome),
        "genome_path": str(genome),
        "genome_sha256": reporter._sha256(genome),
        "organism": {"name": "Synthetic bacterium", "tax_id": 0, "strain": "fixture"},
        "assembly": {"accession": "GCF_000000001.1"},
        "status": "ready",
        "error": None,
    }
    write_json(manifest, {"schema_version": 1, "total": 1, "entries": [entry]})
    runs = tmp_path / "runs"
    run_dir = runs / entry["accession"]
    run_dir.mkdir(parents=True)
    cds, ids = run_dir / "predicted_cds.fna", run_dir / "highly_expressed_ids.txt"
    cds.write_text(">g1\nATGAAATAA\n>g2\nATGCCCTAA\n>g3\nATGGGGTAA\n", encoding="ascii")
    ids.write_text("g1\ng2\n", encoding="ascii")
    settings = {
        "mode": "full",
        "training_set": "madin",
        "temperature": None,
        "genetic_code": 11,
        "confidence": 0.95,
        "threads": 1,
        "reference": "bacteria-swissprot-v1",
    }
    stats = {
        "CUBHE": 0.8,
        "GC": 0.5,
        "GCdiv": 0.0,
        "ConsistencyHE": 0.2,
        "CUB": 0.4,
        "CPB": -0.3,
        "FilteredSequences": 0,
        "nHE": 2,
        "dCUB": -0.5,
        "d": 2.0,
        "LowerCI": 1.0,
        "UpperCI": 3.0,
    }
    python = {
        "status": "complete",
        "r_comparison_status": "complete",
        "settings": settings,
        "statistics": stats,
        "input": str(genome),
        "input_sha256": entry["genome_sha256"],
        "predicted_genes": 3,
        "ribosomal_genes": 2,
        "retained_genes": 3,
        "retained_ribosomal_genes": 2,
        "warnings": ["Synthetic filter warning"],
        "model": "gRodon_model_base_madin",
        "r_matches": True,
    }
    rsettings = {
        key: value for key, value in settings.items() if key not in {"threads", "reference"}
    }
    rsettings.update(temperature="none", fragments=False)
    original = {
        "settings": rsettings,
        "statistics": stats,
        "versions": {"R": "R version fixture", "gRodon": "2.7.3"},
        "input_genes": 3,
        "highly_expressed_before_filtering": 2,
        "warnings": ["Synthetic R warning"],
    }
    parameters = dict(rsettings, temperature=None, cds_path=str(cds), ids_path=str(ids))
    write_json(run_dir / "result.json", python)
    write_json(run_dir / "r_result.json", original)
    write_json(run_dir / "r_parameters.json", parameters)
    write_json(run_dir / "r_comparison.json", {"matches": True})
    (run_dir / "r_comparison.log").write_text("Synthetic R log\n", encoding="ascii")
    run = {
        "status": "complete",
        "input": {"path": str(genome), "sha256": entry["genome_sha256"]},
        "settings": settings,
        "software": {"python": "3.fixture", "pyrodigal": "3.7.1"},
        "dependencies": {"blastp": {"version": "2.17.0"}},
        "reference": {"name": "bacteria-swissprot-v1", "sha256": "fixture"},
    }
    write_json(run_dir / "run.json", run)
    value = {
        "manifest": manifest,
        "runs": runs,
        "entry": entry,
        "genome": genome,
        "run_dir": run_dir,
        "run": run_dir / "run.json",
        "python": run_dir / "result.json",
        "r": run_dir / "r_result.json",
    }
    refresh_artifacts(value)
    return value


def summarize(case, **kwargs):
    return reporter.summarize(case["manifest"], case["runs"], expected_count=1, **kwargs)


def alter(path, callback):
    value = json.loads(path.read_text())
    callback(value)
    write_json(path, value)


def test_complete_audit_checks_all_full_statistics_and_provenance(case):
    report = summarize(case)
    assert report["all_pass"] and report["all_bitwise_equal"]
    assert report["counts"] == {"pass": 1, "mismatch": 0, "error": 0, "not_run": 0}
    row = report["genomes"][0]
    assert row["artifact_verification"]["verified"]
    assert len(row["metric_comparisons"]) == 12
    assert all(check["matches"] for check in row["input_checks"] + row["agreement_checks"])
    assert row["doubling_time"]["python_minutes"] == 120
    assert row["counts"]["prepared_cds"] == 3
    assert row["counts"]["r_nHE"] == 2
    assert report["observed_versions"]["r"]["gRodon"] == ["2.7.3"]
    assert report["blastp_versions"] == ["2.17.0"]
    assert {value["metric"] for value in report["metric_summary"]} == set(reporter.FULL_METRICS)


def test_rounding_is_numerical_match_but_not_exact(case):
    alter(case["r"], lambda value: value["statistics"].update(d=2.0 + 1e-11))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["all_pass"] and not report["all_bitwise_equal"]
    metric = next(item for item in report["metric_summary"] if item["metric"] == "d")
    assert metric["matching_genomes"] == 1
    assert metric["bitwise_equal_genomes"] == 0
    assert metric["max_absolute_difference"] == pytest.approx(1e-11, abs=1e-15)


def test_true_mismatch_is_recomputed_despite_stored_matches_true(case):
    alter(case["r"], lambda value: value["statistics"].update(d=2.1))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["counts"]["mismatch"] == 1
    assert not report["all_pass"]
    assert "d: Outside numerical tolerance" in report["genomes"][0]["mismatches"]


def test_symmetric_max_tolerance_not_numpy_additive_tolerance(case):
    # At scale 0.01 the relative and absolute tolerances are both about 1e-12.
    # A difference of 1.5e-12 would pass atol+rtol*abs(r), but must fail max().
    alter(case["python"], lambda value: value["statistics"].update(GC=0.01))
    alter(case["r"], lambda value: value["statistics"].update(GC=0.01 + 1.5e-12))
    refresh_artifacts(case)
    assert summarize(case)["counts"]["mismatch"] == 1


def test_negative_zero_is_numerically_equal_but_not_bitwise_equal(case):
    alter(case["r"], lambda value: value["statistics"].update(GCdiv=-0.0))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["all_pass"]
    assert not report["all_bitwise_equal"]


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -float("inf"), True, "2.0"])
def test_matching_nonfinite_or_nonnumeric_values_never_pass(case, value):
    for path in (case["r"], case["python"]):
        alter(path, lambda data: data["statistics"].update(d=value))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["counts"]["mismatch"] == 1
    assert not report["genomes"][0]["numerically_matches"]
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("key", ["CUBHE", "LowerCI", "FilteredSequences"])
def test_missing_metric_fails_even_when_missing_from_both(case, key):
    for path in (case["r"], case["python"]):
        alter(path, lambda data: data["statistics"].pop(key))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["counts"]["mismatch"] == 1
    assert not report["genomes"][0]["full_metric_set"]


def test_unknown_extra_statistic_is_not_silently_ignored(case):
    alter(case["r"], lambda value: value["statistics"].update(Extra=2.0))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["counts"]["mismatch"] == 1
    assert not report["genomes"][0]["numeric_key_sets_equal"]
    assert any(row["metric"] == "Extra" for row in report["metric_summary"])


def test_missing_input_genome_is_error(case):
    case["genome"].unlink()
    report = summarize(case)
    assert report["counts"]["error"] == 1
    assert any("genome input" in error for error in report["genomes"][0]["errors"])


def test_missing_run_is_not_run_and_is_included_in_markdown(case):
    case["run_dir"].rename(case["runs"] / "other")
    report = summarize(case)
    assert report["counts"]["not_run"] == 1
    text = reporter.render_markdown(report)
    assert "GCF_000000001.1" in text and "not_run" in text
    assert "INCOMPLETE OR FAILED" in text


def test_input_fetch_error_is_not_hidden_as_not_run(case):
    case["run_dir"].rename(case["runs"] / "other")
    alter(
        case["manifest"],
        lambda value: value["entries"][0].update(
            status="error", genome_path=None, error="download unavailable"
        ),
    )
    report = summarize(case)
    assert report["counts"]["error"] == 1
    assert "download unavailable" in reporter.render_markdown(report)


def test_r_process_failure_is_an_error_with_diagnostics(case):
    alter(case["run"], lambda value: value.update(status="failed", error={"message": "R exited 1"}))
    alter(
        case["python"],
        lambda value: value.update(
            status="prediction_complete_r_failed",
            r_comparison_status="failed",
            r_comparison_error="R exited 1",
        ),
    )
    case["r"].write_text("", encoding="ascii")
    refresh_artifacts(case)
    report = summarize(case)
    assert report["counts"]["error"] == 1
    assert "R exited 1" in reporter.render_markdown(report)


def test_input_hash_mismatch_invalidates_matching_statistics(case):
    case["genome"].write_text(">changed\nATGCCCTAA\n", encoding="ascii")
    report = summarize(case)
    assert report["counts"]["error"] == 1
    assert report["genomes"][0]["numerically_matches"]
    assert not report["all_pass"]


def test_run_hash_cannot_point_to_another_input(case):
    alter(case["run"], lambda value: value["input"].update(sha256="0" * 64))
    report = summarize(case)
    assert report["counts"]["error"] == 1
    assert any("Run genome SHA256" in error for error in report["genomes"][0]["errors"])


def test_modified_artifact_fails_default_audit_and_opt_out_is_explicit(case):
    (case["run_dir"] / "r_comparison.log").write_text("changed", encoding="ascii")
    report = summarize(case)
    assert report["counts"]["error"] == 1
    assert "Artifact SHA256 or size mismatch" in reporter.render_markdown(report)
    unchecked = summarize(case, verify_artifacts=False)
    assert unchecked["all_pass"]
    assert "explicitly disabled" in reporter.render_markdown(unchecked)


@pytest.mark.parametrize("action", ["extra", "missing", "directory"])
def test_exact_artifact_inventory_is_required(case, action):
    if action == "extra":
        (case["run_dir"] / "extra.txt").write_text("extra", encoding="ascii")
    elif action == "missing":
        (case["run_dir"] / "r_comparison.log").unlink()
    else:
        (case["run_dir"] / "unrecorded_directory").mkdir()
    assert summarize(case)["counts"]["error"] == 1


def test_path_traversal_in_artifact_manifest_is_rejected(case):
    alter(
        case["run"],
        lambda value: value["artifacts"].update(
            {"../outside.txt": {"sha256": "0" * 64, "bytes": 0}}
        ),
    )
    report = summarize(case)
    assert report["counts"]["error"] == 1
    assert "Unsafe or malformed artifact entry" in reporter.render_markdown(report)


def test_r_must_read_same_prepared_inputs(case):
    alter(
        case["run_dir"] / "r_parameters.json",
        lambda value: value.update(cds_path=str(case["genome"])),
    )
    refresh_artifacts(case)
    assert summarize(case)["counts"]["error"] == 1


def test_same_gene_counts_and_settings_are_required(case):
    alter(case["r"], lambda value: value.update(input_genes=4))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["counts"]["mismatch"] == 1
    assert "Prepared/R CDS count" in reporter.render_markdown(report)


def test_incomplete_inventory_cannot_claim_all_40_passed(case):
    report = reporter.summarize(case["manifest"], case["runs"])
    assert report["counts"]["pass"] == 1
    assert not report["complete_inventory"] and not report["all_pass"]
    assert "Expected 40 assemblies, found 1" in reporter.render_markdown(report)


def test_all_inventory_rows_included_even_with_no_runs(case):
    entries = []
    for index in range(40):
        entry = deepcopy(case["entry"])
        entry["accession"] = f"GCF_{index + 100:09d}.1"
        entry["assembly"]["accession"] = entry["accession"]
        entry["source_kind"] = "genomic" if index < 3 else "protein"
        entries.append(entry)
    write_json(case["manifest"], {"entries": entries})
    report = reporter.summarize(case["manifest"], case["runs"])
    assert report["complete_inventory"]
    assert report["counts"]["not_run"] == 40
    text = reporter.render_markdown(report)
    assert "3 genomic DNA files and 37 protein-only files" in text
    assert "not back-translated" in text
    for entry in entries:
        assert f"| {entry['accession']} |" in text


@pytest.mark.parametrize("accession", ["..", "../outside", "A/B", "A\\B", "C:escape"])
def test_unsafe_accession_rejected_before_reading_run_paths(case, accession):
    alter(case["manifest"], lambda value: value["entries"][0].update(accession=accession))
    with pytest.raises(ValueError, match="Unsafe input accession"):
        summarize(case)


def test_duplicate_inventory_entries_are_not_silently_collapsed(case):
    alter(case["manifest"], lambda value: value["entries"].append(value["entries"][0]))
    with pytest.raises(ValueError, match="Duplicate input accession"):
        summarize(case)


def test_cli_writes_both_reports_and_exit_status(case, tmp_path):
    output_json, output_md = tmp_path / "report.json", tmp_path / "report.md"
    arguments = [
        "--manifest",
        str(case["manifest"]),
        "--runs-dir",
        str(case["runs"]),
        "--output-json",
        str(output_json),
        "--output-md",
        str(output_md),
        "--expected-count",
        "1",
    ]
    assert reporter.main(arguments) == 0
    assert json.loads(output_json.read_text())["all_pass"]
    assert "1 of 1 assemblies passed" in output_md.read_text(encoding="utf-8")
    case["genome"].unlink()
    assert reporter.main(arguments) == 1
    assert not json.loads(output_json.read_text())["all_pass"]


def test_cli_csv_preserves_raw_numbers_and_all_assembly_rows(case, tmp_path):
    output_csv = tmp_path / "report.csv"
    status = reporter.main(
        [
            "--manifest",
            str(case["manifest"]),
            "--runs-dir",
            str(case["runs"]),
            "--output-json",
            str(tmp_path / "report.json"),
            "--output-md",
            str(tmp_path / "report.md"),
            "--output-csv",
            str(output_csv),
            "--expected-count",
            "1",
        ]
    )
    assert status == 0
    rows = list(csv.DictReader(io.StringIO(output_csv.read_text())))
    assert len(rows) == 1
    assert rows[0]["accession"] == case["entry"]["accession"]
    assert rows[0]["python_doubling_hours"] == "2.0"
    assert rows[0]["r_upper_ci_hours"] == "3.0"
    assert rows[0]["predicted_genes"] == "3"
    assert rows[0]["matches"] == "True"
    assert rows[0]["bitwise_equal_statistics"] == "12"


def test_single_r_warning_is_not_split_into_characters(case):
    alter(case["r"], lambda value: value.update(warnings="One R warning"))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["all_pass"]
    assert report["genomes"][0]["warnings"]["r"] == ["One R warning"]
    assert "r: One R warning" in reporter.render_markdown(report)


@pytest.mark.parametrize("value", [None, 42, ["valid", 42]])
def test_malformed_warning_data_is_reported_as_error(case, value):
    alter(case["r"], lambda data: data.update(warnings=value))
    refresh_artifacts(case)
    report = summarize(case)
    assert report["counts"]["error"] == 1
    assert "warnings must be" in reporter.render_markdown(report)


def test_cli_cannot_overwrite_audited_run_artifacts(case, tmp_path):
    original = case["python"].read_bytes()
    with pytest.raises(SystemExit) as error:
        reporter.main(
            [
                "--manifest",
                str(case["manifest"]),
                "--runs-dir",
                str(case["runs"]),
                "--output-json",
                str(case["python"]),
                "--output-md",
                str(tmp_path / "report.md"),
            ]
        )
    assert error.value.code == 2
    assert case["python"].read_bytes() == original


def test_annotation_header_differences_do_not_claim_different_genomes(case):
    alter(
        case["manifest"],
        lambda value: value["entries"][0].update(
            source_kind="protein",
            download={
                "source_protein_md5_matches": False,
                "source_protein_content_matches": True,
                "source_protein_order_matches": True,
                "source_protein_headers_match": False,
                "source_protein_difference_kind": "annotation_headers_differ_protein_sequences_identical",
            },
        ),
    )
    report = summarize(case)
    assert report["all_pass"]
    assert report["protein_input_audit"]["different_compressed_archives"] == 1
    assert "protein IDs, sequences, and order were identical" in reporter.render_markdown(report)
    assert "not evidence of different assembly genomes" in reporter.render_markdown(report)


def test_public_reproduction_commands_are_portable_and_preserve_archived_outputs(case):
    text = reporter.render_markdown(summarize(case))
    commands = "\n".join(re.findall(r"```powershell\n(.*?)\n```", text, flags=re.DOTALL))
    assert "python -m venv .venv" in commands
    assert "-m grodon setup --tools-dir .tools/blast-2.17.0 --version 2.17.0" in commands
    assert (
        "scripts/fetch_benchmark_genomes.py --manifest validation/test_genomes_inputs_manifest.json"
        in commands
    )
    assert "--input-dir validation/test_genomes_inputs/public-genomes" in commands
    assert "--manifest validation/reproduced_inputs.json" in commands
    assert '--rscript "C:\\path\\to\\R\\bin\\Rscript.exe"' in commands
    assert "--output-json validation/reproduced_comparison.json" in commands
    assert "--output-csv validation/reproduced_comparison.csv" in commands
    assert "--output-md validation/reproduced_comparison.md" in commands
    assert "--output validation/reproduced_growth_reference.json" in commands
    assert "C:\\Users\\" not in commands and ".reference" not in commands
    assert "--output-json validation/test_genomes_comparison.json" not in commands
    assert "--output tests/data/folder_growth_reference.json" not in commands
    assert "Historical absolute paths are provenance" in text
    assert "forty DNA inputs" in text
