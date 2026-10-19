"""Run-independent tests for the optional actual-R comparison boundary."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

import grodon.r_comparison as comparison


@pytest.fixture
def case(tmp_path, monkeypatch):
    source = tmp_path / "inputs with spaces"
    source.mkdir()
    cds = source / "coding sequences.fna"
    ids = source / "highly expressed ids.txt"
    rscript = source / "Rscript executable"
    cds.write_text(">gene_1 30S ribosomal protein\nATGAAATAA\n", encoding="utf-8")
    ids.write_text("gene_1\n", encoding="utf-8")
    rscript.write_text("fake executable", encoding="utf-8")
    monkeypatch.setattr(comparison.shutil, "which", lambda executable: None)
    return {
        "cds_path": cds,
        "ids_path": ids,
        "output_dir": tmp_path / "results with spaces",
        "rscript": rscript,
        "python_statistics": {"d": 2.0, "CUBHE": 0.7, "nHE": 12},
    }


@pytest.fixture
def r_process(monkeypatch):
    state = {
        "result": {
            "statistics": {"d": 2.0, "CUBHE": 0.7, "nHE": 12},
            "versions": {"R": "R version test", "gRodon": "2.7.3"},
            "session_info": ["Test R session"],
        },
        "status": 0,
        "stderr": "R diagnostic\n",
        "calls": [],
    }

    def run(command, **kwargs):
        state["calls"].append((command, kwargs))
        stdout = state.get("stdout", json.dumps(state["result"]))
        return subprocess.CompletedProcess(command, state["status"], stdout, state["stderr"])

    monkeypatch.setattr(comparison.subprocess, "run", run)
    return state


def test_actual_r_command_boundary_and_audit_artifacts(case, r_process):
    report = comparison.compare_with_r(**case)
    assert report["matches"] is True
    assert report["tolerance"] == {"relative": 1e-10, "absolute": 1e-12}
    assert len(report["comparison"]) == 3
    assert report["r_result"]["versions"]["gRodon"] == "2.7.3"
    for path in report["files"].values():
        assert Path(path).is_file()
    saved = json.loads(Path(report["files"]["comparison"]).read_text())
    assert saved == report
    assert "R diagnostic" in Path(report["files"]["log"]).read_text()
    parameters = json.loads(Path(report["files"]["parameters"]).read_text())
    assert parameters["cds_path"] == str(case["cds_path"])
    assert parameters["ids_path"] == str(case["ids_path"])
    assert parameters["fragments"] is False
    assert parameters["confidence"] == 0.95
    assert parameters["temperature"] is None
    command, kwargs = r_process["calls"][0]
    assert command[0] == str(case["rscript"])
    assert command[1] == "--vanilla"
    assert command[-1] == report["files"]["parameters"]
    assert kwargs["shell"] is False
    assert kwargs["timeout"] == 600


def test_numeric_mismatches_return_detailed_report_without_raising(case, r_process):
    r_process["result"]["statistics"]["d"] = 3.0
    report = comparison.compare_with_r(**case)
    assert report["matches"] is False
    difference = next(row for row in report["comparison"] if row["metric"] == "d")
    assert difference["python"] == 2.0
    assert difference["r"] == 3.0
    assert difference["absolute_difference"] == 1.0
    assert difference["relative_difference"] == pytest.approx(1 / 3)
    assert difference["reason"] == "outside numerical tolerance"


def test_all_statistic_keys_are_compared_including_missing_values(case, r_process):
    del r_process["result"]["statistics"]["CUBHE"]
    r_process["result"]["statistics"]["CPB"] = None
    report = comparison.compare_with_r(**case)
    rows = {row["metric"]: row for row in report["comparison"]}
    assert rows["CUBHE"]["reason"] == "missing from R"
    assert rows["CPB"]["reason"] == "missing from Python"
    assert report["matches"] is False


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_statistics_are_equivalent_to_json_null(case, r_process, value):
    case["python_statistics"]["CPB"] = value
    r_process["result"]["statistics"]["CPB"] = None
    report = comparison.compare_with_r(**case)
    assert report["matches"] is True
    assert report["python_statistics"]["CPB"] is None
    assert "NaN" not in Path(report["files"]["comparison"]).read_text()


def test_null_does_not_match_a_finite_number(case, r_process):
    r_process["result"]["statistics"]["d"] = None
    report = comparison.compare_with_r(**case)
    assert not report["matches"]
    row = next(row for row in report["comparison"] if row["metric"] == "d")
    assert row["reason"] == "finite and nonfinite values differ"


def test_tolerances_allow_roundoff_and_preserve_numbers(case, r_process):
    r_process["result"]["statistics"]["d"] += 1e-11
    report = comparison.compare_with_r(**case)
    assert report["matches"] is True
    assert report["r_result"]["statistics"]["d"] != 2.0


@pytest.mark.parametrize(
    "name", ["r_parameters.json", "r_result.json", "r_comparison.json", "r_comparison.log"]
)
def test_preexisting_artifact_is_never_overwritten(case, r_process, name):
    output = case["output_dir"]
    output.mkdir()
    protected = output / name
    protected.write_text("existing work")
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        comparison.compare_with_r(**case)
    assert protected.read_text() == "existing work"
    assert list(output.iterdir()) == [protected]
    assert not r_process["calls"]


def test_unrelated_pipeline_files_can_share_output_directory(case, r_process):
    case["output_dir"].mkdir()
    existing = case["output_dir"] / "python_result.json"
    existing.write_text("preserve")
    comparison.compare_with_r(**case)
    assert existing.read_text() == "preserve"


def test_r_failures_preserve_logs_and_raise(case, r_process):
    r_process["status"] = 1
    r_process["stderr"] = "Error: missing package coRdon"
    with pytest.raises(RuntimeError, match="exit status 1"):
        comparison.compare_with_r(**case)
    log = (case["output_dir"] / "r_comparison.log").read_text()
    assert "missing package coRdon" in log
    assert (case["output_dir"] / "r_parameters.json").is_file()


@pytest.mark.parametrize(
    "stdout",
    [
        "not JSON",
        "[]",
        '{"statistics": {"d":2}}',
        '{"statistics":{"d":[2]},"versions":{"R":"test","gRodon":"2"}}',
        '{"statistics":{"d":true},"versions":{"R":"test","gRodon":"2"}}',
        '{"statistics":{"d":NaN},"versions":{"R":"test","gRodon":"2"}}',
    ],
)
def test_malformed_r_output_is_not_claimed_as_a_comparison(case, r_process, stdout):
    r_process["stdout"] = stdout
    with pytest.raises(RuntimeError, match="malformed comparison"):
        comparison.compare_with_r(**case)
    assert stdout in (case["output_dir"] / "r_comparison.log").read_text()


def test_timeout_has_diagnostic_log_and_no_shell(case, monkeypatch):
    def run(command, **kwargs):
        assert kwargs["shell"] is False
        raise subprocess.TimeoutExpired(
            command, kwargs["timeout"], output=b"partial stdout", stderr=b"partial error"
        )

    monkeypatch.setattr(comparison.subprocess, "run", run)
    with pytest.raises(RuntimeError, match="timed out after 2"):
        comparison.compare_with_r(**case, timeout=2)
    log = (case["output_dir"] / "r_comparison.log").read_text()
    assert "partial stdout" in log and "partial error" in log


def test_invalid_input_is_rejected_before_creating_output(case, r_process):
    case["cds_path"].unlink()
    with pytest.raises(FileNotFoundError, match="input does not exist"):
        comparison.compare_with_r(**case)
    assert not case["output_dir"].exists()


@pytest.mark.parametrize("timeout", [0, -1, float("inf"), True])
def test_invalid_timeout(case, r_process, timeout):
    with pytest.raises(ValueError, match="timeout"):
        comparison.compare_with_r(**case, timeout=timeout)
    assert not r_process["calls"]


def test_arguments_stay_json_data_not_r_code(case, r_process):
    unusual = "full; this is not R code"
    report = comparison.compare_with_r(**case, mode=unusual, temperature=37)
    parameters = json.loads(Path(report["files"]["parameters"]).read_text())
    assert parameters["mode"] == unusual
    assert parameters["temperature"] == 37
    assert all(unusual not in argument for argument in r_process["calls"][0][0])
