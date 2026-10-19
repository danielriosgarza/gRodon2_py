"""The public genome commands route settings and errors without hidden downloads."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from grodon.genome_cli import COMMANDS, register_commands, run_command
from grodon.references import DEFAULT_REFERENCE


@pytest.fixture
def parser():
    parser = argparse.ArgumentParser(prog="grodon")
    register_commands(parser.add_subparsers(dest="command", required=True))
    return parser


@pytest.fixture
def pipeline(monkeypatch):
    fake = SimpleNamespace(
        GenomePipelineConfig=Mock(side_effect=lambda **kwargs: SimpleNamespace(**kwargs)),
        run_genome=Mock(return_value=SimpleNamespace(to_dict=lambda: {"status": "complete"})),
        run_genomes=Mock(return_value={"succeeded": 2, "failed": 0}),
    )
    monkeypatch.setitem(sys.modules, "grodon.pipeline", fake)
    return fake


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_subcommand_help(command, parser, capsys):
    with pytest.raises(SystemExit) as error:
        parser.parse_args([command, "--help"])
    assert error.value.code == 0
    assert f"grodon {command}" in capsys.readouterr().out


@pytest.mark.parametrize("arguments", [["genome", "assembly.fna"], ["genome-batch"], ["setup"]])
def test_required_inputs_and_output_directories(arguments, parser, capsys):
    with pytest.raises(SystemExit) as error:
        parser.parse_args(arguments)
    assert error.value.code == 2
    assert "required" in capsys.readouterr().err


def test_default_single_genome_routes_reproducible_settings(parser, pipeline):
    args = parser.parse_args(["genome", "assembly.fna.gz", "--outdir", "results"])
    assert run_command(args) == {"status": "complete"}
    call = pipeline.run_genome.call_args
    assert call.args == (Path("assembly.fna.gz"), Path("results"))
    assert call.kwargs["rscript"] is None
    config = call.kwargs["config"]
    assert config.reference == DEFAULT_REFERENCE
    assert config.blastp == "blastp"
    assert config.genetic_code == 11
    assert config.gene_mode == "single"
    assert config.closed is config.mask is True
    assert config.min_identity == 70
    assert config.min_query_coverage == config.min_reference_coverage == 0.8
    assert config.max_evalue == 1e-8
    assert config.threads == 1
    assert config.timeout == 600
    assert args.output is None
    pipeline.run_genomes.assert_not_called()


def test_all_single_genome_settings_are_passed_to_config(parser, pipeline):
    args = parser.parse_args(
        [
            "genome",
            "assembly.fna",
            "--outdir",
            "results",
            "--reference",
            "custom.faa",
            "--blastp",
            "bin/blastp",
            "--rscript",
            "bin/Rscript",
            "--mode",
            "partial",
            "--training-set",
            "vs",
            "--temperature",
            "30",
            "--genetic-code",
            "4",
            "--gene-mode",
            "meta",
            "--no-closed",
            "--no-mask",
            "--min-identity",
            "85",
            "--min-query-coverage",
            "0.9",
            "--min-reference-coverage",
            "0.95",
            "--max-evalue",
            "1e-12",
            "--min-ribosomal",
            "20",
            "--confidence",
            "0.9",
            "--threads",
            "2",
            "--timeout",
            "120",
        ]
    )
    run_command(args)
    assert pipeline.GenomePipelineConfig.call_args.kwargs == {
        "reference": "custom.faa",
        "blastp": "bin/blastp",
        "mode": "partial",
        "training_set": "vs",
        "temperature": 30.0,
        "genetic_code": 4,
        "gene_mode": "meta",
        "closed": False,
        "mask": False,
        "min_identity": 85.0,
        "min_query_coverage": 0.9,
        "min_reference_coverage": 0.95,
        "max_evalue": 1e-12,
        "min_ribosomal": 20,
        "confidence": 0.9,
        "threads": 2,
        "timeout": 120.0,
    }
    assert pipeline.run_genome.call_args.kwargs["rscript"] == Path("bin/Rscript")


def test_batch_routes_all_inputs_and_failure_policy(parser, pipeline):
    args = parser.parse_args(
        ["genome-batch", "one.fna", "two.fna.gz", "--outdir", "results", "--continue-on-error"]
    )
    assert run_command(args) == {"succeeded": 2, "failed": 0}
    call = pipeline.run_genomes.call_args
    assert call.args == ([Path("one.fna"), Path("two.fna.gz")], Path("results"))
    assert call.kwargs["continue_on_error"] is True
    pipeline.run_genome.assert_not_called()


def test_doctor_does_not_invoke_installer(parser, monkeypatch):
    report = {"ready": False, "blastp": {"available": False}}
    doctor = Mock(return_value=report)
    installer = Mock(side_effect=AssertionError("diagnostics must not download tools"))
    monkeypatch.setattr("grodon.tools.doctor", doctor)
    monkeypatch.setattr("grodon.tools.install_blast", installer)
    assert run_command(parser.parse_args(["doctor", "--blastp", "custom/blastp"])) == report
    doctor.assert_called_once_with(blastp="custom/blastp")
    installer.assert_not_called()


def test_setup_only_invokes_explicit_install_target(parser, monkeypatch):
    installer = Mock(return_value={"blastp": "tools/bin/blastp"})
    monkeypatch.setattr("grodon.tools.install_blast", installer)
    args = parser.parse_args(["setup", "--tools-dir", "tools", "--version", "2.17.0"])
    assert run_command(args) == {"blastp": "tools/bin/blastp"}
    installer.assert_called_once_with(Path("tools"), version="2.17.0")


def test_references_lists_metadata_without_unpacking_or_download(parser, monkeypatch):
    monkeypatch.setattr(
        "grodon.references.gzip.decompress", Mock(side_effect=AssertionError("unnecessary unpack"))
    )
    monkeypatch.setattr(
        "urllib.request.urlopen", Mock(side_effect=AssertionError("unexpected network"))
    )
    report = run_command(parser.parse_args(["references"]))
    assert report["default"] == DEFAULT_REFERENCE
    profiles = {profile["name"]: profile for profile in report["references"]}
    assert profiles["bacteria-swissprot-v1"]["reference_count"] == 39743
    assert profiles["bacteria-swissprot-v1"]["taxon_count"] == 934
    assert profiles["bacteria-swissprot-v1"]["license"] == "CC-BY-4.0"
    assert profiles["lactococcus-v1"]["reference_count"] == 88


@pytest.mark.parametrize("exception", [ValueError("invalid setting"), FileExistsError("exists")])
def test_pipeline_errors_propagate_to_root_handler(parser, pipeline, exception):
    pipeline.run_genome.side_effect = exception
    args = parser.parse_args(["genome", "assembly.fna", "--outdir", "results"])
    with pytest.raises(type(exception), match=str(exception)):
        run_command(args)


def test_unknown_command_rejected():
    with pytest.raises(ValueError, match="Unknown genome command"):
        run_command(argparse.Namespace(command="other"))


def test_root_cli_registration_and_reference_json_in_fresh_process(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "grodon", "references"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["default"] == DEFAULT_REFERENCE
    assert "Traceback" not in result.stderr


def test_root_cli_prints_clean_pipeline_error(monkeypatch, capsys, pipeline):
    from grodon.cli import main

    pipeline.run_genome.side_effect = RuntimeError("BLASTP is unavailable")
    with pytest.raises(SystemExit) as error:
        main(["genome", "assembly.fna", "--outdir", "results"])
    assert error.value.code == 2
    output = capsys.readouterr()
    assert output.out == ""
    assert "BLASTP is unavailable" in output.err
    assert "Traceback" not in output.err


def test_root_cli_batch_failures_are_visible_to_shell(monkeypatch, capsys):
    from grodon.cli import main

    monkeypatch.setattr(
        "grodon.genome_cli.run_command",
        lambda args: {"status": "completed_with_errors", "failed": 1},
    )
    assert main(["genome-batch", "missing.fna", "--outdir", "unused"]) == 1
    assert json.loads(capsys.readouterr().out)["failed"] == 1


def test_root_cli_doctor_missing_dependencies_exits_nonzero(monkeypatch, capsys):
    from grodon.cli import main

    monkeypatch.setattr("grodon.genome_cli.run_command", lambda args: {"ready": False})
    assert main(["doctor"]) == 1
    assert json.loads(capsys.readouterr().out)["ready"] is False
