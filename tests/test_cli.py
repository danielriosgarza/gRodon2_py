"""Exercise the public command-line interface in fresh Python processes."""

import gzip
import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from grodon.codon import CODONS
from grodon.copiotrophy import predict_ioc


def invoke(*arguments, cwd):
    return subprocess.run(
        [sys.executable, "-m", "grodon", *(str(value) for value in arguments)],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )


@pytest.fixture
def synthetic_fasta(tmp_path):
    rng = np.random.default_rng(90210)
    sequences = ["".join(rng.choice(CODONS, size=120)) for _ in range(8)]
    content = "".join(
        f">g{index} {'30S ribosomal protein S1' if index < 3 else 'synthetic protein'}\n"
        f"{sequence}\n"
        for index, sequence in enumerate(sequences)
    )
    fasta = tmp_path / "synthetic.fna"
    fasta.write_text(content, encoding="utf-8")
    return fasta, sequences


@pytest.fixture
def trait_table(tmp_path):
    traits = pd.DataFrame(
        {
            "dCUB": [-0.1, -0.2, -0.5, -0.3, -0.04, -0.7],
            "nCAZy": [10, 50, 80, 45, 3, 130],
            "nGenes": [1500, 2200, 4000, 3000, 800, 5100],
        }
    )
    path = tmp_path / "traits.csv"
    traits.to_csv(path, index=False)
    return path, traits


@pytest.mark.parametrize("arguments", [("--help",), ("--version",), ("predict", "--help")])
def test_help_and_version(arguments, tmp_path):
    result = invoke(*arguments, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert "grodon" in result.stdout


def test_fasta_statistics_and_gzip_produce_identical_json(synthetic_fasta, tmp_path):
    fasta, _ = synthetic_fasta
    compressed = tmp_path / "synthetic.fna.gz"
    with gzip.open(compressed, "wt") as handle:
        handle.write(fasta.read_text(encoding="utf-8"))
    plain = invoke("statistics", fasta, cwd=tmp_path)
    zipped = invoke("statistics", compressed, cwd=tmp_path)
    assert plain.returncode == zipped.returncode == 0, (plain.stderr, zipped.stderr)
    actual = json.loads(plain.stdout)
    assert actual == json.loads(zipped.stdout)
    assert actual["nHE"] == 3
    assert actual["FilteredSequences"] == 0
    assert np.isfinite(actual["CUBHE"])


def test_prediction_writes_json_metadata_and_keeps_stdout_clean(synthetic_fasta, tmp_path):
    fasta, _ = synthetic_fasta
    output = tmp_path / "prediction.json"
    result = invoke(
        "predict",
        fasta,
        "--mode",
        "full",
        "--training-set",
        "vs",
        "--temperature",
        "30",
        "-o",
        output,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    payload = json.loads(output.read_text())
    assert payload["metadata"] == {
        "model": "gRodon_model_temp",
        "mode": "full",
        "training_set": "vs",
        "background": "all",
    }
    assert {"d", "LowerCI", "UpperCI", "CUBHE", "OGT"}.issubset(payload)
    assert payload["OGT"] == 30
    assert "NaN" not in output.read_text() and "Infinity" not in output.read_text()
    assert "Less than 10" in result.stderr


def test_coverage_aligns_by_id_and_nan_serializes_as_null(synthetic_fasta, tmp_path):
    fasta, sequences = synthetic_fasta
    coverage = tmp_path / "coverage.tsv"
    weights = np.arange(1, 9)
    coverage.write_text(
        "gene\tdepth\n" + "".join(f"g{i}\t{weights[i]}\n" for i in reversed(range(8)))
    )
    result = invoke(
        "statistics",
        fasta,
        "--coverage",
        coverage,
        "--coverage-id-column",
        "gene",
        "--coverage-column",
        "depth",
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    expected_gc = np.average(
        [(seq.count("G") + seq.count("C")) / len(seq) for seq in sequences], weights=weights
    )
    assert payload["GC"] == pytest.approx(expected_gc)
    assert payload["CPB"] is None
    assert "NaN" not in result.stdout


def test_highly_expressed_ids_override_annotations_and_support_comments(synthetic_fasta, tmp_path):
    fasta, _ = synthetic_fasta
    selected = tmp_path / "he.txt"
    selected.write_text("\ufeff  # selected synthetic genes\n g3 \n\ng5\n", encoding="utf-8")
    result = invoke("statistics", fasta, "--he-ids", selected, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["nHE"] == 2


def test_individual_background_seed_is_reproducible(synthetic_fasta, tmp_path):
    fasta, _ = synthetic_fasta
    arguments = ("statistics", fasta, "--background", "individual", "--n-le", "3", "--seed", "42")
    first = invoke(*arguments, cwd=tmp_path)
    second = invoke(*arguments, cwd=tmp_path)
    assert first.returncode == second.returncode == 0, (first.stderr, second.stderr)
    assert first.stdout == second.stdout
    assert json.loads(first.stdout)["CPB"] is None


def test_fit_ioc_serialized_model_roundtrip(trait_table, tmp_path):
    table, _ = trait_table
    model = tmp_path / "custom-model.json"
    scores = tmp_path / "training-scores.json"
    fitted = invoke(
        "fit-ioc", table, "--save-model", model, "-o", scores, "--completeness", "90", cwd=tmp_path
    )
    assert fitted.returncode == 0, fitted.stderr
    assert fitted.stdout == ""
    serialized = json.loads(model.read_text())
    assert serialized["features"] == ["rCAZy", "dCUB", "nGenes"]
    projected = invoke("ioc", table, "--model", model, "--completeness", "90", cwd=tmp_path)
    assert projected.returncode == 0, projected.stderr
    np.testing.assert_allclose(
        json.loads(projected.stdout)["IoC"], json.loads(scores.read_text())["IoC"]
    )


def test_builtin_ioc_matches_python_api(trait_table, tmp_path):
    table, traits = trait_table
    result = invoke("ioc", table, "--model", "pacific", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    np.testing.assert_allclose(
        json.loads(result.stdout)["IoC"], predict_ioc(traits, model="pacific")
    )


@pytest.mark.parametrize(
    ("options", "message"),
    [
        (("--he-pattern", "["), "unterminated"),
        (("--genetic-code", "0"), "genetic code"),
        (("--n-le", "0"), "n_le"),
        (("--he-pattern", "missing annotation"), "No highly expressed"),
    ],
)
def test_invalid_statistics_arguments_are_clean_errors(synthetic_fasta, tmp_path, options, message):
    fasta, _ = synthetic_fasta
    result = invoke("statistics", fasta, *options, cwd=tmp_path)
    assert result.returncode == 2
    assert message.lower() in result.stderr.lower()
    assert "Traceback" not in result.stderr
    assert result.stdout == ""


def test_missing_coverage_and_bad_model_have_clean_errors(synthetic_fasta, trait_table, tmp_path):
    fasta, _ = synthetic_fasta
    coverage = tmp_path / "coverage.tsv"
    coverage.write_text("id\tcoverage\ng0\t1\n")
    missing = invoke("statistics", fasta, "--coverage", coverage, cwd=tmp_path)
    assert missing.returncode == 2
    assert "missing FASTA IDs" in missing.stderr
    malformed = tmp_path / "model.json"
    malformed.write_text("{broken")
    invalid = invoke("ioc", trait_table[0], "--model", malformed, cwd=tmp_path)
    assert invalid.returncode == 2
    assert "Traceback" not in invalid.stderr


def test_fit_model_output_collision_is_rejected_without_overwriting(trait_table, tmp_path):
    output = tmp_path / "same.json"
    output.write_text("existing content")
    result = invoke("fit-ioc", trait_table[0], "--save-model", output, "-o", output, cwd=tmp_path)
    assert result.returncode == 2
    assert "different files" in result.stderr
    assert output.read_text() == "existing content"
