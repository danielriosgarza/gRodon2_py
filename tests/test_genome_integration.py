"""Opt-in reproduction of the previously validated benign Lactococcus assembly.

Set GRODON_TEST_GENOME to CIRM-BIA2698's assembly and GRODON_TEST_BLASTP to a
local BLASTP executable. GRODON_TEST_RSCRIPT optionally checks original R too.
No genome is downloaded and this test does not run when those paths are absent.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from grodon.io import read_fasta
from grodon.pipeline import GenomePipelineConfig, run_genome


@pytest.mark.genome_integration
def test_lactococcus_pipeline_reproduces_validated_result(tmp_path):
    genome = os.environ.get("GRODON_TEST_GENOME")
    blastp = os.environ.get("GRODON_TEST_BLASTP")
    if not genome or not blastp:
        pytest.skip(
            "Set GRODON_TEST_GENOME and GRODON_TEST_BLASTP for the opt-in Lactococcus integration"
        )
    pytest.importorskip("pyrodigal")
    rscript = os.environ.get("GRODON_TEST_RSCRIPT")
    result = run_genome(
        Path(genome),
        tmp_path / "lactococcus_run",
        config=GenomePipelineConfig(reference="lactococcus-v1", blastp=blastp),
        rscript=rscript,
    )
    summary = result.to_dict()
    assert summary["status"] == "complete"
    assert summary["predicted_genes"] == 2832
    assert summary["ribosomal_genes"] == 54
    assert summary["retained_genes"] == 2600
    assert summary["retained_ribosomal_genes"] == 42
    assert summary["statistics"]["FilteredSequences"] == 232
    assert summary["model"] == "gRodon_model_base_madin"
    assert result.prediction.d == pytest.approx(0.120516194130022, rel=1e-10, abs=1e-12)
    assert len(read_fasta(result.output_dir / "predicted_cds.fna")) == 2832
    selected = (result.output_dir / "highly_expressed_ids.txt").read_text().splitlines()
    assert len(selected) == len(set(selected)) == 54
    assert len(summary["scientific_fingerprint"]) == 64
    manifest = json.loads((result.output_dir / "run.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["scientific_fingerprint"] == summary["scientific_fingerprint"]
    for name, artifact in manifest["artifacts"].items():
        data = (result.output_dir / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == artifact["sha256"], name
        assert len(data) == artifact["bytes"], name
    saved = json.loads((result.output_dir / "result.json").read_text())
    assert saved == summary
    assert (result.output_dir / "reference.faa.gz").is_file()
    if rscript:
        assert summary["r_matches"] is True
        assert (result.output_dir / "r_comparison.json").is_file()
