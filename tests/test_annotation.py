"""Synthetic assembly and homology fixtures; external tools are optional."""

from __future__ import annotations

import csv
import gzip
import json
import os
import random
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

import grodon.annotation as annotation
from grodon.annotation import annotate_ribosomal, call_genes
from grodon.io import read_fasta


@pytest.fixture
def fake_pyrodigal(monkeypatch):
    calls: dict[str, Any] = {"train": [], "find": []}

    class Gene:
        begin, end, strand = 1, 9, 1
        partial_begin, partial_end = False, False

        def sequence(self):
            return "ATGGCTTAA"

        def translate(self):
            return "MA*"

    class Finder:
        def __init__(self, **kwargs):
            calls["options"] = kwargs

        def train(self, *sequences, **kwargs):
            calls["train"].append((sequences, kwargs))

        def find_genes(self, sequence):
            calls["find"].append(sequence)
            return [Gene()]

    library = SimpleNamespace(
        GeneFinder=Finder,
        __version__="3.7.1",
        MetagenomicBins=list,
        METAGENOMIC_BINS=[
            SimpleNamespace(training_info=SimpleNamespace(translation_table=code))
            for code in (4, 11, 11)
        ],
    )
    monkeypatch.setattr(annotation, "_load_pyrodigal", lambda: library)
    return calls


def test_gene_calling_trains_all_contigs_and_preserves_safe_ids(tmp_path, fake_pyrodigal):
    genome = tmp_path / "assembly.fna.gz"
    with gzip.open(genome, "wt") as stream:
        stream.write(
            ">good.1 comment\nATGGCTTAA\n>weird|id\nATGRYSWKMBDHVNTAA\n>contig_00002\nACGTACGT\n"
        )
    result = call_genes(genome, tmp_path / "output")
    assert fake_pyrodigal["options"] == {"meta": False, "closed": True, "mask": True}
    assert fake_pyrodigal["train"] == [
        (("ATGGCTTAA", "ATGRYSWKMBDHVNTAA", "ACGTACGT"), {"translation_table": 11})
    ]
    assert result["total_predicted_genes"] == 3
    assert result["counts_by_contig"] == {"good.1": 1, "weird|id": 1, "contig_00002": 1}
    assert read_fasta(result["cds_path"]).ids == (
        "good.1_CDS00001",
        "contig_00002_1_CDS00001",
        "contig_00002_CDS00001",
    )
    assert read_fasta(result["proteins_path"]).sequences == ("MA", "MA", "MA")
    with result["gene_table_path"].open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert rows[1]["contig"] == "weird|id"
    assert rows[1]["strand"] == "1"
    assert rows[1]["length_bp"] == "9"
    for key in ("cds_path", "proteins_path", "gene_table_path", "metadata_path"):
        assert b"\r\n" not in result[key].read_bytes()
    metadata = json.loads(result["metadata_path"].read_text())
    assert metadata["cds_path"] == "predicted_cds.fna"
    assert metadata["input"] == str(genome.resolve())
    assert metadata["pyrodigal_version"] == "3.7.1"
    assert len(metadata["cds_sha256"]) == 64


def test_gene_calling_metagenome_mode_filters_translation_bins(tmp_path, fake_pyrodigal):
    genome = tmp_path / "assembly.fna"
    genome.write_text(">a\nATGGCTTAA\n")
    result = call_genes(genome, tmp_path / "output", meta=True, closed=False)
    assert fake_pyrodigal["train"] == []
    assert len(fake_pyrodigal["options"]["metagenomic_bins"]) == 2
    assert all(
        item.training_info.translation_table == 11
        for item in fake_pyrodigal["options"]["metagenomic_bins"]
    )
    assert result["settings"]["closed"] is False
    with pytest.raises(ValueError, match="no metagenomic training bins"):
        call_genes(genome, tmp_path / "other", meta=True, genetic_code=6)


@pytest.mark.parametrize("sequence", ["MELKQW", "ACGT--ACGT", "ACGT*", "ACGTU"])
def test_gene_calling_rejects_protein_and_gaps_before_loading_dependency(
    tmp_path, monkeypatch, sequence
):
    genome = tmp_path / "assembly.fna"
    genome.write_text(f">contig\n{sequence}\n")
    monkeypatch.setattr(
        annotation,
        "_load_pyrodigal",
        lambda: pytest.fail("Dependency should not be loaded for invalid DNA"),
    )
    with pytest.raises(ValueError, match="DNA IUPAC"):
        call_genes(genome, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_missing_pyrodigal_has_install_instruction(tmp_path, monkeypatch):
    genome = tmp_path / "assembly.fna"
    genome.write_text(">a\nATGGCTTAA\n")

    def missing(name):
        raise ModuleNotFoundError(name)

    monkeypatch.setattr(annotation.importlib, "import_module", missing)
    with pytest.raises(ImportError, match=r"grodon2-py\[genome\]"):
        call_genes(genome, tmp_path / "output")


def test_gene_calling_never_overwrites_previous_output(tmp_path, fake_pyrodigal):
    genome = tmp_path / "assembly.fna"
    genome.write_text(">a\nATGGCTTAA\n")
    result = call_genes(genome, tmp_path / "output")
    before = result["cds_path"].read_bytes()
    with pytest.raises(FileExistsError):
        call_genes(genome, tmp_path / "output")
    assert result["cds_path"].read_bytes() == before


@pytest.mark.parametrize("code", [True, 11.5, "11"])
def test_gene_caller_requires_integer_translation_code(tmp_path, code):
    with pytest.raises(ValueError, match="genetic_code"):
        call_genes(tmp_path / "unused", tmp_path / "output", genetic_code=code)


@pytest.fixture
def proteins(tmp_path):
    source = tmp_path / "proteins.faa"
    reference = tmp_path / "references.faa"
    source.write_text("".join(f">q{i}\nM{'A' * 99}\n" for i in range(1, 4)))
    reference.write_text(
        "".join(f">ref{letter} ribosomal protein {letter}\nM{'K' * 99}\n" for letter in ("A", "B"))
    )
    return source, reference


def _hit(query, reference, *, identity=80, length=90, evalue=1e-12, bitscore=100):
    return f"{query}\t{reference}\t{identity}\t{length}\t100\t100\t1\t{length}\t1\t{length}\t{evalue}\t{bitscore}\n"


def _fake_blast(monkeypatch, content, *, exit_code=0, timeout=False):
    calls = []
    monkeypatch.setattr(annotation, "_resolve_blastp", lambda executable: Path("synthetic-blastp"))

    def run(command, **kwargs):
        calls.append((command, kwargs))
        if command[-1] == "-version":
            return subprocess.CompletedProcess(command, 0, "blastp: 2.16.0+\n", "")
        output = Path(kwargs["cwd"]) / command[command.index("-out") + 1]
        output.write_text(content, encoding="utf-8", newline="\n")
        if timeout:
            raise subprocess.TimeoutExpired(
                command, kwargs["timeout"], output="partial", stderr="timed out"
            )
        return subprocess.CompletedProcess(
            command, exit_code, "", "synthetic failure" if exit_code else ""
        )

    monkeypatch.setattr(annotation.subprocess, "run", run)
    return calls


def test_ribosomal_filters_boundaries_ties_and_query_order(tmp_path, proteins, monkeypatch):
    content = (
        _hit("q2", "refA", identity=69.99)
        + _hit("q1", "refB")
        + _hit("q3", "refA", identity=100, length=79)
        + _hit("q1", "refA")
        + _hit("q2", "refB", identity=70, length=80, evalue=1e-8, bitscore=80)
    )
    calls = _fake_blast(monkeypatch, content)
    result = annotate_ribosomal(*proteins, tmp_path / "output", threads=2, timeout=45)
    assert result["selected_ids"] == ["q1", "q2"]
    assert result["selected_hits"][0]["reference"] == "refA"
    assert result["selected_hits"][0]["description"] == "refA ribosomal protein A"
    assert result["selected_hits"][1]["query_coverage"] == 0.8
    assert result["selected_count"] == 2
    command, settings = calls[1]
    assert settings["shell"] is False
    assert settings["timeout"] == 45
    assert command[command.index("-seg") + 1] == "no"
    assert command[command.index("-comp_based_stats") + 1] == "0"
    assert command[command.index("-max_hsps") + 1] == "1"
    assert command[command.index("-max_target_seqs") + 1] == str(result["reference_count"])
    assert result["settings"]["max_target_seqs"] == result["reference_count"]
    assert command[command.index("-num_threads") + 1] == "2"
    assert command[command.index("-evalue") + 1] == "1e-08"
    assert result["ids_path"].read_text() == "q1\nq2\n"
    metadata = json.loads(result["metadata_path"].read_text())
    assert metadata["hits_path"] == "ribosomal_hits.tsv"
    assert metadata["blastp_version"] == "blastp: 2.16.0+"


def test_ribosomal_tie_selection_independent_of_hit_row_order(tmp_path, proteins, monkeypatch):
    rows = [_hit("q1", "refB"), _hit("q1", "refA")]
    _fake_blast(monkeypatch, "".join(rows))
    first = annotate_ribosomal(*proteins, tmp_path / "first")
    _fake_blast(monkeypatch, "".join(reversed(rows)))
    second = annotate_ribosomal(*proteins, tmp_path / "second")
    assert first["selected_hits"] == second["selected_hits"]


def test_ribosomal_empty_matches_are_a_valid_annotation_result(tmp_path, proteins, monkeypatch):
    _fake_blast(monkeypatch, "")
    result = annotate_ribosomal(*proteins, tmp_path / "output")
    assert result["selected_ids"] == []
    assert result["selected_count"] == 0
    assert result["ids_path"].read_text() == ""
    assert result["selected_hits_path"].read_text().startswith("id\treference\t")


@pytest.mark.parametrize("failure", ["timeout", "nonzero"])
def test_blast_failure_preserves_log_and_removes_partial_hits(
    tmp_path, proteins, monkeypatch, failure
):
    _fake_blast(
        monkeypatch,
        _hit("q1", "refA"),
        timeout=failure == "timeout",
        exit_code=2 if failure == "nonzero" else 0,
    )
    output = tmp_path / "output"
    with pytest.raises(RuntimeError, match="BLAST"):
        annotate_ribosomal(*proteins, output, timeout=0.01)
    assert (output / "blastp.log").is_file()
    assert not (output / "ribosomal_hits.tsv").exists()
    assert not (output / "ribosomal_annotation.json").exists()
    if failure == "timeout":
        assert "partial" in (output / "blastp.log").read_text()


def test_blast_gzip_inputs_are_expanded_and_scratch_removed(tmp_path, proteins, monkeypatch):
    compressed = []
    for path in proteins:
        gzpath = path.with_suffix(path.suffix + ".gz")
        with gzip.open(gzpath, "wt") as stream:
            stream.write(path.read_text())
        compressed.append(gzpath)
    calls = _fake_blast(monkeypatch, _hit("q1", "refA"))
    output = tmp_path / "output"
    result = annotate_ribosomal(*compressed, output)
    assert result["selected_ids"] == ["q1"]
    assert calls[1][0][calls[1][0].index("-query") + 1] == ".blast_query.faa"
    assert not list(output.glob(".blast_*.faa"))


@pytest.mark.parametrize(
    "content",
    [
        "bad\tcolumns\n",
        _hit("missing", "refA"),
        _hit("q1", "refA", identity=float("nan")),
        _hit("q1", "refA", length=101),
    ],
)
def test_malformed_blast_results_fail_instead_of_selecting_genes(
    tmp_path, proteins, monkeypatch, content
):
    _fake_blast(monkeypatch, content)
    with pytest.raises(ValueError, match="BLAST"):
        annotate_ribosomal(*proteins, tmp_path / "output")
    assert not (tmp_path / "output" / "highly_expressed_ids.txt").exists()


@pytest.mark.parametrize(
    "options",
    [
        {"min_identity": -1},
        {"min_query_coverage": 1.1},
        {"min_reference_coverage": float("nan")},
        {"max_evalue": 0},
        {"threads": True},
        {"threads": 1.5},
        {"timeout": 0},
    ],
)
def test_invalid_blast_options_fail_before_io(tmp_path, options):
    with pytest.raises(ValueError):
        annotate_ribosomal(tmp_path / "unused", tmp_path / "unused", tmp_path / "output", **options)


def test_missing_blast_has_clear_local_install_message(tmp_path, proteins, monkeypatch):
    monkeypatch.setattr(annotation.shutil, "which", lambda value: None)
    with pytest.raises(FileNotFoundError, match="Install NCBI BLAST"):
        annotate_ribosomal(*proteins, tmp_path / "output", blastp="not-installed-blastp")


def test_real_pyrodigal_on_artificial_sequence_when_installed(tmp_path):
    pytest.importorskip("pyrodigal")

    rng = random.Random(19)
    codons = ["GCT", "GCC", "GAT", "GAC", "GAA", "GAG", "AAA", "AAG", "TCT", "AGC", "CCT", "CCG"]
    dna = "TAA" * 30 + "AGGAGGAAAAAAAATG" + "".join(rng.choices(codons, k=300)) + "TAA" + "TAA" * 30
    genome = tmp_path / "artificial.fna"
    genome.write_text(">artificial_contig\n" + dna + "\n")
    result = call_genes(genome, tmp_path / "output", meta=True)
    assert result["total_predicted_genes"] > 0
    genes = read_fasta(result["cds_path"])
    assert all(len(sequence) % 3 == 0 for sequence in genes)
    assert all(not sequence.endswith("*") for sequence in read_fasta(result["proteins_path"]))


def test_real_blast_on_artificial_proteins_when_available(tmp_path):
    executable = os.environ.get("GRODON_TEST_BLASTP") or shutil.which("blastp")
    if executable is None:
        pytest.skip("Set GRODON_TEST_BLASTP or install BLAST+ to run the local integration test")
    rng = random.Random(190)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    sequence = "M" + "".join(rng.choices(alphabet, k=119))
    different = "M" + "".join(rng.choices(alphabet, k=119))
    query = tmp_path / "query.faa"
    reference = tmp_path / "reference.faa"
    query.write_text(f">synthetic_match\n{sequence}\n>synthetic_other\n{different}\n")
    reference.write_text(f">synthetic_reference artificial test reference\n{sequence}\n")
    result = annotate_ribosomal(query, reference, tmp_path / "output", blastp=executable)
    assert result["selected_ids"] == ["synthetic_match"]
    assert result["selected_hits"][0]["identity_percent"] == 100
    assert result["selected_hits"][0]["query_coverage"] == 1
    assert "blastp:" in result["blastp_version"]
    assert b"\r\n" not in result["hits_path"].read_bytes()
