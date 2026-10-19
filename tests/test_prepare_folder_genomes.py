"""Offline checks for exact-accession input preparation and archive provenance."""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
from pathlib import Path
from unittest.mock import Mock

import pytest

_SCRIPT = Path(__file__).parents[1] / "scripts" / "prepare_folder_genomes.py"
_SPEC = importlib.util.spec_from_file_location("prepare_folder_genomes", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
preparer = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(preparer)

ACCESSION = "GCF_000000001.1"
REPORT_URL = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{ACCESSION}/dataset_report"
BASE = f"{ACCESSION}_ASM1v1"
DIRECTORY = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/001/{BASE}/"


def cached_file(path, content, url):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    receipt = {
        "url": url,
        "sha256": preparer.digest(path),
        "md5": preparer.digest(path, "md5"),
        "bytes": len(content),
        "retrieved_utc": "2026-01-01T00:00:00+00:00",
    }
    preparer.save_json(path.with_name(path.name + ".receipt.json"), receipt)
    return path


@pytest.fixture
def no_network(monkeypatch):
    network = Mock(side_effect=AssertionError("Unexpected network access"))
    monkeypatch.setattr(preparer.urllib.request, "urlopen", network)
    return network


@pytest.fixture
def assembly_cache(tmp_path):
    cache = tmp_path / "cache"
    folder = cache / ACCESSION
    report = {
        "reports": [
            {
                "accession": ACCESSION,
                "assembly_info": {"assembly_name": "ASM1v1", "assembly_level": "Complete Genome"},
                "organism": {"organism_name": "Synthetic test", "tax_id": 1},
            }
        ]
    }
    report_path = cached_file(
        folder / "dataset_report.json", json.dumps(report).encode(), REPORT_URL
    )
    genome_name, protein_name = BASE + "_genomic.fna.gz", BASE + "_protein.faa.gz"
    genome = cached_file(
        folder / genome_name,
        gzip.compress(b">contig\nATGAAATAA\n", mtime=0),
        DIRECTORY + genome_name,
    )
    protein = cached_file(
        folder / protein_name,
        gzip.compress(b">p1 current label\nMK\n>p2 current label\nMA\n", mtime=0),
        DIRECTORY + protein_name,
    )
    checksums = (
        f"{preparer.digest(genome, 'md5')}  ./{genome_name}\n"
        f"{preparer.digest(protein, 'md5')}  ./{protein_name}\n"
    )
    cached_file(folder / "md5checksums.txt", checksums.encode(), DIRECTORY + "md5checksums.txt")
    return cache, genome, protein, report_path


def test_valid_cache_works_offline_and_detects_tampering(tmp_path, no_network):
    url = "https://ftp.ncbi.nlm.nih.gov/test"
    path = cached_file(tmp_path / "file", b"original", url)
    result = preparer.cached_resource(url, path, fetch=False, limit=100)
    assert result["sha256"] == preparer.digest(path)
    path.write_bytes(b"changed")
    with pytest.raises(preparer.NeedsFetch):
        preparer.cached_resource(url, path, fetch=False, limit=100)
    no_network.assert_not_called()


def test_missing_cache_requires_explicit_fetch(tmp_path, no_network):
    with pytest.raises(preparer.NeedsFetch, match="--fetch"):
        preparer.cached_resource(REPORT_URL, tmp_path / "report", fetch=False, limit=100)
    no_network.assert_not_called()


def test_download_verifies_checksum_and_leaves_no_false_cache(tmp_path, monkeypatch):
    class Response(io.BytesIO):
        headers = {"Content-Length": "4"}

        def geturl(self):
            return REPORT_URL

    monkeypatch.setattr(preparer, "_throttle", lambda: None)
    monkeypatch.setattr(preparer.urllib.request, "urlopen", lambda *a, **kw: Response(b"data"))
    target = tmp_path / "target"
    with pytest.raises(ValueError, match="MD5 mismatch"):
        preparer.cached_resource(REPORT_URL, target, fetch=True, limit=100, expected_md5="0" * 32)
    assert not target.exists()
    assert not target.with_name("target.part").exists()
    assert not target.with_name("target.receipt.json").exists()


def test_exact_assembly_report_is_required(assembly_cache, no_network):
    cache, genome, _, report_path = assembly_cache
    payload = json.loads(report_path.read_text())
    payload["reports"][0]["accession"] = "GCF_000000001.2"
    cached_file(report_path, json.dumps(payload).encode(), REPORT_URL)
    result = preparer.prepare_entry(genome, ACCESSION, "genomic", cache, False)
    assert result["status"] == "error"
    assert "exact assembly" in result["error"]
    assert result["genome_path"] is None


def test_original_genome_is_preserved_and_verified(assembly_cache, tmp_path, no_network):
    cache, genome, _, _ = assembly_cache
    source = tmp_path / "original.fna.gz"
    source.write_bytes(genome.read_bytes())
    original = source.read_bytes()
    result = preparer.prepare_entry(source, ACCESSION, "genomic", cache, False)
    assert result["status"] == "ready"
    assert result["genome_path"] == str(source)
    assert result["download"]["used_original_genome"] is True
    assert result["download"]["genome_md5_matches"] is True
    assert source.read_bytes() == original
    source.write_bytes(gzip.compress(b">different\nATGTAATAA\n", mtime=0))
    result = preparer.prepare_entry(source, ACCESSION, "genomic", cache, False)
    assert result["status"] == "error"
    assert "MD5 mismatch" in result["error"]


def test_protein_header_changes_are_not_sequence_changes(assembly_cache, tmp_path, no_network):
    cache, genome, protein, _ = assembly_cache
    source = tmp_path / "original.faa.gz"
    source.write_bytes(gzip.compress(b">p1 old label\nMK\n>p2 old label\nMA\n", mtime=0))
    original = source.read_bytes()
    result = preparer.prepare_entry(source, ACCESSION, "protein", cache, False)
    details = result["download"]
    assert result["status"] == "ready"
    assert result["genome_path"] == str(genome)
    assert details["official_protein_path"] == str(protein)
    assert details["source_protein_md5_matches"] is False
    assert details["source_protein_content_matches"] is True
    assert details["source_protein_order_matches"] is True
    assert details["source_protein_headers_match"] is False
    assert details["source_protein_original_count"] == details["source_protein_official_count"] == 2
    assert details["source_protein_changed_sequence_count"] == 0
    assert source.read_bytes() == original


def test_protein_changes_and_record_order_are_reported(tmp_path):
    original, official = tmp_path / "old.gz", tmp_path / "new.gz"
    original.write_bytes(gzip.compress(b">p1\nMK\n>p2\nMA\n>p3\nMP\n", mtime=0))
    official.write_bytes(gzip.compress(b">p3\nMP\n>p1\nMM\n>p4\nMT\n", mtime=0))
    result = preparer.protein_comparison(original, official)
    assert result["source_protein_content_matches"] is False
    assert result["source_protein_order_matches"] is False
    assert result["source_protein_shared_id_count"] == 2
    assert result["source_protein_changed_sequence_count"] == 1
    assert result["source_protein_only_original_id_count"] == 1
    assert result["source_protein_only_official_id_count"] == 1


@pytest.mark.parametrize("content", [b"MK\n", b">p1\nMK\n>p1\nMA\n", b">p1\n", b">\nMK\n"])
def test_invalid_protein_records_are_rejected(tmp_path, content):
    source = tmp_path / "invalid.gz"
    source.write_bytes(gzip.compress(content, mtime=0))
    with pytest.raises(ValueError):
        preparer.protein_snapshot(source)


def test_checksum_manifest_rejects_duplicate_paths(tmp_path):
    path = tmp_path / "checksums"
    path.write_text("0" * 32 + "  ./file\n" + "1" * 32 + "  ./file\n")
    with pytest.raises(ValueError, match="Duplicate"):
        preparer.parse_checksums(path)


@pytest.mark.parametrize(
    "filenames",
    [
        ["unversioned_protein.faa.gz"],
        [f"{ACCESSION}_one_protein.faa.gz", f"{ACCESSION}_two_genomic.fna.gz"],
    ],
)
def test_cli_rejects_ambiguous_input_before_network(tmp_path, monkeypatch, no_network, filenames):
    source = tmp_path / "input"
    source.mkdir()
    for filename in filenames:
        (source / filename).write_bytes(b"irrelevant")
    monkeypatch.setattr(
        preparer.sys,
        "argv",
        [
            str(_SCRIPT),
            "--input-dir",
            str(source),
            "--cache-dir",
            str(tmp_path / "cache"),
            "--manifest",
            str(tmp_path / "manifest.json"),
        ],
    )
    with pytest.raises(SystemExit) as raised:
        preparer.main()
    assert raised.value.code == 2
    no_network.assert_not_called()
