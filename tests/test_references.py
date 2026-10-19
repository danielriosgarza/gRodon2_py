"""Reference resolution is deterministic, offline, validated, and auditable."""

from __future__ import annotations

import gzip
import hashlib
import json

import pytest

from grodon.references import (
    DEFAULT_REFERENCE,
    REFERENCE_NAMES,
    resolve_reference,
    sha256_file,
    validate_protein_fasta,
)


@pytest.mark.parametrize("profile", REFERENCE_NAMES)
def test_bundled_profile_integrity_and_cleanup(profile):
    with resolve_reference(profile) as reference:
        assert reference.name == profile
        assert reference.version == "1"
        assert reference.path.is_file()
        assert sha256_file(reference.path) == reference.sha256
        assert reference.count == reference.source_metadata["reference_count"]
        assert len(validate_protein_fasta(reference.path)) == reference.count
        temporary_path = reference.path
    assert not temporary_path.exists()


def test_bundled_lactococcus_is_original_exact_snapshot():
    with resolve_reference("lactococcus-v1") as reference:
        assert reference.count == 88
        assert reference.sha256 == (
            "075848596f8e2f356d8c6ad244cc24ab27ea376a5305b32a59f790e029bb8d80"
        )
        assert "Lactococcus only" in reference.source_metadata["scope"]


def test_default_reference_is_broad_and_has_source_attribution():
    assert DEFAULT_REFERENCE == "bacteria-swissprot-v1"
    with resolve_reference() as reference:
        assert reference.count == 39743
        assert reference.source_metadata["taxon_count"] > 500
        assert reference.source_metadata["database_release"] == "2026_03"
        assert reference.source_metadata["license"] == "CC-BY-4.0"
        assert reference.source_metadata["source_url"].startswith("https://rest.uniprot.org/")


def test_resolution_never_uses_network(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Reference resolution attempted network access")

    monkeypatch.setattr("urllib.request.urlopen", forbidden)
    with resolve_reference("lactococcus-v1") as reference:
        assert reference.count == 88


def test_custom_reference_normalization_and_original_integrity(tmp_path):
    source = tmp_path / "custom.faa"
    original = b">r1 ribosomal protein L1\r\nmkvv*\r\n>r2 ribosomal protein S1\r\nMaX\r\n"
    source.write_bytes(original)
    with resolve_reference(source) as reference:
        assert reference.path != source
        assert reference.name == "custom.faa"
        assert reference.version == "custom"
        assert reference.path.read_bytes() == (
            b">r1 ribosomal protein L1\nMKVV\n>r2 ribosomal protein S1\nMAX\n"
        )
        assert reference.descriptions == {
            "r1": "r1 ribosomal protein L1",
            "r2": "r2 ribosomal protein S1",
        }
        assert reference.source_metadata["source_sha256"] == hashlib.sha256(original).hexdigest()
        assert reference.sha256 == sha256_file(reference.path)
    assert source.read_bytes() == original


def test_custom_gzip_supported(tmp_path):
    source = tmp_path / "custom.faa.gz"
    source.write_bytes(gzip.compress(b">r1 annotation\nMKVL\n", mtime=0))
    with resolve_reference(source) as reference:
        assert reference.path.read_text() == ">r1 annotation\nMKVL\n"
        assert reference.count == 1


@pytest.mark.parametrize(
    ("content", "message"),
    [
        (">r1\nMKL\n>r1 duplicate\nMKL\n", "identifiers must be unique"),
        (">r1\n", "Empty sequence"),
        (">r1\n*\n", "Empty protein sequence"),
        (">r1\nMK*L\n", "Invalid protein sequence"),
        (">r1\nMK-L\n", "Invalid protein sequence"),
        (">r1\nMK1L\n", "Invalid protein sequence"),
        (">r1 bad\tlabel\nMKL\n", "Control character"),
        (">\nMKL\n", "Empty FASTA header"),
        ("MKL\n>r1\nMKL\n", "Sequence before first FASTA header"),
        ("", "FASTA contains no sequences"),
    ],
)
def test_invalid_custom_reference_rejected(tmp_path, content, message):
    source = tmp_path / "invalid.faa"
    source.write_text(content)
    with pytest.raises(ValueError, match=message), resolve_reference(source):
        pass


def test_unknown_profile_or_missing_file_reports_choices():
    with pytest.raises(FileNotFoundError, match="neither a file nor a bundled profile"):
        with resolve_reference("missing-profile-v99"):
            pass


@pytest.mark.parametrize("tamper", ["compressed", "fasta_hash", "count", "asset_name", "name"])
def test_bundled_corruption_rejected(tmp_path, monkeypatch, tamper):
    directory = tmp_path / "data" / "references"
    directory.mkdir(parents=True)
    profile = "lactococcus-v1"
    payload = b">r1 synthetic ribosomal reference\nMKVL\n"
    compressed = gzip.compress(payload, mtime=0)
    manifest = {
        "schema_version": 1,
        "name": profile,
        "version": "1",
        "reference_count": 1,
        "asset_filename": f"{profile}.faa.gz",
        "asset_sha256": hashlib.sha256(compressed).hexdigest(),
        "fasta_sha256": hashlib.sha256(payload).hexdigest(),
    }
    if tamper == "compressed":
        compressed += b"corrupted"
    elif tamper == "fasta_hash":
        manifest["fasta_sha256"] = "bad-hash"
    elif tamper == "count":
        manifest["reference_count"] = 2
    elif tamper == "asset_name":
        manifest["asset_filename"] = "../outside.faa.gz"
    else:
        manifest["name"] = "other-profile"
    (directory / f"{profile}.faa.gz").write_bytes(compressed)
    (directory / f"{profile}.json").write_text(json.dumps(manifest))
    monkeypatch.setattr("grodon.references.resources.files", lambda _: tmp_path)
    with pytest.raises(ValueError), resolve_reference(profile):
        pass
