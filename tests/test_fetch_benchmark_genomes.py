"""Downloader contract tests using synthetic bytes and a mocked HTTPS transport."""

from __future__ import annotations

import gzip
import hashlib
import http.client
import importlib.util
import io
import json
import urllib.request
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

ACCESSION = "GCF_000000001.1"
URL = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{ACCESSION}_synthetic_genomic.fna.gz"
CONTENT = gzip.compress(b">synthetic\nATGAAATAA\n", mtime=0)


@pytest.fixture(scope="module")
def fetcher():
    path = Path(__file__).resolve().parents[1] / "scripts" / "fetch_benchmark_genomes.py"
    spec = importlib.util.spec_from_file_location("benchmark_genome_fetcher", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def entry():
    return {
        "accession": ACCESSION,
        "genome_sha256": hashlib.sha256(CONTENT).hexdigest(),
        "download": {"genome_url": URL},
    }


def write_manifest(tmp_path, entries):
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps({"entries": entries}), encoding="utf-8")
    return path


def target_in(folder, accession=ACCESSION):
    return folder / f"{accession}_benchmark_genomic.fna.gz"


@pytest.fixture
def no_network(fetcher, monkeypatch):
    mock = Mock(side_effect=AssertionError("Network access must be explicitly requested"))
    monkeypatch.setattr(fetcher.urllib.request, "build_opener", mock)
    return mock


class Response(io.BytesIO):
    def __init__(self, content=CONTENT, *, url=URL, length=None):
        super().__init__(content)
        self.url = url
        self.headers = {} if length is None else {"Content-Length": str(length)}

    def geturl(self):
        return self.url


def install_transport(fetcher, monkeypatch, **response_options):
    opened = Mock(side_effect=lambda *args, **kwargs: Response(**response_options))
    build = Mock(return_value=SimpleNamespace(open=opened))
    monkeypatch.setattr(fetcher.urllib.request, "build_opener", build)
    return build, opened


def test_valid_existing_archive_is_verified_without_network(fetcher, tmp_path, entry, no_network):
    manifest = write_manifest(tmp_path, [entry])
    output = tmp_path / "downloads"
    output.mkdir()
    target = target_in(output)
    target.write_bytes(CONTENT)
    result = fetcher.fetch_benchmark_genomes(manifest, output)
    assert result["ready"] == 1 and result["errors"] == 0 and result["downloaded"] == 0
    assert result["results"][0]["path"] == str(target)
    assert target.read_bytes() == CONTENT
    no_network.assert_not_called()


def test_missing_offline_cache_reports_error_without_creating_directory(
    fetcher, tmp_path, entry, no_network
):
    manifest = write_manifest(tmp_path, [entry])
    output = tmp_path / "absent"
    result = fetcher.fetch_benchmark_genomes(manifest, output)
    assert result["ready"] == 0 and result["errors"] == 1
    assert "--fetch" in result["results"][0]["error"]
    assert not output.exists()
    no_network.assert_not_called()


@pytest.mark.parametrize("fetch", [False, True])
def test_existing_mismatched_archive_is_preserved_even_when_fetch_is_enabled(
    fetcher, tmp_path, entry, no_network, fetch
):
    manifest = write_manifest(tmp_path, [entry])
    output = tmp_path / "downloads"
    output.mkdir()
    target = target_in(output)
    target.write_bytes(b"valuable existing bytes")
    result = fetcher.fetch_benchmark_genomes(manifest, output, fetch=fetch)
    assert result["errors"] == 1 and result["downloaded"] == 0
    assert "SHA256 mismatch; file preserved" in result["results"][0]["error"]
    assert target.read_bytes() == b"valuable existing bytes"
    assert list(output.iterdir()) == [target]
    no_network.assert_not_called()


def test_fetch_uses_exact_archived_url_timeout_hash_and_compatible_safe_filename(
    fetcher, monkeypatch, tmp_path, entry
):
    manifest = write_manifest(tmp_path, [entry])
    output = tmp_path / "downloads"
    build, opened = install_transport(fetcher, monkeypatch, length=len(CONTENT))
    result = fetcher.fetch_benchmark_genomes(manifest, output, fetch=True, timeout=7)
    assert result["ready"] == 1 and result["errors"] == 0 and result["downloaded"] == 1
    target = target_in(output)
    assert target.read_bytes() == CONTENT
    assert list(output.iterdir()) == [target]
    request = opened.call_args.args[0]
    assert request.full_url == entry["download"]["genome_url"]
    assert opened.call_args.kwargs == {"timeout": 7.0}
    assert isinstance(build.call_args.args[0], fetcher._NcbiRedirectHandler)
    # Matches prepare_folder_genomes.py's required accession/label/genomic suffix.
    assert target.name == "GCF_000000001.1_benchmark_genomic.fna.gz"


def test_checksum_failure_removes_temporary_archive_and_never_creates_cache(
    fetcher, monkeypatch, tmp_path, entry
):
    install_transport(fetcher, monkeypatch, content=b"wrong archive bytes")
    result = fetcher.fetch_benchmark_genomes(
        write_manifest(tmp_path, [entry]), tmp_path / "downloads", fetch=True
    )
    assert result["ready"] == 0 and result["errors"] == 1
    assert "Downloaded archive SHA256 mismatch" in result["results"][0]["error"]
    assert not list((tmp_path / "downloads").iterdir())


@pytest.mark.parametrize("declared", [True, False])
def test_download_size_bound_checks_both_header_and_actual_stream(
    fetcher, monkeypatch, tmp_path, entry, declared
):
    monkeypatch.setattr(fetcher, "MAX_BYTES", len(CONTENT) - 1)
    install_transport(fetcher, monkeypatch, length=len(CONTENT) if declared else None)
    output = tmp_path / "downloads"
    result = fetcher.fetch_benchmark_genomes(write_manifest(tmp_path, [entry]), output, fetch=True)
    assert result["errors"] == 1
    assert "download limit" in result["results"][0]["error"]
    assert not list(output.iterdir())


def test_cache_publication_cannot_overwrite_file_created_during_download(
    fetcher, monkeypatch, tmp_path, entry
):
    install_transport(fetcher, monkeypatch)
    output = tmp_path / "downloads"

    def competing_writer(source, target):
        target.write_bytes(b"another process created this")
        raise FileExistsError("destination appeared")

    monkeypatch.setattr(fetcher.os, "link", competing_writer)
    result = fetcher.fetch_benchmark_genomes(write_manifest(tmp_path, [entry]), output, fetch=True)
    assert result["errors"] == 1 and result["ready"] == 0
    assert target_in(output).read_bytes() == b"another process created this"
    assert list(output.iterdir()) == [target_in(output)]


@pytest.mark.parametrize(
    "failure", [TimeoutError("timed out"), http.client.IncompleteRead(b"part")]
)
def test_stream_failure_removes_partial_file_and_is_reported(
    fetcher, monkeypatch, tmp_path, entry, failure
):
    response = Response()
    response.read = Mock(side_effect=[CONTENT[:5], failure])
    monkeypatch.setattr(
        fetcher.urllib.request,
        "build_opener",
        lambda *args: SimpleNamespace(open=lambda *args, **kwargs: response),
    )
    output = tmp_path / "downloads"
    result = fetcher.fetch_benchmark_genomes(write_manifest(tmp_path, [entry]), output, fetch=True)
    assert result["errors"] == 1
    assert not list(output.iterdir())


def test_failure_for_one_entry_does_not_hide_other_verified_inputs(
    fetcher, tmp_path, entry, no_network
):
    other = deepcopy(entry)
    other["accession"] = "GCA_000000002.1"
    other["download"]["genome_url"] = URL.replace(ACCESSION, other["accession"])
    output = tmp_path / "downloads"
    output.mkdir()
    target_in(output, other["accession"]).write_bytes(CONTENT)
    result = fetcher.fetch_benchmark_genomes(write_manifest(tmp_path, [entry, other]), output)
    assert (result["requested"], result["ready"], result["errors"]) == (2, 1, 1)
    assert [row["status"] for row in result["results"]] == ["error", "ready"]
    no_network.assert_not_called()


@pytest.mark.parametrize(
    "url",
    [
        URL.replace("https:", "http:"),
        URL.replace("ftp.ncbi.nlm.nih.gov", "example.com"),
        URL.replace("ftp.ncbi.nlm.nih.gov", "ftp.ncbi.nlm.nih.gov.example.com"),
        URL.replace("ftp.ncbi.nlm.nih.gov", "user@ftp.ncbi.nlm.nih.gov"),
        URL.replace("ftp.ncbi.nlm.nih.gov", "ftp.ncbi.nlm.nih.gov:8443"),
        URL.replace("_genomic.fna.gz", "_protein.faa.gz"),
        URL + "#fragment",
        URL + "?redirect=elsewhere",
        URL.replace(ACCESSION, "GCF_000000002.1"),
    ],
)
def test_invalid_remote_source_is_rejected_before_output_or_network(
    fetcher, tmp_path, entry, no_network, url
):
    entry["download"]["genome_url"] = url
    output = tmp_path / "downloads"
    with pytest.raises(ValueError):
        fetcher.fetch_benchmark_genomes(write_manifest(tmp_path, [entry]), output, fetch=True)
    assert not output.exists()
    no_network.assert_not_called()


@pytest.mark.parametrize("accession", ["../GCF_000000001.1", "GCF_000000001", "C:/escape", None])
def test_bad_accession_cannot_escape_output_directory(
    fetcher, tmp_path, entry, no_network, accession
):
    entry["accession"] = accession
    output = tmp_path / "downloads"
    with pytest.raises(ValueError, match="accession"):
        fetcher.fetch_benchmark_genomes(write_manifest(tmp_path, [entry]), output, fetch=True)
    assert not output.exists()
    no_network.assert_not_called()


@pytest.mark.parametrize("checksum", [None, "", "0" * 63, "z" * 64])
def test_missing_or_malformed_checksum_fails_before_output(
    fetcher, tmp_path, entry, no_network, checksum
):
    entry["genome_sha256"] = checksum
    output = tmp_path / "downloads"
    with pytest.raises(ValueError, match="SHA256"):
        fetcher.fetch_benchmark_genomes(write_manifest(tmp_path, [entry]), output, fetch=True)
    assert not output.exists()
    no_network.assert_not_called()


def test_duplicate_accessions_reject_entire_manifest_before_fetch(
    fetcher, tmp_path, entry, no_network
):
    output = tmp_path / "downloads"
    with pytest.raises(ValueError, match="Duplicate"):
        fetcher.fetch_benchmark_genomes(
            write_manifest(tmp_path, [entry, deepcopy(entry)]), output, fetch=True
        )
    assert not output.exists()
    no_network.assert_not_called()


def test_redirect_to_untrusted_destination_is_rejected_before_following(fetcher):
    handler = fetcher._NcbiRedirectHandler()
    request = urllib.request.Request(URL)
    with pytest.raises(ValueError, match="NCBI HTTPS"):
        handler.redirect_request(
            request, None, 302, "Found", {}, "https://example.com/file_genomic.fna.gz"
        )


def test_final_response_url_is_also_checked(fetcher, monkeypatch, tmp_path, entry):
    install_transport(fetcher, monkeypatch, url="https://example.com/file_genomic.fna.gz")
    output = tmp_path / "downloads"
    result = fetcher.fetch_benchmark_genomes(write_manifest(tmp_path, [entry]), output, fetch=True)
    assert result["errors"] == 1
    assert not list(output.iterdir())


@pytest.mark.parametrize("timeout", [0, -1, float("nan"), float("inf"), True, "60"])
def test_invalid_timeout_fails_before_loading_or_writing(fetcher, tmp_path, timeout):
    with pytest.raises(ValueError, match="timeout"):
        fetcher.fetch_benchmark_genomes(tmp_path / "absent.json", tmp_path / "out", timeout=timeout)
    assert not (tmp_path / "out").exists()


def test_cli_returns_nonzero_with_ready_error_totals_for_missing_offline_cache(
    fetcher, tmp_path, entry, no_network, capsys
):
    manifest = write_manifest(tmp_path, [entry])
    code = fetcher.main(["--manifest", str(manifest), "--output-dir", str(tmp_path / "downloads")])
    captured = capsys.readouterr()
    assert code == 1
    assert json.loads(captured.out) == {"requested": 1, "ready": 0, "errors": 1, "downloaded": 0}
    assert "--fetch" in captured.err
    no_network.assert_not_called()


def test_bundled_manifest_exposes_all_forty_public_sources_without_private_files(fetcher):
    manifest = (
        Path(__file__).resolve().parents[1] / "validation" / "test_genomes_inputs_manifest.json"
    )
    entries = fetcher.load_manifest(manifest)
    assert len(entries) == 40
    assert len({entry["accession"] for entry in entries}) == 40
    assert all(set(entry) == {"accession", "url", "sha256"} for entry in entries)
