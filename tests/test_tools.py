"""Dependency diagnostics and bounded, verified tool installation."""

from __future__ import annotations

import hashlib
import io
import json
import subprocess
import tarfile
import types
import urllib.error
from pathlib import Path

import pytest

import grodon.tools as tools

VERSION = "2.17.0"
ROOT = f"ncbi-blast-{VERSION}+"
ARCHIVE = f"{ROOT}-x64-win64.tar.gz"


def archive_bytes(entries=None):
    if entries is None:
        entries = [
            (f"{ROOT}/bin/blastp.exe", b"fake executable"),
            (f"{ROOT}/bin/ncbi-vdb-md.dll", b"fake library"),
            (f"{ROOT}/bin/makeblastdb.exe", b"must not install"),
            (f"{ROOT}/doc/README", b"must not install"),
        ]
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w:gz") as package:
        for name, data in entries:
            if isinstance(data, tarfile.TarInfo):
                package.addfile(data)
            else:
                member = tarfile.TarInfo(name)
                member.size = len(data)
                package.addfile(member, io.BytesIO(data))
    return stream.getvalue()


class FakeResponse(io.BytesIO):
    def __init__(self, data, url, *, content_length=None):
        super().__init__(data)
        self.url = url
        self.headers = {
            "Content-Length": str(len(data) if content_length is None else content_length)
        }
        self.read_sizes = []

    def geturl(self):
        return self.url

    def read(self, size=-1):
        self.read_sizes.append(size)
        return super().read(size)


@pytest.fixture
def fake_network(monkeypatch):
    state = {"archive": archive_bytes(), "checksum": None, "calls": []}

    def open_url(request, timeout):
        assert timeout == 30
        url = request.full_url
        assert url.startswith("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/")
        assert "/LATEST/" not in url
        state["calls"].append(url)
        if url.endswith(".md5"):
            checksum = state["checksum"] or hashlib.md5(state["archive"]).hexdigest()
            body = f"{checksum}  {ARCHIVE}\n".encode()
        else:
            body = state["archive"]
        return FakeResponse(body, url)

    monkeypatch.setattr(tools.urllib.request, "urlopen", open_url)
    monkeypatch.setattr(
        tools,
        "_probe_blastp",
        lambda path: {
            "available": True,
            "path": str(path),
            "version": VERSION,
            "error": None,
        },
    )
    return state


def install(path):
    return tools.install_blast(path, platform="Windows", architecture="AMD64")


def test_installs_only_blastp_and_libraries_with_verifiable_manifest(tmp_path, fake_network):
    destination = tmp_path / "tools"
    result = install(destination)
    assert result["blastp"] == str(destination / "bin" / "blastp.exe")
    assert result["version"] == VERSION
    assert result["platform"] == "x64-win64"
    assert result["reused"] is False
    assert sorted(path.name for path in (destination / "bin").iterdir()) == [
        "blastp.exe",
        "ncbi-vdb-md.dll",
    ]
    manifest = json.loads(Path(result["manifest"]).read_text())
    assert manifest["archive_md5"] == hashlib.md5(fake_network["archive"]).hexdigest()
    assert manifest["archive_sha256"] == hashlib.sha256(fake_network["archive"]).hexdigest()
    for relative, metadata in manifest["files"].items():
        contents = (destination / relative).read_bytes()
        assert metadata == {"size": len(contents), "sha256": hashlib.sha256(contents).hexdigest()}
    assert len(fake_network["calls"]) == 2
    assert not list(tmp_path.glob(".tools.grodon-*"))


def test_verified_reuse_is_offline(tmp_path, fake_network):
    first = install(tmp_path / "tools")
    calls = list(fake_network["calls"])
    second = install(tmp_path / "tools")
    assert first["blastp"] == second["blastp"]
    assert second["reused"] is True
    assert fake_network["calls"] == calls


def test_modified_executable_cannot_be_reused_or_overwritten(tmp_path, fake_network):
    result = install(tmp_path / "tools")
    binary = Path(result["blastp"])
    binary.write_bytes(b"tampered")
    with pytest.raises(ValueError, match="SHA-256"):
        install(tmp_path / "tools")
    assert binary.read_bytes() == b"tampered"
    assert len(fake_network["calls"]) == 2


def test_nonempty_destination_is_untouched(tmp_path, fake_network):
    destination = tmp_path / "tools"
    destination.mkdir()
    preserved = destination / "personal.txt"
    preserved.write_text("keep this")
    with pytest.raises(FileExistsError, match="nonempty"):
        install(destination)
    assert preserved.read_text() == "keep this"
    assert fake_network["calls"] == []


def test_existing_empty_destination_is_allowed(tmp_path, fake_network):
    destination = tmp_path / "tools"
    destination.mkdir()
    assert install(destination)["reused"] is False


def test_destination_file_is_not_overwritten(tmp_path, fake_network):
    destination = tmp_path / "tools"
    destination.write_text("keep")
    with pytest.raises(FileExistsError, match="file"):
        install(destination)
    assert destination.read_text() == "keep"
    assert not fake_network["calls"]


def test_wrong_checksum_never_extracts(tmp_path, fake_network):
    fake_network["checksum"] = "0" * 32
    with pytest.raises(ValueError, match="MD5 verification"):
        install(tmp_path / "tools")
    assert not (tmp_path / "tools").exists()
    assert not list(tmp_path.glob(".tools.grodon-*"))


@pytest.mark.parametrize(
    "name",
    [
        "../outside.txt",
        "/absolute.txt",
        f"{ROOT}/../../outside.txt",
        f"{ROOT}\\bin\\blastp.exe",
        f"{ROOT}/bin/C:outside.dll",
        "other-root/bin/blastp.exe",
        f"{ROOT}/bin/.. /.. /outside.dll",
        f"{ROOT}/bin/NUL.dll",
        f"{ROOT}/bin/con.dll",
        f"{ROOT}/bin/library.dll.",
    ],
)
def test_archive_path_traversal_rejected_even_for_unused_entries(tmp_path, fake_network, name):
    fake_network["archive"] = archive_bytes([(f"{ROOT}/bin/blastp.exe", b"ok"), (name, b"unsafe")])
    with pytest.raises(ValueError, match="archive path|expected directory"):
        install(tmp_path / "tools")
    assert not (tmp_path / "tools").exists()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("kind", [tarfile.SYMTYPE, tarfile.LNKTYPE])
def test_links_cannot_escape_installation(tmp_path, fake_network, kind):
    link = tarfile.TarInfo(f"{ROOT}/bin/evil.dll")
    link.type, link.linkname = kind, "../../../outside"
    fake_network["archive"] = archive_bytes([(f"{ROOT}/bin/blastp.exe", b"ok"), (link.name, link)])
    with pytest.raises(ValueError, match="expected directory"):
        install(tmp_path / "tools")
    assert not (tmp_path / "tools").exists()


def test_internal_library_symlink_becomes_regular_file_copy(tmp_path, fake_network):
    link = tarfile.TarInfo(f"{ROOT}/bin/alias.dll")
    link.type, link.linkname = tarfile.SYMTYPE, "library.dll"
    fake_network["archive"] = archive_bytes(
        [
            (f"{ROOT}/bin/blastp.exe", b"ok"),
            (f"{ROOT}/bin/library.dll", b"library"),
            (link.name, link),
        ]
    )
    install(tmp_path / "tools")
    alias = tmp_path / "tools" / "bin" / "alias.dll"
    assert alias.read_bytes() == b"library"
    assert not alias.is_symlink()


def test_archive_duplicate_and_case_collision_are_rejected(tmp_path, fake_network):
    fake_network["archive"] = archive_bytes(
        [
            (f"{ROOT}/bin/blastp.exe", b"ok"),
            (f"{ROOT}/bin/BLASTP.EXE", b"overwrite"),
        ]
    )
    with pytest.raises(ValueError, match="duplicate or case-colliding"):
        install(tmp_path / "tools")


def test_archive_special_device_rejected(tmp_path, fake_network):
    device = tarfile.TarInfo(f"{ROOT}/device")
    device.type = tarfile.CHRTYPE
    fake_network["archive"] = archive_bytes(
        [(f"{ROOT}/bin/blastp.exe", b"ok"), (device.name, device)]
    )
    with pytest.raises(ValueError, match="special archive"):
        install(tmp_path / "tools")


def test_missing_binary_and_version_mismatch_do_not_install(tmp_path, fake_network, monkeypatch):
    fake_network["archive"] = archive_bytes([(f"{ROOT}/README", b"empty")])
    with pytest.raises(ValueError, match="does not contain"):
        install(tmp_path / "tools")
    fake_network["archive"] = archive_bytes()
    monkeypatch.setattr(
        tools, "_probe_blastp", lambda path: {"available": True, "version": "2.16.0", "error": None}
    )
    with pytest.raises(RuntimeError, match="version check"):
        install(tmp_path / "tools")
    assert not (tmp_path / "tools").exists()


def test_destination_created_during_install_is_not_overwritten(tmp_path, fake_network, monkeypatch):
    destination = tmp_path / "tools"

    def finish_version_check(*args):
        destination.mkdir()
        (destination / "personal.txt").write_text("keep")

    monkeypatch.setattr(tools, "_require_version", finish_version_check)
    with pytest.raises(FileExistsError, match="changed while downloading"):
        install(destination)
    assert (destination / "personal.txt").read_text() == "keep"


def test_extra_files_and_manifest_path_traversal_prevent_reuse(tmp_path, fake_network):
    destination = tmp_path / "tools"
    result = install(destination)
    extra = destination / "personal.txt"
    extra.write_text("keep")
    with pytest.raises(ValueError, match="outside its managed manifest"):
        install(destination)
    extra.unlink()
    manifest_path = Path(result["manifest"])
    manifest = json.loads(manifest_path.read_text())
    manifest["files"]["../outside"] = {"size": 0, "sha256": "0" * 64}
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="verification"):
        install(destination)


def test_extra_empty_directory_is_not_claimed_by_managed_installation(tmp_path, fake_network):
    destination = tmp_path / "tools"
    install(destination)
    (destination / "personal").mkdir()
    with pytest.raises(ValueError, match="directories outside"):
        install(destination)
    assert (destination / "personal").is_dir()


def test_existing_installation_cannot_escape_through_a_windows_junction(
    tmp_path, fake_network, monkeypatch
):
    destination = tmp_path / "tools"
    result = install(destination)
    binary = Path(result["blastp"])
    original_resolve = Path.resolve

    def resolve(path, *args, **kwargs):
        if path == binary:
            return tmp_path / "outside" / "blastp.exe"
        return original_resolve(path, *args, **kwargs)

    monkeypatch.setattr(Path, "resolve", resolve)
    with pytest.raises(ValueError, match="junction"):
        install(destination)
    assert binary.read_bytes() == b"fake executable"


def test_install_lock_is_not_overwritten(tmp_path, fake_network):
    lock = tmp_path / ".tools.grodon-install.lock"
    lock.write_text("another install")
    with pytest.raises(FileExistsError, match="Another installation"):
        install(tmp_path / "tools")
    assert lock.read_text() == "another install"
    assert not fake_network["calls"]


@pytest.mark.parametrize(
    "version", ["latest", "2.17", "../2.17.0", "2.17.0/../../x", "2.17.0+", None]
)
def test_exact_version_required_before_network(tmp_path, fake_network, version):
    with pytest.raises(ValueError, match="exact numeric"):
        tools.install_blast(tmp_path / "tools", version=version)
    assert not fake_network["calls"]


@pytest.mark.parametrize(
    "system,arch,expected,binary",
    [
        ("Windows", "AMD64", "x64-win64", "blastp.exe"),
        ("Linux", "x86_64", "x64-linux", "blastp"),
        ("Linux", "aarch64", "aarch64-linux", "blastp"),
        ("Darwin", "x86_64", "x64-macosx", "blastp"),
        ("Darwin", "arm64", "aarch64-macosx", "blastp"),
    ],
)
def test_supported_platforms(system, arch, expected, binary):
    assert tools._release_platform(system, arch) == (expected, binary)


@pytest.mark.parametrize(
    "system,arch", [("Windows", "ARM64"), ("Linux", "i386"), ("FreeBSD", "x86_64")]
)
def test_unsupported_platform_is_helpful(system, arch):
    with pytest.raises(ValueError, match="Install BLAST separately"):
        tools._release_platform(system, arch)


def test_download_streams_with_timeout_and_checks_length(tmp_path, monkeypatch):
    body = b"a" * (tools._CHUNK_SIZE + 10)
    url = "https://ftp.ncbi.nlm.nih.gov/test"
    response = FakeResponse(body, url)
    monkeypatch.setattr(tools.urllib.request, "urlopen", lambda request, timeout: response)
    md5, sha256 = tools._download(url, tmp_path / "download", max_bytes=len(body))
    assert md5 == hashlib.md5(body).hexdigest()
    assert sha256 == hashlib.sha256(body).hexdigest()
    assert all(0 < size <= tools._CHUNK_SIZE for size in response.read_sizes)
    truncated = FakeResponse(b"too short", url, content_length=50)
    monkeypatch.setattr(tools.urllib.request, "urlopen", lambda request, timeout: truncated)
    with pytest.raises(ValueError, match="truncated"):
        tools._download(url, tmp_path / "truncated", max_bytes=100)


def test_download_rejects_redirects_limits_and_network_errors(tmp_path, monkeypatch):
    url = "https://ftp.ncbi.nlm.nih.gov/test"
    monkeypatch.setattr(
        tools.urllib.request,
        "urlopen",
        lambda request, timeout: FakeResponse(b"x", "http://example.com/file"),
    )
    with pytest.raises(ValueError, match="redirected outside"):
        tools._download(url, tmp_path / "redirect", max_bytes=10)
    monkeypatch.setattr(
        tools.urllib.request, "urlopen", lambda request, timeout: FakeResponse(b"too large", url)
    )
    with pytest.raises(ValueError, match="allowed file size"):
        tools._download(url, tmp_path / "large", max_bytes=2)

    def missing(request, timeout):
        raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(tools.urllib.request, "urlopen", missing)
    with pytest.raises(RuntimeError, match="Could not download.*404"):
        tools._download(url, tmp_path / "missing", max_bytes=10)


def test_doctor_reports_both_dependencies_without_network(monkeypatch):
    monkeypatch.setattr(
        tools.importlib, "import_module", lambda name: types.SimpleNamespace(__version__="3.6.3")
    )
    monkeypatch.setattr(tools.shutil, "which", lambda executable: "/tools/blastp")
    calls = []

    def run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, 0, "blastp: 2.17.0+\nPackage: blast\n", "")

    monkeypatch.setattr(tools.subprocess, "run", run)
    result = tools.doctor("some-blastp")
    assert result["ready"] is True
    assert result["pyrodigal"]["version"] == "3.6.3"
    assert result["blastp"]["version"] == "2.17.0"
    assert calls[0][0][-1] == "-version"
    assert calls[0][1]["shell"] is False
    assert calls[0][1]["timeout"] == 10


def test_doctor_reports_missing_dependencies(monkeypatch):
    def unavailable(name):
        raise ImportError("missing")

    monkeypatch.setattr(tools.importlib, "import_module", unavailable)
    monkeypatch.setattr(tools.shutil, "which", lambda executable: None)
    result = tools.doctor("this-binary-does-not-exist")
    assert result["ready"] is False
    assert result["blastp"]["available"] is False
    assert "grodon setup" in result["blastp"]["error"]
    assert "[genome]" in result["pyrodigal"]["error"]


@pytest.mark.parametrize("failure", ["timeout", "error", "wrong-program"])
def test_doctor_reports_unusable_executable(monkeypatch, failure):
    monkeypatch.setattr(tools.shutil, "which", lambda executable: "/tools/blastp")

    def run(command, **kwargs):
        if failure == "timeout":
            raise subprocess.TimeoutExpired(command, 10)
        if failure == "error":
            raise OSError("cannot start")
        return subprocess.CompletedProcess(command, 0, "some other program", "")

    monkeypatch.setattr(tools.subprocess, "run", run)
    assert tools._probe_blastp("blastp")["available"] is False
