"""Explicit, isolated installation and diagnostics for optional genome tools.

Prediction never calls this installer. Downloads use NCBI's versioned HTTPS
release directory, verify its published MD5, and record SHA-256 for the archive
and every installed file. MD5 provides the publisher's integrity check; HTTPS
provides transport authentication. No installer modifies PATH or system files.
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import os
import platform as host_platform
import posixpath
import re
import shutil
import subprocess
import tarfile
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path, PurePosixPath
from typing import Any

DEFAULT_BLAST_VERSION = "2.17.0"
_DOWNLOAD_BASE = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+"
_MANIFEST = "grodon-blast-manifest.json"
_CHUNK_SIZE = 1024 * 1024
_MAX_ARCHIVE_SIZE = 2 * 1024**3
_MAX_INSTALLED_SIZE = 1024**3
_TIMEOUT = 30
_DOWNLOAD_DEADLINE = 900
_VERSION = re.compile(r"\d+\.\d+\.\d+\Z")
_BLAST_VERSION = re.compile(r"(?m)^blastp:\s*(\d+\.\d+\.\d+)\+")
_LIBRARY = re.compile(r".+(?:\.dll|\.dylib|\.so(?:\.[0-9]+)*)\Z", re.IGNORECASE)
_WINDOWS_DEVICES = {"CON", "PRN", "AUX", "NUL"} | {
    f"{prefix}{number}" for prefix in ("COM", "LPT") for number in "123456789¹²³"
}


def _probe_blastp(blastp: str | os.PathLike[str]) -> dict[str, Any]:
    requested = os.fspath(blastp)
    resolved = shutil.which(requested)
    if resolved is None:
        candidate = Path(requested)
        if candidate.is_file():
            resolved = str(candidate.resolve())
    report: dict[str, Any] = {"available": False, "path": resolved, "version": None, "error": None}
    if resolved is None:
        report["error"] = (
            "BLASTP was not found. Run grodon setup --tools-dir DIR or pass --blastp PATH."
        )
        return report
    report["path"] = str(Path(resolved).resolve())
    try:
        process = subprocess.run(
            [report["path"], "-version"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=10,
            check=False,
            shell=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        report["error"] = f"BLASTP could not run: {exc}"
        return report
    match = _BLAST_VERSION.search(process.stdout)
    if process.returncode != 0 or match is None:
        detail = (process.stderr or process.stdout).strip()[:1000]
        report["error"] = (
            f"The executable did not report a working BLASTP version (exit {process.returncode}): {detail}"
        )
        return report
    report.update(available=True, version=match.group(1))
    return report


def doctor(blastp: str | os.PathLike[str] = "blastp") -> dict[str, Any]:
    """Report optional Pyrodigal and BLASTP availability without downloading.

    ``ready`` indicates that both optional genome-pipeline dependencies work.
    Ordinary CDS-based gRodon predictions require neither dependency.
    """
    pyrodigal: dict[str, Any] = {"available": False, "version": None, "error": None}
    try:
        module = importlib.import_module("pyrodigal")
        version = getattr(module, "__version__", None)
        if version is None:
            version = importlib.metadata.version("pyrodigal")
        pyrodigal.update(available=True, version=str(version))
    except Exception as exc:
        pyrodigal["error"] = f"Pyrodigal is unavailable: {exc}. Install grodon2-py[genome]."
    blast_report = _probe_blastp(blastp)
    return {
        "pyrodigal": pyrodigal,
        "blastp": blast_report,
        "ready": bool(pyrodigal["available"] and blast_report["available"]),
    }


def _release_platform(platform: str | None, architecture: str | None) -> tuple[str, str]:
    system = (host_platform.system() if platform is None else platform).lower()
    machine = (host_platform.machine() if architecture is None else architecture).lower()
    systems = {
        "windows": "windows",
        "win32": "windows",
        "linux": "linux",
        "darwin": "macos",
        "macos": "macos",
    }
    machines = {
        "amd64": "x64",
        "x86_64": "x64",
        "x64": "x64",
        "aarch64": "aarch64",
        "arm64": "aarch64",
    }
    system, machine = systems.get(system, system), machines.get(machine, machine)
    builds = {
        ("windows", "x64"): "x64-win64",
        ("linux", "x64"): "x64-linux",
        ("linux", "aarch64"): "aarch64-linux",
        ("macos", "x64"): "x64-macosx",
        ("macos", "aarch64"): "aarch64-macosx",
    }
    try:
        return builds[(system, machine)], "blastp.exe" if system == "windows" else "blastp"
    except KeyError as exc:
        raise ValueError(
            f"No managed NCBI BLAST build for {system}/{machine}. Supported builds are "
            "Windows x64, Linux x64/ARM64, and macOS x64/ARM64. Install BLAST separately and pass --blastp PATH."
        ) from exc


def _download(url: str, destination: Path, *, max_bytes: int) -> tuple[str, str]:
    """Stream a bounded NCBI HTTPS response, returning MD5 and SHA-256."""
    request = urllib.request.Request(
        url, headers={"User-Agent": "grodon2-py/optional-tool-installer"}
    )
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    start = time.monotonic()
    try:
        with urllib.request.urlopen(request, timeout=_TIMEOUT) as response:
            final_url = urllib.parse.urlparse(response.geturl())
            if final_url.scheme != "https" or final_url.hostname != "ftp.ncbi.nlm.nih.gov":
                raise ValueError("Refusing a download redirected outside NCBI HTTPS")
            content_length = response.headers.get("Content-Length")
            if content_length is not None and int(content_length) > max_bytes:
                raise ValueError("NCBI download exceeds the allowed file size")
            written = 0
            with destination.open("xb") as stream:
                while chunk := response.read(min(_CHUNK_SIZE, max_bytes - written + 1)):
                    written += len(chunk)
                    if written > max_bytes:
                        raise ValueError("NCBI download exceeds the allowed file size")
                    if time.monotonic() - start > _DOWNLOAD_DEADLINE:
                        raise TimeoutError("NCBI download exceeded its 15-minute time limit")
                    stream.write(chunk)
                    md5.update(chunk)
                    sha256.update(chunk)
            if content_length is not None and written != int(content_length):
                raise ValueError("NCBI download was truncated (Content-Length mismatch)")
    except (urllib.error.URLError, OSError) as exc:
        raise RuntimeError(f"Could not download {url}: {exc}") from exc
    return md5.hexdigest(), sha256.hexdigest()


def _archive_path(name: str, root: str) -> PurePosixPath:
    if not name or "\\" in name or ":" in name or "\x00" in name:
        raise ValueError(f"Unsafe archive path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or not path.parts or path.parts[0] != root:
        raise ValueError(f"Archive entry leaves its expected directory: {name!r}")
    if any(
        part.endswith((".", " ")) or part.split(".", 1)[0].upper() in _WINDOWS_DEVICES
        for part in path.parts
    ):
        raise ValueError(f"Unsafe archive path under Windows filename rules: {name!r}")
    return path


def _link_target(member: tarfile.TarInfo, root: str) -> PurePosixPath:
    link = member.linkname
    if not link or "\\" in link or ":" in link or "\x00" in link or link.startswith("/"):
        raise ValueError(f"Unsafe archive link: {member.name!r} -> {link!r}")
    joined = posixpath.join(posixpath.dirname(member.name), link) if member.issym() else link
    return _archive_path(posixpath.normpath(joined), root)


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(_CHUNK_SIZE):
            digest.update(block)
    return digest.hexdigest()


def _extract_blast(
    archive: Path, installation: Path, version: str, executable: str
) -> dict[str, Any]:
    """Copy selected regular files; safe internal library links become copies."""
    root = f"ncbi-blast-{version}+"
    installed: dict[str, Any] = {}
    with tarfile.open(archive, mode="r:gz") as package:
        members = package.getmembers()
        by_name: dict[str, tarfile.TarInfo] = {}
        case_names: set[str] = set()
        for member in members:
            path = _archive_path(member.name, root)
            normalized = path.as_posix()
            if normalized.casefold() in case_names:
                raise ValueError(
                    f"Archive contains duplicate or case-colliding paths: {member.name!r}"
                )
            case_names.add(normalized.casefold())
            by_name[normalized] = member
            if member.issym() or member.islnk():
                _link_target(member, root)
            elif not member.isfile() and not member.isdir():
                raise ValueError(f"Unsupported special archive entry: {member.name!r}")

        selected = []
        for name, member in by_name.items():
            relative = PurePosixPath(name).relative_to(root)
            if len(relative.parts) < 2 or relative.parts[0] not in {"bin", "lib"}:
                continue
            if relative.as_posix() == f"bin/{executable}" or _LIBRARY.fullmatch(relative.name):
                if member.isdir():
                    raise ValueError(f"Expected a file, found an archive directory: {name}")
                selected.append((relative, member))
        if not any(path.as_posix() == f"bin/{executable}" for path, _ in selected):
            raise ValueError(f"The NCBI archive does not contain bin/{executable}")

        total = 0
        for relative, member in selected:
            target = member
            visited: set[str] = set()
            while target.issym() or target.islnk():
                link = _link_target(target, root).as_posix()
                if link in visited or link not in by_name:
                    raise ValueError(f"Missing or cyclic archive link target: {link!r}")
                visited.add(link)
                target = by_name[link]
            if not target.isfile() or target.size < 0:
                raise ValueError(f"Archive link does not resolve to a regular file: {member.name}")
            total += target.size
            if total > _MAX_INSTALLED_SIZE:
                raise ValueError("Selected BLAST files exceed the allowed installed size")
            destination = installation.joinpath(*relative.parts)
            destination.parent.mkdir(parents=True, exist_ok=True)
            source = package.extractfile(target)
            if source is None:
                raise ValueError(f"Could not read archive file: {target.name}")
            with source, destination.open("xb") as stream:
                shutil.copyfileobj(source, stream, length=_CHUNK_SIZE)
            if destination.stat().st_size != target.size:
                raise ValueError(f"Truncated archive file: {target.name}")
            destination.chmod(0o755)
            installed[relative.as_posix()] = {
                "size": target.size,
                "sha256": _hash_file(destination),
            }
    return installed


def _manifest_path(name: str) -> PurePosixPath:
    path = _archive_path("installation/" + name, "installation").relative_to("installation")
    if len(path.parts) < 2 or path.parts[0] not in {"bin", "lib"}:
        raise ValueError("Invalid file path in managed BLAST manifest")
    return path


def _verify_existing(
    destination: Path, version: str, build: str, url: str, executable: str
) -> dict[str, Any]:
    manifest_file = destination / _MANIFEST
    if (
        manifest_file.is_symlink()
        or not manifest_file.is_file()
        or manifest_file.stat().st_size > 1024**2
    ):
        raise FileExistsError(
            "Tool directory is nonempty and has no valid managed BLAST manifest; choose an empty directory"
        )
    try:
        manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
        if (
            manifest["schema_version"] != 1
            or manifest["version"] != version
            or manifest["platform"] != build
            or manifest["archive_url"] != url
            or not re.fullmatch(r"[0-9a-f]{32}", manifest["archive_md5"])
            or not re.fullmatch(r"[0-9a-f]{64}", manifest["archive_sha256"])
        ):
            raise ValueError("Manifest does not match the requested BLAST release")
        files = manifest["files"]
        if not isinstance(files, dict) or f"bin/{executable}" not in files:
            raise ValueError("Manifest does not include the required BLASTP executable")
        for name, metadata in files.items():
            relative = _manifest_path(name)
            path = destination.joinpath(*relative.parts)
            if not path.resolve().is_relative_to(destination):
                raise ValueError(
                    "Managed BLAST files must not leave the installation through a junction"
                )
            if any(
                parent.is_symlink()
                for parent in (path, *path.parents)
                if parent != destination.parent
            ):
                raise ValueError("Managed BLAST files must not be symbolic links")
            if (
                not path.is_file()
                or path.stat().st_size != metadata["size"]
                or _hash_file(path) != metadata["sha256"]
            ):
                raise ValueError(f"Managed BLAST file failed its SHA-256 check: {name}")
        entries = []
        pending = [destination]
        visited = {destination}
        while pending:
            for path in pending.pop().iterdir():
                canonical = path.resolve()
                if path.is_symlink() or not canonical.is_relative_to(destination):
                    raise ValueError(
                        "Managed BLAST directories must not contain symbolic links or external junctions"
                    )
                entries.append(path)
                if path.is_dir():
                    if canonical in visited:
                        raise ValueError("Managed BLAST directories contain a directory cycle")
                    visited.add(canonical)
                    pending.append(path)
        expected_directories = {
            parent.as_posix()
            for name in files
            for parent in PurePosixPath(name).parents
            if parent != PurePosixPath(".")
        }
        actual_directories = {
            path.relative_to(destination).as_posix() for path in entries if path.is_dir()
        }
        actual = {path.relative_to(destination).as_posix() for path in entries if not path.is_dir()}
        if actual != set(files) | {_MANIFEST}:
            raise ValueError("Tool directory contains files outside its managed manifest")
        if actual_directories != expected_directories:
            raise ValueError("Tool directory contains directories outside its managed manifest")
    except (KeyError, TypeError, ValueError, OSError) as exc:
        raise ValueError(
            f"Existing BLAST installation failed verification; nothing was overwritten: {exc}"
        ) from exc
    _require_version(destination / "bin" / executable, version)
    return manifest


def _require_version(executable: Path, version: str) -> None:
    report = _probe_blastp(executable)
    if not report["available"] or report["version"] != version:
        raise RuntimeError(
            f"Installed BLASTP did not pass its version check: {report['error'] or report['version']}"
        )


def _installation_result(
    destination: Path, executable: str, manifest: dict[str, Any], reused: bool
) -> dict[str, Any]:
    return {
        "blastp": str(destination / "bin" / executable),
        "version": manifest["version"],
        "platform": manifest["platform"],
        "manifest": str(destination / _MANIFEST),
        "archive_url": manifest["archive_url"],
        "archive_md5": manifest["archive_md5"],
        "archive_sha256": manifest["archive_sha256"],
        "reused": reused,
    }


def install_blast(
    destination: str | os.PathLike[str],
    *,
    version: str = DEFAULT_BLAST_VERSION,
    platform: str | None = None,
    architecture: str | None = None,
) -> dict[str, Any]:
    """Install a pinned official BLASTP release in a dedicated directory.

    The directory must be absent, empty, or contain this installer's matching,
    intact manifest. Reuse verifies every installed file and requires no network.
    ``platform`` and ``architecture`` normally come from this machine; explicit
    overrides support managed environments and tests. Unsupported combinations
    fail before downloading. No administrator access or PATH changes are used.
    """
    if not isinstance(version, str) or not _VERSION.fullmatch(version):
        raise ValueError("version must be an exact numeric release such as '2.17.0'")
    build, executable = _release_platform(platform, architecture)
    requested = Path(destination).expanduser()
    if requested.is_symlink():
        raise ValueError("The installation directory must not be a symbolic link")
    target = requested.resolve()
    archive_name = f"ncbi-blast-{version}+-{build}.tar.gz"
    url = f"{_DOWNLOAD_BASE}/{version}/{archive_name}"
    if target.exists():
        if not target.is_dir():
            raise FileExistsError("The installation destination already exists as a file")
        if any(target.iterdir()):
            manifest = _verify_existing(target, version, build, url, executable)
            return _installation_result(target, executable, manifest, True)
    target.parent.mkdir(parents=True, exist_ok=True)
    lock = target.parent / f".{target.name}.grodon-install.lock"
    try:
        lock_stream = lock.open("x", encoding="ascii")
    except FileExistsError as exc:
        raise FileExistsError(
            f"Another installation may be using this directory (lock {lock}); nothing was overwritten"
        ) from exc
    try:
        with lock_stream:
            lock_stream.write(str(os.getpid()))
        with tempfile.TemporaryDirectory(
            prefix=f".{target.name}.grodon-", dir=target.parent
        ) as temporary:
            staging = Path(temporary)
            checksum_path, archive_path = staging / "release.md5", staging / "release.tar.gz"
            _download(url + ".md5", checksum_path, max_bytes=4096)
            checksum_text = checksum_path.read_text(encoding="ascii").strip()
            match = re.fullmatch(r"([0-9a-fA-F]{32})\s+\*?([^\s]+)", checksum_text)
            if match is None or match.group(2) != archive_name:
                raise ValueError("NCBI checksum file does not identify the requested archive")
            expected_md5 = match.group(1).lower()
            archive_md5, archive_sha256 = _download(url, archive_path, max_bytes=_MAX_ARCHIVE_SIZE)
            if archive_md5 != expected_md5:
                raise ValueError("NCBI archive MD5 verification failed; no tools were installed")
            installation = staging / "installation"
            installation.mkdir()
            files = _extract_blast(archive_path, installation, version, executable)
            _require_version(installation / "bin" / executable, version)
            manifest = {
                "schema_version": 1,
                "version": version,
                "platform": build,
                "archive_url": url,
                "checksum_url": url + ".md5",
                "archive_md5": archive_md5,
                "archive_sha256": archive_sha256,
                "files": files,
            }
            (installation / _MANIFEST).write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            if target.is_symlink() or (
                target.exists() and (not target.is_dir() or any(target.iterdir()))
            ):
                raise FileExistsError(
                    "The installation destination changed while downloading; nothing was overwritten"
                )
            if target.exists():
                target.rmdir()  # Only the explicitly requested, verified-empty directory.
            installation.rename(target)
        return _installation_result(target, executable, manifest, False)
    finally:
        lock.unlink(missing_ok=True)
