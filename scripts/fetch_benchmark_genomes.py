"""Fetch the public benchmark inputs using archived NCBI URLs and SHA256 hashes.

Without --fetch this utility only verifies existing files. It downloads genomic
FASTA archives, never proteins or inferred sequences. Existing files with a
different checksum are preserved and reported as errors. No model code is run.
"""

from __future__ import annotations

import argparse
import hashlib
import http.client
import json
import math
import os
import re
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

MAX_BYTES = 100_000_000
CHUNK_BYTES = 1024 * 1024


def _validate_url(url: Any) -> str:
    if not isinstance(url, str) or not url:
        raise ValueError("A genomic FASTA download URL is required")
    parsed = urllib.parse.urlsplit(url)
    host = parsed.hostname or ""
    if (
        parsed.scheme != "https"
        or not (host == "ncbi.nlm.nih.gov" or host.endswith(".ncbi.nlm.nih.gov"))
        or parsed.username is not None
        or parsed.password is not None
        or parsed.port not in (None, 443)
        or parsed.query
        or parsed.fragment
        or not parsed.path.endswith("_genomic.fna.gz")
    ):
        raise ValueError(f"Expected an NCBI HTTPS genomic FASTA archive URL: {url!r}")
    return url


class _NcbiRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        # Check the next destination before urllib sends any redirected request.
        _validate_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def load_manifest(path: Path) -> list[dict[str, str]]:
    """Validate all identities and locations before network or filesystem writes."""
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    entries = value.get("entries") if isinstance(value, dict) else None
    if not isinstance(entries, list) or not entries:
        raise ValueError("The manifest must contain a nonempty entries list")
    selected: list[dict[str, str]] = []
    seen: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Every manifest entry must be an object")
        accession = entry.get("accession")
        if not isinstance(accession, str) or not re.fullmatch(r"GC[AF]_\d{9}\.\d+", accession):
            raise ValueError(f"Invalid assembly accession: {accession!r}")
        if accession in seen:
            raise ValueError(f"Duplicate assembly accession: {accession}")
        seen.add(accession)
        checksum = entry.get("genome_sha256")
        if not isinstance(checksum, str) or not re.fullmatch(r"[a-fA-F0-9]{64}", checksum):
            raise ValueError(f"Missing or malformed genome SHA256: {accession}")
        download = entry.get("download")
        url = _validate_url(download.get("genome_url") if isinstance(download, dict) else None)
        if not urllib.parse.urlsplit(url).path.rsplit("/", 1)[-1].startswith(accession + "_"):
            raise ValueError(f"URL filename does not identify the exact assembly: {accession}")
        selected.append({"accession": accession, "url": url, "sha256": checksum.lower()})
    return selected


def _digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _fetch_one(
    entry: dict[str, str], folder: Path, *, fetch: bool, timeout: float
) -> dict[str, Any]:
    target = folder / f"{entry['accession']}_benchmark_genomic.fna.gz"
    if target.is_symlink():
        raise ValueError(f"Refusing a symlink cache entry: {target}")
    if target.exists():
        if not target.is_file():
            raise ValueError(f"Cache target is not a regular file: {target}")
        if target.stat().st_size > MAX_BYTES or _digest(target) != entry["sha256"]:
            raise ValueError(f"Existing cache SHA256 mismatch; file preserved: {target}")
        return {"status": "ready", "path": str(target), "downloaded": False}
    if not fetch:
        raise FileNotFoundError(f"Missing verified cache: {target.name}; use --fetch to download")
    folder.mkdir(parents=True, exist_ok=True)
    opener = urllib.request.build_opener(_NcbiRedirectHandler())
    request = urllib.request.Request(
        entry["url"], headers={"User-Agent": "grodon2-py/public-benchmark-fetcher"}
    )
    temporary: Path | None = None
    try:
        with opener.open(request, timeout=timeout) as response:
            _validate_url(response.geturl())
            declared = response.headers.get("Content-Length")
            if declared is not None and (int(declared) < 0 or int(declared) > MAX_BYTES):
                raise ValueError(f"Archive exceeds the {MAX_BYTES}-byte download limit")
            with tempfile.NamedTemporaryFile(
                mode="wb", prefix=target.name + ".", suffix=".part", dir=folder, delete=False
            ) as stream:
                temporary = Path(stream.name)
                checksum = hashlib.sha256()
                size = 0
                while chunk := response.read(CHUNK_BYTES):
                    size += len(chunk)
                    if size > MAX_BYTES:
                        raise ValueError(f"Archive exceeds the {MAX_BYTES}-byte download limit")
                    checksum.update(chunk)
                    stream.write(chunk)
            if checksum.hexdigest() != entry["sha256"]:
                raise ValueError(f"Downloaded archive SHA256 mismatch: {entry['accession']}")
        # A hard link publishes the complete verified file atomically and fails
        # if another process created the destination. Unlike replace/rename,
        # this cannot overwrite a cache file that appeared during the download.
        assert temporary is not None
        os.link(temporary, target)
        return {"status": "ready", "path": str(target), "downloaded": True}
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def fetch_benchmark_genomes(
    manifest: Path, output_dir: Path, *, fetch: bool = False, timeout: float = 60
) -> dict[str, Any]:
    """Verify/download the manifest's archives and retain every per-entry outcome."""
    if (
        isinstance(timeout, bool)
        or not isinstance(timeout, (int, float))
        or not math.isfinite(timeout)
        or timeout <= 0
    ):
        raise ValueError("timeout must be a positive finite number")
    if not isinstance(fetch, bool):
        raise ValueError("fetch must be a boolean")
    entries = load_manifest(manifest)
    output = output_dir.resolve()
    if output.exists() and not output.is_dir():
        raise ValueError(f"Output directory is not a directory: {output}")
    results: list[dict[str, Any]] = []
    for entry in entries:
        record: dict[str, Any] = {"accession": entry["accession"]}
        try:
            record.update(_fetch_one(entry, output, fetch=fetch, timeout=float(timeout)))
        except (OSError, ValueError, http.client.HTTPException) as exc:
            record.update(status="error", error=f"{type(exc).__name__}: {exc}")
        results.append(record)
    return {
        "requested": len(entries),
        "ready": sum(record["status"] == "ready" for record in results),
        "errors": sum(record["status"] == "error" for record in results),
        "downloaded": sum(record.get("downloaded") is True for record in results),
        "results": results,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--fetch", action="store_true", help="explicitly allow NCBI HTTPS downloads"
    )
    parser.add_argument("--timeout", type=float, default=60, help="network timeout in seconds")
    arguments = parser.parse_args(argv)
    try:
        result = fetch_benchmark_genomes(
            arguments.manifest,
            arguments.output_dir,
            fetch=arguments.fetch,
            timeout=arguments.timeout,
        )
    except (OSError, ValueError) as exc:
        parser.exit(2, f"benchmark input error: {exc}\n")
    for record in result["results"]:
        if record["status"] == "error":
            print(f"{record['accession']}: {record['error']}", file=sys.stderr)
    print(json.dumps({key: value for key, value in result.items() if key != "results"}))
    return 1 if result["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
