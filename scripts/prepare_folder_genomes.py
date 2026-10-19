"""Inventory versioned assemblies and explicitly fetch missing genomic FASTA.

Local input files are never changed. Protein-only inputs are resolved to the
same assembly accession/version through official NCBI reports and checksums.
Network access requires --fetch; verified cached downloads work offline.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

_FILENAME = re.compile(r"^(GC[AF]_\d{9}\.\d+)_(.+)_(genomic\.fna|protein\.faa)\.gz$")
_ALLOWED_HOSTS = {"api.ncbi.nlm.nih.gov", "ftp.ncbi.nlm.nih.gov"}
_REQUEST_LOCK = threading.Lock()
_NEXT_REQUEST = 0.0


class NeedsFetch(Exception):
    """A required cache resource is unavailable and network access is disabled."""


def digest(path: Path, algorithm: str = "sha256") -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def protein_snapshot(path: Path) -> dict[str, Any]:
    """Read protein identities and sequence content, independent of gzip encoding."""
    records: dict[str, str] = {}
    headers: dict[str, str] = {}
    checksum = hashlib.sha256()
    identifier: str | None = None
    sequence: list[str] = []
    size = 0

    def finish() -> None:
        if identifier is not None:
            protein = "".join(sequence)
            if not protein:
                raise ValueError(f"Empty protein record in {path.name}: {identifier}")
            records[identifier] = protein

    with gzip.open(path, "rb") as handle:
        for raw_line in handle:
            size += len(raw_line)
            if size > 150_000_000:
                raise ValueError(f"Decompressed protein FASTA exceeds limit: {path.name}")
            checksum.update(raw_line)
            line = raw_line.decode("utf-8-sig").strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty protein FASTA header: {path.name}")
                identifier = header.split()[0]
                if identifier in headers:
                    raise ValueError(f"Duplicate protein identifier: {identifier}")
                headers[identifier] = header
                sequence = []
            else:
                if identifier is None:
                    raise ValueError(f"Sequence before FASTA header: {path.name}")
                sequence.append("".join(line.split()).upper())
    finish()
    if not records:
        raise ValueError(f"Empty protein FASTA: {path.name}")
    return {"records": records, "headers": headers, "sha256": checksum.hexdigest()}


def protein_comparison(source: Path, remote: Path) -> dict[str, Any]:
    original = protein_snapshot(source)
    official = original if source == remote else protein_snapshot(remote)
    old_records, new_records = original["records"], official["records"]
    old_ids, new_ids = set(old_records), set(new_records)
    content_matches = old_records == new_records
    headers_match = original["headers"] == official["headers"]
    bytes_match = original["sha256"] == official["sha256"]
    if source == remote:
        difference = "identical_compressed_archive"
    elif bytes_match:
        difference = "gzip_archive_bytes_only"
    elif content_matches and headers_match:
        difference = "record_order_or_text_formatting_only"
    elif content_matches:
        difference = "annotation_headers_differ_protein_sequences_identical"
    else:
        difference = "protein_records_differ"
    return {
        "source_protein_content_matches": content_matches,
        "source_protein_order_matches": list(old_records) == list(new_records),
        "source_protein_headers_match": headers_match,
        "source_protein_decompressed_bytes_match": bytes_match,
        "source_protein_original_count": len(old_records),
        "source_protein_official_count": len(new_records),
        "source_protein_shared_id_count": len(old_ids & new_ids),
        "source_protein_changed_sequence_count": sum(
            old_records[key] != new_records[key] for key in old_ids & new_ids
        ),
        "source_protein_only_original_id_count": len(old_ids - new_ids),
        "source_protein_only_official_id_count": len(new_ids - old_ids),
        "source_protein_original_decompressed_sha256": original["sha256"],
        "source_protein_official_decompressed_sha256": official["sha256"],
        "source_protein_difference_kind": difference,
    }


def save_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _throttle() -> None:
    global _NEXT_REQUEST
    with _REQUEST_LOCK:
        delay = max(0.0, _NEXT_REQUEST - time.monotonic())
        if delay:
            time.sleep(delay)
        _NEXT_REQUEST = time.monotonic() + 0.35


def cached_resource(
    url: str, target: Path, *, fetch: bool, limit: int, expected_md5: str | None = None
) -> dict[str, Any]:
    """Use a verified cache receipt, or retrieve bounded bytes from NCBI HTTPS."""
    receipt_path = target.with_name(target.name + ".receipt.json")
    if target.is_file() and receipt_path.is_file():
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            valid = receipt["url"] == url and receipt["sha256"] == digest(target)
            valid = valid and target.stat().st_size <= limit
            if expected_md5 is not None:
                valid = valid and digest(target, "md5") == expected_md5
            if valid:
                return receipt
        except (OSError, ValueError, KeyError, TypeError):
            pass
    if not fetch:
        raise NeedsFetch(f"Missing or invalid cache: {target.name}; rerun with --fetch")
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or parsed.hostname not in _ALLOWED_HOSTS:
        raise ValueError(f"Unsupported remote source: {url}")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".part")
    for attempt in range(3):
        try:
            _throttle()
            request = urllib.request.Request(
                url, headers={"User-Agent": "grodon2-py/reproducible-assembly-preparer"}
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                final = urllib.parse.urlparse(response.geturl())
                if final.scheme != "https" or final.hostname not in _ALLOWED_HOSTS:
                    raise ValueError("NCBI download redirected outside the permitted HTTPS hosts")
                size = 0
                with temporary.open("wb") as stream:
                    while block := response.read(128 * 1024):
                        size += len(block)
                        if size > limit:
                            raise ValueError(f"Download exceeds {limit} bytes: {url}")
                        stream.write(block)
                expected_size = response.headers.get("Content-Length")
                if expected_size is not None and size != int(expected_size):
                    raise ValueError(f"Truncated download: {url}")
                receipt = {
                    "url": url,
                    "final_url": response.geturl(),
                    "bytes": size,
                    "sha256": digest(temporary),
                    "md5": digest(temporary, "md5"),
                    "retrieved_utc": datetime.now(UTC).isoformat(),
                    "http_date": response.headers.get("Date"),
                    "http_last_modified": response.headers.get("Last-Modified"),
                }
            if expected_md5 is not None and receipt["md5"] != expected_md5:
                raise ValueError(f"Official MD5 mismatch for {target.name}")
            temporary.replace(target)
            save_json(receipt_path, receipt)
            return receipt
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            temporary.unlink(missing_ok=True)
            if isinstance(exc, urllib.error.HTTPError) and exc.code not in {
                429,
                500,
                502,
                503,
                504,
            }:
                raise
            if attempt == 2:
                raise
            time.sleep(1.0 + attempt * 2.0)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    raise RuntimeError("Unreachable download state")


def parse_checksums(path: Path) -> dict[str, str]:
    result = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-fA-F]{32})\s+\*?(.+)", line.strip())
        if match is None:
            raise ValueError(f"Invalid NCBI checksum line {number}")
        name = match[2].removeprefix("./")
        if name in result:
            raise ValueError(f"Duplicate NCBI checksum path: {name}")
        result[name] = match[1].lower()
    return result


def prepare_entry(source: Path, accession: str, kind: str, cache: Path, fetch: bool) -> dict:
    entry: dict[str, Any] = {
        "accession": accession,
        "source_file": str(source),
        "source_kind": kind,
        "source_sha256": None,
        "genome_path": None,
        "genome_sha256": None,
        "organism": None,
        "assembly": None,
        "download": {},
        "status": "error",
        "error": None,
    }
    try:
        entry["source_sha256"] = digest(source)
        folder = cache / accession
        report_path = folder / "dataset_report.json"
        report_url = (
            f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report"
        )
        report_receipt = cached_resource(report_url, report_path, fetch=fetch, limit=4_000_000)
        payload = json.loads(report_path.read_text(encoding="utf-8"))
        reports = payload.get("reports", [])
        if len(reports) != 1 or reports[0].get("accession") != accession:
            raise ValueError(f"NCBI report does not uniquely identify exact assembly {accession}")
        report = reports[0]
        info, organism = report["assembly_info"], report["organism"]
        assembly_name = info["assembly_name"]
        if (
            not isinstance(assembly_name, str)
            or not assembly_name
            or any(c in assembly_name for c in "/\\\x00")
        ):
            raise ValueError("Invalid assembly name in NCBI report")
        if not organism.get("organism_name") or not isinstance(organism.get("tax_id"), int):
            raise ValueError("NCBI report has no organism name/taxonomy identifier")
        entry["organism"] = {
            "name": organism["organism_name"],
            "tax_id": organism["tax_id"],
            "strain": organism.get("infraspecific_names", {}).get("strain"),
        }
        entry["assembly"] = {
            "accession": report["accession"],
            "name": assembly_name,
            "level": info.get("assembly_level"),
            "status": info.get("assembly_status"),
            "current_accession": report.get("current_accession"),
            "source_database": report.get("source_database"),
            "annotation": report.get("annotation_info"),
            "taxonomy_check": report.get("average_nucleotide_identity", {}).get(
                "taxonomy_check_status"
            ),
        }
        digits = accession.split("_")[1].split(".")[0]
        basename = f"{accession}_{assembly_name.replace(' ', '_')}"
        directory = (
            f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[:3]}/"
            f"{digits[:3]}/{digits[3:6]}/{digits[6:]}/{urllib.parse.quote(basename, safe='')}/"
        )
        checksums_path = folder / "md5checksums.txt"
        checksums_url = directory + "md5checksums.txt"
        checksums_receipt = cached_resource(
            checksums_url, checksums_path, fetch=fetch, limit=2_000_000
        )
        checksums = parse_checksums(checksums_path)
        genome_name, protein_name = basename + "_genomic.fna.gz", basename + "_protein.faa.gz"
        genome_url = directory + urllib.parse.quote(genome_name, safe="")
        if genome_name not in checksums:
            raise ValueError("Official genomic FASTA is absent from NCBI checksum manifest")
        source_md5 = digest(source, "md5")
        entry["download"] = {
            "dataset_report_url": report_url,
            "dataset_report_path": str(report_path),
            "dataset_report_sha256": report_receipt["sha256"],
            "checksums_url": checksums_url,
            "checksums_path": str(checksums_path),
            "checksums_sha256": checksums_receipt["sha256"],
            "genome_url": genome_url,
            "genome_md5_expected": checksums[genome_name],
            "source_md5_actual": source_md5,
            "protein_md5_expected": checksums.get(protein_name),
            "source_protein_md5_matches": None,
            "protein_url": directory + urllib.parse.quote(protein_name, safe=""),
            "metadata_retrieved_utc": report_receipt["retrieved_utc"],
        }
        if kind == "protein":
            if protein_name not in checksums:
                raise ValueError("Official protein FASTA is absent from NCBI checksum manifest")
            entry["download"]["source_protein_md5_matches"] = source_md5 == checksums[protein_name]
            official_protein = source
            if not entry["download"]["source_protein_md5_matches"]:
                official_protein = folder / protein_name
                protein_receipt = cached_resource(
                    entry["download"]["protein_url"],
                    official_protein,
                    fetch=fetch,
                    limit=50_000_000,
                    expected_md5=checksums[protein_name],
                )
                entry["download"].update(
                    official_protein_path=str(official_protein),
                    official_protein_sha256=protein_receipt["sha256"],
                    official_protein_md5_actual=protein_receipt["md5"],
                    official_protein_retrieved_utc=protein_receipt["retrieved_utc"],
                )
            entry["download"].update(protein_comparison(source, official_protein))
            genome = folder / genome_name
            receipt = cached_resource(
                genome_url,
                genome,
                fetch=fetch,
                limit=50_000_000,
                expected_md5=checksums[genome_name],
            )
            entry["download"]["genome_retrieved_utc"] = receipt["retrieved_utc"]
            entry["download"]["used_original_genome"] = False
        else:
            genome = source
            entry["download"]["used_original_genome"] = True
        actual_md5 = digest(genome, "md5")
        entry["download"].update(
            genome_md5_actual=actual_md5,
            genome_md5_matches=actual_md5 == checksums[genome_name],
        )
        if actual_md5 != checksums[genome_name]:
            raise ValueError(f"Official MD5 mismatch for genomic input {genome.name}")
        if digest(source) != entry["source_sha256"]:
            raise RuntimeError("Source file changed during preparation; retry from a stable input")
        entry.update(genome_path=str(genome), genome_sha256=digest(genome), status="ready")
    except NeedsFetch as exc:
        entry.update(status="needs_fetch", error=str(exc))
    except Exception as exc:
        entry.update(status="error", error=f"{type(exc).__name__}: {exc}")
    return entry


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true", help="Explicitly permit NCBI downloads")
    parser.add_argument("--workers", type=int, choices=range(1, 5), default=3)
    args = parser.parse_args()
    source_dir, cache = args.input_dir.resolve(strict=True), args.cache_dir.resolve()
    if not source_dir.is_dir() or cache == source_dir or cache.is_relative_to(source_dir):
        parser.error("Use an existing input directory and a separate cache outside it")
    if args.manifest.resolve().is_relative_to(source_dir):
        parser.error("The output manifest must not be written inside the input directory")
    inputs, seen, ignored = [], set(), []
    for source in sorted(source_dir.iterdir()):
        if not source.is_file():
            continue
        match = _FILENAME.fullmatch(source.name)
        if match is None:
            if source.name.endswith((".fna.gz", ".faa.gz")):
                parser.error(
                    f"Sequence filename lacks a strict versioned assembly accession: {source.name}"
                )
            ignored.append(source.name)
            continue
        accession = match[1]
        if accession in seen:
            parser.error(f"More than one input represents assembly {accession}")
        seen.add(accession)
        inputs.append((source, accession, "genomic" if match[3] == "genomic.fna" else "protein"))
    if not inputs:
        parser.error("No versioned genomic/protein FASTA inputs were found")
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "input_dir": str(source_dir),
        "cache_dir": str(cache),
        "fetch_requested": args.fetch,
        "total": len(inputs),
        "ready": 0,
        "errors": 0,
        "needs_fetch": 0,
        "ignored_files": ignored,
        "entries": [],
        "scope": "Exact-accession genomic inputs for R/Python growth-estimate comparison; source files remain unchanged",
    }
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        jobs = {
            executor.submit(prepare_entry, source, accession, kind, cache, args.fetch): accession
            for source, accession, kind in inputs
        }
        for job in as_completed(jobs):
            entry = job.result()
            manifest["entries"].append(entry)
            manifest["entries"].sort(key=lambda item: item["accession"])
            manifest["ready"] = sum(item["status"] == "ready" for item in manifest["entries"])
            manifest["errors"] = sum(item["status"] == "error" for item in manifest["entries"])
            manifest["needs_fetch"] = sum(
                item["status"] == "needs_fetch" for item in manifest["entries"]
            )
            save_json(args.manifest, manifest)
            print(
                f"{entry['accession']}: {entry['status']}"
                + (f" ({entry['error']})" if entry["error"] else ""),
                file=sys.stderr,
                flush=True,
            )
    print(json.dumps({key: manifest[key] for key in ("total", "ready", "errors", "needs_fetch")}))
    return 0 if manifest["ready"] == manifest["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
