"""Run the prepared folder test case against native Python and original R.

Preparation is separate and explicitly network-enabled. This runner is offline,
uses the package's real genome pipeline, and never changes model parameters or
silently retries with looser annotation thresholds to obtain passing results.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from datetime import UTC, datetime
from numbers import Integral
from pathlib import Path
from typing import Any

from grodon.pipeline import GenomePipelineConfig, run_genome
from grodon.references import sha256_file
from grodon.tools import doctor

FULL_METRICS = frozenset(
    {
        "CUBHE",
        "GC",
        "GCdiv",
        "ConsistencyHE",
        "CUB",
        "CPB",
        "FilteredSequences",
        "nHE",
        "dCUB",
        "d",
        "LowerCI",
        "UpperCI",
    }
)


def write_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
    temporary.replace(path)


def _resolved_program(value: str) -> str:
    found = shutil.which(value)
    if found is None:
        raise FileNotFoundError(f"Executable not found: {value}")
    return str(Path(found).resolve())


def _check_suite_resume(folder: Path, rscript: str) -> str:
    """Require recorded R provenance before reusing an existing suite directory."""
    path = folder / "suite.json"
    if not path.is_file():
        raise ValueError(
            f"Cannot resume an existing directory without suite.json: {folder}; "
            "use a new --runs-dir"
        )
    previous = json.loads(path.read_text(encoding="utf-8"))
    recorded = previous.get("rscript") if isinstance(previous, dict) else None
    if not isinstance(recorded, str) or not recorded or not Path(recorded).is_absolute():
        raise ValueError(f"Cannot resume a suite without an absolute recorded R executable: {path}")
    if Path(recorded).resolve() != Path(rscript).resolve():
        raise ValueError(
            f"Cannot resume with a different R executable: recorded {recorded!r}, "
            f"requested {rscript!r}; use a new --runs-dir"
        )
    return recorded


def _check_entry(entry: dict[str, Any]) -> None:
    accession = entry.get("accession", "")
    if not re.fullmatch(r"GC[AF]_\d{9}\.\d+", accession):
        raise ValueError(f"Invalid accession: {accession!r}")
    if entry.get("status") != "ready":
        raise ValueError(f"Input is not ready: {accession}: {entry.get('error')}")
    for path_key, hash_key in (("source_file", "source_sha256"), ("genome_path", "genome_sha256")):
        if sha256_file(entry[path_key]) != entry[hash_key]:
            raise ValueError(f"Input checksum mismatch: {entry[path_key]}")


def _check_resume(
    folder: Path, entry: dict[str, Any], config: GenomePipelineConfig
) -> dict[str, Any]:
    run = json.loads((folder / "run.json").read_text(encoding="utf-8"))
    result = json.loads((folder / "result.json").read_text(encoding="utf-8"))
    if (
        run["status"] != "complete"
        or result.get("status") != "complete"
        or result.get("r_comparison_status") != "complete"
    ):
        raise ValueError(f"Cannot reuse an incomplete run: {folder}; use a new --runs-dir")
    if run["input"]["sha256"] != entry["genome_sha256"] or run["settings"] != asdict(config):
        raise ValueError(f"Cannot reuse a run with different inputs/settings: {folder}")
    expected = run["artifacts"]
    required = {
        "result.json",
        "r_result.json",
        "r_comparison.json",
        "r_parameters.json",
        "predicted_cds.fna",
        "highly_expressed_ids.txt",
    }
    if not required <= set(expected):
        raise ValueError(f"Missing comparison artifacts: {folder}")
    actual = {path.name for path in folder.iterdir() if path.is_file() and path.name != "run.json"}
    if actual != set(expected):
        raise ValueError(f"Artifact inventory changed: {folder}")
    for name, values in expected.items():
        target = (folder / name).resolve()
        if not target.is_relative_to(folder.resolve()):
            raise ValueError(f"Unsafe artifact path in manifest: {name}")
        if target.stat().st_size != values["bytes"] or sha256_file(target) != values["sha256"]:
            raise ValueError(f"Artifact checksum mismatch: {target}")
    return result


def _one(
    entry: dict[str, Any], folder: str, settings: dict[str, Any], rscript: str, resume: bool
) -> dict[str, Any]:
    started = time.perf_counter()
    record: dict[str, Any] = {"accession": entry["accession"], "output_dir": folder}
    try:
        _check_entry(entry)
        config = GenomePipelineConfig(**settings)
        output = Path(folder)
        reused = output.exists() and resume
        result = (
            _check_resume(output, entry, config)
            if reused
            else run_genome(entry["genome_path"], output, config=config, rscript=rscript).to_dict()
        )
        statistics = result.get("statistics")
        finite = (
            isinstance(statistics, dict)
            and bool(statistics)
            and all(
                isinstance(value, (int, float))
                and not isinstance(value, bool)
                and math.isfinite(value)
                for value in statistics.values()
            )
        )
        full_metrics = isinstance(statistics, dict) and set(statistics) == FULL_METRICS
        complete = (
            result.get("status") == "complete" and result.get("r_comparison_status") == "complete"
        )
        record.update(
            status="pass"
            if complete and full_metrics and finite and result.get("r_matches") is True
            else "mismatch",
            reused=reused,
            pipeline_complete=complete,
            full_metric_set=full_metrics,
            all_python_statistics_finite=finite,
            r_matches=result.get("r_matches"),
            doubling_time_hours=result["doubling_time_hours"],
            scientific_fingerprint=result["scientific_fingerprint"],
        )
    except Exception as exc:
        record.update(status="error", error=f"{type(exc).__name__}: {exc}")
    record["elapsed_seconds"] = time.perf_counter() - started
    return record


def run_suite(
    manifest_path: Path,
    runs_dir: Path,
    *,
    blastp: str,
    rscript: str,
    workers: int = 4,
    timeout: float = 1800,
    resume: bool = False,
) -> dict[str, Any]:
    """Execute every manifest entry and preserve explicit per-input outcomes."""
    if isinstance(workers, bool) or not isinstance(workers, Integral) or workers < 1:
        raise ValueError("workers must be positive")
    manifest_path = manifest_path.resolve(strict=True)
    inventory = json.loads(manifest_path.read_text(encoding="utf-8"))
    entries = inventory["entries"]
    identifiers = [row["accession"] for row in entries]
    if not entries or len(identifiers) != len(set(identifiers)):
        raise ValueError("The inventory must contain nonempty, unique assembly accessions")
    for entry in entries:
        if not re.fullmatch(r"GC[AF]_\d{9}\.\d+", entry["accession"]):
            raise ValueError(f"Invalid accession: {entry['accession']!r}")
    executable = _resolved_program(blastp)
    r_executable = _resolved_program(rscript)
    health = doctor(executable)
    if not health["ready"]:
        raise RuntimeError(f"Genome dependencies are unavailable: {health}")
    config = GenomePipelineConfig(blastp=executable, timeout=timeout)
    output = runs_dir.resolve()
    if output.exists() and not resume:
        raise FileExistsError(f"Run directory exists; use --resume to verify/reuse: {output}")
    if output.exists():
        r_executable = _check_suite_resume(output, r_executable)
    output.mkdir(parents=True, exist_ok=True)
    # Spawned Python and BLAST/R processes each use one numeric-library thread;
    # concurrency is bounded at the assembly level and recorded in suite.json.
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[variable] = "1"
    suite: dict[str, Any] = {
        "schema_version": 1,
        "status": "running",
        "started_at": datetime.now(UTC).isoformat(),
        "manifest": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "settings": asdict(config),
        "rscript": r_executable,
        "workers": workers,
        "numeric_library_threads": 1,
        "requested": len(entries),
        "finished": 0,
        "passed": 0,
        "mismatched": 0,
        "errors": 0,
        "results": [],
    }
    write_json(output / "suite.json", suite)
    write_json(output / "input_manifest.json", inventory)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                _one, entry, str(output / entry["accession"]), asdict(config), r_executable, resume
            ): entry["accession"]
            for entry in entries
        }
        for future in as_completed(futures):
            accession = futures[future]
            try:
                record = future.result()
            except Exception as exc:
                record = {"accession": accession, "status": "error", "error": str(exc)}
            suite["results"].append(record)
            suite["results"].sort(key=lambda row: row["accession"])
            suite["finished"] += 1
            counter = {"pass": "passed", "mismatch": "mismatched", "error": "errors"}[
                record["status"]
            ]
            suite[counter] += 1
            write_json(output / "suite.json", suite)
            print(
                f"[{suite['finished']}/{suite['requested']}] {accession}: {record['status']}"
                + (f" ({record['error']})" if "error" in record else ""),
                file=sys.stderr,
                flush=True,
            )
    suite["status"] = (
        "complete" if suite["passed"] == suite["requested"] else "completed_with_failures"
    )
    suite["finished_at"] = datetime.now(UTC).isoformat()
    write_json(output / "suite.json", suite)
    return suite


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--runs-dir", type=Path, required=True)
    parser.add_argument("--blastp", required=True)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    try:
        result = run_suite(
            args.manifest,
            args.runs_dir,
            blastp=args.blastp,
            rscript=args.rscript,
            workers=args.workers,
            timeout=args.timeout,
            resume=args.resume,
        )
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        parser.exit(2, f"validation error: {exc}\n")
    print(
        json.dumps(
            {key: result[key] for key in ("status", "requested", "passed", "mismatched", "errors")},
            indent=2,
        )
    )
    return 0 if result["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
