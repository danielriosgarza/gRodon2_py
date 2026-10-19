"""Optional comparison with an explicitly supplied original R installation."""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
from collections.abc import Mapping
from contextlib import ExitStack
from importlib.resources import as_file, files
from numbers import Real
from pathlib import Path
from typing import Any

_RTOL = 1e-10
_ATOL = 1e-12


def _numeric_statistics(values: Any, source: str) -> dict[str, float | int | None]:
    if not isinstance(values, Mapping) or not values:
        raise ValueError(f"{source} statistics must be a nonempty object of numeric scalars")
    normalized: dict[str, float | int | None] = {}
    for key, value in values.items():
        if not isinstance(key, str):
            raise ValueError(f"{source} statistic names must be strings")
        if value is None:
            normalized[key] = None
        elif isinstance(value, bool) or not isinstance(value, Real):
            raise ValueError(f"{source} statistic {key!r} is not a numeric scalar")
        elif not math.isfinite(value):
            normalized[key] = None
        else:
            normalized[key] = int(value) if isinstance(value, int) else float(value)
    return normalized


def _comparison_rows(
    python: Mapping[str, float | int | None], original: Mapping[str, float | int | None]
) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(set(python) | set(original)):
        left, right = python.get(key), original.get(key)
        absolute = relative = None
        reason: str | None
        if key not in python or key not in original:
            matches = False
            reason = "missing from Python" if key not in python else "missing from R"
        elif left is None or right is None:
            matches = left is None and right is None
            reason = None if matches else "finite and nonfinite values differ"
        else:
            difference = abs(left - right)
            absolute = difference if math.isfinite(difference) else None
            if right != 0:
                ratio = difference / abs(right)
                relative = ratio if math.isfinite(ratio) else None
            elif difference == 0:
                relative = 0.0
            matches = math.isclose(left, right, rel_tol=_RTOL, abs_tol=_ATOL)
            reason = None if matches else "outside numerical tolerance"
        rows.append(
            {
                "metric": key,
                "python": left,
                "r": right,
                "absolute_difference": absolute,
                "relative_difference": relative,
                "matches": matches,
                "reason": reason,
            }
        )
    return rows


def _text(value: str | bytes | None) -> str:
    return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else value or ""


def _invalid_constant(value: str) -> None:
    raise ValueError(f"Invalid JSON numeric constant: {value}")


def compare_with_r(
    cds_path: str | os.PathLike[str],
    ids_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    *,
    rscript: str | os.PathLike[str],
    mode: str = "full",
    training_set: str = "madin",
    temperature: float | None = None,
    genetic_code: int = 11,
    python_statistics: Mapping[str, float | int],
    timeout: float = 600,
) -> dict[str, Any]:
    """Compare all numeric statistics with actual ``gRodon::predictGrowth``.

    Both implementations must use the same CDS and HE ID files, a 240 bp length
    threshold (``fragments=False``), and 95% confidence intervals. R must already
    contain gRodon, its dependencies, and jsonlite. This function installs
    nothing and performs no network access. Stochastic modes and documented R
    dispatch bugs may produce mismatches; mismatches are reported, not raised.

    Four audit files are created exclusively: ``r_parameters.json``,
    ``r_result.json``, ``r_comparison.json`` and ``r_comparison.log``. Process or
    malformed-result errors raise RuntimeError and retain parameters and logs.
    Empty reserved result files may remain after such an unsuccessful run.
    """
    python = _numeric_statistics(python_statistics, "Python")
    if (
        isinstance(timeout, bool)
        or not isinstance(timeout, Real)
        or not math.isfinite(timeout)
        or timeout <= 0
    ):
        raise ValueError("R comparison timeout must be a positive finite number")
    if temperature is not None and (
        isinstance(temperature, bool)
        or not isinstance(temperature, Real)
        or not math.isfinite(temperature)
    ):
        raise ValueError("temperature must be a finite number or None")
    if isinstance(genetic_code, bool) or not isinstance(genetic_code, int) or genetic_code <= 0:
        raise ValueError("genetic_code must be a positive integer")
    cds, ids = Path(cds_path).resolve(), Path(ids_path).resolve()
    for path in (cds, ids):
        if not path.is_file():
            raise FileNotFoundError(f"R comparison input does not exist: {path}")
    requested = os.fspath(rscript)
    resolved = shutil.which(requested)
    if resolved is None and Path(requested).is_file():
        resolved = str(Path(requested).resolve())
    if resolved is None:
        raise FileNotFoundError(f"Rscript executable not found: {requested}")
    destination = Path(output_dir).resolve()
    paths = {
        "parameters": destination / "r_parameters.json",
        "result": destination / "r_result.json",
        "comparison": destination / "r_comparison.json",
        "log": destination / "r_comparison.log",
    }
    for path in paths.values():
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Refusing to overwrite an R comparison artifact: {path}")
    destination.mkdir(parents=True, exist_ok=True)
    parameters = {
        "cds_path": str(cds),
        "ids_path": str(ids),
        "mode": mode,
        "training_set": "vs" if training_set == "original" else training_set,
        "temperature": None if temperature is None else float(temperature),
        "genetic_code": genetic_code,
        "fragments": False,
        "confidence": 0.95,
    }
    script_resource = files("grodon").joinpath("data", "compare_genome.R")
    with ExitStack() as stack:
        streams = {
            name: stack.enter_context(path.open("x", encoding="utf-8"))
            for name, path in paths.items()
        }
        json.dump(parameters, streams["parameters"], indent=2, allow_nan=False)
        streams["parameters"].write("\n")
        streams["parameters"].flush()
        script = stack.enter_context(as_file(script_resource))
        command = [
            str(Path(resolved).resolve()),
            "--vanilla",
            str(script),
            str(paths["parameters"]),
        ]
        log = streams["log"]
        log.write("Command arguments: " + json.dumps(command) + "\n")
        log.flush()
        try:
            process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=False,
                timeout=float(timeout),
                shell=False,
            )
        except subprocess.TimeoutExpired as exc:
            log.write(
                "Timed out.\nSTDOUT:\n" + _text(exc.stdout) + "\nSTDERR:\n" + _text(exc.stderr)
            )
            raise RuntimeError(
                f"R comparison timed out after {timeout} seconds; see {paths['log']}"
            ) from exc
        except OSError as exc:
            log.write(f"Could not start Rscript: {exc}\n")
            raise RuntimeError(f"Could not run Rscript; see {paths['log']}: {exc}") from exc
        log.write(
            f"Exit status: {process.returncode}\nSTDOUT:\n{_text(process.stdout)}\nSTDERR:\n{_text(process.stderr)}"
        )
        log.flush()
        if process.returncode != 0:
            raise RuntimeError(
                f"Original R gRodon comparison failed with exit status {process.returncode}; see {paths['log']}"
            )
        try:
            original = json.loads(_text(process.stdout), parse_constant=_invalid_constant)
            if not isinstance(original, dict) or not isinstance(original.get("versions"), dict):
                raise ValueError("R result must contain statistics and package versions")
            if not {"R", "gRodon"} <= set(original["versions"]):
                raise ValueError("R result is missing R/gRodon version information")
            original["statistics"] = _numeric_statistics(original.get("statistics"), "R")
        except (ValueError, TypeError) as exc:
            raise RuntimeError(
                f"R returned a malformed comparison result; see {paths['log']}: {exc}"
            ) from exc
        json.dump(original, streams["result"], indent=2, allow_nan=False)
        streams["result"].write("\n")
        comparison = _comparison_rows(python, original["statistics"])
        report = {
            "matches": all(row["matches"] for row in comparison),
            "tolerance": {"relative": _RTOL, "absolute": _ATOL},
            "comparison": comparison,
            "python_statistics": python,
            "r_result": original,
            "files": {name: str(path) for name, path in paths.items()},
        }
        json.dump(report, streams["comparison"], indent=2, allow_nan=False)
        streams["comparison"].write("\n")
    return report
