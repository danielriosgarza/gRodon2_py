"""Export the original R linear models without executing R or refitting them.

Usage: python scripts/export_growth_models.py ../gRodon2/R/sysdata.rda
Requires the optional development dependency ``rdata``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import warnings
from pathlib import Path

import numpy as np
import rdata


def named_fields(obj):
    names = rdata.conversion.convert(obj.attributes)["names"]
    return dict(zip(names, obj.value, strict=True))


def top_level_objects(path: Path):
    node = rdata.parser.parse_file(path).object
    objects = {}
    while node.info.type.name != "NILVALUE":
        tag = node.tag
        while tag.referenced_object is not None:
            tag = tag.referenced_object
        objects[tag.value.value.decode("utf-8")] = node.value[0]
        node = node.value[1]
    return objects


def lambda_name(name: str) -> str:
    if name.endswith("_t150"):
        return "lambda_t150"
    if name.endswith("_t240"):
        return "lambda_t240"
    if name.endswith("_AOANOB"):
        return "lambda_milc_AOANOB"
    if "_euk" in name:
        return "lambda_milc_euk_i" if name.endswith("_i") else "lambda_milc_euk"
    if "_newmeta" in name:
        return "lambda_newmeta_i" if name.endswith("_i") else "lambda_newmeta"
    if "_madin" in name:
        return "lambda_milc_madin_i" if name.endswith("_i") else "lambda_milc_madin"
    return "lambda_milc"


def export(source: Path, destination: Path) -> None:
    objects = top_level_objects(source)
    exported = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        for name, obj in objects.items():
            if not name.startswith("gRodon_model_"):
                continue
            fields = named_fields(obj)
            coefficients = rdata.conversion.convert(fields["coefficients"])
            predictors = coefficients.coords[coefficients.dims[0]].values.tolist()
            coefficients = np.asarray(coefficients, dtype=float)
            qr = rdata.conversion.convert(fields["qr"])
            rank = int(np.asarray(rdata.conversion.convert(fields["rank"])).item())
            if rank != len(coefficients):
                raise ValueError(f"Cannot export rank-deficient model {name}")
            r_matrix = np.triu(np.asarray(qr["qr"], dtype=float)[:rank, :rank])
            inverse_r = np.linalg.solve(r_matrix, np.eye(rank))
            residuals = np.asarray(rdata.conversion.convert(fields["residuals"]), dtype=float)
            df_residual = int(np.asarray(rdata.conversion.convert(fields["df.residual"])).item())
            if "weights" in fields:
                weights = np.asarray(rdata.conversion.convert(fields["weights"]), dtype=float)
                residual_variance = float(np.dot(weights * residuals, residuals) / df_residual)
            else:
                residual_variance = float(np.dot(residuals, residuals) / df_residual)
            pivot = np.asarray(qr["pivot"], dtype=int) - 1
            covariance = np.empty((rank, rank))
            covariance[np.ix_(pivot, pivot)] = inverse_r @ inverse_r.T * residual_variance
            transform = lambda_name(name)
            model_lambda = float(np.asarray(rdata.conversion.convert(objects[transform])).item())
            training_set = (
                "eukaryote"
                if "_euk" in name
                else (
                    "AOA_NOB"
                    if "_AOANOB" in name
                    else ("vs" if transform == "lambda_milc" else "madin")
                )
            )
            exported[name] = {
                "name": name,
                "predictors": predictors,
                "coefficients": coefficients.tolist(),
                "covariance": covariance.tolist(),
                "df_residual": df_residual,
                "sigma": float(np.sqrt(residual_variance)),
                "boxcox_lambda": model_lambda,
                "lambda_source": transform,
                "training_set": training_set,
                "description": "Original gRodon 2.7.3 fitted linear model; covariance from stored R QR decomposition.",
            }
    try:
        revision = subprocess.check_output(
            [
                "git",
                "-c",
                f"safe.directory={source.parent.parent.as_posix()}",
                "-C",
                str(source.parent.parent),
                "rev-parse",
                "HEAD",
            ],
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = None
    bundle = {
        "schema_version": 1,
        "source": {
            "package": "gRodon",
            "version": "2.7.3",
            "file": "R/sysdata.rda",
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "git_revision": revision,
            "export_method": "Original coefficients and R QR covariance; no model refitting",
        },
        "models": exported,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(bundle, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Exported {len(exported)} models to {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, default=Path("src/grodon/data/growth_models.json"))
    args = parser.parse_args()
    export(args.source.resolve(), args.output)
