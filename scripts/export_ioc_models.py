"""Export the five original R PCA models to portable JSON, without refitting.

Usage: python scripts/export_ioc_models.py ../gRodon2/R/sysdata.rda
Requires the optional ``rdata`` development dependency.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import rdata
from export_growth_models import top_level_objects


def export(source: Path, destination: Path) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        source_models = rdata.conversion.convert(top_level_objects(source)["IoC_models"])
    names = {
        "Soil": "soil",
        "Permafrost": "permafrost",
        "Pacific": "pacific",
        "Human": "human",
        "All": "all_habitat",
    }
    models = {}
    for original_name, original in source_models.items():
        name = names[str(original_name)]
        rotation = np.asarray(original["rotation"], dtype=float)
        scores = np.asarray(original["x"], dtype=float)
        fields = {
            key: np.asarray(original[key], dtype=float).tolist()
            for key in ("center", "scale", "rotation", "sdev")
        }
        fields.update(
            features=["rCAZy", "dCUB", "nGenes"],
            min_score=float(np.min(scores[:, 0] * -np.sign(rotation[1, 0]))),
            n_observations=len(scores),
            name=name,
        )
        models[name] = fields
    bundle = {
        "schema_version": 1,
        "source": "gRodon 2.7.3 R/sysdata.rda, object IoC_models",
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "models": models,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(bundle, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Exported {len(models)} original PCA models to {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, default=Path("src/grodon/data/ioc_models.json"))
    args = parser.parse_args()
    export(args.source.resolve(), args.output)
