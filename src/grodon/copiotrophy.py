"""Fit and apply the PCA-based Index of Copiotrophy (IoC).

The implementation follows gRodon 2.7.3: scale the CAZyme fraction, dCUB,
and completeness-corrected gene count; orient PC1 toward decreasing dCUB;
then offset it so that the smallest *training* score is one. Predictions on
new observations can be less than one and are intentionally not clipped.
"""

from __future__ import annotations

import json
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

_FEATURES = ("rCAZy", "dCUB", "nGenes")
_BUILTINS = ("permafrost", "pacific", "soil", "human", "all_habitat")


@dataclass(frozen=True)
class IoCModel:
    """Portable PCA model and its training-set origin.

    ``rotation`` has features in rows and principal components in columns,
    matching R's ``prcomp`` convention. Arrays are read-only. Training scores
    are optional because only the minimum oriented PC1 is needed to predict.
    """

    center: npt.NDArray[np.float64]
    scale: npt.NDArray[np.float64]
    rotation: npt.NDArray[np.float64]
    sdev: npt.NDArray[np.float64]
    min_score: float
    n_observations: int
    name: str = "custom"
    scores: npt.NDArray[np.float64] | None = None

    def __post_init__(self) -> None:
        for field_name in ("center", "scale", "rotation", "sdev", "scores"):
            value = getattr(self, field_name)
            if value is not None:
                array = np.array(value, dtype=np.float64, copy=True)
                if not np.isfinite(array).all():
                    raise ValueError(f"{field_name} must contain only finite numbers")
                array.setflags(write=False)
                object.__setattr__(self, field_name, array)
        if self.center.shape != (3,) or self.scale.shape != (3,):
            raise ValueError("center and scale must each contain three features")
        if self.rotation.ndim != 2 or self.rotation.shape[0] != 3:
            raise ValueError("rotation must have three rows")
        if not 1 <= self.rotation.shape[1] <= 3:
            raise ValueError("rotation must have one to three components")
        if self.sdev.shape != (self.rotation.shape[1],):
            raise ValueError("sdev must contain one value per component")
        if np.any(self.scale <= 0) or np.any(self.sdev < 0):
            raise ValueError("scale must be positive and sdev nonnegative")
        if not np.isfinite(self.min_score) or self.n_observations < 2:
            raise ValueError("model needs a finite minimum and at least two observations")
        if self.rotation[1, 0] == 0:
            raise ValueError("PC1 has zero dCUB loading; copiotrophy direction is undefined")
        if self.scores is not None and self.scores.shape != (
            self.n_observations,
            self.rotation.shape[1],
        ):
            raise ValueError("scores shape does not match the fitted PCA")

    @property
    def direction(self) -> float:
        """Multiplier orienting PC1 so decreasing dCUB increases IoC."""
        return -float(np.sign(self.rotation[1, 0]))

    @property
    def variance_explained(self) -> npt.NDArray[np.float64]:
        """Fraction of total variance represented by each component."""
        variance = self.sdev**2
        return variance / variance.sum()

    @property
    def x(self) -> npt.NDArray[np.float64] | None:
        """Training PCA scores, when available (R ``prcomp$x``)."""
        return self.scores

    def transform(self, features: npt.ArrayLike) -> npt.NDArray[np.float64]:
        """Project already prepared [rCAZy, dCUB, nGenes] feature rows."""
        values = np.asarray(features, dtype=float)
        if values.ndim != 2 or values.shape[1] != 3:
            raise ValueError("features must be a two-dimensional array with three columns")
        if not np.isfinite(values).all():
            raise ValueError("features must contain only finite numbers")
        return ((values - self.center) / self.scale) @ self.rotation

    def to_dict(self, *, include_scores: bool = False) -> dict[str, Any]:
        """Return JSON-serializable model data without runtime dependencies."""
        result: dict[str, Any] = {
            "features": list(_FEATURES),
            "center": self.center.tolist(),
            "scale": self.scale.tolist(),
            "rotation": self.rotation.tolist(),
            "sdev": self.sdev.tolist(),
            "min_score": self.min_score,
            "n_observations": self.n_observations,
            "name": self.name,
        }
        if include_scores and self.scores is not None:
            result["scores"] = self.scores.tolist()
        return result

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> IoCModel:
        """Reconstruct a model produced by :meth:`to_dict`."""
        if tuple(data.get("features", _FEATURES)) != _FEATURES:
            raise ValueError("IoC feature order must be rCAZy, dCUB, nGenes")
        return cls(**{key: value for key, value in data.items() if key != "features"})


@dataclass(frozen=True)
class IoCFit:
    """Fitted index values and model, with R-compatible dictionary key access."""

    values: npt.NDArray[np.float64]
    model: IoCModel

    @property
    def IoC(self) -> npt.NDArray[np.float64]:  # noqa: N802
        """The R-compatible name for the fitted index values."""
        return self.values

    def __getitem__(self, key: str) -> Any:
        if key == "IoC":
            return self.values
        if key == "model":
            return self.model
        raise KeyError(key)

    def __iter__(self) -> Iterator[str]:
        return iter(("IoC", "model"))

    def __len__(self) -> int:
        return 2

    def keys(self) -> tuple[str, str]:
        """Return the R result field names; enables ``dict(result)``."""
        return ("IoC", "model")


def _features(
    x: Any,
    dcub: npt.ArrayLike | None,
    n_cazy: npt.ArrayLike | None,
    n_genes: npt.ArrayLike | None,
    completeness: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    if x is not None:
        if any(value is not None for value in (dcub, n_cazy, n_genes)):
            raise ValueError("provide either x or the three feature vectors, not both")
        try:
            vectors = [
                np.atleast_1d(np.asarray(x[key], dtype=float))
                for key in ("dCUB", "nCAZy", "nGenes")
            ]
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("x must contain numeric dCUB, nCAZy, and nGenes columns") from exc
    else:
        if any(value is None for value in (dcub, n_cazy, n_genes)):
            raise ValueError("provide x or all of dcub, n_cazy, and n_genes")
        vectors = [
            np.atleast_1d(np.asarray(value, dtype=float)) for value in (dcub, n_cazy, n_genes)
        ]
    if any(value.ndim != 1 for value in vectors):
        raise ValueError("feature columns must be one-dimensional")
    if len({len(value) for value in vectors}) != 1:
        raise ValueError("dCUB, nCAZy, and nGenes must have the same length")
    if not len(vectors[0]):
        raise ValueError("at least one observation is required")
    if not all(np.isfinite(value).all() for value in vectors):
        raise ValueError("feature columns must contain only finite numbers")
    if np.any(vectors[1] < 0) or np.any(vectors[2] <= 0):
        raise ValueError("nCAZy must be nonnegative and nGenes must be positive")
    comp = np.asarray(completeness, dtype=float)
    if comp.ndim > 1 or (comp.ndim == 1 and len(comp) not in (1, len(vectors[0]))):
        raise ValueError("completeness must be scalar or have one value per observation")
    if not np.isfinite(comp).all() or np.any(comp <= 0) or np.any(comp > 100):
        raise ValueError("completeness must be finite percentages in (0, 100]")
    # As in R, calculate the CAZyme fraction before correcting the gene count.
    return np.column_stack((vectors[1] / vectors[2], vectors[0], vectors[2] / (comp / 100)))


def fit_ioc(
    x: pd.DataFrame | Mapping[str, Any] | None = None,
    *,
    dcub: npt.ArrayLike | None = None,
    n_cazy: npt.ArrayLike | None = None,
    n_genes: npt.ArrayLike | None = None,
    completeness: npt.ArrayLike = 100,
) -> IoCFit:
    """Fit an IoC to a table with dCUB, nCAZy, and nGenes columns.

    Alternatively, supply the three keyword vectors. ``completeness`` is a
    scalar percentage or one percentage per row. Scaling uses sample standard
    deviations (``ddof=1``), exactly as R's ``prcomp(scale=TRUE)``.
    """
    values = _features(x, dcub, n_cazy, n_genes, completeness)
    if len(values) < 2:
        raise ValueError("at least two observations are needed to fit an IoC")
    center = values.mean(axis=0)
    scale = values.std(axis=0, ddof=1)
    if np.any(scale == 0):
        constant = ", ".join(np.asarray(_FEATURES)[scale == 0])
        raise ValueError(f"cannot scale constant feature(s): {constant}")
    standardized = (values - center) / scale
    _, singular_values, vt = np.linalg.svd(standardized, full_matrices=False)
    rotation = vt.T
    scores = standardized @ rotation
    direction = -float(np.sign(rotation[1, 0]))
    oriented = scores[:, 0] * direction
    model = IoCModel(
        center=center,
        scale=scale,
        rotation=rotation,
        sdev=singular_values / np.sqrt(len(values) - 1),
        min_score=float(np.min(oriented)),
        n_observations=len(values),
        scores=scores,
    )
    return IoCFit(values=oriented - model.min_score + 1, model=model)


def available_ioc_models() -> tuple[str, ...]:
    """Return the names of the five models distributed with gRodon 2.7.3."""
    return _BUILTINS


@lru_cache(maxsize=5)
def load_ioc_model(name: str) -> IoCModel:
    """Load an original R PCA model from losslessly exported numeric data."""
    if name not in _BUILTINS:
        raise ValueError(f"unknown IoC model {name!r}; choose one of {', '.join(_BUILTINS)}")
    resource = files("grodon").joinpath("data", "ioc_models.json")
    data = json.loads(resource.read_text(encoding="utf-8"))
    return IoCModel.from_dict(data["models"][name])


def predict_ioc(
    x: pd.DataFrame | Mapping[str, Any] | None = None,
    *,
    model: str | IoCModel | IoCFit | Mapping[str, Any],
    dcub: npt.ArrayLike | None = None,
    n_cazy: npt.ArrayLike | None = None,
    n_genes: npt.ArrayLike | None = None,
    completeness: npt.ArrayLike = 100,
) -> npt.NDArray[np.float64]:
    """Place observations on a fitted or bundled IoC, keeping its origin."""
    if isinstance(model, str):
        model = load_ioc_model(model)
    elif isinstance(model, IoCFit):
        model = model.model
    elif isinstance(model, Mapping):
        model = model["model"] if "model" in model else IoCModel.from_dict(model)
    if not isinstance(model, IoCModel):
        raise TypeError("model must be an IoCModel, fit_ioc result, or a built-in model name")
    values = _features(x, dcub, n_cazy, n_genes, completeness)
    return model.transform(values)[:, 0] * model.direction - model.min_score + 1


def fitIoC(x=None, dCUB=None, nCAZy=None, nGenes=None, completeness=100):  # noqa: N802, N803
    """Compatibility wrapper using the original R argument names."""
    return fit_ioc(x, dcub=dCUB, n_cazy=nCAZy, n_genes=nGenes, completeness=completeness)


def predictIoC(  # noqa: N802
    x=None,
    dCUB=None,
    nCAZy=None,
    nGenes=None,
    completeness=100,
    model=None,  # noqa: N803
):
    """Compatibility wrapper using the original R argument names."""
    return predict_ioc(
        x, model=model, dcub=dCUB, n_cazy=nCAZy, n_genes=nGenes, completeness=completeness
    )
