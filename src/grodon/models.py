"""Portable fitted models and R-compatible confidence intervals.

The shipped parameters come directly from gRodon 2.7.3. Predictions use the
confidence interval for the regression mean, as in R's ``predict.lm`` with
``interval='confidence'``; they are not observation prediction intervals.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache, lru_cache
from importlib.resources import files
from typing import Any

import numpy as np
from scipy import special, stats


def boxcox_transform(x: Any, lam: float, *, back_transform: bool = False):
    """Apply Box-Cox, including its continuous logarithmic limit at lambda=0.

    The inverse returns NaN outside its positive-real domain. This matches the
    bundled R models and avoids spurious powers of negative numbers for custom
    integer exponents. Callers can distinguish undefined extrapolations from
    very large valid estimates.
    """
    if not np.isfinite(lam):
        raise ValueError("Box-Cox lambda must be finite")
    values = np.asarray(x, dtype=float)
    if not back_transform and (np.any(values <= 0) or not np.all(np.isfinite(values))):
        raise ValueError("Box-Cox input must be finite and strictly positive")
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        result = special.inv_boxcox(values, lam) if back_transform else special.boxcox(values, lam)
    return float(result) if values.ndim == 0 else result


@dataclass(frozen=True)
class LinearModel:
    """A fitted linear model of Box-Cox transformed doubling time.

    ``covariance`` is the covariance of fitted coefficients, already scaled by
    residual variance. Predictor order includes ``'(Intercept)'`` explicitly.
    """

    name: str
    predictors: tuple[str, ...]
    coefficients: tuple[float, ...]
    covariance: tuple[tuple[float, ...], ...]
    df_residual: int
    sigma: float
    boxcox_lambda: float
    training_set: str = ""
    description: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "predictors", tuple(self.predictors))
        object.__setattr__(self, "coefficients", tuple(float(x) for x in self.coefficients))
        object.__setattr__(
            self, "covariance", tuple(tuple(float(x) for x in row) for row in self.covariance)
        )
        size = len(self.predictors)
        covariance = np.asarray(self.covariance)
        if not size or len(set(self.predictors)) != size or len(self.coefficients) != size:
            raise ValueError(
                "Model predictors and coefficients must have equal, nonzero lengths and unique names"
            )
        if covariance.shape != (size, size) or not np.all(np.isfinite(covariance)):
            raise ValueError("Coefficient covariance must be a finite square matrix")
        if not np.allclose(covariance, covariance.T, atol=1e-12):
            raise ValueError("Coefficient covariance must be symmetric")
        if np.min(np.linalg.eigvalsh(covariance)) < -1e-10 * max(
            1.0, float(np.linalg.norm(covariance))
        ):
            raise ValueError("Coefficient covariance must be positive semidefinite")
        if not np.all(np.isfinite(self.coefficients)) or not np.isfinite(self.boxcox_lambda):
            raise ValueError("Model coefficients and Box-Cox lambda must be finite")
        if (
            isinstance(self.df_residual, (bool, np.bool_))
            or not np.isfinite(self.df_residual)
            or self.df_residual <= 0
            or int(self.df_residual) != self.df_residual
        ):
            raise ValueError("Residual degrees of freedom must be a positive integer")
        object.__setattr__(self, "df_residual", int(self.df_residual))
        if not np.isfinite(self.sigma) or self.sigma < 0:
            raise ValueError("Residual sigma must be finite and nonnegative")

    def design_row(self, statistics: Mapping[str, float]) -> np.ndarray:
        """Build a numeric design row in the original fitted coefficient order."""
        missing = [
            name for name in self.predictors if name != "(Intercept)" and name not in statistics
        ]
        if missing:
            raise ValueError(f"Missing statistics for {self.name}: {', '.join(missing)}")
        row = np.array(
            [1.0 if name == "(Intercept)" else statistics[name] for name in self.predictors],
            dtype=float,
        )
        if not np.all(np.isfinite(row)):
            raise ValueError(f"Statistics used by model {self.name} must be finite")
        return row

    def predict_transformed(
        self, statistics: Mapping[str, float], confidence: float = 0.95
    ) -> tuple[float, float, float]:
        """Return fit, lower and upper confidence limits on the Box-Cox scale."""
        if not 0 < confidence < 1:
            raise ValueError("confidence must lie strictly between 0 and 1")
        row = self.design_row(statistics)
        fit = float(row @ np.asarray(self.coefficients))
        variance = float(row @ np.asarray(self.covariance) @ row)
        margin = float(
            stats.t.ppf((1 + confidence) / 2, self.df_residual) * np.sqrt(max(variance, 0.0))
        )
        return fit, fit - margin, fit + margin

    def predict(
        self, statistics: Mapping[str, float], confidence: float = 0.95
    ) -> tuple[float, float, float]:
        """Return doubling time and its confidence limits in hours."""
        transformed = self.predict_transformed(statistics, confidence)
        result = boxcox_transform(transformed, self.boxcox_lambda, back_transform=True)
        return float(result[0]), float(result[1]), float(result[2])


@lru_cache(maxsize=1)
def _model_bundle() -> dict[str, Any]:
    resource = files("grodon").joinpath("data", "growth_models.json")
    bundle = json.loads(resource.read_text(encoding="utf-8"))
    if bundle.get("schema_version") != 1:
        raise RuntimeError("Unsupported bundled growth model schema")
    return bundle


def available_models() -> tuple[str, ...]:
    """Names of all shipped R models, including individual-background variants."""
    return tuple(sorted(_model_bundle()["models"]))


@cache
def load_model(name: str) -> LinearModel:
    """Load one immutable fitted model by its original R object name."""
    try:
        fields = dict(_model_bundle()["models"][name])
    except KeyError as exc:
        raise ValueError(f"Unknown model {name!r}; use available_models() to list models") from exc
    fields.pop("lambda_source", None)
    return LinearModel(**fields)


def model_provenance() -> dict[str, Any]:
    """Return source version, revision and checksum for the bundled parameters."""
    return dict(_model_bundle()["source"])
