"""Growth prediction from coding sequences or precomputed codon statistics."""

from __future__ import annotations

import warnings
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any

import numpy as np

from .models import LinearModel, load_model

PREDICTION_MODES = (
    "full",
    "partial",
    "metagenome_v1",
    "metagenome_v2",
    "metagenome_150bp",
    "metagenome_240bp",
    "eukaryote",
    "metagenome_euk",
    "meta_testing",
    "meta_nogc_testing",
)
_EUK_MODES = frozenset(("eukaryote", "metagenome_euk"))
_SHORT_MODES = frozenset(("metagenome_150bp", "metagenome_240bp"))
_META_MODES = frozenset(PREDICTION_MODES) - {"full", "partial", "eukaryote"}


class GrowthPredictionWarning(UserWarning):
    """A prediction setting or estimate warrants the caller's attention."""


@dataclass(frozen=True, eq=False)
class PredictionResult(Mapping[str, float | int]):
    """Immutable result with R-compatible keys and Python-friendly attributes.

    ``dict(result)`` yields the original statistics plus ``d``, ``LowerCI`` and
    ``UpperCI``. Model selection metadata is exposed as attributes, outside the
    numerical mapping, so the original result schema stays compatible.
    """

    statistics: Mapping[str, float | int]
    model_name: str
    mode: str
    training_set: str
    background: str
    confidence: float = 0.95

    def __post_init__(self) -> None:
        object.__setattr__(self, "statistics", MappingProxyType(dict(self.statistics)))

    def __getitem__(self, key: str) -> float | int:
        return self.statistics[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self.statistics)

    def __len__(self) -> int:
        return len(self.statistics)

    @property
    def doubling_time(self) -> float:
        """Estimated doubling time in hours."""
        return float(self["d"])

    @property
    def d(self) -> float:
        return self.doubling_time

    @property
    def lower_ci(self) -> float:
        return float(self["LowerCI"])

    @property
    def upper_ci(self) -> float:
        return float(self["UpperCI"])

    def to_dict(self) -> dict[str, float | int]:
        """Copy the numeric result to a mutable dictionary."""
        return dict(self.statistics)


def _temperature(value: Any) -> float | None:
    if value is None or (isinstance(value, str) and value.lower() == "none"):
        return None
    if isinstance(value, (bool, np.bool_)):
        raise ValueError("temperature must be a numeric temperature in Celsius or None")
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("temperature must be a numeric temperature in Celsius or None") from exc
    if np.isnan(numeric):
        warnings.warn(
            "Missing temperature; using the model without temperature.",
            GrowthPredictionWarning,
            stacklevel=3,
        )
        return None
    if not np.isfinite(numeric):
        raise ValueError("temperature must be finite")
    return numeric


def _settings(mode: str, training_set: str, background: str) -> tuple[str, str]:
    if mode not in PREDICTION_MODES:
        raise ValueError(f"Invalid mode {mode!r}; choose one of {', '.join(PREDICTION_MODES)}")
    if background not in {"all", "individual"}:
        raise ValueError("background must be 'all' or 'individual'")
    aliases = {"vs": "vs", "original": "vs", "madin": "madin", "aoa_nob": "AOA_NOB"}
    if mode in _EUK_MODES:
        training_set = "eukaryote"
    else:
        try:
            training_set = aliases[training_set.lower()]
        except (KeyError, AttributeError) as exc:
            raise ValueError(
                "training_set must be 'vs' (or 'original'), 'madin', or 'AOA_NOB'"
            ) from exc
    if mode in _SHORT_MODES or mode == "metagenome_v2":
        if training_set != "madin":
            warnings.warn(
                f"Training set automatically set to 'madin' for {mode} mode.",
                GrowthPredictionWarning,
                stacklevel=3,
            )
            training_set = "madin"
    if mode in _SHORT_MODES:
        background = "individual"
    if background == "individual" and mode in {"full", "partial"}:
        raise ValueError(f"Mode {mode!r} has no model fitted with an individual gene background")
    if background == "individual" and mode == "metagenome_v1" and training_set != "madin":
        raise ValueError("Individual-background metagenome_v1 models require training_set='madin'")
    return training_set, background


def _select_model(
    statistics: Mapping[str, float],
    mode: str,
    temperature: float | None,
    training_set: str,
    background: str,
) -> LinearModel:
    suffix = "_temp" if temperature is not None else ""
    if mode in _SHORT_MODES:
        stem = "temp" if temperature is not None else "base"
        name = f"gRodon_model_{stem}_t{150 if mode == 'metagenome_150bp' else 240}"
    elif mode in _EUK_MODES:
        if mode == "metagenome_euk" and background == "all":
            if "ConsistencyHE" not in statistics or not np.isfinite(statistics["ConsistencyHE"]):
                raise ValueError("metagenome_euk requires finite ConsistencyHE statistics")
            if statistics["ConsistencyHE"] >= 0.6:
                raise ValueError(
                    "ConsistencyHE >= 0.6 requires recomputing individual-background statistics; use predict_growth or supply background='individual'"
                )
        stem = "temp" if temperature is not None else "base"
        name = f"gRodon_model_{stem}_euk" + ("_i" if background == "individual" else "")
    elif mode in {"meta_testing", "meta_nogc_testing"}:
        stem = "newmeta_nogc" if mode == "meta_nogc_testing" else "newmeta"
        name = f"gRodon_model_{stem}{suffix}" + ("_i" if background == "individual" else "")
    elif mode == "metagenome_v2":
        if background == "individual":
            name = f"gRodon_model_newmeta{suffix}_i"
        else:
            if "ConsistencyHE" not in statistics or not np.isfinite(statistics["ConsistencyHE"]):
                raise ValueError("metagenome_v2 requires finite ConsistencyHE statistics")
            if statistics["ConsistencyHE"] >= 0.6:
                raise ValueError(
                    "ConsistencyHE >= 0.6 requires recomputing individual-background statistics; use predict_growth or supply background='individual'"
                )
            name = f"gRodon_model_meta{suffix}_madin"
    else:
        if mode == "full":
            stem = "temp" if temperature is not None else "base"
        else:
            stem = ("partial" if mode == "partial" else "meta") + suffix
        dataset = {"vs": "", "madin": "_madin", "AOA_NOB": "_AOANOB"}[training_set]
        name = f"gRodon_model_{stem}{dataset}" + ("_i" if background == "individual" else "")
    return load_model(name)


def predict_from_statistics(
    statistics: Mapping[str, float | int],
    *,
    mode: str = "full",
    temperature: float | str | None = None,
    training_set: str = "madin",
    background: str = "all",
    confidence: float = 0.95,
) -> PredictionResult:
    """Predict from previously calculated statistics without reading sequences.

    Statistics must match the chosen background and short-read preprocessing.
    For metagenome_v2/euk with pooled ConsistencyHE >= 0.6, call
    :func:`predict_growth` to recalculate the individual-background statistics,
    or supply statistics already calculated with ``background='individual'``.
    """
    training_set, background = _settings(mode, training_set, background)
    ogt = _temperature(temperature)
    values = dict(statistics)
    if "CUB" in values and "CUBHE" in values:
        values["dCUB"] = (
            (values["CUB"] - values["CUBHE"]) / values["CUB"] if values["CUB"] else float("nan")
        )
    if ogt is not None:
        values["OGT"] = ogt
    model = _select_model(values, mode, ogt, training_set, background)
    fit, lower, upper = model.predict(values, confidence)
    values.update(d=fit, LowerCI=lower, UpperCI=upper)
    if not np.isfinite(fit):
        warnings.warn(
            "Estimated doubling time is outside the inverse Box-Cox domain; the fitted model is extrapolating beyond a finite estimate.",
            GrowthPredictionWarning,
            stacklevel=2,
        )
    elif fit > (40 if mode in _EUK_MODES else 5):
        threshold = 40 if mode in _EUK_MODES else 5
        warnings.warn(
            f"Estimated doubling time >{threshold} hours. The codon usage signal saturates around this range; consider reporting '>{threshold} hours'.",
            GrowthPredictionWarning,
            stacklevel=2,
        )
    return PredictionResult(values, model.name, mode, model.training_set, background, confidence)


def predict_growth(
    genes: Any,
    highly_expressed: Any,
    *,
    mode: str = "full",
    temperature: float | str | None = None,
    training_set: str = "madin",
    depth_of_coverage: Any = None,
    fragments: bool | None = False,
    exclude_short: bool = False,
    genetic_code: int | str | None = None,
    n_le: int = 100,
    background: str = "all",
    bg: str | None = None,
    random_state: Any = None,
    confidence: float = 0.95,
) -> PredictionResult:
    """Estimate doubling time in hours from in-frame coding sequences.

    Parameters follow gRodon 2.7.3. Use ``None`` for R's ``NA`` fragments value.
    ``random_state`` fixes individual-background sampling for reproducibility.
    All published models, their original coefficients and confidence intervals
    are included. Explicitly selected original and AOA/NOB training sets are
    honored; upstream model dispatch and inverse-transform bugs are corrected.
    """
    from .statistics import get_codon_statistics

    if isinstance(genes, (str, bytes)):
        raise TypeError("genes must be a collection of sequences; use read_fasta() to read a file")
    if bg is not None:
        if background != "all" and background != bg:
            raise ValueError("background and bg disagree; provide only one")
        background = bg
    training_set, background = _settings(mode, training_set, background)
    ogt = _temperature(temperature)
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie strictly between 0 and 1")
    mask = np.asarray(list(highly_expressed))
    if mask.ndim != 1 or mask.dtype.kind != "b":
        raise ValueError("highly_expressed must be a one-dimensional sequence of booleans")
    if mask.sum() < 10:
        warnings.warn(
            "Less than 10 highly expressed genes provided; performance may suffer.",
            GrowthPredictionWarning,
            stacklevel=2,
        )
    if depth_of_coverage is not None and mode not in _META_MODES:
        warnings.warn(
            "Ignoring depth_of_coverage because the mode is not a metagenome mode.",
            GrowthPredictionWarning,
            stacklevel=2,
        )
        depth_of_coverage = None
    if depth_of_coverage is None and mode in {"metagenome_v1", "metagenome_v2", "metagenome_euk"}:
        warnings.warn(
            "Provide depth_of_coverage for a more realistic average community growth rate.",
            GrowthPredictionWarning,
            stacklevel=2,
        )
    if ogt is None and mode in _EUK_MODES:
        warnings.warn(
            "An optimal growth temperature improves eukaryotic growth predictions.",
            GrowthPredictionWarning,
            stacklevel=2,
        )
    if genetic_code is None:
        genetic_code = 1 if mode in _EUK_MODES else 11
    else:
        warnings.warn(
            "Models were trained with genetic code 1 for eukaryotes and 11 for prokaryotes.",
            GrowthPredictionWarning,
            stacklevel=2,
        )
    materialized = genes if isinstance(genes, Mapping) else list(genes)
    common = {
        "fragments": fragments,
        "depth_of_coverage": depth_of_coverage,
        "genetic_code": genetic_code,
        "n_le": n_le,
        "random_state": random_state,
    }
    if mode in _SHORT_MODES:
        common.update(
            fragments=None,
            exclude_short=exclude_short,
            trimlen=150 if mode == "metagenome_150bp" else 240,
            trimside="start",
            all_metrics=True,
        )
    values = get_codon_statistics(materialized, mask, background=background, **common)
    if (
        mode in {"metagenome_v2", "metagenome_euk"}
        and background == "all"
        and values["ConsistencyHE"] >= 0.6
    ):
        background = "individual"
        values = get_codon_statistics(materialized, mask, background=background, **common)
    return predict_from_statistics(
        values,
        mode=mode,
        temperature=ogt,
        training_set=training_set,
        background=background,
        confidence=confidence,
    )


def predictGrowth(
    genes: Any,
    highly_expressed: Any,
    mode: str = "full",
    temperature: float | str | None = "none",
    training_set: str = "madin",
    depth_of_coverage: Any = None,
    fragments: bool | None = False,
    exclude_short: bool = False,
    genetic_code: int | str | None = None,
    n_le: int = 100,
    bg: str = "all",
    **kwargs: Any,
) -> PredictionResult:
    """R-compatible positional interface; see :func:`predict_growth`."""
    return predict_growth(
        genes,
        highly_expressed,
        mode=mode,
        temperature=temperature,
        training_set=training_set,
        depth_of_coverage=depth_of_coverage,
        fragments=fragments,
        exclude_short=exclude_short,
        genetic_code=genetic_code,
        n_le=n_le,
        bg=bg,
        **kwargs,
    )
