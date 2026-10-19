"""Train gRodon's published regression families on measured doubling times.

These utilities are independent of R. Named results replace R's unnamed
lists; the camelCase compatibility wrappers retain the original list order.
The original ``fitGCModels`` accidentally returns its new temperature model
twice. Here the legacy temperature model is returned in the legacy slot.
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.interpolate import make_interp_spline
from scipy.special import boxcox

from .io import GeneSequences, highly_expressed_mask, read_fasta
from .models import LinearModel

_RIBOSOMAL = re.compile(r"^(?!.*(?:methyl|hydroxy)).*0S ribosomal protein", re.IGNORECASE)


def _design(
    data: pd.DataFrame | Mapping[str, Any], predictors: Sequence[str]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    try:
        frame = pd.DataFrame(data)
    except (TypeError, ValueError) as exc:
        raise ValueError("training data must be a table with numeric columns") from exc
    required = ["d", *predictors]
    missing = [name for name in required if name not in frame]
    if missing:
        raise ValueError(f"missing training columns: {', '.join(missing)}")
    try:
        values = frame[required].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("training columns must be numeric") from exc
    if not np.isfinite(values).all():
        raise ValueError(
            "training data must contain only finite numbers; remove missing values first"
        )
    response = values[:, 0]
    if np.any(response <= 0):
        raise ValueError("doubling times in column d must be positive")
    design = np.column_stack((np.ones(len(values)), values[:, 1:]))
    if len(values) <= design.shape[1]:
        raise ValueError("training needs more observations than fitted coefficients")
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("training predictors are rank deficient (constant or collinear columns)")
    return design, response


def estimate_boxcox_lambda(
    data: pd.DataFrame | Mapping[str, Any],
    predictors: Sequence[str],
    *,
    lambda_grid: npt.ArrayLike | None = None,
    interpolate: bool = True,
) -> float:
    """Choose the Box-Cox exponent by regression profile likelihood.

    Defaults reproduce ``MASS::boxcox`` as used by gRodon's fitting helpers:
    evaluate 41 exponents from -2 to 2, interpolate the profile with an FMM
    cubic spline at 100 positions, then select the first maximum. Set
    ``interpolate=False`` to select directly from ``lambda_grid``.
    """
    design, response = _design(data, predictors)
    grid = (
        np.linspace(-2.0, 2.0, 41) if lambda_grid is None else np.asarray(lambda_grid, dtype=float)
    )
    if grid.ndim != 1 or not len(grid) or not np.isfinite(grid).all():
        raise ValueError("lambda_grid must be a nonempty, finite, one-dimensional vector")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("lambda_grid must be strictly increasing")
    log_scaled = np.log(response) - np.log(response).mean()
    with np.errstate(over="ignore", under="ignore"):
        scaled = np.exp(log_scaled)
    if not np.isfinite(scaled).all() or np.any(scaled <= 0):
        raise ValueError("doubling times span too wide a range for Box-Cox scaling")
    q, _ = np.linalg.qr(design, mode="reduced")
    likelihood: npt.NDArray[np.float64] = np.empty(len(grid), dtype=float)
    for index, exponent in enumerate(grid):
        transformed = boxcox(scaled, exponent)
        if not np.isfinite(transformed).all():
            raise ValueError(
                "Box-Cox search overflowed; narrow lambda_grid or check doubling-time range"
            )
        with np.errstate(over="ignore", invalid="ignore"):
            residuals = transformed - q @ (q.T @ transformed)
            residual_sum = float(residuals @ residuals)
        if not np.isfinite(residual_sum):
            raise ValueError(
                "Box-Cox residuals overflowed; narrow lambda_grid or check doubling-time range"
            )
        if residual_sum <= np.finfo(float).tiny:
            raise ValueError("Box-Cox likelihood is undefined for an exact regression fit")
        likelihood[index] = -len(response) / 2 * np.log(residual_sum)
    if interpolate and len(grid) < 100:
        if len(grid) < 4:
            raise ValueError("interpolation requires at least four lambda_grid values")

        def third_derivative(x: np.ndarray, y: np.ndarray) -> float:
            differences = y.copy()
            for order in range(1, 4):
                differences = np.diff(differences) / (x[order:] - x[:-order])
            return float(6 * differences[0])

        # FMM boundary conditions match third derivatives of the endpoint
        # cubics through the first/last four observations. scipy constructs
        # the spline directly from these mathematical boundary conditions.
        boundary = (
            [(3, third_derivative(grid[:4], likelihood[:4]))],
            [(3, third_derivative(grid[-4:], likelihood[-4:]))],
        )
        spline = make_interp_spline(grid, likelihood, bc_type=boundary)
        grid = np.linspace(grid[0], grid[-1], 100)
        likelihood = spline(grid)
    return float(grid[np.argmax(likelihood)])


def fit_linear_model(
    data: pd.DataFrame | Mapping[str, Any],
    predictors: Sequence[str],
    *,
    boxcox_lambda: float,
    name: str = "custom",
) -> LinearModel:
    """Fit an intercept and numeric predictors to transformed doubling times.

    The coefficient covariance and residual degrees of freedom support the
    same t-based confidence intervals as R's ``predict.lm``. Unlike R's
    silent row omission, missing or nonfinite training values raise an error.
    """
    if not np.isfinite(boxcox_lambda):
        raise ValueError("boxcox_lambda must be finite")
    design, response = _design(data, predictors)
    transformed = boxcox(response, boxcox_lambda)
    if not np.isfinite(transformed).all():
        raise ValueError("Box-Cox transformed doubling times overflowed; check scale and exponent")
    coefficients, _, _, _ = np.linalg.lstsq(design, transformed, rcond=None)
    residuals = transformed - design @ coefficients
    df_residual = len(response) - design.shape[1]
    variance = float(residuals @ residuals / df_residual)
    inverse_design = np.linalg.pinv(design)
    covariance = (inverse_design @ inverse_design.T) * variance
    return LinearModel(
        name=name,
        predictors=("(Intercept)", *predictors),
        coefficients=tuple(float(value) for value in coefficients),
        covariance=tuple(tuple(float(value) for value in row) for row in covariance),
        df_residual=df_residual,
        sigma=float(np.sqrt(variance)),
        boxcox_lambda=float(boxcox_lambda),
        training_set="user supplied",
    )


def fit_models(
    stat_data: pd.DataFrame | Mapping[str, Any],
    stat_data_extremo: pd.DataFrame | Mapping[str, Any],
    *,
    lambda_value: float | None = None,
    eukaryotic_lambda: float | None = None,
) -> dict[str, LinearModel | float]:
    """Fit full, partial, metagenome-v1, and eukaryotic model families.

    Both tables contain measured hours in ``d`` and codon statistics in
    ``CUBHE``, ``ConsistencyHE``, and ``CPB``. The temperature table additionally
    has ``OGT`` in degrees Celsius. Lambda is selected on ``stat_data`` only.
    """
    full = ("CUBHE", "ConsistencyHE", "CPB")
    lam = estimate_boxcox_lambda(stat_data, full) if lambda_value is None else float(lambda_value)
    euk_lam = (
        estimate_boxcox_lambda(stat_data, ("CUBHE",))
        if eukaryotic_lambda is None
        else float(eukaryotic_lambda)
    )
    result: dict[str, LinearModel | float] = {}
    for stem, predictors, exponent in (
        ("base", full, lam),
        ("partial", full[:2], lam),
        ("meta", full[:1], lam),
        ("euk", full[:1], euk_lam),
    ):
        name = f"gRodon_model_{stem}"
        result[name] = fit_linear_model(stat_data, predictors, boxcox_lambda=exponent, name=name)
        temp_name = "gRodon_model_temp" if stem == "base" else name + "_temp"
        result[temp_name] = fit_linear_model(
            stat_data_extremo, (*predictors, "OGT"), boxcox_lambda=exponent, name=temp_name
        )
    result["lambda_milc"] = lam
    result["lambda_milc_euk"] = euk_lam
    return result


def fit_gc_models(
    stat_data: pd.DataFrame | Mapping[str, Any],
    stat_data_extremo: pd.DataFrame | Mapping[str, Any],
    *,
    lambda_value: float | None = None,
    legacy_lambda: float | None = None,
) -> dict[str, LinearModel | float]:
    """Fit GC-corrected, uncorrected, and legacy metagenome model families.

    ``dCUB`` is recomputed as ``(CUB - CUBHE) / CUB`` as in the R package.
    The legacy temperature model is returned correctly, fixing an upstream
    typo in the R helper's seventh return element.
    """
    tables = []
    for original in (stat_data, stat_data_extremo):
        table = pd.DataFrame(original).copy()
        if "CUB" not in table or "CUBHE" not in table:
            raise ValueError("GC model training requires CUB and CUBHE columns")
        if np.any(table["CUB"].to_numpy(dtype=float) == 0):
            raise ValueError("CUB must be nonzero to calculate dCUB")
        table["dCUB"] = (table["CUB"] - table["CUBHE"]) / table["CUB"]
        tables.append(table)
    normal, extremo = tables
    lam = (
        estimate_boxcox_lambda(normal, ("dCUB", "GCdiv"))
        if lambda_value is None
        else float(lambda_value)
    )
    old_lam = (
        estimate_boxcox_lambda(normal, ("CUBHE",))
        if legacy_lambda is None
        else float(legacy_lambda)
    )
    result: dict[str, LinearModel | float] = {}
    for stem, predictors, exponent in (
        ("meta_model", ("dCUB", "GCdiv"), lam),
        ("meta_nogc_model", ("dCUB",), lam),
        ("meta_old_model", ("CUBHE",), old_lam),
    ):
        for suffix, data, terms in (
            ("base", normal, predictors),
            ("temp", extremo, (*predictors, "OGT")),
        ):
            name = f"{stem}_{suffix}"
            result[name] = fit_linear_model(data, terms, boxcox_lambda=exponent, name=name)
    result["lambda_meta"] = lam
    result["lambda_oldmeta"] = old_lam
    return result


def fit_trim_models(
    stat_data: pd.DataFrame | Mapping[str, Any],
    stat_data_extremo: pd.DataFrame | Mapping[str, Any],
    *,
    lambda_value: float | None = None,
) -> dict[str, LinearModel | float]:
    """Fit trimmed-gene models with MILC, ENCprime, B, SCUO, and MCB."""
    terms = ("MILC", "ENCprime", "B", "SCUO", "MCB")
    lam = estimate_boxcox_lambda(stat_data, terms) if lambda_value is None else float(lambda_value)
    return {
        "gRodon_model_trim": fit_linear_model(
            stat_data, terms, boxcox_lambda=lam, name="gRodon_model_trim"
        ),
        "gRodon_model_trim_temp": fit_linear_model(
            stat_data_extremo, (*terms, "OGT"), boxcox_lambda=lam, name="gRodon_model_trim_temp"
        ),
        "lambda_cub": lam,
    }


def get_statistics(
    gene_file: str | Path | None = None,
    *,
    genes: Mapping[str, Any] | Sequence[Any] | None = None,
    genetic_code: int | str = 11,
    bg: str = "all",
    fragments: bool | None = False,
    exclude_short: bool = False,
    trimlen: int | None = None,
    trimside: str = "start",
    all_metrics: bool = False,
    random_state: Any = None,
) -> dict[str, Any]:
    """Calculate codon statistics using ribosomal annotations for HE genes.

    Supply an annotated FASTA file (optionally gzip-compressed), a mapping
    of descriptions to sequences, or SeqRecord objects. The annotation match
    follows R: ``0S ribosomal protein``, excluding methyl/hydroxy annotations.
    """
    from .statistics import get_codon_statistics

    if gene_file is not None:
        if genes is not None:
            raise ValueError("provide gene_file or genes, not both")
        path = Path(gene_file)
        records = read_fasta(path)
        names = list(records.names)
        sequences = list(records)
        label = path.name
    elif isinstance(genes, GeneSequences):
        names = list(genes.names)
        sequences = list(genes)
        label = "genome provided directly"
    elif isinstance(genes, Mapping):
        names = [str(name) for name in genes]
        sequences = [str(getattr(value, "seq", value)) for value in genes.values()]
        label = "genome provided directly"
    elif genes is not None:
        named_records = list(genes)
        names = [
            str(getattr(record, "description", getattr(record, "id", "")))
            for record in named_records
        ]
        sequences = [str(getattr(record, "seq", record)) for record in named_records]
        label = "genome provided directly"
    else:
        raise ValueError("provide gene_file or annotated genes")
    highly_expressed = highly_expressed_mask(
        GeneSequences(tuple(sequences), tuple(names)), pattern=_RIBOSOMAL.pattern
    )
    if not any(highly_expressed):
        raise ValueError("no highly expressed genes: expected '0S ribosomal protein' annotations")
    if bg not in ("all", "individual"):
        raise ValueError("bg must be 'all' or 'individual'")
    result: dict[str, Any] = dict(
        get_codon_statistics(
            sequences,
            highly_expressed,
            genetic_code=genetic_code,
            background=bg,
            fragments=fragments,
            exclude_short=exclude_short,
            trimlen=trimlen,
            trimside=trimside,
            all_metrics=all_metrics,
            random_state=random_state,
        )
    )
    result["File"] = label
    return result


def get_statistics_batch(
    directory: str | Path,
    *,
    workers: int = 1,
    pattern: str = "*",
    errors: str = "raise",
    **kwargs: Any,
) -> pd.DataFrame:
    """Calculate statistics for sorted files in a directory, optionally in parallel.

    ``pattern`` can restrict inputs, for example ``'*.fna.gz'``. By default any
    failed file raises an exception. ``errors='skip'`` warns and omits it.
    Randomized calculations receive independent deterministic streams when
    an integer ``random_state`` is supplied, independent of worker count.
    """
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError("workers must be a positive integer")
    if errors not in ("raise", "skip"):
        raise ValueError("errors must be 'raise' or 'skip'")
    folder = Path(directory)
    if not folder.is_dir():
        raise NotADirectoryError(folder)
    paths = sorted(path for path in folder.glob(pattern) if path.is_file())
    seed = kwargs.pop("random_state", None)
    seeds: list[int | None]
    if seed is None:
        seeds = [None] * len(paths)
    elif isinstance(seed, (int, np.integer)) and not isinstance(seed, bool):
        seeds = [
            int(child.generate_state(1)[0])
            for child in np.random.SeedSequence(int(seed)).spawn(len(paths))
        ]
    else:
        raise ValueError("batch random_state must be an integer or None")

    def run(item: tuple[Path, int | None]) -> dict[str, Any] | None:
        path, state = item
        try:
            return get_statistics(path, random_state=state, **kwargs)
        except Exception as exc:
            if errors == "raise":
                raise
            warnings.warn(f"skipping {path.name}: {exc}", RuntimeWarning, stacklevel=2)
            return None

    if workers == 1:
        results = [run(item) for item in zip(paths, seeds, strict=True)]
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(run, zip(paths, seeds, strict=True)))
    return pd.DataFrame([result for result in results if result is not None])


def fitModels(stat_data, stat_data_extremo):  # noqa: N802
    """R-compatible ordered result for :func:`fit_models`."""
    result = fit_models(stat_data, stat_data_extremo)
    order = (
        "gRodon_model_base",
        "gRodon_model_temp",
        "gRodon_model_partial",
        "gRodon_model_partial_temp",
        "gRodon_model_meta",
        "gRodon_model_meta_temp",
        "lambda_milc",
        "gRodon_model_euk",
        "gRodon_model_euk_temp",
        "lambda_milc_euk",
    )
    return [result[key] for key in order]


def fitGCModels(stat_data, stat_data_extremo):  # noqa: N802
    """R-compatible ordering with the legacy-temperature return bug corrected."""
    result = fit_gc_models(stat_data, stat_data_extremo)
    order = (
        "meta_model_base",
        "meta_model_temp",
        "meta_nogc_model_base",
        "meta_nogc_model_temp",
        "lambda_meta",
        "meta_old_model_base",
        "meta_old_model_temp",
        "lambda_oldmeta",
    )
    return [result[key] for key in order]


def fitTrimModels(stat_data, stat_data_extremo):  # noqa: N802
    """R-compatible ordered result for :func:`fit_trim_models`."""
    result = fit_trim_models(stat_data, stat_data_extremo)
    return [result[key] for key in ("gRodon_model_trim", "gRodon_model_trim_temp", "lambda_cub")]


def getStatistics(  # noqa: N802
    gene_file=None,
    genes=None,
    genetic_code="11",
    bg="all",
    fragments=False,
    exclude_short=False,
    trimlen=None,
    trimside="start",
    all_metrics=False,
    **kwargs,
):
    """R argument-order compatibility for :func:`get_statistics`."""
    return get_statistics(
        gene_file,
        genes=genes,
        genetic_code=genetic_code,
        bg=bg,
        fragments=fragments,
        exclude_short=exclude_short,
        trimlen=trimlen,
        trimside=trimside,
        all_metrics=all_metrics,
        **kwargs,
    )


def getStatisticsBatch(  # noqa: N802
    directory,
    genetic_code="11",
    mc_cores=1,
    bg="all",
    fragments=False,
    exclude_short=False,
    trimlen=None,
    trimside="start",
    all_metrics=False,
    **kwargs,
):
    """R argument-order compatibility; ``mc_cores`` replaces R's ``mc.cores``."""
    workers = kwargs.pop("workers", kwargs.pop("mc.cores", mc_cores))
    return get_statistics_batch(
        directory,
        workers=workers,
        genetic_code=genetic_code,
        bg=bg,
        fragments=fragments,
        exclude_short=exclude_short,
        trimlen=trimlen,
        trimside=trimside,
        all_metrics=all_metrics,
        **kwargs,
    )


getStatisticsBatchWindows = getStatisticsBatch
