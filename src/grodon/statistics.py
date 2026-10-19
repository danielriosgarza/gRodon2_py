"""Sequence filtering and genome-level codon statistics."""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .codon import (
    Metric,
    _genetic_code,
    codon_counts,
    codon_pair_bias,
    codon_usage,
    normalize_genes,
)


class SequenceFilterWarning(UserWarning):
    """Some input sequences were excluded from a codon-usage calculation."""


@dataclass(frozen=True)
class FilteredSequences:
    """Sequences and aligned metadata surviving a quality filter."""

    genes: tuple[str, ...]
    highly_expressed: NDArray[np.bool_]
    depth_of_coverage: NDArray[np.float64] | None
    filtered: int
    indices: NDArray[np.int64]


def weighted_median(values: ArrayLike, weights: ArrayLike) -> float:
    """Interpolated weighted median, matching matrixStats::weightedMedian.

    Interpolation uses cumulative-weight midpoints. Zero-weight observations are
    excluded. Negative, nonfinite or all-zero weights are rejected.
    """
    x, w = np.asarray(values, dtype=float), np.asarray(weights, dtype=float)
    if x.ndim != 1 or w.shape != x.shape:
        raise ValueError("values and weights must be equally sized one-dimensional arrays")
    if not np.all(np.isfinite(w)) or np.any(w < 0) or not np.any(w > 0):
        raise ValueError("weights must be finite, nonnegative and have a positive total")
    positive = w > 0
    x, w = x[positive], w[positive]
    if np.any(np.isnan(x)):
        return float("nan")
    order = np.argsort(x, kind="stable")
    x, w = x[order], w[order]
    w = w / w.max()
    midpoints = (np.cumsum(w) - w / 2) / w.sum()
    return float(np.interp(0.5, midpoints, x))


def filter_sequences(
    genes: Iterable[object],
    highly_expressed: Sequence[bool] | NDArray[np.bool_],
    *,
    length_threshold: int = 240,
    depth_of_coverage: ArrayLike | None = None,
) -> FilteredSequences:
    """Keep unambiguous DNA longer than a threshold and divisible by three."""
    sequences = normalize_genes(genes)
    he = np.asarray(highly_expressed)
    if he.ndim != 1 or len(he) != len(sequences):
        raise ValueError("highly_expressed must have one boolean per gene")
    if he.dtype.kind != "b":
        raise TypeError("highly_expressed must contain booleans")
    if not np.any(he):
        raise ValueError("No highly expressed genes were specified")
    if (
        isinstance(length_threshold, bool)
        or not isinstance(length_threshold, (int, np.integer))
        or length_threshold < 0
    ):
        raise ValueError("length_threshold must be a nonnegative integer")
    depth = None
    if depth_of_coverage is not None:
        depth = np.asarray(depth_of_coverage, dtype=float)
        if depth.ndim != 1 or len(depth) != len(sequences):
            raise ValueError("depth_of_coverage must have one value per gene")
        if not np.all(np.isfinite(depth)) or np.any(depth < 0):
            raise ValueError("depth_of_coverage must be finite and nonnegative")
    valid = np.array(
        [
            len(gene) > length_threshold and len(gene) % 3 == 0 and set(gene) <= set("ACGT")
            for gene in sequences
        ],
        dtype=bool,
    )
    filtered = int((~valid).sum())
    if filtered:
        warnings.warn(
            f"Filtered {filtered} of {len(sequences)} genes: lengths must exceed "
            f"{length_threshold} bp, be divisible by three, and contain only A/C/G/T",
            SequenceFilterWarning,
            stacklevel=2,
        )
    if not np.any(he[valid]):
        raise ValueError("No highly expressed genes remain after filtering")
    if depth is not None:
        depth = depth[valid]
        if not np.any(depth > 0) or not np.any(depth[he[valid]] > 0):
            raise ValueError("Retained genes and highly expressed genes need positive coverage")
    return FilteredSequences(
        genes=tuple(gene for gene, keep in zip(sequences, valid, strict=True) if keep),
        highly_expressed=he[valid].copy(),
        depth_of_coverage=depth,
        filtered=filtered,
        indices=np.flatnonzero(valid),
    )


def _trim_sequences(
    genes: Sequence[str], trimlen: int, side: str, rng: np.random.Generator
) -> tuple[str, ...]:
    result = []
    for gene in genes:
        if len(gene) <= trimlen:
            result.append(gene)
            continue
        if side == "start":
            start = 0
        elif side == "end":
            start = len(gene) - trimlen
        else:
            start = 3 * int(rng.integers((len(gene) - trimlen) // 3 + 1))
        result.append(gene[start : start + trimlen])
    return tuple(result)


def _individual_statistics(
    genes: Sequence[str],
    he: NDArray[np.bool_],
    weights: NDArray[np.float64] | None,
    n_le: int,
    metrics: Sequence[Metric],
    rng: np.random.Generator,
) -> dict[str, tuple[float, float]]:
    low_indices = np.flatnonzero(~he)
    if not len(low_indices):
        raise ValueError("Individual background requires at least one non-highly-expressed gene")
    if weights is not None and not np.any(weights[low_indices] > 0):
        raise ValueError(
            "Non-highly-expressed genes need positive coverage for individual background"
        )
    if weights is not None:
        # Zero-weight genes cannot contribute to the weighted low-expression
        # estimate. Excluding them avoids randomly sampling an all-zero subset.
        low_indices = low_indices[weights[low_indices] > 0]
    sampled = rng.choice(low_indices, size=n_le, replace=True)
    high_indices = np.flatnonzero(he)
    # The reference includes the original sequence and 100 nucleotide shuffles,
    # exactly as singleMILC() does. Work on base integers, avoiding string copies.
    values = {method: np.empty(len(high_indices) + n_le) for method in metrics}
    for row, index in enumerate(np.concatenate([high_indices, sampled])):
        gene = genes[int(index)]
        bases = np.frombuffer(gene.encode("ascii"), dtype=np.uint8)
        base_map = np.zeros(256, dtype=np.int64)
        base_map[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)
        encoded = base_map[bases]
        original = codon_counts([gene])
        pooled = original[0].copy()
        for _ in range(100):
            shuffled = rng.permutation(encoded).reshape(-1, 3)
            indices = shuffled @ np.array([16, 4, 1])
            pooled += np.bincount(indices, minlength=64)
        for method in metrics:
            values[method][row] = codon_usage(
                original,
                method=method,
                reference_counts=pooled,
                # R singleMILC() uses coRdon's default code 1, independently of
                # the public genetic_code argument. Preserve calibration parity.
                genetic_code=1,
            )[0]
    result: dict[str, tuple[float, float]] = {}
    for method, metric_values in values.items():
        high_values, low_values = (
            metric_values[: len(high_indices)],
            metric_values[len(high_indices) :],
        )
        if weights is None:
            result[method] = (float(np.median(low_values)), float(np.median(high_values)))
        else:
            result[method] = (
                weighted_median(low_values, weights[sampled]),
                weighted_median(high_values, weights[high_indices]),
            )
    return result


def get_codon_statistics(
    genes: Iterable[object],
    highly_expressed: Sequence[bool] | NDArray[np.bool_],
    *,
    fragments: bool | None = False,
    exclude_short: bool = False,
    depth_of_coverage: ArrayLike | None = None,
    genetic_code: int | str = 11,
    n_le: int = 100,
    trimlen: int | None = None,
    trimside: str = "start",
    all_metrics: bool = False,
    background: str = "all",
    random_state: int | np.random.Generator | None = None,
) -> dict[str, float | int]:
    """Calculate gRodon's genome or community codon statistics.

    ``fragments=False`` requires >240 bp, ``True`` requires >120 bp, and
    ``None`` disables the minimum-length filter. If ``exclude_short`` and
    ``trimlen`` are supplied with ``fragments=None``, genes shorter than trimlen
    are removed. Trimming is in-frame and keeps metadata aligned.

    ``background='individual'`` compares each gene with 100 nucleotide shuffles
    and samples ``n_le`` low-expression genes with replacement. A fixed
    ``random_state`` reproduces the sampling in Python (not R's random stream).
    Coverage weights affect aggregation and GC, never reference codon counts.
    Zero-coverage low-expression genes are excluded from individual-background
    sampling, preventing an otherwise random all-zero sampled weight total.
    For compatibility with the fitted R models, individual-background metrics
    always use genetic code 1, as upstream ``singleMILC`` does. The requested
    ``genetic_code`` still controls the highly expressed consistency statistic.

    Returns the R field names, including ``CUBHE``, ``CUB``, ``ConsistencyHE``,
    ``GC``, ``GCdiv``, ``CPB``, ``FilteredSequences`` and ``nHE``. With
    ``all_metrics=True``, returns MILC, ENCprime, B, SCUO, MCB and nHE instead.
    """
    if fragments is not None and not isinstance(fragments, (bool, np.bool_)):
        raise TypeError("fragments must be True, False or None")
    if background not in ("all", "individual"):
        raise ValueError("background must be 'all' or 'individual'")
    if trimside not in ("start", "end", "random"):
        raise ValueError("trimside must be 'start', 'end' or 'random'")
    if trimlen is not None and (
        isinstance(trimlen, bool)
        or not isinstance(trimlen, (int, np.integer))
        or trimlen <= 0
        or trimlen % 3
    ):
        raise ValueError("trimlen must be a positive integer divisible by three")
    if isinstance(n_le, bool) or not isinstance(n_le, (int, np.integer)) or n_le <= 0:
        raise ValueError("n_le must be a positive integer")
    _genetic_code(genetic_code)
    rng = np.random.default_rng(random_state)
    threshold = (120 if bool(fragments) else 240) if fragments is not None else 1
    if fragments is None and exclude_short and trimlen is not None:
        threshold = trimlen - 1
    filtered = filter_sequences(
        genes, highly_expressed, length_threshold=threshold, depth_of_coverage=depth_of_coverage
    )
    sequences, he, depth = filtered.genes, filtered.highly_expressed, filtered.depth_of_coverage
    if trimlen is not None:
        # R puts trimmable genes before shorter retained genes. Preserve that
        # order because codon-pair bias includes inter-gene junction pairs.
        order = np.array(
            [i for i, gene in enumerate(sequences) if len(gene) >= trimlen]
            + [i for i, gene in enumerate(sequences) if len(gene) < trimlen],
            dtype=int,
        )
        sequences = tuple(sequences[i] for i in order)
        he = he[order]
        depth = depth[order] if depth is not None else None
        sequences = _trim_sequences(sequences, trimlen, trimside, rng)
    counts = codon_counts(sequences)
    metrics: tuple[Metric, ...] = (
        ("MILC", "ENCprime", "B", "SCUO", "MCB") if all_metrics else ("MILC",)
    )
    aggregates = {}
    if background == "individual":
        aggregates = _individual_statistics(sequences, he, depth, n_le, metrics, rng)
    else:
        for method in metrics:
            values = codon_usage(counts, method=method, genetic_code=genetic_code)
            if depth is None:
                aggregates[method] = (float(np.median(values)), float(np.median(values[he])))
            else:
                aggregates[method] = (
                    weighted_median(values, depth),
                    weighted_median(values[he], depth[he]),
                )
    if all_metrics:
        return {**{method: aggregates[method][1] for method in metrics}, "nHE": int(he.sum())}
    consistency = codon_usage(counts[he], reference_counts=counts[he], genetic_code=genetic_code)
    gc_counts = np.array([gene.count("G") + gene.count("C") for gene in sequences], dtype=float)
    lengths = np.array([len(gene) for gene in sequences], dtype=float)
    if depth is None:
        gc = float(np.sum(gc_counts) / np.sum(lengths))
        consistency_he = float(np.mean(consistency))
    else:
        # Rescaling weights leaves the estimates unchanged and prevents overflow
        # when a caller supplies large finite coverage units.
        gc_weights = depth / np.max(depth)
        he_weights = depth[he] / np.max(depth[he])
        gc = float(np.dot(gc_counts, gc_weights) / np.dot(lengths, gc_weights))
        consistency_he = float(np.average(consistency, weights=he_weights))
    return {
        "CUBHE": aggregates["MILC"][1],
        "GC": gc,
        "GCdiv": abs(0.5 - gc),
        "ConsistencyHE": consistency_he,
        "CUB": aggregates["MILC"][0],
        "CPB": codon_pair_bias(sequences, genetic_code=genetic_code)
        if depth is None and background == "all"
        else float("nan"),
        "FilteredSequences": filtered.filtered,
        "nHE": int(he.sum()),
    }


def get_codon_statistics_i(
    genes: Iterable[object], highly_expressed: Sequence[bool], **kwargs: Any
) -> dict[str, float | int]:
    """Convenience alias for statistics with ``background='individual'``."""
    kwargs["background"] = "individual"
    return get_codon_statistics(genes, highly_expressed, **kwargs)
