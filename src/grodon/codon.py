"""Codon counts and compositional statistics used by gRodon.

The definitions follow coRdon's codonUsage.R and functions.R, including its
alternative-start convention: alternative start codons belong to the methionine
family at *every* position. This convention matters for the fitted R models.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
from functools import lru_cache
from itertools import product
from typing import Any, Literal

import numpy as np
from Bio.Data import CodonTable
from numpy.typing import ArrayLike, NDArray

CODONS: tuple[str, ...] = tuple("".join(x) for x in product("ACGT", repeat=3))
_BASE_INDEX = np.full(256, -1, dtype=np.int16)
_BASE_INDEX[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)
Metric = Literal["MILC", "ENCprime", "B", "SCUO", "MCB"]
_METRICS: tuple[Metric, ...] = ("MILC", "ENCprime", "B", "SCUO", "MCB")


def normalize_genes(genes: Iterable[object]) -> list[str]:
    """Normalize sequences, name-to-sequence mappings or SeqRecords to DNA."""
    if isinstance(genes, (str, bytes)):
        raise TypeError("genes must be an iterable of sequences, not a single string")
    if isinstance(genes, Mapping):
        genes = genes.values()
    return [str(getattr(gene, "seq", gene)).upper() for gene in genes]


def _codon_indices(gene: str) -> NDArray[np.int64]:
    try:
        bases = _BASE_INDEX[np.frombuffer(gene.encode("ascii"), dtype=np.uint8)]
    except UnicodeEncodeError as exc:
        raise ValueError("DNA sequences may only contain A, C, G and T") from exc
    if len(bases) % 3:
        raise ValueError("DNA sequence lengths must be divisible by three")
    if np.any(bases < 0):
        raise ValueError("DNA sequences may only contain A, C, G and T")
    triplets = bases.reshape(-1, 3)
    return (triplets @ np.array([16, 4, 1])).astype(np.int64)


def codon_counts(genes: Iterable[object]) -> NDArray[np.int64]:
    """Count in-frame triplets; columns follow the alphabetical ``CODONS`` order.

    Unlike the high-level statistics API, this function does not discard invalid
    sequences. It raises an informative error instead.
    """
    sequences = normalize_genes(genes)
    counts = np.zeros((len(sequences), 64), dtype=np.int64)
    for row, gene in enumerate(sequences):
        counts[row] = np.bincount(_codon_indices(gene), minlength=64)
    return counts


@lru_cache(maxsize=64, typed=True)
def _genetic_code(genetic_code: int | str) -> CodonTable.CodonTable:
    if isinstance(genetic_code, (bool, np.bool_)) or not isinstance(
        genetic_code, (int, np.integer, str)
    ):
        raise ValueError("genetic_code must be an integer NCBI code or a genetic-code name")
    try:
        code = int(genetic_code)
    except (TypeError, ValueError):
        try:
            return CodonTable.unambiguous_dna_by_name[str(genetic_code)]
        except KeyError as exc:
            raise ValueError(f"Unknown genetic code: {genetic_code!r}") from exc
    try:
        return CodonTable.unambiguous_dna_by_id[code]
    except KeyError as exc:
        raise ValueError(f"Unknown genetic code: {genetic_code!r}") from exc


@lru_cache(maxsize=128)
def _families(
    genetic_code: int | str, alternative_starts: bool, remove_stops: bool
) -> tuple[NDArray[np.int64], tuple[NDArray[np.int64], ...], NDArray[np.int64]]:
    code = _genetic_code(genetic_code)
    amino_acids = [code.forward_table.get(codon, "*") for codon in CODONS]
    if alternative_starts:
        amino_acids = [
            "M" if codon in code.start_codons else aa
            for codon, aa in zip(CODONS, amino_acids, strict=True)
        ]
    kept = np.array(
        [i for i, aa in enumerate(amino_acids) if not remove_stops or aa != "*"],
        dtype=np.int64,
    )
    labels = np.array(amino_acids)[kept]
    families = tuple(np.flatnonzero(labels == aa) for aa in sorted(set(labels)))
    degeneracies = np.array([len(group) for group in families], dtype=np.int64)
    return kept, families, degeneracies


def _validate_counts(counts: ArrayLike) -> NDArray[np.float64]:
    result = np.asarray(counts, dtype=float)
    if result.ndim == 1:
        result = result[np.newaxis, :]
    if result.ndim != 2 or result.shape[1] != 64 or result.shape[0] == 0:
        raise ValueError("codon counts must have shape (number of genes, 64)")
    if not np.all(np.isfinite(result)) or np.any(result < 0):
        raise ValueError("codon counts must be finite and nonnegative")
    if np.any(result.sum(axis=1) == 0):
        raise ValueError("each gene must contain at least one codon")
    return result


def _sum_families(
    values: NDArray[np.float64], families: tuple[NDArray[np.int64], ...]
) -> NDArray[np.float64]:
    return np.column_stack([values[:, group].sum(axis=1) for group in families])


def _nansum(values: NDArray[np.float64], axis: int) -> NDArray[np.float64]:
    # R's na.rm=TRUE drops NaN, but preserves positive/negative infinity.
    return np.where(np.isnan(values), 0.0, values).sum(axis=axis)


def _effective_codons(
    homozygosities: NDArray[np.float64], degeneracies: NDArray[np.int64]
) -> NDArray[np.float64]:
    """Match coRdon's redundancy-class averaging and 61-codon upper cap."""
    classes = np.unique(degeneracies[degeneracies != 1])
    terms: dict[int, NDArray[np.float64]] = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        for redundancy in classes:
            values = homozygosities[:, degeneracies == redundancy]
            present = ~np.isnan(values) & (values != 0)
            average = _nansum(values, axis=1) / present.sum(axis=1)
            if redundancy != 3:
                average[np.isnan(average)] = 1 / redundancy
            terms[int(redundancy)] = values.shape[1] / average
        if 3 in terms:
            absent = np.isnan(terms[3])
            if np.any(absent):
                if 2 in terms and 4 in terms:
                    # This reproduces the upstream class-3 fallback, including
                    # its use of reciprocal homozygosities rather than F values.
                    if np.all(terms[2][absent] > 0) and np.all(terms[4][absent] > 0):
                        terms[3][absent] = (
                            terms[2][absent] / np.sum(degeneracies == 2)
                            + terms[4][absent] / np.sum(degeneracies == 4)
                        ) / 2
                    else:
                        terms[3][:] = 1 / np.sum(degeneracies == 3)
                else:
                    # coRdon errors on unusual codes missing classes 2 or 4.
                    # An absent class uses the same unbiased default as others.
                    terms[3][absent] = 3 * np.sum(degeneracies == 3)
    total = np.full(len(homozygosities), np.sum(degeneracies == 1), dtype=float)
    for term in terms.values():
        total += term
    return np.minimum(total, 61.0)


def codon_usage(
    counts: ArrayLike,
    method: Metric = "MILC",
    *,
    reference_counts: ArrayLike | None = None,
    genetic_code: int | str = 11,
    alternative_starts: bool = True,
    remove_stops: bool | None = None,
) -> NDArray[np.float64]:
    """Compute one codon-usage statistic per gene.

    ``reference_counts`` is a 64-column count matrix or an aggregate 64-vector.
    Its pooled synonymous-codon frequencies form the background; omitted means
    all input genes. ENCprime removes stops by default; other metrics keep them.
    SCUO does not depend on a reference. Undefined values are represented by NaN.
    """
    if method not in _METRICS:
        raise ValueError(f"method must be one of {_METRICS}")
    _genetic_code(genetic_code)
    full_counts = _validate_counts(counts)
    if remove_stops is None:
        remove_stops = method == "ENCprime"
    kept, families, deg = _families(genetic_code, alternative_starts, remove_stops)
    observed = full_counts[:, kept]
    totals = observed.sum(axis=1)
    aa_counts = _sum_families(observed, families)
    frequencies = np.zeros_like(observed)
    with np.errstate(divide="ignore", invalid="ignore"):
        for family, indices in enumerate(families):
            frequencies[:, indices] = observed[:, indices] / aa_counts[:, [family]]
        if method == "SCUO":
            entropy = _sum_families(
                np.where(frequencies > 0, -frequencies * np.log10(frequencies), 0),
                families,
            )
            maximum = np.log10(deg)
            orderliness = (maximum - entropy) / maximum
            abundance = aa_counts / aa_counts[:, deg > 1].sum(axis=1, keepdims=True)
            return _nansum(orderliness * abundance, axis=1)

        reference = full_counts if reference_counts is None else _validate_counts(reference_counts)
        pooled = reference[:, kept].sum(axis=0)
        expected = np.zeros_like(pooled)
        for indices in families:
            expected[indices] = pooled[indices] / pooled[indices].sum()

        if method == "MILC":
            terms = 2 * observed * np.log(frequencies / expected)
            correction = (aa_counts > 0) @ (deg - 1) / totals - 0.5
            return _nansum(terms, axis=1) / totals - correction
        if method == "B":
            distances = _sum_families(np.abs(frequencies - expected), families)
            return _nansum(aa_counts / totals[:, None] * distances, axis=1)

        chi = (frequencies - expected) ** 2 / expected
        chi[np.isnan(chi)] = 0
        if method == "MCB":
            chi[observed <= 0] = 0
            distances = _sum_families(chi, families)
            terms = distances[:, deg > 1] * np.log10(aa_counts[:, deg > 1])
            return _nansum(terms, axis=1) / (aa_counts[:, deg > 1] > 0).sum(axis=1)
        chi_square = aa_counts * _sum_families(chi, families)
        homozygosity = (chi_square + aa_counts - deg) / (deg * (aa_counts - 1))
        homozygosity[aa_counts < 5] = np.nan
        return _effective_codons(homozygosity, deg)


def milc(counts: ArrayLike, **kwargs: Any) -> NDArray[np.float64]:
    """Measure Independent of Length and Composition (MILC)."""
    return codon_usage(counts, "MILC", **kwargs)


def enc_prime(counts: ArrayLike, **kwargs: Any) -> NDArray[np.float64]:
    """Effective number of codons corrected for background composition."""
    return codon_usage(counts, "ENCprime", **kwargs)


def karlin_b(counts: ArrayLike, **kwargs: Any) -> NDArray[np.float64]:
    """Karlin's B statistic."""
    return codon_usage(counts, "B", **kwargs)


def scuo(counts: ArrayLike, **kwargs: Any) -> NDArray[np.float64]:
    """Synonymous codon usage orderliness."""
    return codon_usage(counts, "SCUO", **kwargs)


def mcb(counts: ArrayLike, **kwargs: Any) -> NDArray[np.float64]:
    """Maximum-likelihood codon bias."""
    return codon_usage(counts, "MCB", **kwargs)


def codon_pair_bias(genes: Iterable[object], *, genetic_code: int | str = 11) -> float:
    """Compute gRodon's codon-pair bias, preserving the R model convention.

    Genes are concatenated, so junction pairs are included. Scores are summed
    over distinct pairs, then divided by the number of distinct pairs minus one;
    they are not weighted by pair frequency. Degenerate inputs return NaN.
    """
    sequences = normalize_genes(genes)
    code = _genetic_code(genetic_code)
    indices = [_codon_indices(gene) for gene in sequences]
    if not indices:
        return float("nan")
    genome = np.concatenate(indices)
    if len(genome) < 2:
        return float("nan")
    codon_n = np.bincount(genome, minlength=64)
    pair_keys, pair_n = np.unique(genome[:-1] * 64 + genome[1:], return_counts=True)
    if len(pair_keys) <= 1:
        return float("nan")
    # Translating each isolated codon in R treats alternative start codons as
    # methionine. Translating the whole concatenated genome only does so at its
    # first codon. Retain this asymmetry for exact compatibility.
    ordinary_aa = np.array([code.forward_table.get(c, "*") for c in CODONS])
    pair_aa = np.array(
        ["M" if c in code.start_codons else aa for c, aa in zip(CODONS, ordinary_aa, strict=True)]
    )
    translated = ordinary_aa[genome].tolist()
    if CODONS[int(genome[0])] in code.start_codons:
        translated[0] = "M"
    aa_n = Counter(translated)
    first, second = pair_keys // 64, pair_keys % 64
    aa_pairs = list(zip(pair_aa[first], pair_aa[second], strict=True))
    aa_pair_n: Counter[tuple[str, str]] = Counter()
    for pair, count in zip(aa_pairs, pair_n, strict=True):
        aa_pair_n[pair] += int(count)
    scores = np.empty(len(pair_keys), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        for i, (a1, a2) in enumerate(aa_pairs):
            if not aa_n[a1] or not aa_n[a2]:
                return float("nan")
            expected = aa_pair_n[(a1, a2)] * (
                float(codon_n[first[i]]) * codon_n[second[i]] / (aa_n[a1] * aa_n[a2])
            )
            scores[i] = np.log(pair_n[i] / expected)
    return float(scores.sum() / (len(scores) - 1))
