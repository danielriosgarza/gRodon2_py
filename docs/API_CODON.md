# Low-level codon API

[API reference](API.md) · [Growth models and fitting](API_FITTING.md) · [Copiotrophy](API_COPIOTROPHY.md)

These functions are imported from `grodon.codon`. They expose numerical building
blocks; for filtering, highly expressed gene selection, and growth prediction,
use the [high-level API](API.md). Their results follow the coRdon conventions
used by the bundled R models. They do not call R or require genome-pipeline tools.

## Constants and input conventions

`CODONS` is a tuple of 64 strings in the order produced by
`itertools.product("ACGT", repeat=3)`: `AAA`, `AAC`, `AAG`, `AAT`, …, `TTT`.
Every count vector and matrix below uses this column order, including stop
codons. `Metric` is the typing alias
`Literal["MILC", "ENCprime", "B", "SCUO", "MCB"]`.

Sequence inputs are iterables of DNA strings, Biopython sequence objects or
`SeqRecord` objects, or mappings whose values are sequences. Mapping insertion
order is retained. Pass a list containing one sequence for a single-gene input;
a bare string or bytes object raises `TypeError`.

## `normalize_genes`

```text
normalize_genes(genes: Iterable[object]) -> list[str]
```

Converts each item, or its `.seq` attribute when present, to uppercase text.
For a mapping, uses the values in insertion order. This helper does not validate
DNA characters or frame, remove whitespace, translate proteins, or replace RNA
`U` with DNA `T`. Subsequent counting functions perform DNA validation.

## `codon_counts`

```text
codon_counts(genes: Iterable[object]) -> numpy.ndarray
```

Returns an `int64` array of shape `(number_of_genes, 64)`. It counts every
non-overlapping triplet starting at position zero; start and stop codons are
retained. Input case is normalized. No genes are silently removed.

Raises `ValueError` for non-ACGT characters or a sequence length not divisible
by three. An empty iterable produces an array of shape `(0, 64)`; an empty
sequence produces a zero row. Such empty inputs are not accepted by
`codon_usage`. Unlike the high-level growth API, there is no minimum gene length.

```python
from grodon.codon import CODONS, codon_counts

counts = codon_counts(["ATGGCTGCTTAA", "atggccgcctag"])
assert counts.shape == (2, 64)
assert counts[0, CODONS.index("GCT")] == 2
assert counts[1, CODONS.index("GCC")] == 2
assert counts.sum(axis=1).tolist() == [4, 4]
```

These short artificial sequences illustrate counting only; they would be
excluded by the default growth-prediction length filter.

## `codon_usage`

```text
codon_usage(
    counts,
    method="MILC",
    *,
    reference_counts=None,
    genetic_code=11,
    alternative_starts=True,
    remove_stops=None,
) -> numpy.ndarray
```

Returns one `float64` value per input row, shape `(n_genes,)`. A single
64-element vector is promoted to a one-row matrix, so its return shape is `(1,)`.

| Parameter | Type and default | Meaning |
| --- | --- | --- |
| `counts` | Array-like, required | A `(n, 64)` count matrix or one 64-element vector, in `CODONS` order. Counts may be floating-point, but must be finite and nonnegative; every row must have positive total count. |
| `method` | `Metric`, `"MILC"` | Exact, case-sensitive metric name listed below. |
| `reference_counts` | Array-like or `None` | A count matrix or aggregate 64-vector. Rows are pooled to define synonymous-codon frequencies. `None` pools `counts`. Unused for SCUO. |
| `genetic_code` | `int` or `str`, `11` | An NCBI translation-table number, its numeric string, or a DNA code name recognized by Biopython. Booleans and unknown codes are rejected. |
| `alternative_starts` | `bool`, `True` | Treat every codon listed as an alternative start in the selected code as methionine at every position, following coRdon's model convention. This is not restricted to the first triplet of a gene. |
| `remove_stops` | `bool` or `None` | With `None`, excludes codons translated as stop for ENCprime and retains them for the other metrics. An explicit boolean overrides that choice. Dual-coding codons with a sense assignment in the selected table are retained. |

| Method | Quantity | Wrapper |
| --- | --- | --- |
| `MILC` | Measure Independent of Length and Composition | `milc` |
| `ENCprime` | Effective number of codons corrected for background composition | `enc_prime` |
| `B` | Karlin's codon-bias statistic | `karlin_b` |
| `SCUO` | Synonymous codon usage orderliness; independent of reference counts | `scuo` |
| `MCB` | Maximum-likelihood codon bias | `mcb` |

Raises `ValueError` for an unknown method or genetic code, an invalid count
shape, an empty matrix, negative/nonfinite counts, or a zero-total row. The same
count validation applies to reference counts when the chosen metric uses them.
Undefined statistics can be `NaN`; zero reference frequencies can also produce
infinity. These outputs are not replaced with zero or a sentinel number.

## Metric wrappers

```text
milc(counts, **kwargs) -> numpy.ndarray
enc_prime(counts, **kwargs) -> numpy.ndarray
karlin_b(counts, **kwargs) -> numpy.ndarray
scuo(counts, **kwargs) -> numpy.ndarray
mcb(counts, **kwargs) -> numpy.ndarray
```

Each wrapper selects the corresponding `codon_usage` method. `kwargs` accepts
`reference_counts`, `genetic_code`, `alternative_starts`, and `remove_stops`.
Return shapes, validation, and undefined-value behavior are the same as above.
Do not pass a second `method` value through a wrapper.

```python
import numpy as np
from grodon.codon import codon_counts, codon_usage, milc

# Artificial, in-frame sequences for demonstrating the numerical interface.
counts = codon_counts([
    "ATGGCTGCCGCGGCA" * 25,
    "ATGGCCGCCGCTGCA" * 25,
    "ATGGCGGCGGCTGCA" * 25,
])
reference = counts.sum(axis=0)
scores = milc(counts, reference_counts=reference, genetic_code=11)
np.testing.assert_allclose(
    scores,
    codon_usage(counts, method="MILC", reference_counts=reference),
)
assert scores.shape == (3,)
```

## `codon_pair_bias`

```text
codon_pair_bias(genes: Iterable[object], *, genetic_code=11) -> float
```

Computes one aggregate codon-pair bias value. `genes` uses the sequence-input
conventions above. `genetic_code` accepts the same values as `codon_usage`.
Invalid sequence characters, frames, or genetic codes raise `ValueError`.

To preserve the original R model calculation, the function:

1. Concatenates genes in their supplied order, including pairs crossing gene
   boundaries. Changing gene order can therefore change the result.
2. Sums scores over distinct observed codon pairs and divides by the number of
   distinct pairs minus one. It does not weight that sum by pair frequency.
3. Retains the upstream alternative-start convention: isolated codons used for
   pair labels treat alternative starts as methionine, whereas translation of
   the concatenated sequence applies that rule only at its first codon.

Returns `NaN` for degenerate inputs, including no genes, fewer than two total
codons, or only one distinct observed codon pair. Check `math.isfinite(value)`
before passing the result to a fitted regression. There is no length filtering,
highly expressed subset, or coverage weighting in this low-level function.

See [compatibility notes](COMPATIBILITY.md) for the R conventions and the scope
of numerical-equivalence testing.
