"""Run a complete workflow on generated, nonbiological test sequences.

Install the project first, then run ``python examples/synthetic_workflow.py``.
The numbers demonstrate the software interface, not biological conclusions.
"""

import json
import warnings

import numpy as np

from grodon import GrowthPredictionWarning, fit_ioc, predict_growth, predict_ioc
from grodon.codon import CODONS


def main() -> None:
    rng = np.random.default_rng(2026)
    codons = np.array([codon for codon in CODONS if codon not in {"TAA", "TAG", "TGA"}])
    genes = ["ATG" + "".join(rng.choice(codons, size=149)) for _ in range(160)]
    he = [index < 20 for index in range(len(genes))]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", GrowthPredictionWarning)
        result = predict_growth(genes, he)
    print("Synthetic example; estimates are not biological observations.")
    print(json.dumps(dict(result), indent=2))

    fitted = fit_ioc(
        dcub=[-0.8, -0.4, -0.6, -0.2, -0.9],
        n_cazy=[90, 40, 100, 30, 160],
        n_genes=[3500, 2200, 4200, 1900, 5000],
    )
    reapplied = predict_ioc(
        dcub=[-0.8, -0.4, -0.6, -0.2, -0.9],
        n_cazy=[90, 40, 100, 30, 160],
        n_genes=[3500, 2200, 4200, 1900, 5000],
        model=fitted.model,
    )
    np.testing.assert_allclose(reapplied, fitted.values)
    print("IoC:", reapplied.tolist())


if __name__ == "__main__":
    main()
