"""Command registration and dispatch for reproducible bacterial-genome analysis."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from .references import DEFAULT_REFERENCE, reference_catalog
from .tools import DEFAULT_BLAST_VERSION

COMMANDS = frozenset({"genome", "genome-batch", "doctor", "setup", "references"})
_CONFIG_ARGUMENTS = (
    "reference",
    "blastp",
    "mode",
    "training_set",
    "temperature",
    "genetic_code",
    "gene_mode",
    "closed",
    "mask",
    "min_identity",
    "min_query_coverage",
    "min_reference_coverage",
    "max_evalue",
    "min_ribosomal",
    "confidence",
    "threads",
    "timeout",
)


def _genome_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="New output directory; existing files are preserved",
    )
    parser.add_argument(
        "--reference",
        default=DEFAULT_REFERENCE,
        help=f"Bundled profile or ribosomal protein FASTA (default: {DEFAULT_REFERENCE})",
    )
    parser.add_argument("--blastp", default="blastp", help="BLASTP executable name or full path")
    parser.add_argument(
        "--rscript",
        type=Path,
        help="Optionally compare to the original R package using this Rscript executable",
    )
    parser.add_argument("--mode", choices=("full", "partial"), default="full")
    parser.add_argument("--training-set", choices=("madin", "vs", "AOA_NOB"), default="madin")
    parser.add_argument("--temperature", type=float, help="Optimal growth temperature in Celsius")
    parser.add_argument(
        "--genetic-code", type=int, default=11, help="Translation table (default: 11)"
    )
    parser.add_argument(
        "--gene-mode",
        choices=("single", "meta"),
        default="single",
        help="Single-genome training or pretrained metagenomic gene calling (default: single)",
    )
    parser.add_argument(
        "--closed",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Prevent genes crossing contig ends (default: enabled)",
    )
    parser.add_argument(
        "--mask",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Mask ambiguous sequence gaps during gene calling (default: enabled)",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=70.0,
        help="Minimum amino-acid identity, percent (0–100)",
    )
    parser.add_argument(
        "--min-query-coverage",
        type=float,
        default=0.8,
        help="Minimum aligned fraction of each query protein (0–1; default: 0.8)",
    )
    parser.add_argument(
        "--min-reference-coverage",
        type=float,
        default=0.8,
        help="Minimum aligned fraction of each reference protein (0–1; default: 0.8)",
    )
    parser.add_argument("--max-evalue", type=float, default=1e-8)
    parser.add_argument(
        "--min-ribosomal",
        type=int,
        default=10,
        help="Minimum retained ribosomal genes required for prediction (default: 10)",
    )
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--threads", type=int, default=1, help="BLASTP CPU threads (default: 1)")
    parser.add_argument(
        "--timeout",
        type=float,
        default=600.0,
        help="External command timeout, seconds (default: 600)",
    )
    parser.set_defaults(output=None)


def register_commands(subparsers: argparse._SubParsersAction) -> None:
    """Add genome preparation, batch analysis, diagnostics, and setup commands."""
    genome = subparsers.add_parser(
        "genome",
        help="Predict growth directly from a bacterial genome FASTA",
        description=(
            "Call genes, identify ribosomal proteins, and predict growth from a bacterial "
            "genome. Uses local tools and pinned reference snapshots; no LLM or network "
            "access is used during analysis. From the package checkout, install optional "
            "dependencies with pip install '.[genome]' and grodon setup."
        ),
    )
    genome.add_argument("genome", type=Path, help="Genome FASTA or gzip-compressed FASTA")
    _genome_arguments(genome)

    batch = subparsers.add_parser(
        "genome-batch",
        help="Analyze several bacterial genomes into separate output directories",
        description="Run the same reproducible genome pipeline for each input assembly.",
    )
    batch.add_argument(
        "genomes", nargs="+", type=Path, help="Genome FASTA files, optionally gzipped"
    )
    _genome_arguments(batch)
    batch.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue with other genomes after an input fails; retain failures in the summary",
    )

    doctor = subparsers.add_parser(
        "doctor", help="Check optional gene-calling and BLASTP dependencies without downloading"
    )
    doctor.add_argument("--blastp", default="blastp", help="BLASTP executable name or full path")
    doctor.set_defaults(output=None)

    setup = subparsers.add_parser(
        "setup",
        help="Explicitly download and verify BLASTP into a dedicated tools directory",
        description=(
            "Download BLASTP from the official NCBI release, verify its checksum, and "
            "install into the supplied directory. This command does not modify PATH."
        ),
    )
    setup.add_argument(
        "--tools-dir", type=Path, required=True, help="Dedicated local tools directory"
    )
    setup.add_argument(
        "--version", default=DEFAULT_BLAST_VERSION, help="Pinned NCBI BLAST+ release"
    )
    setup.set_defaults(output=None)

    references = subparsers.add_parser(
        "references", help="List bundled offline ribosomal profiles and their provenance"
    )
    references.set_defaults(output=None)


def run_command(args: argparse.Namespace) -> dict[str, Any]:
    """Execute a registered command and return its machine-readable summary.

    Exceptions propagate to the root CLI's consistent error handler. Pipeline
    imports are delayed so diagnostics/help/reference listings work without the
    optional gene-calling dependency.
    """
    if args.command == "references":
        return {
            "default": DEFAULT_REFERENCE,
            "references": reference_catalog(),
            "custom_reference": "Pass a curated ribosomal protein FASTA using --reference PATH",
        }
    if args.command == "doctor":
        from .tools import doctor

        return doctor(blastp=args.blastp)
    if args.command == "setup":
        from .tools import install_blast

        return install_blast(args.tools_dir, version=args.version)
    if args.command not in {"genome", "genome-batch"}:
        raise ValueError(f"Unknown genome command: {args.command!r}")

    from .pipeline import GenomePipelineConfig, run_genome, run_genomes

    config = GenomePipelineConfig(**{field: getattr(args, field) for field in _CONFIG_ARGUMENTS})
    if args.command == "genome":
        return run_genome(args.genome, args.outdir, config=config, rscript=args.rscript).to_dict()
    return run_genomes(
        args.genomes,
        args.outdir,
        config=config,
        rscript=args.rscript,
        continue_on_error=args.continue_on_error,
    )
