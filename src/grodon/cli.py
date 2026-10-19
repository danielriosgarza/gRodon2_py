"""The grodon command-line interface."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .io import DEFAULT_HE_PATTERN, highly_expressed_mask, read_coverage, read_fasta


def _json_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    if hasattr(value, "tolist"):
        return _json_value(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _write_json(value: Any, output: Path | None) -> None:
    text = json.dumps(_json_value(value), indent=2, allow_nan=False) + "\n"
    if output is None:
        sys.stdout.write(text)
    else:
        output.write_text(text, encoding="utf-8")


def _parser() -> argparse.ArgumentParser:
    from . import __version__

    parser = argparse.ArgumentParser(prog="grodon", description="gRodon microbial trait inference")
    parser.add_argument("--version", action="version", version=f"grodon {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    from .genome_cli import register_commands

    register_commands(commands)
    for command in ("predict", "statistics"):
        sub = commands.add_parser(command, help=f"Calculate {command} from coding-sequence FASTA")
        sub.add_argument("fasta", type=Path)
        selection = sub.add_mutually_exclusive_group()
        selection.add_argument(
            "--he-pattern",
            default=DEFAULT_HE_PATTERN,
            help="Case-insensitive regex matching highly expressed annotations",
        )
        selection.add_argument(
            "--he-ids", type=Path, help="Text file containing one FASTA ID per line"
        )
        sub.add_argument("--coverage", type=Path, help="TSV with id and coverage columns")
        sub.add_argument("--coverage-id-column", default="id")
        sub.add_argument("--coverage-column", default="coverage")
        sub.add_argument("--genetic-code", type=int)
        sub.add_argument(
            "--fragments",
            choices=("false", "true", "none"),
            default="false",
            help="Length filter: >240 bp, >120 bp, or no minimum (in-frame still required)",
        )
        sub.add_argument("--exclude-short", action="store_true")
        sub.add_argument("--background", choices=("all", "individual"), default="all")
        sub.add_argument("--n-le", type=int, default=100)
        sub.add_argument("--seed", type=int, default=None)
        sub.add_argument("-o", "--output", type=Path, help="Output JSON (default stdout)")
        if command == "predict":
            sub.add_argument(
                "--mode",
                default="full",
                choices=(
                    "full",
                    "partial",
                    "metagenome_v1",
                    "metagenome_v2",
                    "eukaryote",
                    "metagenome_euk",
                    "metagenome_150bp",
                    "metagenome_240bp",
                ),
            )
            sub.add_argument("--temperature", type=float)
            sub.add_argument("--training-set", choices=("madin", "vs", "AOA_NOB"), default="madin")
            sub.add_argument("--confidence", type=float, default=0.95)
        else:
            sub.add_argument("--all-metrics", action="store_true")
            sub.add_argument("--trimlen", type=int)
            sub.add_argument("--trimside", choices=("start", "end", "random"), default="start")
    for command in ("ioc", "fit-ioc"):
        sub = commands.add_parser(command, help="Predict or fit the Index of Copiotrophy")
        sub.add_argument("table", type=Path, help="CSV with dCUB, nCAZy, nGenes columns")
        sub.add_argument("--completeness", type=float, default=100.0)
        sub.add_argument("-o", "--output", type=Path)
        if command == "ioc":
            sub.add_argument(
                "--model",
                required=True,
                help="permafrost, pacific, soil, human, all_habitat, or model JSON",
            )
        else:
            sub.add_argument("--save-model", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run a command; diagnostics go to stderr and machine-readable output to stdout."""
    parser = _parser()
    args = parser.parse_args(argv)
    value: dict[str, Any]
    try:
        from .genome_cli import COMMANDS, run_command

        if args.command in COMMANDS:
            value = run_command(args)
        elif args.command in ("predict", "statistics"):
            genes = read_fasta(args.fasta)
            ids = None
            if args.he_ids:
                ids = [
                    line.strip()
                    for line in args.he_ids.read_text(encoding="utf-8-sig").splitlines()
                    if line.strip() and not line.strip().startswith("#")
                ]
            he = highly_expressed_mask(genes, pattern=args.he_pattern, ids=ids)
            depth = (
                None
                if args.coverage is None
                else read_coverage(
                    args.coverage,
                    genes,
                    id_column=args.coverage_id_column,
                    coverage_column=args.coverage_column,
                )
            )
            kwargs = dict(
                depth_of_coverage=depth,
                fragments={"false": False, "true": True, "none": None}[args.fragments],
                exclude_short=args.exclude_short,
                background=args.background,
                n_le=args.n_le,
                random_state=args.seed,
            )
            if args.command == "predict":
                from .prediction import predict_growth

                result = predict_growth(
                    genes,
                    he,
                    mode=args.mode,
                    temperature=args.temperature,
                    training_set=args.training_set,
                    genetic_code=args.genetic_code,
                    confidence=args.confidence,
                    **kwargs,
                )
                value = dict(result)
                value["metadata"] = {
                    "model": result.model_name,
                    "mode": result.mode,
                    "training_set": result.training_set,
                    "background": result.background,
                }
            else:
                from .statistics import get_codon_statistics

                value = get_codon_statistics(
                    genes,
                    he,
                    genetic_code=11 if args.genetic_code is None else args.genetic_code,
                    all_metrics=args.all_metrics,
                    trimlen=args.trimlen,
                    trimside=args.trimside,
                    **kwargs,
                )
        else:
            import pandas as pd

            from .copiotrophy import IoCModel, fit_ioc, predict_ioc

            table = pd.read_csv(args.table)
            if args.command == "fit-ioc":
                if args.output is not None and args.output.resolve() == args.save_model.resolve():
                    raise ValueError("--output and --save-model must name different files")
                fitted = fit_ioc(table, completeness=args.completeness)
                value = {"IoC": fitted.values}
                _write_json(fitted.model.to_dict(), args.save_model)
            else:
                model = args.model
                if Path(model).is_file():
                    model = IoCModel.from_dict(json.loads(Path(model).read_text(encoding="utf-8")))
                value = {"IoC": predict_ioc(table, model=model, completeness=args.completeness)}
        _write_json(value, args.output)
        if args.command == "genome-batch" and value.get("failed", 0):
            return 1
        if args.command == "doctor" and not value.get("ready", False):
            return 1
    except (ValueError, TypeError, OSError, KeyError, EOFError, RuntimeError, re.error) as exc:
        parser.exit(2, f"grodon: error: {exc}\n")
    return 0
