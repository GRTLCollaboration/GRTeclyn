#!/usr/bin/env python3
"""Entry point for metric guesser batch validation.

Run from grteclyn-wrapper directory::

    uv run python tests/validate_metric_guesser.py
    uv run python tests/validate_metric_guesser.py --output-dir ./validation_out
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from grteclyn_wrapper.initial_data.validate_guesser import run_validation


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate the constraint-satisfying metric guesser on synthetic candidates.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for candidate generation (default: 42).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_out"),
        help="Directory for guesser_validation.csv and summary JSON.",
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="Print summary only; do not write CSV/JSON files.",
    )
    args = parser.parse_args()

    output_dir = None if args.no_write else args.output_dir
    run_validation(seed=args.seed, output_dir=output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
