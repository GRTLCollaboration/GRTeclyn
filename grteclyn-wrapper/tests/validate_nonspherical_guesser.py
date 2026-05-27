#!/usr/bin/env python3
"""One-step non-spherical metric validation along 1D radial rays.

Run from grteclyn-wrapper:

    uv run python tests/validate_nonspherical_guesser.py --output-dir validation_out
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from grteclyn_wrapper.nonspherical_guesser import run_nonspherical_validation


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate axisymmetric non-spherical chi proposals along radial rays.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_out"),
        help="Directory for nonspherical ray-validation CSV/JSON.",
    )
    parser.add_argument(
        "--contrast-threshold",
        type=float,
        default=0.5,
        help="Reject if max angular chi contrast exceeds this value.",
    )
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    output_dir = None if args.no_write else args.output_dir
    run_nonspherical_validation(
        seed=args.seed,
        output_dir=output_dir,
        contrast_threshold=args.contrast_threshold,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
