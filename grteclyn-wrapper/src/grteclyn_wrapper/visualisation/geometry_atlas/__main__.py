"""CLI: visualise top geometry-atlas elites as midplane field panels.

Example::

    python -m grteclyn_wrapper.visualisation.geometry_atlas \\
        runs/geometry_atlas/geometry_atlas_topologies_n64_L32_e260_20260723 \\
        --top 5
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .plot_elites import default_output_dir, visualise_top_elites


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Visualise top geometry-atlas MAP-Elites elites: midplane slices of "
            "α, |β|, shift quiver, |γ−I|, required ρ, and an intrinsic "
            "spatial-geometry embedding-style surface."
        )
    )
    parser.add_argument(
        "run_dir",
        type=Path,
        help="Campaign directory with archive.json + elites/ (or evals/).",
    )
    parser.add_argument("--top", type=int, default=5, help="Number of top elites (default 5).")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: <run_dir>/figures/geometry_atlas_elites).",
    )
    parser.add_argument("--n", type=int, default=None, help="Render resolution override.")
    parser.add_argument("--L", type=float, default=None, help="Box length override.")
    parser.add_argument(
        "--quiver-stride",
        type=int,
        default=4,
        help="Subsample stride for shift quiver (default 4).",
    )
    parser.add_argument(
        "--no-gallery",
        action="store_true",
        help="Skip the compact multi-elite gallery figure.",
    )
    args = parser.parse_args()

    run_dir = args.run_dir.expanduser().resolve()
    out_dir = (
        args.out_dir.expanduser().resolve()
        if args.out_dir is not None
        else default_output_dir(run_dir)
    )
    paths = visualise_top_elites(
        run_dir,
        top_n=int(args.top),
        out_dir=out_dir,
        n=args.n,
        L=args.L,
        quiver_stride=int(args.quiver_stride),
        gallery=not bool(args.no_gallery),
    )
    for p in paths:
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
