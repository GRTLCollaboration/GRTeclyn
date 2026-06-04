#!/usr/bin/env python3
"""CLI launcher: write a Teo rotating wormhole .gridinit for GRTeclyn.

All physics lives in ``grteclyn_wrapper.initial_data.teo``; this script only
parses arguments, calls the generator, and writes the file plus a manifest.

Run variants are driven here (spin / source), not by separate params files:

    # geometry-first weak-spin control (effective source in the file)
    make_teo_wormhole_gridinit.py --spin 0.05 --output runs/teo/weak_spin.gridinit
    # matched a=0 baseline for junk-radiation subtraction
    make_teo_wormhole_gridinit.py --spin 0.0  --output runs/teo/a0.gridinit
    # zeroed source for vacuum relaxation control
    make_teo_wormhole_gridinit.py --spin 0.05 --source zero \
        --output runs/teo/weak_spin_zero_source.gridinit
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

WRAPPER_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WRAPPER_ROOT / "src"))

from grteclyn_wrapper.grtresna.io import write_gridinit
from grteclyn_wrapper.initial_data.teo import (
    TeoWormholeConfig,
    build_grid,
    constraint_residuals,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=Path("runs/teo_wormhole/initial_data.gridinit"))
    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--ny", type=int, default=64)
    parser.add_argument("--nz", type=int, default=64)
    parser.add_argument("--lx", type=float, default=64.0)
    parser.add_argument("--ly", type=float, default=64.0)
    parser.add_argument("--lz", type=float, default=64.0)
    parser.add_argument("--center", type=float, nargs=3, default=(32.0, 32.0, 0.0), metavar=("X", "Y", "Z"))
    parser.add_argument("--origin", type=float, nargs=3, default=(0.0, 0.0, 0.0), metavar=("X", "Y", "Z"))
    parser.add_argument("--throat-radius", type=float, default=0.5)
    parser.add_argument(
        "--spin",
        type=float,
        default=0.05,
        help="Dimensionless Teo spin a_hat; ergoregion threshold is near |a_hat|=0.5.",
    )
    parser.add_argument(
        "--source",
        choices=("einstein", "zero"),
        default="einstein",
        help="Populate teo_* fields from the numerical Einstein tensor or zero them.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Also evaluate ADM Hamiltonian/momentum constraint residuals (self-consistency).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cfg = TeoWormholeConfig(
        nx=args.nx,
        ny=args.ny,
        nz=args.nz,
        Lx=args.lx,
        Ly=args.ly,
        Lz=args.lz,
        center=tuple(args.center),
        origin=tuple(args.origin),
        throat_radius=args.throat_radius,
        spin=args.spin,
        source=args.source,
    )

    grid = build_grid(cfg)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_gridinit(grid.data, grid.comp_names, grid.dx_xyz, grid.origin, args.output)

    manifest = {
        "gridinit": str(args.output),
        "config": asdict(cfg),
        "component_names": grid.comp_names,
        "dx": grid.dx_xyz.tolist(),
        "origin": grid.origin.tolist(),
        "metrics": grid.metrics,
    }
    if args.check:
        manifest["constraint_residuals"] = constraint_residuals(cfg).as_dict()

    manifest_path = args.output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    summary = {"gridinit": str(args.output), "manifest": str(manifest_path), **grid.metrics}
    if args.check:
        summary["constraint_residuals"] = manifest["constraint_residuals"]
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
