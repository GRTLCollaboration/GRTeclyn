#!/usr/bin/env python3
"""Generate GRTresna initial data for rotating exotic wormhole candidates."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

WRAPPER_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(WRAPPER_ROOT / "src"))

from grteclyn_wrapper.grtresna.solver import (
    GRTresnaConfig,
    parse_convergence,
    solve,
    write_grtresna_params,
)


LOG = logging.getLogger("make_rotating_wormhole_id")


def _parse_float_list(raw: str) -> list[float]:
    values: list[float] = []
    for part in raw.replace(",", " ").split():
        values.append(float(part))
    if not values:
        raise argparse.ArgumentTypeError("expected at least one floating point value")
    return values


def _case_name(omega: float, amp: float, width: float) -> str:
    sign = "p" if omega >= 0.0 else "m"
    omega_token = f"{abs(omega):.4g}".replace(".", "p")
    amp_token = f"{amp:.4g}".replace(".", "p")
    width_token = f"{width:.4g}".replace(".", "p")
    return f"rotwh_omega_{sign}{omega_token}_amp_{amp_token}_w_{width_token}"


def build_config(args: argparse.Namespace, omega: float) -> GRTresnaConfig:
    """Build one exotic rotating GRTresna configuration."""
    return GRTresnaConfig(
        mpi_ranks=args.ranks,
        max_NL_iterations=args.iterations,
        timeout=args.timeout,
        N=(args.nx, args.ny, args.nz),
        L=args.length,
        max_level=args.max_level,
        block_factor=args.block_factor,
        max_grid_size=args.max_grid_size,
        refine_threshold=args.refine_threshold,
        regrid_radius=args.regrid_radius,
        coefficient_average_type="arithmetic",
        bh1_bare_mass=args.bh_mass,
        bh1_spin=(0.0, 0.0, 0.0),
        bh2_bare_mass=0.0,
        dphi=0.0,
        dpi=0.0,
        scalar_mass=args.scalar_mass,
        lumps=[
            {
                "amp": args.amp,
                "width": args.width,
                "center": (0.0, 0.0, 0.0),
                "velocity": (args.velocity_x, args.velocity_y, args.velocity_z),
                "omega": omega,
                "mode": args.mode,
                "exotic": 1,
            }
        ],
        maximal_slicing=True,
        psi_relaxation=args.psi_relaxation,
        psi_floor=args.psi_floor,
        maximal_jacobian_cap=args.maximal_jacobian_cap,
        gridinit_nx=args.gridinit_nx,
        gridinit_ny=args.gridinit_ny,
        gridinit_nz=args.gridinit_nz,
        target_center=(args.target_center_x, args.target_center_y, args.target_center_z),
        cleanup=not args.keep_hdf5,
    )


def _write_case_manifest(
    case_dir: Path,
    *,
    cfg: GRTresnaConfig,
    gridinit: Path | None,
    convergence: dict[str, float] | None,
    accepted: bool,
) -> None:
    payload = {
        "accepted": accepted,
        "gridinit": str(gridinit) if gridinit is not None else None,
        "convergence": convergence,
        "config": asdict(cfg),
    }
    (case_dir / "manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _within_thresholds(
    convergence: dict[str, float] | None, max_ham_pct: float, max_mom_pct: float
) -> bool:
    if convergence is None:
        return False
    return (
        convergence["ham_pct"] <= max_ham_pct
        and convergence["mom_pct"] <= max_mom_pct
    )


def generate_cases(args: argparse.Namespace, omegas: Iterable[float]) -> int:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    failures = 0

    for omega in omegas:
        cfg = build_config(args, omega)
        case_dir = args.out_dir / _case_name(omega, args.amp, args.width)
        case_dir.mkdir(parents=True, exist_ok=True)

        params_path = case_dir / "params.txt"
        write_grtresna_params(cfg, params_path)
        LOG.info("Wrote %s", params_path)

        if args.dry_run:
            _write_case_manifest(
                case_dir,
                cfg=cfg,
                gridinit=None,
                convergence=None,
                accepted=False,
            )
            continue

        gridinit_path = case_dir / "initial_data.gridinit"
        gridinit = solve(cfg, work_dir=case_dir, gridinit_path=gridinit_path)
        convergence = parse_convergence(case_dir)
        accepted = _within_thresholds(
            convergence, args.max_ham_pct, args.max_mom_pct
        )
        _write_case_manifest(
            case_dir,
            cfg=cfg,
            gridinit=gridinit,
            convergence=convergence,
            accepted=accepted,
        )

        if not accepted:
            failures += 1
            LOG.error(
                "Rejected %s: convergence=%s thresholds=(Ham %.3g%%, Mom %.3g%%)",
                case_dir.name,
                convergence,
                args.max_ham_pct,
                args.max_mom_pct,
            )
        else:
            LOG.info("Accepted %s: %s", case_dir.name, convergence)

    return failures


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate constraint-satisfying rotating exotic initial data via "
            "GRTresna. Includes omega=0 as the non-rotating control unless "
            "--no-control is passed."
        )
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("runs/rotating_wormhole_id"),
        help="Directory that receives one subdirectory per generated case.",
    )
    parser.add_argument(
        "--omegas",
        type=_parse_float_list,
        default=[0.05, 0.1, 0.2],
        help="Comma/space separated rotation rates to solve.",
    )
    parser.add_argument("--no-control", action="store_true", help="Do not add omega=0.")
    parser.add_argument("--amp", type=float, default=0.2, help="Exotic lump amplitude.")
    parser.add_argument("--width", type=float, default=8.0, help="Exotic lump width.")
    parser.add_argument("--mode", type=int, default=2, help="Azimuthal mode for rotation.")
    parser.add_argument("--velocity-x", type=float, default=0.0)
    parser.add_argument("--velocity-y", type=float, default=0.0)
    parser.add_argument("--velocity-z", type=float, default=0.0)
    parser.add_argument("--scalar-mass", type=float, default=0.1)
    parser.add_argument("--bh-mass", type=float, default=0.0)
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--ny", type=int, default=64)
    parser.add_argument("--nz", type=int, default=32)
    parser.add_argument("--length", type=float, default=128.0)
    parser.add_argument("--max-level", type=int, default=3)
    parser.add_argument("--block-factor", type=int, default=16)
    parser.add_argument("--max-grid-size", type=int, default=16)
    parser.add_argument("--refine-threshold", type=float, default=0.35)
    parser.add_argument("--regrid-radius", type=float, default=0.0)
    parser.add_argument("--psi-relaxation", type=float, default=0.6)
    parser.add_argument("--psi-floor", type=float, default=0.1)
    parser.add_argument("--maximal-jacobian-cap", type=float, default=25.0)
    parser.add_argument("--gridinit-nx", type=int, default=64)
    parser.add_argument("--gridinit-ny", type=int, default=64)
    parser.add_argument("--gridinit-nz", type=int, default=64)
    parser.add_argument("--target-center-x", type=float, default=32.0)
    parser.add_argument("--target-center-y", type=float, default=32.0)
    parser.add_argument("--target-center-z", type=float, default=0.0)
    parser.add_argument("--max-ham-pct", type=float, default=1.0)
    parser.add_argument("--max-mom-pct", type=float, default=1.0)
    parser.add_argument("--dry-run", action="store_true", help="Only write GRTresna params.")
    parser.add_argument("--keep-hdf5", action="store_true", help="Keep Chombo HDF5 output.")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    omegas = list(args.omegas)
    if not args.no_control and 0.0 not in omegas:
        omegas.insert(0, 0.0)

    return 1 if generate_cases(args, omegas) else 0


if __name__ == "__main__":
    raise SystemExit(main())
