#!/usr/bin/env python3
"""Replay a GRTresna search/QD eval at higher resolution and longer stop_time."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from grteclyn_wrapper.core.config import resolve_example, resolve_executable
from grteclyn_wrapper.core.evaluation import evaluate_overrides
from grteclyn_wrapper.grtresna.domain import GRTresnaDomainConfig
from grteclyn_wrapper.grtresna.io import read_gridinit
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig


def _load_overrides(source_eval: Path) -> dict:
    meta_path = source_eval / "metadata.json"
    if not meta_path.is_file():
        raise FileNotFoundError(f"missing metadata.json in {source_eval}")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    overrides = dict(meta.get("overrides") or {})
    if not overrides:
        raise ValueError(f"no overrides in {meta_path}")
    return overrides


def _promotion_overrides(
    base: dict,
    *,
    n_full: int,
    l_full: float,
    stop_time: float,
    plot_interval: int,
    max_level: int,
    regrid_threshold: float,
) -> dict:
    overrides = dict(base)
    half = 0.5 * l_full
    overrides.update(
        {
            "L_full": l_full,
            "N_full": int(n_full),
            "center": f"{half:g} {half:g} {half:g}",
            "stop_time": float(stop_time),
            "plot_interval": int(plot_interval),
            "checkpoint_interval": -1,
            "max_level": int(max_level),
            "regrid_threshold": float(regrid_threshold),
            "max_box_size": 32,
            "min_box_size": 16,
            "max_spatial_derivative_order": 4,
            "nan_check": 1,
        }
    )
    return overrides


def _parse_center(center: str | tuple[float, float, float]) -> tuple[float, float, float]:
    if isinstance(center, str):
        parts = center.split()
        return (float(parts[0]), float(parts[1]), float(parts[2]))
    return (float(center[0]), float(center[1]), float(center[2]))


def _validate_gridinit_alignment(
    gridinit_path: Path,
    *,
    evolution_center: tuple[float, float, float],
    tol: float = 1.0,
) -> None:
    """Reject gridinits whose baked origin does not match the evolution center."""
    gi = read_gridinit(gridinit_path)
    extent = gi.dx_xyz * np.array([gi.data.shape[2], gi.data.shape[1], gi.data.shape[0]])
    implied_center = tuple(float(gi.origin[i] + 0.5 * extent[i]) for i in range(3))
    expected_origin = tuple(
        float(evolution_center[i] - 0.5 * extent[i]) for i in range(3)
    )
    origin_delta = tuple(float(gi.origin[i] - expected_origin[i]) for i in range(3))
    center_delta = tuple(float(implied_center[i] - evolution_center[i]) for i in range(3))
    if any(abs(delta) > tol for delta in center_delta):
        raise ValueError(
            "gridinit is not aligned with the requested evolution center: "
            f"gridinit implies matter center {implied_center}, "
            f"evolution center {evolution_center} "
            f"(delta {center_delta}); origin {tuple(gi.origin)} "
            f"expected {expected_origin} (delta {origin_delta}). "
            "Re-solve GRTresna for the larger box instead of reusing an old .gridinit."
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source_eval",
        type=Path,
        help="Source eval dir, e.g. runs/grtresna_qd/qd_.../eval_000057",
    )
    parser.add_argument("--name", required=True, help="Output episode name")
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=None,
        help="Parent output dir (default: runs/grtresna_promote)",
    )
    parser.add_argument("--gpu", default="0", help="CUDA device id(s)")
    parser.add_argument("--n-full", type=int, default=96, help="Promoted base resolution")
    parser.add_argument(
        "--l-full",
        type=float,
        default=None,
        help="Evolution box width in code units (default: same as --n-full)",
    )
    parser.add_argument("--max-level", type=int, default=3, help="GPU AMR max_level")
    parser.add_argument(
        "--regrid-threshold",
        type=float,
        default=0.02,
        help="AMR regrid threshold (lower => more refinement)",
    )
    parser.add_argument("--stop-time", type=float, default=16.0)
    parser.add_argument(
        "--plot-interval",
        type=int,
        default=48,
        help="Plot cadence for long runs (~33 frames at t=16)",
    )
    parser.add_argument("--ftl-L", type=float, default=8.0)
    parser.add_argument("--grtresna-ranks", type=int, default=8)
    parser.add_argument("--grtresna-iterations", type=int, default=30)
    parser.add_argument("--grtresna-max-level", type=int, default=3)
    parser.add_argument("--grtresna-refine-threshold", type=float, default=0.5)
    parser.add_argument("--grtresna-regrid-radius", type=float, default=0.0)
    parser.add_argument("--grtresna-jacobian-cap", type=float, default=25.0)
    parser.add_argument("--grtresna-domain-l", type=float, default=128.0)
    parser.add_argument("--consumer-keep-last", type=int, default=2)
    parser.add_argument(
        "--gridinit",
        type=Path,
        default=None,
        help="Reuse an existing initial_data.gridinit and skip the GRTresna solve",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    source_eval = args.source_eval.expanduser().resolve()
    repo_root = source_eval
    while repo_root != repo_root.parent:
        if (repo_root / "Examples" / "RadialRecipe").is_dir():
            break
        repo_root = repo_root.parent
    else:
        print("Could not locate GRTeclyn repo root from source eval path.", file=sys.stderr)
        return 2

    runs_dir = (
        args.runs_dir.expanduser().resolve()
        if args.runs_dir is not None
        else repo_root / "runs" / "grtresna_promote"
    )
    runs_dir.mkdir(parents=True, exist_ok=True)

    base_overrides = _load_overrides(source_eval)
    n = int(args.n_full)
    l_full = float(args.l_full) if args.l_full is not None else float(n)
    overrides = _promotion_overrides(
        base_overrides,
        n_full=n,
        l_full=l_full,
        stop_time=args.stop_time,
        plot_interval=args.plot_interval,
        max_level=args.max_level,
        regrid_threshold=args.regrid_threshold,
    )

    domain = GRTresnaDomainConfig(
        full_z=True,
        l_full=l_full,
        n_full=n,
        grtresna_l=float(args.grtresna_domain_l),
        grtresna_nx=n,
        grtresna_ny=n,
        grtresna_nz=n,
        gridinit_nx=n,
        gridinit_ny=n,
        gridinit_nz=n,
    )
    overrides = {**overrides, **domain.evolution_overrides()}

    use_grtresna = args.gridinit is None
    grtresna_config = None
    if use_grtresna:
        grtresna_config = GRTresnaConfig(
            mpi_ranks=args.grtresna_ranks,
            max_NL_iterations=args.grtresna_iterations,
            max_level=args.grtresna_max_level,
            refine_threshold=args.grtresna_refine_threshold,
            regrid_radius=args.grtresna_regrid_radius,
            coefficient_average_type="harmonic",
            psi_relaxation=1.0,
            psi_floor=-1.0,
            maximal_jacobian_cap=args.grtresna_jacobian_cap,
            bh1_bare_mass=0.0,
            bh1_spin=(0.0, 0.0, 0.0),
            bh2_bare_mass=0.0,
            dphi=0.0,
            dpi=0.0,
            cleanup=True,
        )
        grtresna_config = domain.apply_to_solver(grtresna_config)
    else:
        gridinit = args.gridinit.expanduser().resolve()
        if not gridinit.is_file():
            raise FileNotFoundError(f"gridinit not found: {gridinit}")
        evolution_center = _parse_center(str(overrides["center"]))
        _validate_gridinit_alignment(
            gridinit,
            evolution_center=evolution_center,
        )
        overrides["recipe_initial_data_file"] = str(gridinit)

    example = resolve_example("RadialRecipe")
    template = example.template
    executable = None
    if not args.dry_run:
        executable = resolve_executable(
            None,
            example=example,
            mpi_ranks=1,
            comp="gnu",
            cuda=True,
            debug=False,
        )

    mode = "grtresna+gpu" if use_grtresna else "gpu-only"
    print(
        f"[replay] {source_eval.name} -> {runs_dir / args.name} "
        f"(L={l_full:g}, N={n}, t={args.stop_time}, GPU={args.gpu}, mode={mode})",
        flush=True,
    )
    result = evaluate_overrides(
        overrides,
        out_dir=runs_dir,
        name=args.name,
        example=example,
        template=template,
        executable=executable,
        constrained=False,
        phantom=False,
        use_preflight=False,
        cuda_devices=args.gpu,
        check_params=True,
        dry_run=args.dry_run,
        target_stop_time=args.stop_time,
        objective_mode="ftl_first",
        ftl_L=args.ftl_L,
        consume_plotfiles=True,
        consumer_radii=(4.0, 8.0),
        consumer_keep_last=args.consumer_keep_last,
        grtresna=use_grtresna,
        grtresna_base=grtresna_config,
        grtresna_solved_ftl_gate=False,
    )
    print(
        json.dumps(
            {
                "episode": result.episode_path,
                "score": result.score,
                "exit_code": result.exit_code,
                "reason": result.reason,
                "notes": result.notes,
            },
            indent=2,
        ),
        flush=True,
    )
    if result.exit_code not in (None, 0):
        return int(result.exit_code or 1)
    if result.reason:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
