#!/usr/bin/env python3
"""Replay a GRTresna search/QD eval at higher resolution and longer stop_time."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from grteclyn_wrapper.core.config import resolve_example, resolve_executable
from grteclyn_wrapper.core.evaluation import evaluate_overrides
from grteclyn_wrapper.core.params import regrid_intervals_for_max_level
from grteclyn_wrapper.grtresna.domain import GRTresnaDomainConfig
from grteclyn_wrapper.grtresna.matter.wiring import (
    EVOLUTION_MATTER_KEYS,
    GRTRESNA_COMPLEX_SCALAR_MODEL,
    evolution_overrides_from_metadata,
    plot_vars_for_bicomplex_scalar,
    plot_vars_for_complex_scalar,
    plot_vars_for_independent_scalars,
    read_matter_metadata,
)
from grteclyn_wrapper.grtresna.matter.models import GRTRESNA_BICOMPLEX_SCALAR_MODEL
from grteclyn_wrapper.grtresna.profiles.qball_couplings import QBallCouplings
from grteclyn_wrapper.grtresna.io import read_gridinit
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig
from grteclyn_wrapper.search.grtresna_convergence_gate import GRTresnaConvergenceConfig
from grteclyn_wrapper.search.optimize.candidates import _clamp_trajectory_speed

# Match scripts/campaigns/lib/search_common.sh (QD stage-0 defaults).
_QD_PLOT_INTERVAL = 320
_BOSON_FRAMES_FIELDS = (
    "scalar_activity phi Pi phi_lump0 Pi_lump0 phi_lump1 "
    "chi chi_minus_1 local_speed shift1 rho_req"
)
_SCALAR_FRAMES_FIELDS = (
    "lump_activity scalar_activity phi_lump_sum Pi_lump_sum "
    "chi chi_minus_1 local_speed shift1 rho_req"
)


def _apply_replay_consumer_env(overrides: Mapping[str, Any]) -> None:
    """Apply QD-matching frame extraction env when not already set by the caller."""
    os.environ.setdefault("GRTECLYN_FRAMES", "1")
    if os.environ.get("GRTECLYN_FRAMES_FIELDS", "").strip():
        return

    model = str(
        overrides.get("grtresna_matter_model")
        or overrides.get("recipe_matter_model", "")
    ).lower()
    sector = str(overrides.get("grtresna_matter_sector", "")).lower()
    if (
        "bicomplex" in model
        or "complex_scalar" in model
        or sector == "boson_star"
    ):
        os.environ["GRTECLYN_FRAMES_FIELDS"] = _BOSON_FRAMES_FIELDS
        os.environ.setdefault("GRTECLYN_PROJECTION_FIELDS", "scalar_activity phi")
    else:
        os.environ["GRTECLYN_FRAMES_FIELDS"] = _SCALAR_FRAMES_FIELDS
        os.environ.setdefault(
            "GRTECLYN_PROJECTION_FIELDS", "lump_activity scalar_activity"
        )
    os.environ.setdefault("GRTECLYN_FRAMES_ZOOM", "none")
    os.environ.setdefault("GRTECLYN_PROJECTION_AXES", "x y z")
    os.environ.setdefault("GRTECLYN_PROJECTION_METHOD", "mip")


def _load_overrides(source_eval: Path) -> dict:
    meta_path = source_eval / "metadata.json"
    if not meta_path.is_file():
        raise FileNotFoundError(f"missing metadata.json in {source_eval}")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    overrides = dict(meta.get("overrides") or {})
    if not overrides:
        raise ValueError(f"no overrides in {meta_path}")
    return overrides


def _parse_params_value(raw: str) -> int | float | str:
    value = raw.strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        return value[1:-1]
    if " " in value:
        return value
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _load_matter_replay_overrides(source_eval: Path) -> dict[str, int | float | str]:
    """Restore evolution matter params for gridinit-only replay.

    Prefer ``initial_data.matter.json`` (full metadata including scalar_mu);
    fall back to scanning ``params.txt`` for ``EVOLUTION_MATTER_KEYS``.
    """
    matter_json = source_eval / "initial_data.matter.json"
    if matter_json.is_file():
        meta = read_matter_metadata(matter_json)
        overrides = evolution_overrides_from_metadata(meta)
        return {
            key: value
            for key, value in overrides.items()
            if key in EVOLUTION_MATTER_KEYS
        }

    params_path = source_eval / "params.txt"
    if not params_path.is_file():
        return {}

    replay_keys = set(EVOLUTION_MATTER_KEYS)
    overrides: dict[str, int | float | str] = {}
    for line in params_path.read_text(encoding="utf-8").splitlines():
        body = line.split("#", 1)[0].strip()
        if "=" not in body:
            continue
        key, raw_value = (part.strip() for part in body.split("=", 1))
        if key in replay_keys:
            overrides[key] = _parse_params_value(raw_value)
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
    promoted_max_level = int(max_level)
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
            "max_level": promoted_max_level,
            # CMA-ES/QD metadata may carry a shorter list from a lower max_level;
            # AMReX aborts if regrid_interval has fewer entries than max_level+1.
            "regrid_interval": regrid_intervals_for_max_level(promoted_max_level),
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
        default=_QD_PLOT_INTERVAL,
        help="Plotfile cadence in steps (QD default: 320 → ~6 dumps at t=16, dt≈0.01).",
    )
    parser.add_argument("--ftl-L", type=float, default=8.0)
    parser.add_argument("--grtresna-ranks", type=int, default=8)
    parser.add_argument("--grtresna-iterations", type=int, default=30)
    parser.add_argument("--grtresna-max-level", type=int, default=3)
    parser.add_argument("--grtresna-refine-threshold", type=float, default=0.5)
    parser.add_argument("--grtresna-regrid-radius", type=float, default=0.0)
    parser.add_argument("--grtresna-jacobian-cap", type=float, default=25.0)
    parser.add_argument(
        "--grtresna-domain-l",
        type=float,
        default=None,
        help="GRTresna solve box width (default: same as --l-full)",
    )
    parser.add_argument("--grtresna-timeout", type=int, default=3600)
    parser.add_argument("--grtresna-max-ham-pct", type=float, default=5.0)
    parser.add_argument("--grtresna-max-mom-pct", type=float, default=5.0)
    parser.add_argument("--consumer-keep-last", type=int, default=2)
    parser.add_argument(
        "--evolving-geodesic",
        action="store_true",
        help="Enable end-of-run 4D evolving null-geodesic trace in scoring.",
    )
    parser.add_argument(
        "--gridinit",
        type=Path,
        default=None,
        help="Reuse an existing initial_data.gridinit and skip the GRTresna solve",
    )
    parser.add_argument(
        "--objective-mode",
        default="general_ftl",
        choices=[
            "weighted",
            "ftl_first",
            "robust_ftl",
            "general_ftl",
            "critical_collapse",
        ],
        help="Scoring objective (default: general_ftl for HQ v20 replays).",
    )
    parser.add_argument(
        "--extra-override",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Extra key=value overrides applied after all other overrides "
        "(e.g. --extra-override lapse_power=0.0). May be repeated.",
    )
    parser.add_argument(
        "--qball-preset",
        choices=("standard", "stiff", "compact"),
        default=None,
        help="Apply QBallCouplings preset overrides (standard=λ160/μ5333 ω0.4, "
        "stiff=λ640/μ85333 ω0.4, compact=λ640/μ85333 ω0.8 — thick-wall, "
        "box-fitting localized soliton; use with --qball-ode-profile).",
    )
    parser.add_argument(
        "--qball-equilibrium-amplitude",
        action="store_true",
        help="Cap trajectory well_depth to flat-space Q-ball core amplitude √(3λ/4μ).",
    )
    parser.add_argument(
        "--qball-ode-profile",
        action="store_true",
        help="Seed lumps from the flat-space Q-ball radial ODE (profile 3) instead "
        "of sech. GRTresna C++ now loads the same tabulated phi0(r) (qball_profile_path) "
        "so the solve and the gridinit repaint are consistent.",
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
    grtresna_domain_l = (
        float(args.grtresna_domain_l) if args.grtresna_domain_l is not None else l_full
    )
    overrides = _promotion_overrides(
        base_overrides,
        n_full=n,
        l_full=l_full,
        stop_time=args.stop_time,
        plot_interval=args.plot_interval,
        max_level=args.max_level,
        regrid_threshold=args.regrid_threshold,
    )

    # Apply --extra-override KEY=VALUE pairs last (highest priority).
    if args.qball_preset is not None:
        preset = {
            "standard": QBallCouplings.standard,
            "stiff": QBallCouplings.stiff,
            "compact": QBallCouplings.compact,
        }[args.qball_preset]()
        overrides.update(
            preset.seed_overrides(
                equilibrium_amplitude=args.qball_equilibrium_amplitude,
                ode_profile=args.qball_ode_profile,
            )
        )
    elif args.qball_equilibrium_amplitude or args.qball_ode_profile:
        if args.qball_equilibrium_amplitude:
            overrides["grtresna_qball_equilibrium_amplitude"] = 1
        if args.qball_ode_profile:
            overrides["grtresna_qball_ode_profile"] = 1

    for token in args.extra_override:
        if "=" not in token:
            print(f"[replay] ignoring malformed --extra-override {token!r}", file=sys.stderr)
            continue
        key, raw = token.split("=", 1)
        overrides[key.strip()] = _parse_params_value(raw)

    # Enforce the sub-luminal / adiabatic trajectory-speed cap on replay too:
    # historical elites (e.g. eval 122) carry v_t = R0*|omega_rot| up to ~6c,
    # which no soliton can follow.  Clamp omega_rot per lump (default 0.3c,
    # override via --extra-override trajectory_v_max=...) before the genome
    # reaches the GRTresna seed and the GRTeclyn co-moving trap.
    _clamp_trajectory_speed(overrides)

    domain = GRTresnaDomainConfig(
        full_z=True,
        l_full=l_full,
        n_full=n,
        grtresna_l=grtresna_domain_l,
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
    grtresna_convergence_config = GRTresnaConvergenceConfig(
        max_ham_pct=args.grtresna_max_ham_pct,
        max_mom_pct=args.grtresna_max_mom_pct,
    )
    if use_grtresna:
        grtresna_config = GRTresnaConfig(
            mpi_ranks=args.grtresna_ranks,
            max_NL_iterations=args.grtresna_iterations,
            timeout=args.grtresna_timeout,
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
        overrides.update(_load_matter_replay_overrides(source_eval))
        matter_model = overrides.get("recipe_matter_model")
        if matter_model == GRTRESNA_BICOMPLEX_SCALAR_MODEL:
            # Force the full canonical+phantom channel set so HQ frames can show
            # both fields even if the source eval was produced before the
            # phantom plot_vars existed.
            overrides["amr.plot_vars"] = plot_vars_for_bicomplex_scalar()
        elif matter_model == GRTRESNA_COMPLEX_SCALAR_MODEL and "amr.plot_vars" not in overrides:
            overrides["amr.plot_vars"] = plot_vars_for_complex_scalar()
        num_fields = overrides.get("recipe_num_scalar_fields")
        if num_fields and "amr.plot_vars" not in overrides:
            overrides["amr.plot_vars"] = plot_vars_for_independent_scalars(int(num_fields))
        overrides["recipe_initial_data_file"] = str(gridinit)

    _apply_replay_consumer_env(overrides)

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
        objective_mode=args.objective_mode,
        ftl_L=args.ftl_L,
        consume_plotfiles=True,
        consumer_radii=(4.0, 8.0),
        consumer_keep_last=args.consumer_keep_last,
        consumer_evolving_geodesic=args.evolving_geodesic,
        grtresna=use_grtresna,
        grtresna_base=grtresna_config,
        grtresna_solved_ftl_gate=False,
        grtresna_convergence_config=grtresna_convergence_config,
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
