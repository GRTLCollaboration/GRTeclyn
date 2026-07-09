#!/usr/bin/env python3
"""Project a geometry-first RadialRecipe motif through GRTresna.

This script is intentionally separate from the matter-first ``optimize`` and
``qd`` launchers.  Fitted matter is never pushed directly into GRTeclyn; the
workflow is:

  geometry-first episode -> motif.json -> fitted_matter.json -> GRTresna solve
  -> motif preservation check -> optional post-load GRTeclyn gate -> evolution
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from grteclyn_wrapper.core.config import (
    DEFAULT_RADIAL_RECIPE_TEMPLATE,
    resolve_example,
    resolve_executable,
)
from grteclyn_wrapper.core.episode import create_episode, update_metadata, write_json
from grteclyn_wrapper.core.evaluation import evaluate_overrides
from grteclyn_wrapper.core.params import write_params
from grteclyn_wrapper.grtresna.fit.motif import (
    build_grtresna_config_from_fitted,
    fit_matter_from_motif,
    fit_momentum_from_motif,
    read_fitted_matter_json,
    write_fitted_matter_json,
    write_momentum_target_json,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, parse_convergence, solve
from grteclyn_wrapper.initial_data.motif import (
    extract_motif_from_episode,
    motif_from_alcubierre,
    read_motif_json,
    write_motif_json,
)
from grteclyn_wrapper.projection.iterate import IterationConfig, run_iterate
from grteclyn_wrapper.projection.mismatch import warp_factory_cross_check
from grteclyn_wrapper.projection.motif_preservation import (
    compare_motif_preservation,
    write_preservation_report,
)
from grteclyn_wrapper.projection.postload_gate import (
    PostLoadGateConfig,
    run_postload_gate,
    write_gate_result,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source_episode",
        type=Path,
        nargs="?",
        default=None,
        help="Geometry-first episode directory (eval_XXXXXX or smoke run). "
             "Optional when --target alcubierre is used.",
    )
    parser.add_argument(
        "--target",
        choices=("episode", "alcubierre"),
        default="episode",
        help="Motif source: 'episode' (default) extracts from source_episode; "
             "'alcubierre' builds a warp-bubble motif from --warp-* params.",
    )
    parser.add_argument(
        "--warp-velocity",
        type=float,
        default=0.5,
        help="Alcubierre bubble velocity v (fraction of c)",
    )
    parser.add_argument(
        "--warp-bubble-radius",
        type=float,
        default=2.0,
        help="Alcubierre bubble radius R",
    )
    parser.add_argument(
        "--warp-sigma",
        type=float,
        default=2.0,
        help="Alcubierre bubble wall thickness sigma",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory for projection artifacts",
    )
    parser.add_argument(
        "--mode",
        choices=("fit-only", "solve-only", "solve-and-evolve"),
        default="fit-only",
        help="fit-only writes motif/fitted matter; solve-only runs GRTresna; solve-and-evolve also replays in GRTeclyn",
    )
    parser.add_argument("--max-lumps", type=int, default=3)
    parser.add_argument("--ftl-L", type=float, default=None)
    parser.add_argument("--phantom", action="store_true", default=True)
    parser.add_argument("--no-phantom", dest="phantom", action="store_false")
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--grtresna-L", type=float, default=128.0)
    parser.add_argument("--gridinit-n", type=int, default=64)
    parser.add_argument("--cuda-device", type=str, default=None)
    parser.add_argument("--stop-time", type=float, default=2.0)
    parser.add_argument("--plot-interval", type=int, default=10)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--motif-json",
        type=Path,
        default=None,
        help="Reuse an existing motif.json (solve-only / solve-and-evolve)",
    )
    parser.add_argument(
        "--fitted-matter-json",
        type=Path,
        default=None,
        help="Reuse an existing fitted_matter.json",
    )
    parser.add_argument(
        "--skip-postload-gate",
        action="store_true",
        help="Skip the short GRTeclyn post-load constraint check",
    )

    # Iteration loop flags (geometry-first matter adjustment)
    parser.add_argument(
        "--iterate",
        type=int,
        default=0,
        metavar="N",
        help="Max CMA-ES evals for iterative matter adjustment (0 = one-shot, default)",
    )
    parser.add_argument("--iterate-popsize", type=int, default=8)
    parser.add_argument("--iterate-sigma0", type=float, default=0.2)
    parser.add_argument(
        "--max-concurrent-grtresna",
        type=int,
        default=6,
        help="Max parallel GRTresna solves during iteration",
    )
    parser.add_argument(
        "--promote-retention-min",
        type=float,
        default=0.5,
        help="Minimum retention_score to promote iterated candidate to evolution",
    )
    return parser


def _default_grtresna_base(args: argparse.Namespace) -> GRTresnaConfig:
    n = args.gridinit_n
    L = args.grtresna_L
    if getattr(args, "target", "episode") == "alcubierre":
        # Full-box domain: the Alcubierre toroidal ring extends to ±R in y/z
        # around the transport axis, so octant symmetry doesn't apply.
        # Centre the physics at (L/2, L/2, L/2).
        return GRTresnaConfig(
            mpi_ranks=args.mpi_ranks,
            N=(n, n, n),
            L=L,
            gridinit_nx=n,
            gridinit_ny=n,
            gridinit_nz=n,
            scalar_mass=0.0,
            dphi=0.0,
            dpi=0.0,
            bh1_bare_mass=0.0,
            bh1_spin=(0.0, 0.0, 0.0),
            max_NL_iterations=200,
            nl_stall_tolerance=0.005,
            lo_boundary=(0, 0, 0),
            target_center=(L / 2.0, L / 2.0, L / 2.0),
        )
    return GRTresnaConfig(
        mpi_ranks=args.mpi_ranks,
        N=(n, n, n // 2),
        L=L,
        gridinit_nx=n,
        gridinit_ny=n,
        gridinit_nz=n,
        scalar_mass=0.0,
        dphi=0.0,
        dpi=0.0,
        bh1_bare_mass=0.0,
        bh1_spin=(0.0, 0.0, 0.0),
        # Exotic matter needs more iterations and tighter stall tolerance
        max_NL_iterations=200,
        nl_stall_tolerance=0.005,
    )


def _write_projection_report(
    out_dir: Path,
    *,
    motif_path: Path,
    fitted_path: Path,
    preservation_path: Path | None,
    gate_path: Path | None,
    convergence: dict[str, float] | None,
    iteration: dict[str, Any] | None = None,
) -> Path:
    report = {
        "motif_json": str(motif_path),
        "fitted_matter_json": str(fitted_path),
        "preservation_report": str(preservation_path) if preservation_path else None,
        "postload_gate": str(gate_path) if gate_path else None,
        "grtresna_convergence": convergence,
        "iteration": iteration,
    }
    report_path = out_dir / "projection_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report_path


def run_projection(args: argparse.Namespace) -> int:
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    motif_path = out_dir / "motif.json"
    fitted_path = out_dir / "fitted_matter.json"
    momentum_path = out_dir / "momentum_target.json"
    gridinit_path = out_dir / "initial_data.gridinit"
    preservation_path = out_dir / "projection_report_preservation.json"
    gate_path = out_dir / "postload_gate.json"

    if args.motif_json is not None:
        motif = read_motif_json(args.motif_json)
    elif args.target == "alcubierre":
        motif = motif_from_alcubierre(
            velocity=args.warp_velocity,
            bubble_radius=args.warp_bubble_radius,
            sigma=args.warp_sigma,
            ftl_L=args.ftl_L,
        )
        write_motif_json(motif, motif_path)
        # Cross-check reconstructed rho against the exact T = G/8pi
        cross_check = warp_factory_cross_check(motif, L=args.ftl_L)
        cross_check_path = out_dir / "warp_factory_crosscheck.json"
        cross_check_path.write_text(
            json.dumps(cross_check.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"Warp-factory cross-check: rho_l2={cross_check.rho_l2:.6f}, "
              f"NEC violation={cross_check.nec_violation_fraction:.3f}, "
              f"exotic_energy={cross_check.exotic_energy:.4f}")
        print(f"Wrote {cross_check_path}")
    else:
        if args.source_episode is None:
            print("Error: source_episode is required when --target episode (default)",
                  file=sys.stderr)
            return 1
        motif = extract_motif_from_episode(
            args.source_episode,
            phantom=args.phantom,
            ftl_L=args.ftl_L,
        )
        write_motif_json(motif, motif_path)

    if args.fitted_matter_json is not None:
        fitted = read_fitted_matter_json(args.fitted_matter_json)
    else:
        fitted = fit_matter_from_motif(motif, max_lumps=args.max_lumps)
        fitted = fit_momentum_from_motif(motif, fitted)
        # For Alcubierre target, shift lump centres from physics origin to
        # the GRTresna domain centre (L/2, L/2, L/2) so they sit inside [0,L]³.
        if args.target == "alcubierre":
            cx = args.grtresna_L / 2.0
            shifted_lumps = []
            for lump in fitted.lumps:
                updated = dict(lump)
                c = updated["center"]
                updated["center"] = (c[0] + cx, c[1] + cx, c[2] + cx)
                shifted_lumps.append(updated)
            from dataclasses import replace as _dc_replace
            fitted = _dc_replace(fitted, lumps=tuple(shifted_lumps))
        write_fitted_matter_json(fitted, fitted_path)
        write_momentum_target_json(fitted.momentum_target, momentum_path)

    if args.mode == "fit-only":
        _write_projection_report(
            out_dir,
            motif_path=motif_path if motif_path.exists() else args.motif_json or motif_path,
            fitted_path=fitted_path if fitted_path.exists() else args.fitted_matter_json or fitted_path,
            preservation_path=None,
            gate_path=None,
            convergence=None,
        )
        print(f"Wrote {motif_path}")
        print(f"Wrote {fitted_path}")
        print(f"Wrote {momentum_path}")
        return 0

    cfg = build_grtresna_config_from_fitted(
        fitted,
        base=_default_grtresna_base(args),
    )
    if args.dry_run:
        from grteclyn_wrapper.grtresna.solver import write_grtresna_params

        write_grtresna_params(cfg, out_dir / "grtresna_params.txt")
        _write_projection_report(
            out_dir,
            motif_path=motif_path,
            fitted_path=fitted_path,
            preservation_path=None,
            gate_path=None,
            convergence=None,
        )
        print(f"Dry-run wrote {out_dir / 'grtresna_params.txt'}")
        return 0

    # --- Iterative matter adjustment (when --iterate N > 0) ---------------
    iterate_n = getattr(args, "iterate", 0)
    if iterate_n > 0:
        iter_cfg = IterationConfig(
            max_evals=iterate_n,
            popsize=args.iterate_popsize,
            sigma0=args.iterate_sigma0,
            max_concurrent_grtresna=args.max_concurrent_grtresna,
            mpi_ranks=args.mpi_ranks,
            retention_threshold=args.promote_retention_min,
            L=args.ftl_L,
            gridinit_n=args.gridinit_n,
            grtresna_L=args.grtresna_L,
        )
        iter_result = run_iterate(
            motif,
            fitted,
            out_dir=out_dir / "iterate",
            config=iter_cfg,
            base_grtresna_cfg=_default_grtresna_base(args),
        )
        print(
            f"Iteration complete: {iter_result.total_evals} evals, "
            f"best_fitness={iter_result.best_fitness:.6f}, "
            f"retention={iter_result.best_preservation_score:.3f}, "
            f"converged={iter_result.converged}"
        )
        # Use the best iterated output as the final gridinit / fitted matter
        if iter_result.best_gridinit_path is not None and iter_result.best_gridinit_path.exists():
            import shutil
            shutil.copy2(iter_result.best_gridinit_path, gridinit_path)
            write_fitted_matter_json(iter_result.best_fitted_matter, fitted_path)
            fitted = iter_result.best_fitted_matter
        else:
            print("Iteration produced no valid gridinit; falling back to one-shot solve",
                  file=sys.stderr)
            solve(cfg, work_dir=out_dir / "grtresna", gridinit_path=gridinit_path)

        convergence = parse_convergence(out_dir / "iterate" / "best" / "grtresna")
        if convergence is None:
            convergence = parse_convergence(out_dir / "grtresna")
    else:
        # --- One-shot solve (original behavior) ----------------------------
        solve(
            cfg,
            work_dir=out_dir / "grtresna",
            gridinit_path=gridinit_path,
        )
        convergence = parse_convergence(out_dir / "grtresna")

    preservation = compare_motif_preservation(
        motif,
        gridinit_path,
        fitted_matter_path=fitted_path,
        ftl_L=args.ftl_L,
    )
    write_preservation_report(preservation, preservation_path)

    gate_result = None
    if not args.skip_postload_gate:
        gate_result = run_postload_gate(
            gridinit_path,
            out_dir=out_dir / "postload_gate_run",
            cuda_devices=args.cuda_device,
            dry_run=False,
        )
        write_gate_result(gate_result, gate_path)

    iteration_info = None
    if iterate_n > 0:
        iteration_info = {
            "max_evals": iterate_n,
            "total_evals": iter_result.total_evals,
            "generations": iter_result.generations,
            "best_fitness": iter_result.best_fitness,
            "best_preservation_score": iter_result.best_preservation_score,
            "converged": iter_result.converged,
            "fitness_history": iter_result.fitness_history,
            "notes": iter_result.notes,
        }

    report_path = _write_projection_report(
        out_dir,
        motif_path=motif_path,
        fitted_path=fitted_path,
        preservation_path=preservation_path,
        gate_path=gate_path if gate_result is not None else None,
        convergence=convergence,
        iteration=iteration_info,
    )

    print(f"Wrote {gridinit_path}")
    print(f"Wrote {preservation_path} passed={preservation.passed}")
    if gate_result is not None:
        print(f"Wrote {gate_path} passed={gate_result.passed}")
    print(f"Wrote {report_path}")

    if not preservation.passed:
        print("Projection rejected: motif not preserved", file=sys.stderr)
        return 2
    if gate_result is not None and not gate_result.passed:
        print("Projection rejected: post-load constraints failed", file=sys.stderr)
        return 3

    if args.mode == "solve-and-evolve":
        example = resolve_example("RadialRecipe")
        executable = resolve_executable(example=example)
        evaluation = evaluate_overrides(
            {
                "recipe_initial_data_file": str(gridinit_path),
                "stop_time": args.stop_time,
                "plot_interval": args.plot_interval,
            },
            out_dir=out_dir,
            name="projection_evolve",
            example=example,
            template=DEFAULT_RADIAL_RECIPE_TEMPLATE,
            executable=executable,
            constrained=False,
            phantom=False,
            use_preflight=False,
            cuda_devices=args.cuda_device,
            dry_run=False,
            target_stop_time=args.stop_time,
            grtresna=False,
        )
        write_json(out_dir / "projection_evolve_score.json", {
            "score": evaluation.score,
            "components": evaluation.components,
            "episode_path": evaluation.episode_path,
            "exit_code": evaluation.exit_code,
        })
        if evaluation.exit_code not in (0, None):
            return int(evaluation.exit_code or 1)

    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    return run_projection(args)


if __name__ == "__main__":
    raise SystemExit(main())
