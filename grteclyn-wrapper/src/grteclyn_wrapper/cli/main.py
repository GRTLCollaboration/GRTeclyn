"""CLI entry-point dispatch."""

from __future__ import annotations

from ..initial_data.candidates import resolve_initial_data_overrides
from ..initial_data.validate_guesser import run_validation
from .args import collect_overrides, parse_score_weights
from .commands.atlas import run_atlas_command
from .commands.geometry_atlas import run_geometry_atlas_command
from .commands.optimize import run_optimize_command
from .commands.pareto import run_pareto_command
from .commands.qd import run_qd_command
from .commands.warpfactory import run_warpfactory_command
from .episode import run_single, run_sweep
from .parser import build_parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.phantom_default = False
    args.initial_data_source = None
    args.score_weights = parse_score_weights(args.score_weight) or None
    overrides = collect_overrides(args.set)
    if (
        getattr(args, "seed_name", None)
        or getattr(args, "candidate_id", None)
        or getattr(args, "nonspherical_id", None)
    ):
        base_overrides, phantom_default, source = resolve_initial_data_overrides(
            seed_name=getattr(args, "seed_name", None),
            candidate_id=getattr(args, "candidate_id", None),
            nonspherical_id=getattr(args, "nonspherical_id", None),
            validation_seed=getattr(args, "validation_seed", 42),
        )
        overrides = {**base_overrides, **overrides}
        args.phantom_default = phantom_default
        args.initial_data_source = source
    if args.command == "optimize":
        return run_optimize_command(args, overrides)
    if args.command == "qd":
        return run_qd_command(args, overrides)
    if args.command == "pareto":
        return run_pareto_command(args)
    if args.command == "warpfactory":
        return run_warpfactory_command(args)
    if args.command == "validate":
        output_dir = None if args.no_write else args.output_dir
        run_validation(seed=args.seed, output_dir=output_dir)
        return 0
    if args.command == "atlas":
        return run_atlas_command(args, overrides)
    if args.command == "geometry_atlas":
        return run_geometry_atlas_command(args)
    if args.command == "sweep":
        return run_sweep(args, overrides)
    return run_single(args, overrides)
