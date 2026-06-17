"""CMA-ES optimize command handler."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ...core.config import resolve_example, resolve_executable
from ...initial_data.seeds import get_seed
from ...search.optimize import run_optimize
from ..grtresna_args import (
    ensure_radial_recipe_for_grtresna,
    grtresna_postload_gate_enabled,
    grtresna_solved_ftl_gate_enabled,
    postload_gate_config_from_args,
    solved_ftl_gate_config_from_args,
)
from ..grtresna_context import build_grtresna_search_context


def run_optimize_command(args: argparse.Namespace, base_overrides: dict) -> int:
    ensure_radial_recipe_for_grtresna(args, label="optimize")
    example = resolve_example(args.example)
    executable = None
    if not args.dry_run:
        executable = resolve_executable(
            args.executable,
            example=example,
            mpi_ranks=args.mpi_ranks,
            comp=args.comp,
            cuda=not args.no_cuda,
            debug=args.debug,
        )

    template = Path(args.template).expanduser().resolve() if args.template else None
    ctx = build_grtresna_search_context(args, base_overrides)

    x0 = None
    if getattr(args, "seed_name", None):
        seed_obj = get_seed(args.seed_name)
        x0 = []
        for dim in ctx.search_space:
            seed_value = float(seed_obj.overrides.get(dim.param_key, dim.center))
            x0.append(max(dim.lower, min(dim.upper, seed_value)))

    gpu_ids = getattr(args, "gpu_ids", None)

    use_constrained = getattr(args, "constrained", False)
    use_phantom = getattr(args, "phantom", False)
    use_preflight = getattr(args, "preflight", False)
    if args.command == "optimize":
        use_constrained = True if not hasattr(args, "_no_constrained") else use_constrained
        use_phantom = not getattr(args, "no_phantom", False)
        use_preflight = True if not hasattr(args, "_no_preflight") else use_preflight

    result = run_optimize(
        search_space=ctx.search_space,
        runs_dir=Path(args.runs_dir).expanduser().resolve(),
        executable=executable,
        max_generations=args.max_generations,
        population_size=args.population_size,
        sigma0=args.sigma0,
        seed=args.seed,
        base_overrides=ctx.base_overrides,
        template=template,
        example=example,
        name=args.name,
        dry_run=args.dry_run,
        constrained=use_constrained,
        phantom=use_phantom,
        use_preflight=use_preflight,
        cuda_devices=args.cuda_devices,
        gpu_ids=gpu_ids,
        check_params=not args.skip_check_params,
        x0=x0,
        consume_plotfiles=not getattr(args, "no_consume_plotfiles", False),
        consumer_radii=getattr(args, "consumer_radii", [4.0, 8.0]),
        consumer_keep_last=getattr(args, "consumer_keep_last", 1),
        score_weights=getattr(args, "score_weights", None),
        objective_mode=getattr(args, "objective_mode", "weighted"),
        ftl_L=getattr(args, "ftl_L", None),
        surrogate=getattr(args, "surrogate", False),
        surrogate_keep_fraction=getattr(args, "surrogate_keep_fraction", 0.5),
        grtresna=ctx.use_grtresna,
        grtresna_config=ctx.grtresna_config,
        grtresna_solved_ftl_gate=grtresna_solved_ftl_gate_enabled(
            args,
            use_grtresna=ctx.use_grtresna,
            objective_mode=getattr(args, "objective_mode", "weighted"),
        ),
        solved_ftl_gate_config=solved_ftl_gate_config_from_args(args),
        grtresna_convergence_config=ctx.grtresna_convergence_config,
        grtresna_postload_gate=grtresna_postload_gate_enabled(args, use_grtresna=ctx.use_grtresna),
        postload_gate_config=(
            postload_gate_config_from_args(args) if ctx.use_grtresna else None
        ),
        warm_start_trajectories=[
            Path(p).expanduser().resolve()
            for p in getattr(args, "warm_start_trajectory", [])
        ],
        warm_start_top_k=getattr(args, "warm_start_top_k", 8),
        warm_start_jitter=getattr(args, "warm_start_jitter", 0.08),
        random_injection_fraction=getattr(args, "random_injection_fraction", 0.0),
        exotic_injection_fraction=getattr(args, "exotic_injection_fraction", 0.0),
        keep_top_eval_dirs=getattr(args, "keep_top_eval_dirs", 0),
        ftl_retention_enabled=getattr(args, "ftl_retention", False),
    )
    print(json.dumps({
        "best_score": result.best_score,
        "best_episode": result.best_episode,
        "best_params": result.best_params,
        "generations": result.generations,
        "evaluations": result.evaluations,
    }, indent=2))
    return 0
