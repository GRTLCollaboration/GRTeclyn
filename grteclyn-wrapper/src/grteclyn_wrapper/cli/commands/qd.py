"""MAP-Elites qd command handler."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from ...core.config import resolve_example, resolve_executable
from ...search.qd_search import run_qd_search
from ..grtresna_args import (
    ensure_radial_recipe_for_grtresna,
    grtresna_postload_gate_enabled,
    grtresna_solved_ftl_gate_enabled,
    postload_gate_config_from_args,
)
from ..grtresna_context import build_grtresna_search_context
from ..pipeline_args import pipeline_settings_from_args


def load_seed_overrides(eval_dirs: list[str] | None) -> list[dict[str, Any]]:
    """Read ``metadata.json`` overrides from prior eval dirs to warm-start QD."""
    seeds: list[dict[str, Any]] = []
    for raw in eval_dirs or []:
        meta_path = Path(raw).expanduser() / "metadata.json"
        if not meta_path.is_file():
            print(f"[qd] seed skip (no metadata.json): {raw}")
            continue
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        overrides = dict(meta.get("overrides") or {})
        if overrides:
            seeds.append(overrides)
            print(f"[qd] seed loaded: {raw} ({len(overrides)} overrides)")
    return seeds


def _load_motif_if_needed(args: argparse.Namespace) -> Any:
    """Load a geometry motif when running in geometry_first mode."""
    motif_path = getattr(args, "motif_json", None)
    if not motif_path:
        return None
    from ...initial_data.motif import read_motif_json
    motif = read_motif_json(Path(motif_path).expanduser().resolve())
    print(f"[qd] geometry_first motif loaded: {motif_path}")
    return motif


def run_qd_command(args: argparse.Namespace, base_overrides: dict) -> int:
    ensure_radial_recipe_for_grtresna(args, label="qd")
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
    pipeline_cfg = pipeline_settings_from_args(args)

    archive = run_qd_search(
        runs_dir=Path(args.runs_dir).expanduser().resolve(),
        executable=executable,
        iterations=args.iterations,
        batch_size=args.batch_size,
        bins=args.bins,
        init_random=args.init_random,
        seed=args.seed,
        base_overrides=ctx.base_overrides,
        search_space=ctx.search_space,
        template=template,
        example=example,
        name=args.name,
        dry_run=args.dry_run,
        constrained=not ctx.use_grtresna,
        phantom=getattr(args, "phantom", False) or getattr(args, "phantom_default", False),
        use_preflight=(False if ctx.use_grtresna else getattr(args, "preflight", False)),
        cuda_devices=args.cuda_devices,
        gpu_ids=getattr(args, "gpu_ids", None),
        check_params=not args.skip_check_params,
        score_weights=getattr(args, "score_weights", None),
        objective_mode=getattr(args, "objective_mode", "weighted"),
        ftl_L=getattr(args, "ftl_L", None),
        consume_plotfiles=getattr(args, "consume_plotfiles", True),
        consumer_radii=getattr(args, "consumer_radii", [4.0, 8.0]),
        consumer_keep_last=getattr(args, "consumer_keep_last", 1),
        descriptor_mode=getattr(args, "descriptor_mode", "legacy"),
        grtresna=ctx.use_grtresna,
        grtresna_config=ctx.grtresna_config,
        grtresna_solved_ftl_gate=grtresna_solved_ftl_gate_enabled(
            args,
            use_grtresna=ctx.use_grtresna,
            objective_mode=getattr(args, "objective_mode", "weighted"),
        ),
        solved_ftl_gate_config=ctx.solved_ftl_gate_config,
        grtresna_convergence_config=ctx.grtresna_convergence_config,
        grtresna_postload_gate=grtresna_postload_gate_enabled(args, use_grtresna=ctx.use_grtresna),
        postload_gate_config=(
            postload_gate_config_from_args(args) if ctx.use_grtresna else None
        ),
        resume=getattr(args, "resume", False),
        target_evals=getattr(args, "target_evals", None),
        seed_overrides=load_seed_overrides(getattr(args, "seed_eval_dirs", None)),
        keep_top_eval_dirs=getattr(args, "keep_top_eval_dirs", 0),
        ftl_retention_enabled=getattr(args, "ftl_retention", True),
        slots_per_gpu=pipeline_cfg["slots_per_gpu"],
        cluster_cpu_fraction=pipeline_cfg["cluster_cpu_fraction"],
        pipeline_share=pipeline_cfg["pipeline_share"],
        max_concurrent_grtresna_override=pipeline_cfg["max_concurrent_grtresna"],
        reserve_cores=pipeline_cfg["reserve_cores"],
        grtresna_mpi_ranks=getattr(args, "grtresna_ranks", 8),
        remove_partial=pipeline_cfg["remove_partial"],
        use_pipeline=pipeline_cfg["use_pipeline"],
        pre_gpu_learning=getattr(args, "pre_gpu_learning", None),
        near_miss_pool_size=getattr(args, "near_miss_pool_size", 32),
        motif=_load_motif_if_needed(args),
        geometry_first_L=getattr(args, "ftl_L", None),
    )
    best = archive.best
    print(json.dumps({
        "num_elites": len(archive.cells),
        "coverage": archive.coverage,
        "best_score": best.score if best else None,
        "best_episode": best.episode if best else None,
    }, indent=2))
    return 0
