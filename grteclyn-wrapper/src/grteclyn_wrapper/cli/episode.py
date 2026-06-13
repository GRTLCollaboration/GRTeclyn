"""Single-episode and sweep command handlers."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any, Mapping

from ..core.config import resolve_example, resolve_executable
from ..core.episode import create_episode, update_metadata, write_json
from ..core.params import write_params
from ..core.runner import run_episode
from ..initial_data.constrained_recipe import constrained_overrides
from ..metrics import dataclass_to_dict, read_episode_metrics
from ..metrics.score import score_episode

SWEEP_RANGES = {
    "wormhole_phi_perturbation_amplitude": (-0.04, 0.04),
    "wormhole_support_strength": (0.2, 1.0),
    "wormhole_phi_perturbation_width": (0.25, 1.0),
}


def finalize_score(
    episode_dir: Path,
    target_stop_time: float | None,
    *,
    score_weights: Mapping[str, float] | None = None,
    ftl_L: float | None = None,
) -> int:
    metrics = read_episode_metrics(episode_dir, ftl_L=ftl_L)
    score = score_episode(metrics, target_stop_time=target_stop_time, weights=score_weights)
    write_json(
        episode_dir / "score.json",
        {
            "score": dataclass_to_dict(score),
            "metrics": dataclass_to_dict(metrics),
        },
    )
    print(json.dumps({"episode": str(episode_dir), "score": score.total}, indent=2))
    return 0


def run_single(args: argparse.Namespace, overrides: dict[str, Any]) -> int:
    example = resolve_example(args.example)
    phantom = getattr(args, "phantom", False) or getattr(args, "phantom_default", False)
    if getattr(args, "constrained", False) and example.name == "RadialRecipe":
        constrained_overrides(overrides, phantom=phantom)
    if phantom and example.name == "RadialRecipe":
        overrides.setdefault("recipe_exotic_matter", 1)
    runs_dir = Path(args.runs_dir).expanduser().resolve()
    metadata: dict[str, Any] = {
        "mode": args.command,
        "example": example.name,
        "overrides": overrides,
    }
    if getattr(args, "initial_data_source", None):
        metadata["initial_data_source"] = args.initial_data_source
    episode = create_episode(
        runs_dir,
        name=args.name,
        metadata=metadata,
    )
    template = Path(args.template).expanduser().resolve() if args.template else example.template
    write_params(
        template,
        episode.params_path,
        episode_dir=episode.path,
        example=example,
        overrides=overrides,
    )

    if args.dry_run:
        update_metadata(episode, {"dry_run": True})
        print(f"Wrote dry-run episode: {episode.path}")
        return 0

    executable = resolve_executable(
        args.executable,
        example=example,
        mpi_ranks=args.mpi_ranks,
        comp=args.comp,
        cuda=not args.no_cuda,
        debug=args.debug,
    )
    result = run_episode(
        episode,
        executable,
        check_params=not args.skip_check_params,
        cuda_devices=args.cuda_devices,
        consume_plotfiles=args.consume_plotfiles,
        consumer_radii=args.consumer_radii,
        consumer_delete=args.consumer_delete,
        consumer_keep_last=getattr(args, "consumer_keep_last", 1),
    )
    finalize_score(
        episode.path,
        target_stop_time=overrides.get("stop_time"),
        score_weights=getattr(args, "score_weights", None),
        ftl_L=getattr(args, "ftl_L", None),
    )
    return result.returncode


def _sample_overrides(base: dict[str, Any], rng: random.Random) -> dict[str, Any]:
    overrides = dict(base)
    for key, (lo, hi) in SWEEP_RANGES.items():
        overrides.setdefault(key, rng.uniform(lo, hi))
    return overrides


def run_sweep(args: argparse.Namespace, base_overrides: dict[str, Any]) -> int:
    rng = random.Random(args.seed)
    status = 0
    for index in range(args.count):
        overrides = _sample_overrides(base_overrides, rng)
        name = args.name or f"sweep_{index + 1:06d}"
        per_args = argparse.Namespace(**vars(args))
        per_args.name = name
        per_args.command = "sweep"
        rc = run_single(per_args, overrides)
        status = status or rc
        if rc != 0 and args.stop_on_failure:
            break
    return status
