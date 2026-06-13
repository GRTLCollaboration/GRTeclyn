"""Random sampling driver for the Spacetime Failure Atlas."""

from __future__ import annotations

import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from ...core.config import DEFAULT_TEMPLATE, ExecutableConfig, ExampleConfig, resolve_example
from ...core.episode import create_episode, update_metadata, write_json
from ...core.params import write_params
from ...core.runner import run_episode
from ...initial_data.constrained_recipe import constrained_overrides
from ...initial_data.preflight import preflight_check
from ...metrics import read_episode_metrics
from ...metrics.score import score_episode
from .config import (
    DEFAULT_LOW_RES_OVERRIDES,
    AtlasPaths,
    atlas_ranges_for_example,
    make_atlas_dir,
    sample_overrides,
)
from .io import append_csv, append_jsonl
from .records import build_record, classify_episode, summarize_records, write_score


def run_atlas(
    *,
    runs_dir: Path,
    executable: ExecutableConfig | None,
    count: int,
    seed: int | None,
    base_overrides: Mapping[str, Any] | None = None,
    template: Path = DEFAULT_TEMPLATE,
    example: ExampleConfig | str = "SupportedWormholeCollapse",
    name: str | None = None,
    dry_run: bool = False,
    stop_on_failure: bool = False,
    cuda_devices: str | None = None,
    check_params: bool = True,
    constrained: bool = False,
    phantom: bool = False,
    preflight: bool = False,
) -> tuple[AtlasPaths, list[dict[str, Any]], dict[str, Any]]:
    example_cfg = example if isinstance(example, ExampleConfig) else resolve_example(example)
    ranges = atlas_ranges_for_example(example_cfg)
    paths = make_atlas_dir(runs_dir, name=name)
    rng = random.Random(seed)
    records: list[dict[str, Any]] = []

    write_json(
        paths.root / "metadata.json",
        {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "count": count,
            "seed": seed,
            "example": example_cfg.name,
            "dry_run": dry_run,
            "base_overrides": dict(base_overrides or {}),
            "ranges": ranges,
            "low_res_defaults": DEFAULT_LOW_RES_OVERRIDES,
        },
    )

    for index in range(1, count + 1):
        overrides = sample_overrides(rng, base=base_overrides, ranges=ranges)
        if constrained and example_cfg.name == "RadialRecipe":
            constrained_overrides(overrides, phantom=phantom)

        if preflight and example_cfg.name == "RadialRecipe":
            pf = preflight_check(overrides, phantom=phantom)
            if not pf.passed:
                target_stop_time = float(overrides["stop_time"]) if "stop_time" in overrides else None
                episode = create_episode(
                    paths.root,
                    name=f"episode_{index:06d}",
                    metadata={
                        "mode": "atlas", "example": example_cfg.name,
                        "atlas_index": index, "overrides": overrides,
                        "preflight_rejected": True, "preflight_reason": pf.reason,
                    },
                )
                write_params(template, episode.params_path, episode_dir=episode.path, example=example_cfg, overrides=overrides)
                metrics = read_episode_metrics(episode.path)
                score = score_episode(metrics, target_stop_time=target_stop_time)
                labels = ["preflight_rejected"]
                write_score(episode, metrics, score, labels)
                record = build_record(
                    index=index, episode=episode, overrides=overrides,
                    exit_code=None, metrics=metrics, score=score,
                    labels=labels, target_stop_time=target_stop_time,
                )
                append_jsonl(paths.jsonl, record)
                append_csv(paths.csv, record)
                records.append(record)
                continue

        target_stop_time = float(overrides["stop_time"]) if "stop_time" in overrides else None
        episode = create_episode(
            paths.root,
            name=f"episode_{index:06d}",
            metadata={"mode": "atlas", "example": example_cfg.name, "atlas_index": index, "overrides": overrides},
        )
        write_params(template, episode.params_path, episode_dir=episode.path, example=example_cfg, overrides=overrides)

        exit_code: int | None = None
        if dry_run:
            update_metadata(episode, {"dry_run": True})
            metrics = read_episode_metrics(episode.path)
        else:
            if executable is None:
                raise ValueError("executable is required unless dry_run=True")
            try:
                result = run_episode(
                    episode,
                    executable,
                    check_params=check_params,
                    cuda_devices=cuda_devices,
                )
                exit_code = result.returncode
            except Exception as exc:
                exit_code = 1
                update_metadata(episode, {"simulation_error": repr(exc), "simulation_exit_code": exit_code})
            metrics = read_episode_metrics(episode.path)

        score = score_episode(metrics, target_stop_time=target_stop_time)
        labels = classify_episode(metrics, score, exit_code=exit_code, target_stop_time=target_stop_time)
        write_score(episode, metrics, score, labels)

        record = build_record(
            index=index,
            episode=episode,
            overrides=overrides,
            exit_code=exit_code,
            metrics=metrics,
            score=score,
            labels=labels,
            target_stop_time=target_stop_time,
        )
        append_jsonl(paths.jsonl, record)
        append_csv(paths.csv, record)
        records.append(record)

        if exit_code not in (None, 0) and stop_on_failure:
            break

    summary = summarize_records(records)
    write_json(paths.summary, summary)
    return paths, records, summary
