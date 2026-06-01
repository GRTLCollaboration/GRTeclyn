from __future__ import annotations

import csv
import json
import random
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..core.config import DEFAULT_TEMPLATE, ExecutableConfig, ExampleConfig, resolve_example
from ..core.episode import Episode, create_episode, update_metadata, write_json
from ..core.params import write_params
from ..core.runner import run_episode
from ..initial_data.constrained_recipe import constrained_overrides
from ..initial_data.preflight import PreflightResult, preflight_check
from ..metrics.episode_metrics import EpisodeMetrics, dataclass_to_dict, read_episode_metrics
from ..metrics.score import Score, score_episode


DEFAULT_RANGES: dict[str, tuple[float, float]] = {
    "wormhole_phi_perturbation_amplitude": (-0.04, 0.04),
    "wormhole_support_strength": (0.2, 1.0),
    "wormhole_phi_perturbation_width": (0.25, 1.0),
    "wormhole_phi_monopole_amplitude": (-0.02, 0.02),
}

DEFAULT_RECIPE_RANGES: dict[str, tuple[float, float]] = {
    "recipe_chi_coeff_0": (-0.5, 0.1),
    "recipe_chi_coeff_1": (-0.2, 0.2),
    "recipe_phi_coeff_0": (-0.1, 0.1),
    "recipe_basis_width": (0.5, 2.0),
}


def atlas_ranges_for_example(example: ExampleConfig) -> dict[str, tuple[float, float]]:
    if example.name == "RadialRecipe":
        return DEFAULT_RECIPE_RANGES
    return DEFAULT_RANGES

DEFAULT_LOW_RES_OVERRIDES: dict[str, Any] = {
    "N_full": 32,
    "max_level": 0,
    "stop_time": 0.04,
    "plot_interval": 1000,
    "checkpoint_interval": 1000,
}

CSV_FIELDS = [
    "index",
    "episode",
    "exit_code",
    "score_total",
    "labels",
    "final_time",
    "target_stop_time",
    "min_lapse",
    "min_chi",
    "max_abs_k",
    "max_horizon_radius",
    "min_theta_plus",
    "max_hamiltonian_l2",
    "max_momentum_l2",
    "min_rho_required",
    "integral_negative_rho",
    "overrides_json",
]


@dataclass(frozen=True)
class AtlasThresholds:
    constraint_blowup: float = 1.0e-1
    lapse_collapse: float = 1.0e-3
    horizon_radius: float = 1.0e-8
    trivial_geometry: float = 1.0e-3


@dataclass(frozen=True)
class AtlasPaths:
    root: Path
    jsonl: Path
    csv: Path
    summary: Path


def make_atlas_dir(runs_dir: Path, name: str | None = None) -> AtlasPaths:
    if name is None:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"atlas_{timestamp}"
    root = (runs_dir / name).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=False)
    return AtlasPaths(
        root=root,
        jsonl=root / "atlas.jsonl",
        csv=root / "atlas.csv",
        summary=root / "summary.json",
    )


def sample_overrides(
    rng: random.Random,
    *,
    base: Mapping[str, Any] | None = None,
    ranges: Mapping[str, tuple[float, float]] | None = None,
) -> dict[str, Any]:
    overrides = dict(DEFAULT_LOW_RES_OVERRIDES)
    if base:
        overrides.update(base)

    for key, (lo, hi) in (ranges or DEFAULT_RANGES).items():
        overrides.setdefault(key, rng.uniform(lo, hi))
    return overrides


def classify_episode(
    metrics: EpisodeMetrics,
    score: Score,
    *,
    exit_code: int | None,
    target_stop_time: float | None,
    thresholds: AtlasThresholds = AtlasThresholds(),
) -> list[str]:
    labels: list[str] = []

    if exit_code is not None and exit_code != 0:
        labels.append("solver_failed")

    if metrics.collapse is None and metrics.constraints is None:
        labels.append("missing_diagnostics")
        return labels

    constraints = metrics.constraints
    if constraints is not None:
        max_constraint = max(
            constraints.max_hamiltonian_l2 or 0.0,
            constraints.max_momentum_l2 or 0.0,
        )
        if max_constraint >= thresholds.constraint_blowup:
            labels.append("constraint_blowup")

    collapse = metrics.collapse
    if collapse is not None:
        if collapse.min_lapse is not None and collapse.min_lapse <= thresholds.lapse_collapse:
            labels.append("lapse_collapse")
        if collapse.max_horizon_radius is not None and collapse.max_horizon_radius > thresholds.horizon_radius:
            labels.append("horizon_formed")
        if abs(score.components.get("nontrivial_geometry", 0.0)) <= thresholds.trivial_geometry:
            labels.append("trivial_geometry")

    if (
        exit_code == 0
        and "missing_diagnostics" not in labels
        and (target_stop_time is None or score.components.get("survival", 0.0) >= 1.0)
    ):
        labels.append("completed")

    return labels or ["completed"]


def build_record(
    *,
    index: int,
    episode: Episode,
    overrides: Mapping[str, Any],
    exit_code: int | None,
    metrics: EpisodeMetrics,
    score: Score,
    labels: Sequence[str],
    target_stop_time: float | None,
) -> dict[str, Any]:
    collapse = metrics.collapse
    constraints = metrics.constraints
    final_time = None
    if collapse and collapse.final_time is not None:
        final_time = collapse.final_time
    elif constraints and constraints.final_time is not None:
        final_time = constraints.final_time

    return {
        "index": index,
        "episode": str(episode.path),
        "exit_code": exit_code,
        "overrides": dict(overrides),
        "metrics": dataclass_to_dict(metrics),
        "score": dataclass_to_dict(score),
        "labels": list(labels),
        "final_time": final_time,
        "target_stop_time": target_stop_time,
    }


def append_jsonl(path: Path, record: Mapping[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def append_csv(path: Path, record: Mapping[str, Any]) -> None:
    row = flatten_record(record)
    should_write_header = not path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        if should_write_header:
            writer.writeheader()
        writer.writerow(row)


def flatten_record(record: Mapping[str, Any]) -> dict[str, Any]:
    metrics = record.get("metrics") or {}
    collapse = metrics.get("collapse") or {}
    constraints = metrics.get("constraints") or {}
    score = record.get("score") or {}
    return {
        "index": record.get("index"),
        "episode": record.get("episode"),
        "exit_code": record.get("exit_code"),
        "score_total": score.get("total"),
        "labels": ",".join(record.get("labels") or []),
        "final_time": record.get("final_time"),
        "target_stop_time": record.get("target_stop_time"),
        "min_lapse": collapse.get("min_lapse"),
        "min_chi": collapse.get("min_chi"),
        "max_abs_k": collapse.get("max_abs_k"),
        "max_horizon_radius": collapse.get("max_horizon_radius"),
        "min_theta_plus": collapse.get("min_theta_plus"),
        "max_hamiltonian_l2": constraints.get("max_hamiltonian_l2"),
        "max_momentum_l2": constraints.get("max_momentum_l2"),
        "min_rho_required": constraints.get("min_rho_required"),
        "integral_negative_rho": constraints.get("integral_negative_rho"),
        "overrides_json": json.dumps(record.get("overrides") or {}, sort_keys=True),
    }


def summarize_records(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    label_counts: Counter[str] = Counter()
    best_record = None
    for record in records:
        label_counts.update(record.get("labels") or [])
        if best_record is None:
            best_record = record
            continue
        score = (record.get("score") or {}).get("total")
        best_score = (best_record.get("score") or {}).get("total")
        if score is not None and (best_score is None or score > best_score):
            best_record = record

    return {
        "count": len(records),
        "label_counts": dict(sorted(label_counts.items())),
        "best": {
            "episode": best_record.get("episode") if best_record else None,
            "score_total": (best_record.get("score") or {}).get("total") if best_record else None,
            "labels": best_record.get("labels") if best_record else [],
        },
    }


def write_score(episode: Episode, metrics: EpisodeMetrics, score: Score, labels: Sequence[str]) -> None:
    write_json(
        episode.score_path,
        {
            "labels": list(labels),
            "score": dataclass_to_dict(score),
            "metrics": dataclass_to_dict(metrics),
        },
    )


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
