"""Single-candidate and parallel GPU evaluation for CMA-ES."""

from __future__ import annotations

import threading
from pathlib import Path
from typing import Any, Mapping, Sequence

from ...core.config import ExecutableConfig, ExampleConfig
from ...core.episode import create_episode, update_metadata, write_json
from ...core.params import write_params
from ...core.runner import run_episode
from ...initial_data.constrained_recipe import constrained_overrides
from ...initial_data.preflight import preflight_check
from ...metrics import dataclass_to_dict, read_episode_metrics
from ...metrics.score import domain_half_width_for_episode, score_episode
from ..grtresna_convergence_gate import (
    GRTRESNA_REJECTION_BASE_FITNESS,
    GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
    GRTresnaConvergenceConfig,
    convergence_rejection_fitness as _grtresna_rejection_fitness,
    convergence_rejection_reason as _grtresna_convergence_rejection_reason,
)
from ..trajectory_log import format_eval_log_line, infer_trajectory_status
from .candidates import _vector_to_overrides
from .config import build_grtresna_config, parse_convergence_safe
from .dimension import SearchDimension

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]

ObjectiveMode = str

SOLVED_FTL_REJECTION_BASE_FITNESS = 75.0
SOLVED_FTL_REJECTION_MAX_EXTRA_FITNESS = 20.0

def _assign_gpu(index: int, gpu_ids: Sequence[int]) -> str:
    """Round-robin GPU assignment for parallel evaluation."""
    return str(gpu_ids[index % len(gpu_ids)])


def _collect_training(
    trajectory: Sequence[Mapping[str, Any]],
    dims: Sequence[SearchDimension],
):
    """Build (X, y) arrays of evaluated candidates for surrogate fitting."""
    xs: list[list[float]] = []
    ys: list[float] = []
    for rec in trajectory:
        if rec.get("preflight_rejected") or rec.get("surrogate_predicted"):
            continue
        overrides = rec.get("overrides")
        score = rec.get("score")
        if not isinstance(overrides, dict) or score is None:
            continue
        try:
            xs.append([float(overrides[d.param_key]) for d in dims])
            ys.append(float(score))
        except (KeyError, TypeError, ValueError):
            continue
    if not xs:
        return None, None
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def _track_trajectory(
    trajectory: list[dict[str, Any]],
    record: dict[str, Any],
) -> None:
    record.setdefault("status", infer_trajectory_status(record))
    trajectory.append(record)
    print(format_eval_log_line(record), flush=True)


def _objective(
    x: Sequence[float],
    *,
    dims: Sequence[SearchDimension],
    base_overrides: Mapping[str, Any],
    opt_dir: Path,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    eval_counter: list[int],
    constrained: bool,
    phantom: bool,
    use_preflight: bool,
    cuda_devices: str | None,
    check_params: bool,
    dry_run: bool,
    trajectory: list[dict[str, Any]],
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    objective_mode: ObjectiveMode = "weighted",
    ftl_L: float | None = None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
    consumer_keep_last: int = 1,
    grtresna: bool = False,
    grtresna_base: "GRTresnaConfig | None" = None,
    grtresna_solved_ftl_gate: bool = False,
    solved_ftl_gate_config: Any | None = None,
    grtresna_convergence_config: GRTresnaConvergenceConfig | None = None,
    grtresna_postload_gate: bool = False,
    postload_gate_config: Any | None = None,
) -> float:
    """Evaluate one candidate.  Returns negative score (CMA-ES minimizes)."""
    eval_counter[0] += 1
    idx = eval_counter[0]

    overrides = _vector_to_overrides(x, dims, base_overrides)

    # In GRTresna mode the constraint solve replaces the 1D radial recipe, so
    # the recipe-specific constrained/preflight steps are skipped.
    if constrained and not grtresna:
        constrained_overrides(overrides, phantom=phantom)

    if use_preflight and not grtresna:
        pf = preflight_check(overrides, phantom=phantom)
        if not pf.passed:
            record = {
                "eval": idx,
                "score": 0.0,
                "preflight_rejected": True,
                "reason": pf.reason,
                "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
            }
            _track_trajectory(trajectory, record)
            return 100.0

    episode = create_episode(
        opt_dir,
        name=f"eval_{idx:06d}",
        metadata={
            "mode": "optimize",
            "example": example.name,
            "eval_index": idx,
            "overrides": overrides,
        },
    )

    # GRTeclyn params get everything except the grtresna_* search keys (which
    # drive the upstream solver, not GRTeclyn's params.txt).
    gte_overrides = {
        k: v for k, v in overrides.items() if not str(k).startswith("grtresna_")
    }

    if grtresna and not dry_run:
        from ...grtresna.matter_wiring import merge_evolution_overrides
        from ...grtresna.solver import solve as grtresna_solve
        from ..grtresna_evaluation_gates import (
            GRTresnaPreEvolutionGateConfig,
            apply_grtresna_pre_evolution_gates,
        )

        cfg = build_grtresna_config(overrides, grtresna_base)
        try:
            gridinit = grtresna_solve(
                cfg,
                work_dir=episode.path / "grtresna",
                gridinit_path=episode.path / "initial_data.gridinit",
            )
            gte_overrides["recipe_initial_data_file"] = str(gridinit)
            gte_overrides.update(merge_evolution_overrides(gte_overrides, cfg))
            convergence = parse_convergence_safe(episode.path / "grtresna")
            update_metadata(episode, {"grtresna_convergence": convergence})

            gate_rejection = apply_grtresna_pre_evolution_gates(
                episode=episode,
                convergence=convergence,
                gridinit_path=Path(gridinit),
                gte_overrides=gte_overrides,
                cuda_devices=cuda_devices,
                config=GRTresnaPreEvolutionGateConfig(
                    convergence_config=grtresna_convergence_config,
                    solved_ftl_enabled=grtresna_solved_ftl_gate,
                    solved_ftl_config=solved_ftl_gate_config,
                    postload_enabled=grtresna_postload_gate,
                    postload_config=postload_gate_config,
                    ftl_L=ftl_L,
                ),
            )
            if gate_rejection is not None:
                update_metadata(
                    episode,
                    {
                        **gate_rejection.metadata,
                        "simulation_exit_code": None,
                    },
                )
                record = {
                    "eval": idx,
                    "episode": str(episode.path),
                    "score": -gate_rejection.fitness,
                    "fitness": gate_rejection.fitness,
                    "reason": gate_rejection.reason,
                    "components": {
                        gate_rejection.component_key: -gate_rejection.fitness
                    },
                    "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
                    **gate_rejection.metadata,
                }
                if gate_rejection.metrics:
                    record.update(gate_rejection.metrics)
                _track_trajectory(trajectory, record)
                return gate_rejection.fitness
        except Exception as exc:  # solver failure -> penalise, keep searching
            fitness = GRTRESNA_REJECTION_BASE_FITNESS + GRTRESNA_REJECTION_MAX_EXTRA_FITNESS
            update_metadata(episode, {"grtresna_error": repr(exc)})
            record = {
                "eval": idx,
                "episode": str(episode.path),
                "score": -fitness,
                "fitness": fitness,
                "grtresna_failed": True,
                "reason": repr(exc),
                "components": {"grtresna_rejection": -fitness},
                "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
            }
            _track_trajectory(trajectory, record)
            return fitness

    write_params(
        template, episode.params_path,
        episode_dir=episode.path, example=example, overrides=gte_overrides,
    )

    exit_code: int | None = None
    if dry_run:
        update_metadata(episode, {"dry_run": True})
    else:
        if executable is None:
            raise ValueError("executable is required unless dry_run=True")
        try:
            result = run_episode(
                episode, executable,
                check_params=check_params,
                cuda_devices=cuda_devices,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_delete=True,
                consumer_keep_last=consumer_keep_last,
            )
            exit_code = result.returncode
        except Exception as exc:
            exit_code = 1
            update_metadata(episode, {
                "simulation_error": repr(exc),
                "simulation_exit_code": exit_code,
            })

    metrics = read_episode_metrics(episode.path, ftl_L=ftl_L)
    score = score_episode(
        metrics,
        target_stop_time=target_stop_time,
        weights=score_weights,
        objective_mode=objective_mode,
        domain_half_width=domain_half_width_for_episode(episode.path, gte_overrides),
    )

    write_json(episode.score_path, {
        "score": dataclass_to_dict(score),
        "metrics": dataclass_to_dict(metrics),
    })

    record = {
        "eval": idx,
        "episode": str(episode.path),
        "exit_code": exit_code,
        "score": score.total,
        "components": score.components,
        "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
    }
    _track_trajectory(trajectory, record)
    return -score.total
def _evaluate_generation_parallel(
    solutions: list,
    *,
    dims: Sequence[SearchDimension],
    base_overrides: Mapping[str, Any],
    opt_dir: Path,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    eval_counter: list[int],
    constrained: bool,
    phantom: bool,
    use_preflight: bool,
    gpu_ids: Sequence[int],
    check_params: bool,
    trajectory: list[dict[str, Any]],
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    objective_mode: ObjectiveMode,
    ftl_L: float | None,
    consume_plotfiles: bool,
    consumer_radii: Sequence[float],
    consumer_keep_last: int = 1,
    grtresna: bool = False,
    grtresna_config: "GRTresnaConfig | None" = None,
    grtresna_solved_ftl_gate: bool = False,
    solved_ftl_gate_config: Any | None = None,
    grtresna_convergence_config: GRTresnaConvergenceConfig | None = None,
    grtresna_postload_gate: bool = False,
    postload_gate_config: Any | None = None,
) -> list[float]:
    """Evaluate an entire CMA-ES generation in parallel across GPUs.

    Each solution is assigned a GPU round-robin, and all are launched
    concurrently via threads (GIL released during subprocess waits).
    """
    import threading

    fitnesses: list[float | None] = [None] * len(solutions)
    lock = threading.Lock()

    def _eval_one(idx_in_gen: int, sol) -> None:
        gpu = _assign_gpu(idx_in_gen, gpu_ids)
        f = _objective(
            sol,
            dims=dims,
            base_overrides=base_overrides,
            opt_dir=opt_dir,
            example=example,
            template=template,
            executable=executable,
            eval_counter=eval_counter,
            constrained=constrained,
            phantom=phantom,
            use_preflight=use_preflight,
            cuda_devices=gpu,
            check_params=check_params,
            dry_run=False,
            trajectory=trajectory,
            target_stop_time=target_stop_time,
            score_weights=score_weights,
            objective_mode=objective_mode,
            ftl_L=ftl_L,
            consume_plotfiles=consume_plotfiles,
            consumer_radii=consumer_radii,
            consumer_keep_last=consumer_keep_last,
            grtresna=grtresna,
            grtresna_base=grtresna_config,
            grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
            solved_ftl_gate_config=solved_ftl_gate_config,
            grtresna_convergence_config=grtresna_convergence_config,
            grtresna_postload_gate=grtresna_postload_gate,
            postload_gate_config=postload_gate_config,
        )
        with lock:
            fitnesses[idx_in_gen] = f

    threads = []
    for i, sol in enumerate(solutions):
        t = threading.Thread(target=_eval_one, args=(i, sol))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()

    return [f if f is not None else 100.0 for f in fitnesses]
