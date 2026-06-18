"""Single-candidate and parallel GPU evaluation for CMA-ES."""

from __future__ import annotations

import threading
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping, Sequence

from ...core.config import ExecutableConfig, ExampleConfig
from ...core.episode import create_episode, update_metadata, write_json
from ...core.evaluation import Evaluation, _resolve_consumer_evolving_geodesic
from ...core.runner import run_episode
from ...initial_data.constrained_recipe import constrained_overrides
from ...initial_data.preflight import preflight_check
from ...metrics import dataclass_to_dict, read_episode_metrics
from ...metrics.score import domain_half_width_for_episode, score_episode
from ..eval_pipeline import EvalJob, EvalPipeline
from ..grtresna_convergence_gate import (
    GRTresnaConvergenceConfig,
    convergence_rejection_fitness as _grtresna_rejection_fitness,
    convergence_rejection_reason as _grtresna_convergence_rejection_reason,
)
from ..ftl_peak_metrics import peak_fields_for_descriptor_details
from ..trajectory_log import format_eval_log_line, format_trajectory_line, infer_trajectory_status, trajectory_flags_from_evaluation
from .candidates import _vector_to_overrides
from .dimension import SearchDimension

if TYPE_CHECKING:
    from ..gpu_pool import CpuAdmissionController, GpuPool

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
    *,
    trajectory_lock: threading.Lock | None = None,
    live_trajectory_path: Path | None = None,
) -> None:
    record.setdefault("status", infer_trajectory_status(record))
    if trajectory_lock is not None:
        with trajectory_lock:
            trajectory.append(record)
            if live_trajectory_path is not None:
                with live_trajectory_path.open("a", encoding="utf-8") as fh:
                    fh.write(format_trajectory_line(record))
    else:
        trajectory.append(record)
        if live_trajectory_path is not None:
            with live_trajectory_path.open("a", encoding="utf-8") as fh:
                fh.write(format_trajectory_line(record))
    print(format_eval_log_line(record), flush=True)


def _fitness_from_evaluation(res: Evaluation) -> float:
    if res.preflight_rejected:
        return 100.0
    return -res.score


def _trajectory_record_from_evaluation(
    res: Evaluation,
    *,
    idx: int,
    dims: Sequence[SearchDimension],
    overrides: Mapping[str, Any],
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "eval": idx,
        "score": res.score,
        "components": dict(res.components),
        "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
        **trajectory_flags_from_evaluation(res),
    }
    if res.episode_path:
        record["episode"] = res.episode_path
    if res.metrics:
        record.update(res.metrics)
    if not res.preflight_rejected and res.reason is None:
        record["descriptor_details"] = peak_fields_for_descriptor_details(
            res.metrics.get("ftl_timeseries") if res.metrics else None,
            components=res.components,
            metrics=res.metrics,
        )
    record["status"] = infer_trajectory_status(record)
    return record


def _run_grtresna_evaluation(
    *,
    overrides: Mapping[str, Any],
    opt_dir: Path,
    name: str,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    check_params: bool,
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    objective_mode: ObjectiveMode,
    ftl_L: float | None,
    consume_plotfiles: bool,
    consumer_radii: Sequence[float],
    consumer_keep_last: int,
    grtresna_base: "GRTresnaConfig | None",
    grtresna_solved_ftl_gate: bool,
    solved_ftl_gate_config: Any | None,
    grtresna_convergence_config: GRTresnaConvergenceConfig | None,
    grtresna_postload_gate: bool,
    postload_gate_config: Any | None,
    cuda_devices: str | None,
    cpu_admission: "CpuAdmissionController | None",
    gpu_pool: "GpuPool | None",
) -> Evaluation:
    """GRTresna eval via the same CPU/GPU session path as MAP-Elites QD."""
    from ...core.evaluation import _run_cpu_grtresna_gates, _run_gpu_session

    def _cpu_phase() -> Evaluation | Any:
        if cpu_admission is not None:
            with cpu_admission.admit():
                return _run_cpu_grtresna_gates(
                    overrides=overrides,
                    out_dir=opt_dir,
                    name=name,
                    example=example,
                    grtresna_base=grtresna_base,
                    grtresna_convergence_config=grtresna_convergence_config,
                    grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
                    solved_ftl_gate_config=solved_ftl_gate_config,
                    ftl_L=ftl_L,
                    dry_run=False,
                )
        return _run_cpu_grtresna_gates(
            overrides=overrides,
            out_dir=opt_dir,
            name=name,
            example=example,
            grtresna_base=grtresna_base,
            grtresna_convergence_config=grtresna_convergence_config,
            grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
            solved_ftl_gate_config=solved_ftl_gate_config,
            ftl_L=ftl_L,
            dry_run=False,
        )

    def _gpu_session(cpu_result, leased_gpu: str | None) -> Evaluation:
        return _run_gpu_session(
            cpu_result=cpu_result,
            example=example,
            template=template,
            executable=executable,
            cuda_devices=leased_gpu,
            check_params=check_params,
            dry_run=False,
            target_stop_time=target_stop_time,
            score_weights=score_weights,
            objective_mode=objective_mode,
            ftl_L=ftl_L,
            consume_plotfiles=consume_plotfiles,
            consumer_radii=consumer_radii,
            consumer_keep_last=consumer_keep_last,
            consumer_ftl_timeseries=True,
            consumer_evolving_geodesic=None,
            grtresna_postload_gate=grtresna_postload_gate,
            postload_gate_config=postload_gate_config,
        )

    cpu_phase = _cpu_phase()
    if isinstance(cpu_phase, Evaluation):
        return cpu_phase

    if gpu_pool is not None:
        with gpu_pool.lease() as leased_gpu:
            return _gpu_session(cpu_phase, leased_gpu)
    return _gpu_session(cpu_phase, cuda_devices)


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
    counter_lock: threading.Lock | None = None,
    trajectory_lock: threading.Lock | None = None,
    live_trajectory_path: Path | None = None,
    eval_id: int | None = None,
    gpu_pool: "GpuPool | None" = None,
    cpu_admission: "CpuAdmissionController | None" = None,
    postload_cuda_devices: str | None = None,
) -> float:
    """Evaluate one candidate.  Returns negative score (CMA-ES minimizes)."""
    if eval_id is not None:
        idx = eval_id
    elif counter_lock is not None:
        with counter_lock:
            eval_counter[0] += 1
            idx = eval_counter[0]
    else:
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
            _track_trajectory(
                trajectory, record,
                trajectory_lock=trajectory_lock,
                live_trajectory_path=live_trajectory_path,
            )
            return 100.0

    if grtresna and not dry_run:
        res = _run_grtresna_evaluation(
            overrides=overrides,
            opt_dir=opt_dir,
            name=f"eval_{idx:06d}",
            example=example,
            template=template,
            executable=executable,
            check_params=check_params,
            target_stop_time=target_stop_time,
            score_weights=score_weights,
            objective_mode=objective_mode,
            ftl_L=ftl_L,
            consume_plotfiles=consume_plotfiles,
            consumer_radii=consumer_radii,
            consumer_keep_last=consumer_keep_last,
            grtresna_base=grtresna_base,
            grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
            solved_ftl_gate_config=solved_ftl_gate_config,
            grtresna_convergence_config=grtresna_convergence_config,
            grtresna_postload_gate=grtresna_postload_gate,
            postload_gate_config=postload_gate_config,
            cuda_devices=cuda_devices,
            cpu_admission=cpu_admission,
            gpu_pool=gpu_pool,
        )
        record = _trajectory_record_from_evaluation(
            res, idx=idx, dims=dims, overrides=overrides,
        )
        _track_trajectory(
            trajectory, record,
            trajectory_lock=trajectory_lock,
            live_trajectory_path=live_trajectory_path,
        )
        return _fitness_from_evaluation(res)

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

    from ...core.params import write_params

    write_params(
        template, episode.params_path,
        episode_dir=episode.path, example=example, overrides=gte_overrides,
    )

    exit_code: int | None = None
    evolving_enabled = _resolve_consumer_evolving_geodesic(None)
    if dry_run:
        update_metadata(episode, {"dry_run": True})
    else:
        if executable is None:
            raise ValueError("executable is required unless dry_run=True")

        def _run_gpu_evolution(leased_gpu: str | None) -> None:
            nonlocal exit_code
            try:
                extra_env = {"GRTECLYN_EVOLVING_GEODESIC": "1"} if evolving_enabled else None
                result = run_episode(
                    episode, executable,
                    check_params=check_params,
                    cuda_devices=leased_gpu if leased_gpu is not None else cuda_devices,
                    extra_env=extra_env,
                    consume_plotfiles=consume_plotfiles,
                    consumer_radii=consumer_radii,
                    consumer_delete=True,
                    consumer_keep_last=consumer_keep_last,
                    consumer_ftl_timeseries=consume_plotfiles,
                    consumer_ftl_L=ftl_L,
                    consumer_incremental_score=consume_plotfiles,
                    consumer_objective_mode=objective_mode,
                    consumer_target_stop_time=target_stop_time,
                    consumer_score_weights=score_weights,
                    consumer_evolving_geodesic=evolving_enabled,
                )
                exit_code = result.returncode
            except Exception as exc:
                exit_code = 1
                update_metadata(episode, {
                    "simulation_error": repr(exc),
                    "simulation_exit_code": exit_code,
                })

        try:
            if gpu_pool is not None:
                with gpu_pool.lease() as leased_gpu:
                    _run_gpu_evolution(leased_gpu)
            else:
                _run_gpu_evolution(cuda_devices)
        except Exception as exc:
            exit_code = 1
            update_metadata(episode, {
                "simulation_error": repr(exc),
                "simulation_exit_code": exit_code,
            })

    metrics = read_episode_metrics(
        episode.path,
        ftl_L=ftl_L,
        evolving_geodesic=evolving_enabled,
        objective_mode=objective_mode,
    )
    score = score_episode(
        metrics,
        target_stop_time=target_stop_time,
        weights=score_weights,
        objective_mode=objective_mode,
        domain_half_width=domain_half_width_for_episode(episode.path, gte_overrides),
    )

    metrics_dict = dataclass_to_dict(metrics)
    write_json(episode.score_path, {
        "score": dataclass_to_dict(score),
        "metrics": metrics_dict,
    })

    # Time-resolved FTL peak fields (f_geo/f_op/speed/superluminal/lifetime).
    # Mirrors the QD trajectory enrichment so the CMA-ES path can rank pruned
    # evals and feed the FTL champion board (ftl_retention).
    record = {
        "eval": idx,
        "episode": str(episode.path),
        "exit_code": exit_code,
        "score": score.total,
        "components": score.components,
        "descriptor_details": peak_fields_for_descriptor_details(
            metrics_dict.get("ftl_timeseries") if isinstance(metrics_dict, dict) else None,
            components=score.components,
            metrics=metrics_dict if isinstance(metrics_dict, dict) else None,
        ),
        "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
    }
    _track_trajectory(
        trajectory, record,
        trajectory_lock=trajectory_lock,
        live_trajectory_path=live_trajectory_path,
    )
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
    live_trajectory_path: Path | None = None,
) -> list[float]:
    """Evaluate an entire CMA-ES generation in parallel across GPUs.

    Each solution is assigned a GPU round-robin, and all are launched
    concurrently via threads (GIL released during subprocess waits).
    """
    import threading

    fitnesses: list[float | None] = [None] * len(solutions)
    lock = threading.Lock()
    counter_lock = threading.Lock()
    trajectory_lock = threading.Lock()

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
            counter_lock=counter_lock,
            trajectory_lock=trajectory_lock,
            live_trajectory_path=live_trajectory_path,
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


def _evaluate_generation_pipelined(
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
    gpu_pool: "GpuPool",
    cpu_admission: "CpuAdmissionController",
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
    metadata_path: Path | None = None,
    live_trajectory_path: Path | None = None,
) -> list[float]:
    """Evaluate a CMA-ES generation via EvalPipeline (GpuPool + CPU admission)."""
    fitnesses: list[float | None] = [None] * len(solutions)
    counter_lock = threading.Lock()
    trajectory_lock = threading.Lock()
    postload_cuda = str(gpu_ids[0]) if gpu_ids else None
    common = dict(
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
        counter_lock=counter_lock,
        trajectory_lock=trajectory_lock,
        live_trajectory_path=live_trajectory_path,
        gpu_pool=gpu_pool,
        cpu_admission=cpu_admission,
        postload_cuda_devices=postload_cuda,
    )

    def evaluate_fn(job: EvalJob[tuple[int, list]]) -> Evaluation:
        gen_idx, sol = job.payload
        fitness = _objective(sol, eval_id=job.eval_id, cuda_devices=None, **common)
        fitnesses[gen_idx] = fitness
        return Evaluation(
            score=-float(fitness),
            components={},
            notes=[],
            episode_path=None,
            exit_code=0,
            preflight_rejected=False,
            reason=None,
            metrics={},
        )

    def on_complete(_job: EvalJob[tuple[int, list]], _res: Evaluation) -> None:
        return

    pipeline = EvalPipeline(
        gpu_pool=gpu_pool,
        cpu_admission=cpu_admission,
        evaluate_fn=evaluate_fn,
        on_complete=on_complete,
        counter_lock=counter_lock,
        eval_counter=eval_counter,
        metadata_path=metadata_path,
    )
    for i, sol in enumerate(solutions):
        pipeline.submit((i, list(sol)))
    pipeline.wait_until_idle()
    pipeline.shutdown()
    return [f if f is not None else 100.0 for f in fitnesses]
