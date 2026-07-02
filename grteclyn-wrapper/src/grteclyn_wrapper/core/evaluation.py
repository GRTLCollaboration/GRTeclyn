"""Shared single-candidate evaluation used by the optimizer and QD search.

Wraps the constrained-recipe / pre-flight / episode-run / score pipeline into
one reusable call so that CMA-ES, MAP-Elites, and any future driver evaluate
candidates identically.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping, Sequence

if TYPE_CHECKING:
    from ..search.scoring_pool import ScoringPool

from .artifact_cleanup import cleanup_episode_artifacts
from .config import ExampleConfig, ExecutableConfig
from .episode import Episode, create_episode, update_metadata, write_json
from .params import write_params
from .runner import run_episode
from ..initial_data.constrained_recipe import constrained_overrides
from ..initial_data.preflight import preflight_check
from ..metrics import dataclass_to_dict, read_episode_metrics
from ..metrics.score import domain_half_width_for_episode, score_episode
from ..search.gpu_pool import GpuPool
from ..search.grtresna_evaluation_gates import (
    GRTresnaGateRejection,
    GRTresnaPreEvolutionGateConfig,
    reject_postload_gate,
    run_cpu_gates,
)


@dataclass
class Evaluation:
    score: float
    components: dict[str, float]
    notes: list[str]
    episode_path: str | None
    exit_code: int | None
    preflight_rejected: bool
    reason: str | None
    metrics: dict[str, Any]


@dataclass(frozen=True)
class CpuGateResult:
    episode: Episode
    gte_overrides: dict[str, Any]
    gridinit_path: Path
    gate_config: GRTresnaPreEvolutionGateConfig


def _resolve_consumer_evolving_geodesic(value: bool | None) -> bool:
    if value is not None:
        return value
    return os.environ.get("GRTECLYN_EVOLVING_GEODESIC", "").strip().lower() in {
        "1",
        "on",
        "yes",
        "true",
    }


def _gate_rejection_to_evaluation(
    gate_rejection: GRTresnaGateRejection,
    *,
    episode: Episode,
) -> Evaluation:
    update_metadata(
        episode,
        {
            **gate_rejection.metadata,
            "simulation_exit_code": None,
        },
    )
    return Evaluation(
        score=-gate_rejection.fitness,
        components={gate_rejection.component_key: -gate_rejection.fitness},
        notes=list(gate_rejection.notes) or [gate_rejection.reason],
        episode_path=str(episode.path),
        exit_code=None,
        preflight_rejected=False,
        reason=gate_rejection.reason,
        metrics=gate_rejection.metrics or {},
    )


def _run_cpu_grtresna_gates(
    *,
    overrides: Mapping[str, Any],
    out_dir: Path,
    name: str,
    example: ExampleConfig,
    grtresna_base: Any | None,
    grtresna_convergence_config: Any | None,
    grtresna_solved_ftl_gate: bool,
    solved_ftl_gate_config: Any | None,
    ftl_L: float | None,
    dry_run: bool,
) -> CpuGateResult | Evaluation:
    from ..grtresna.matter.wiring import merge_evolution_overrides
    from ..grtresna.solver import solve as grtresna_solve
    from ..search.grtresna_convergence_gate import (
        GRTRESNA_REJECTION_BASE_FITNESS,
        GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
    )
    from ..search.optimize import build_grtresna_config, parse_convergence_safe

    episode = create_episode(
        out_dir,
        name=name,
        metadata={"mode": "evaluate", "example": example.name, "overrides": dict(overrides)},
    )
    gte_overrides = {
        k: v for k, v in overrides.items() if not str(k).startswith("grtresna_")
    }
    gate_config = GRTresnaPreEvolutionGateConfig(
        convergence_config=grtresna_convergence_config,
        solved_ftl_enabled=grtresna_solved_ftl_gate,
        solved_ftl_config=solved_ftl_gate_config,
        postload_enabled=False,
        ftl_L=ftl_L,
    )

    try:
        cfg = build_grtresna_config(overrides, grtresna_base)
        if dry_run:
            gridinit = episode.path / "initial_data.gridinit"
            gridinit.write_bytes(b"")
            return CpuGateResult(
                episode=episode,
                gte_overrides=gte_overrides,
                gridinit_path=gridinit,
                gate_config=gate_config,
            )

        gridinit = grtresna_solve(
            cfg,
            work_dir=episode.path / "grtresna",
            gridinit_path=episode.path / "initial_data.gridinit",
        )
        gte_overrides["recipe_initial_data_file"] = str(gridinit)
        gte_overrides.update(merge_evolution_overrides(gte_overrides, cfg))
        convergence = parse_convergence_safe(episode.path / "grtresna")
        update_metadata(episode, {"grtresna_convergence": convergence})

        rejection = run_cpu_gates(
            episode=episode,
            convergence=convergence,
            gridinit_path=Path(gridinit),
            config=gate_config,
        )
        if rejection is not None:
            return _gate_rejection_to_evaluation(rejection, episode=episode)

        return CpuGateResult(
            episode=episode,
            gte_overrides=gte_overrides,
            gridinit_path=Path(gridinit),
            gate_config=gate_config,
        )
    except Exception as exc:  # noqa: BLE001 - record and continue
        fitness = GRTRESNA_REJECTION_BASE_FITNESS + GRTRESNA_REJECTION_MAX_EXTRA_FITNESS
        update_metadata(episode, {
            "grtresna_error": repr(exc),
            "simulation_exit_code": None,
        })
        return Evaluation(
            score=-fitness,
            components={"grtresna_rejection": -fitness},
            notes=[repr(exc)],
            episode_path=str(episode.path),
            exit_code=None,
            preflight_rejected=False,
            reason=repr(exc),
            metrics={},
        )


@dataclass(frozen=True)
class ScoringRequest:
    """Everything the CPU scoring stage needs, decoupled from the GPU lease.

    All fields are plain data so this can be shipped to a separate scoring
    process (``ScoringPool``), which lets scoring run in parallel across evals
    without the process-wide ``PLOTFILE_READ_LOCK`` (yt is not thread-safe) and
    without holding a GPU slot while the CPU-bound yt analysis runs.
    """

    episode_path: str
    exit_code: int | None
    dry_run: bool
    target_stop_time: float | None
    score_weights: dict[str, float] | None
    objective_mode: str
    ftl_L: float | None
    evolving_enabled: bool
    gte_overrides: dict[str, Any]
    cleanup_plotfiles: bool


def _run_gpu_binary_phase(
    *,
    cpu_result: CpuGateResult,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    cuda_devices: str | None,
    check_params: bool,
    dry_run: bool,
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    objective_mode: str,
    ftl_L: float | None,
    consume_plotfiles: bool,
    consumer_radii: Sequence[float],
    consumer_keep_last: int,
    consumer_ftl_timeseries: bool,
    consumer_evolving_geodesic: bool | None,
    grtresna_postload_gate: bool,
    postload_gate_config: Any | None,
    cleanup_plotfiles: bool,
) -> ScoringRequest | Evaluation:
    """GPU-bound phase: postload gate + write params + run the evolution binary.

    Runs while a GPU slot is held. Returns a :class:`ScoringRequest` describing
    the CPU-only scoring work to do next, or an :class:`Evaluation` if the run
    was rejected before scoring (postload gate).
    """
    episode = cpu_result.episode
    gte_overrides = dict(cpu_result.gte_overrides)
    if objective_mode == "critical_collapse":
        from ..grtresna.matter.splash import apply_splash_overrides

        gte_overrides = apply_splash_overrides(gte_overrides)

    if grtresna_postload_gate:
        rejection = reject_postload_gate(
            cpu_result.gridinit_path,
            episode=episode,
            gte_overrides=gte_overrides,
            cuda_devices=cuda_devices,
            config=postload_gate_config,
            dry_run=dry_run,
        )
        if rejection is not None:
            return _gate_rejection_to_evaluation(rejection, episode=episode)

    write_params(
        template, episode.params_path,
        episode_dir=episode.path, example=example, overrides=gte_overrides,
    )

    exit_code: int | None = None
    evolving_enabled = _resolve_consumer_evolving_geodesic(consumer_evolving_geodesic)
    if dry_run:
        update_metadata(episode, {"dry_run": True})
    else:
        if executable is None:
            raise ValueError("executable is required unless dry_run=True")
        try:
            extra_env = None
            if evolving_enabled:
                extra_env = {"GRTECLYN_EVOLVING_GEODESIC": "1"}
            central_timeseries_enabled = os.environ.get("GRTECLYN_CENTRAL_TIMESERIES", "").strip().lower() in {
                "1",
                "on",
                "yes",
                "true",
            }
            consumer_incremental = consumer_ftl_timeseries or (
                objective_mode == "critical_collapse" and central_timeseries_enabled
            )
            result = run_episode(
                episode, executable,
                check_params=check_params,
                cuda_devices=cuda_devices,
                extra_env=extra_env,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_delete=True,
                consumer_keep_last=consumer_keep_last,
                consumer_ftl_timeseries=consumer_ftl_timeseries,
                consumer_ftl_L=ftl_L,
                consumer_incremental_score=consumer_incremental,
                consumer_objective_mode=objective_mode,
                consumer_target_stop_time=target_stop_time,
                consumer_score_weights=score_weights,
                consumer_evolving_geodesic=evolving_enabled,
            )
            exit_code = result.returncode
            # Early termination via .stop_sim is intentional (e.g.
            # dispersion_complete); the SIGTERM exit code should not mark
            # the evaluation as gpu_failed.
            if exit_code != 0 and (episode.path / ".stop_sim").is_file():
                exit_code = 0
        except Exception as exc:  # noqa: BLE001 - record and continue
            exit_code = 1
            update_metadata(episode, {"simulation_error": repr(exc), "simulation_exit_code": exit_code})

    return ScoringRequest(
        episode_path=str(episode.path),
        exit_code=exit_code,
        dry_run=dry_run,
        target_stop_time=target_stop_time,
        score_weights=dict(score_weights) if score_weights is not None else None,
        objective_mode=objective_mode,
        ftl_L=ftl_L,
        evolving_enabled=evolving_enabled,
        gte_overrides=gte_overrides,
        cleanup_plotfiles=cleanup_plotfiles,
    )


def _score_episode_phase(request: ScoringRequest) -> Evaluation:
    """CPU-only scoring phase: read metrics, score, persist, clean up.

    Pure function of an on-disk episode directory so it can run inline or in a
    separate process. Holds no GPU slot.
    """
    episode = Episode(path=Path(request.episode_path))

    metrics = read_episode_metrics(
        episode.path,
        ftl_L=request.ftl_L,
        evolving_geodesic=request.evolving_enabled,
        objective_mode=request.objective_mode,
    )
    score = score_episode(
        metrics,
        target_stop_time=request.target_stop_time,
        weights=request.score_weights,
        objective_mode=request.objective_mode,
        domain_half_width=domain_half_width_for_episode(
            episode.path, request.gte_overrides
        ),
    )

    write_json(episode.score_path, {
        "score": dataclass_to_dict(score),
        "metrics": dataclass_to_dict(metrics),
    })

    if request.cleanup_plotfiles and episode.score_path.is_file():
        cleanup_episode_artifacts(episode.path, tier="plotfiles_only")

    return Evaluation(
        score=score.total,
        components=dict(score.components),
        notes=[*score.notes, *(["dry_run"] if request.dry_run else [])],
        episode_path=str(episode.path),
        exit_code=request.exit_code,
        preflight_rejected=False,
        reason=None,
        metrics={k: dataclass_to_dict(getattr(metrics, k)) for k in metrics.__dataclass_fields__},
    )


def _run_gpu_session(
    *,
    cpu_result: CpuGateResult,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    cuda_devices: str | None = None,
    check_params: bool,
    dry_run: bool,
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    objective_mode: str,
    ftl_L: float | None,
    consume_plotfiles: bool,
    consumer_radii: Sequence[float],
    consumer_keep_last: int,
    consumer_ftl_timeseries: bool,
    consumer_evolving_geodesic: bool | None,
    grtresna_postload_gate: bool,
    postload_gate_config: Any | None,
    cleanup_plotfiles: bool = True,
    gpu_pool: GpuPool | None = None,
    scoring_pool: "ScoringPool | None" = None,
) -> Evaluation:
    """Run the GPU evolution then score the episode.

    When ``gpu_pool`` is given the GPU slot is leased only for the evolution
    binary and released before the CPU-bound scoring stage. When ``scoring_pool``
    is given, scoring runs in a separate process so it neither holds a GPU slot
    nor contends on the process-wide yt read lock.
    """
    binary_kwargs = dict(
        cpu_result=cpu_result,
        example=example,
        template=template,
        executable=executable,
        check_params=check_params,
        dry_run=dry_run,
        target_stop_time=target_stop_time,
        score_weights=score_weights,
        objective_mode=objective_mode,
        ftl_L=ftl_L,
        consume_plotfiles=consume_plotfiles,
        consumer_radii=consumer_radii,
        consumer_keep_last=consumer_keep_last,
        consumer_ftl_timeseries=consumer_ftl_timeseries,
        consumer_evolving_geodesic=consumer_evolving_geodesic,
        grtresna_postload_gate=grtresna_postload_gate,
        postload_gate_config=postload_gate_config,
        cleanup_plotfiles=cleanup_plotfiles,
    )

    if gpu_pool is not None:
        with gpu_pool.lease() as leased_gpu:
            phase = _run_gpu_binary_phase(cuda_devices=leased_gpu, **binary_kwargs)
    else:
        phase = _run_gpu_binary_phase(cuda_devices=cuda_devices, **binary_kwargs)

    if isinstance(phase, Evaluation):
        return phase

    if scoring_pool is not None:
        return scoring_pool.score(phase)
    return _score_episode_phase(phase)


def evaluate_overrides(
    overrides: Mapping[str, Any],
    *,
    out_dir: Path,
    name: str,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    constrained: bool = True,
    phantom: bool = True,
    use_preflight: bool = True,
    cuda_devices: str | None = None,
    check_params: bool = True,
    dry_run: bool = False,
    target_stop_time: float | None = None,
    score_weights: Mapping[str, float] | None = None,
    objective_mode: str = "weighted",
    ftl_L: float | None = None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
    consumer_keep_last: int = 1,
    consumer_ftl_timeseries: bool = True,
    consumer_evolving_geodesic: bool | None = None,
    grtresna: bool = False,
    grtresna_base: Any | None = None,
    grtresna_solved_ftl_gate: bool = False,
    solved_ftl_gate_config: Any | None = None,
    grtresna_convergence_config: Any | None = None,
    grtresna_postload_gate: bool = False,
    postload_gate_config: Any | None = None,
    gpu_pool: GpuPool | None = None,
    scoring_pool: "ScoringPool | None" = None,
    cleanup_plotfiles: bool = True,
) -> Evaluation:
    overrides = dict(overrides)

    if constrained and not grtresna:
        constrained_overrides(overrides, phantom=phantom)

    if phantom and not grtresna:
        overrides.setdefault("recipe_exotic_matter", 1)

    if use_preflight and not grtresna:
        pf = preflight_check(overrides, phantom=phantom)
        if not pf.passed:
            return Evaluation(
                score=0.0,
                components={},
                notes=[f"preflight_rejected: {pf.reason}"],
                episode_path=None,
                exit_code=None,
                preflight_rejected=True,
                reason=pf.reason,
                metrics={},
            )

    gte_overrides = {
        k: v for k, v in overrides.items() if not str(k).startswith("grtresna_")
    }

    if grtresna:
        cpu_phase = _run_cpu_grtresna_gates(
            overrides=overrides,
            out_dir=out_dir,
            name=name,
            example=example,
            grtresna_base=grtresna_base,
            grtresna_convergence_config=grtresna_convergence_config,
            grtresna_solved_ftl_gate=grtresna_solved_ftl_gate,
            solved_ftl_gate_config=solved_ftl_gate_config,
            ftl_L=ftl_L,
            dry_run=dry_run,
        )
        if isinstance(cpu_phase, Evaluation):
            return cpu_phase
        cpu_result = cpu_phase
        episode = cpu_result.episode
    else:
        episode = create_episode(
            out_dir,
            name=name,
            metadata={"mode": "evaluate", "example": example.name, "overrides": overrides},
        )
        cpu_result = CpuGateResult(
            episode=episode,
            gte_overrides=gte_overrides,
            gridinit_path=episode.path / "initial_data.gridinit",
            gate_config=GRTresnaPreEvolutionGateConfig(),
        )

    gpu_kwargs = dict(
        cpu_result=cpu_result,
        example=example,
        template=template,
        executable=executable,
        check_params=check_params,
        dry_run=dry_run,
        target_stop_time=target_stop_time,
        score_weights=score_weights,
        objective_mode=objective_mode,
        ftl_L=ftl_L,
        consume_plotfiles=consume_plotfiles,
        consumer_radii=consumer_radii,
        consumer_keep_last=consumer_keep_last,
        consumer_ftl_timeseries=consumer_ftl_timeseries,
        consumer_evolving_geodesic=consumer_evolving_geodesic,
        grtresna_postload_gate=grtresna_postload_gate,
        postload_gate_config=postload_gate_config,
        cleanup_plotfiles=cleanup_plotfiles,
    )

    # _run_gpu_session leases the GPU only for the evolution binary and releases
    # it before the (optionally out-of-process) scoring stage.
    return _run_gpu_session(
        cuda_devices=cuda_devices,
        gpu_pool=gpu_pool,
        scoring_pool=scoring_pool,
        **gpu_kwargs,
    )
