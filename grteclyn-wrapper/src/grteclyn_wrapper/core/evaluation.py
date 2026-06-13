"""Shared single-candidate evaluation used by the optimizer and QD search.

Wraps the constrained-recipe / pre-flight / episode-run / score pipeline into
one reusable call so that CMA-ES, MAP-Elites, and any future driver evaluate
candidates identically.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from .config import ExampleConfig, ExecutableConfig
from .episode import create_episode, update_metadata, write_json
from .params import write_params
from .runner import run_episode
from ..initial_data.constrained_recipe import constrained_overrides
from ..initial_data.preflight import preflight_check
from ..metrics import dataclass_to_dict, read_episode_metrics
from ..metrics.score import domain_half_width_from_overrides, score_episode


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
    grtresna: bool = False,
    grtresna_base: Any | None = None,
    grtresna_solved_ftl_gate: bool = False,
    solved_ftl_gate_config: Any | None = None,
    grtresna_convergence_config: Any | None = None,
    grtresna_postload_gate: bool = False,
    postload_gate_config: Any | None = None,
) -> Evaluation:
    overrides = dict(overrides)

    if constrained and not grtresna:
        constrained_overrides(overrides, phantom=phantom)

    # The constrained recipe solves for phi assuming phantom (rho <= 0) coupling
    # when phantom=True.  Tell the C++ level to evolve the matching exotic matter
    # so the geometry is faithfully sourced (otherwise a canonical, rho >= 0
    # field is evolved and the Hamiltonian constraint is violated at t=0).
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
    episode = create_episode(
        out_dir,
        name=name,
        metadata={"mode": "evaluate", "example": example.name, "overrides": overrides},
    )

    if grtresna:
        from ..grtresna.matter_wiring import merge_evolution_overrides
        from ..grtresna.solver import solve as grtresna_solve
        from ..search.grtresna_convergence_gate import (
            GRTRESNA_REJECTION_BASE_FITNESS,
            GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
        )
        from ..search.grtresna_evaluation_gates import (
            GRTresnaPreEvolutionGateConfig,
            apply_grtresna_pre_evolution_gates,
        )
        from ..search.optimize import build_grtresna_config, parse_convergence_safe

        try:
            cfg = build_grtresna_config(overrides, grtresna_base)
            if not dry_run:
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
                    return Evaluation(
                        score=-gate_rejection.fitness,
                        components={
                            gate_rejection.component_key: -gate_rejection.fitness
                        },
                        notes=list(gate_rejection.notes) or [gate_rejection.reason],
                        episode_path=str(episode.path),
                        exit_code=None,
                        preflight_rejected=False,
                        reason=gate_rejection.reason,
                        metrics=gate_rejection.metrics or {},
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
                consumer_ftl_timeseries=consumer_ftl_timeseries,
                consumer_ftl_L=ftl_L,
            )
            exit_code = result.returncode
        except Exception as exc:  # noqa: BLE001 - record and continue
            exit_code = 1
            update_metadata(episode, {"simulation_error": repr(exc), "simulation_exit_code": exit_code})

    metrics = read_episode_metrics(episode.path, ftl_L=ftl_L)
    score = score_episode(
        metrics,
        target_stop_time=target_stop_time,
        weights=score_weights,
        objective_mode=objective_mode,
        domain_half_width=domain_half_width_from_overrides(gte_overrides),
    )

    write_json(episode.score_path, {
        "score": dataclass_to_dict(score),
        "metrics": dataclass_to_dict(metrics),
    })

    return Evaluation(
        score=score.total,
        components=dict(score.components),
        notes=list(score.notes),
        episode_path=str(episode.path),
        exit_code=exit_code,
        preflight_rejected=False,
        reason=None,
        metrics={k: dataclass_to_dict(getattr(metrics, k)) for k in metrics.__dataclass_fields__},
    )
