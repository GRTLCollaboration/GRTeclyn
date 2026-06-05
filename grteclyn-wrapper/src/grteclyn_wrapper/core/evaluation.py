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
from ..metrics.episode_metrics import dataclass_to_dict, read_episode_metrics
from ..metrics.score import score_episode


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
    grtresna: bool = False,
    grtresna_base: Any | None = None,
    grtresna_solved_ftl_gate: bool = False,
    solved_ftl_gate_config: Any | None = None,
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

    if grtresna and not dry_run:
        from ..grtresna.solver import solve as grtresna_solve
        from ..metrics.ftl_solved_geometry import compute_solved_geometry_ftl
        from ..search.optimize import (
            GRTRESNA_REJECTION_BASE_FITNESS,
            GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
            SOLVED_FTL_REJECTION_BASE_FITNESS,
            SOLVED_FTL_REJECTION_MAX_EXTRA_FITNESS,
            _grtresna_convergence_rejection_reason,
            _grtresna_rejection_fitness,
            build_grtresna_config,
            parse_convergence_safe,
        )
        from ..search.solved_ftl_gate import (
            solved_ftl_has_signal,
            solved_geometry_ftl_to_dict,
            solved_geometry_rejection_fitness,
        )

        try:
            cfg = build_grtresna_config(overrides, grtresna_base)
            gridinit = grtresna_solve(
                cfg,
                work_dir=episode.path / "grtresna",
                gridinit_path=episode.path / "initial_data.gridinit",
            )
            gte_overrides["recipe_initial_data_file"] = gridinit
            convergence = parse_convergence_safe(episode.path / "grtresna")
            update_metadata(episode, {"grtresna_convergence": convergence})

            rejection_reason = _grtresna_convergence_rejection_reason(convergence)
            if rejection_reason is not None:
                fitness = _grtresna_rejection_fitness(convergence)
                update_metadata(episode, {
                    "grtresna_rejected": True,
                    "grtresna_rejection_reason": rejection_reason,
                    "grtresna_rejection_fitness": fitness,
                    "simulation_exit_code": None,
                })
                return Evaluation(
                    score=-fitness,
                    components={"grtresna_rejection": -fitness},
                    notes=[rejection_reason],
                    episode_path=str(episode.path),
                    exit_code=None,
                    preflight_rejected=False,
                    reason=rejection_reason,
                    metrics={},
                )

            solved_ftl = compute_solved_geometry_ftl(
                episode.path / "initial_data.gridinit", L=ftl_L,
            )
            if solved_ftl is not None:
                update_metadata(
                    episode,
                    {
                        "solved_geometry_ftl": solved_geometry_ftl_to_dict(
                            solved_ftl,
                            config=solved_ftl_gate_config,
                        )
                    },
                )
            if (
                grtresna_solved_ftl_gate
                and not solved_ftl_has_signal(solved_ftl, config=solved_ftl_gate_config)
            ):
                fitness = solved_geometry_rejection_fitness(
                    solved_ftl,
                    config=solved_ftl_gate_config,
                    base=SOLVED_FTL_REJECTION_BASE_FITNESS,
                    max_extra=SOLVED_FTL_REJECTION_MAX_EXTRA_FITNESS,
                )
                solved_dict = (
                    solved_geometry_ftl_to_dict(
                        solved_ftl,
                        config=solved_ftl_gate_config,
                    )
                    if solved_ftl is not None
                    else None
                )
                update_metadata(episode, {
                    "solved_ftl_rejected": True,
                    "solved_ftl_rejection_fitness": fitness,
                    "simulation_exit_code": None,
                })
                return Evaluation(
                    score=-fitness,
                    components={"solved_ftl_rejection": -fitness},
                    notes=["solved geometry has no FTL/precursor signal"],
                    episode_path=str(episode.path),
                    exit_code=None,
                    preflight_rejected=False,
                    reason="solved_ftl_rejected",
                    metrics={"solved_geometry_ftl": solved_dict},
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
