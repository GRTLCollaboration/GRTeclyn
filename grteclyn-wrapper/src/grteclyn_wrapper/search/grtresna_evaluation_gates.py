"""Pre-evolution gates for GRTresna-in-the-loop search evaluations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from ..core.episode import Episode, update_metadata
from ..projection.postload_gate import PostLoadGateConfig, run_postload_gate
from ..metrics.ftl_solved_geometry import compute_solved_geometry_ftl
from .grtresna_convergence_gate import (
    GRTRESNA_REJECTION_BASE_FITNESS,
    GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
    GRTresnaConvergenceConfig,
    convergence_rejection_fitness,
    convergence_rejection_reason,
)
from .solved_ftl_gate import (
    SolvedFtlGateConfig,
    solved_ftl_has_signal,
    solved_geometry_ftl_to_dict,
    solved_geometry_rejection_fitness,
)

POSTLOAD_REJECTION_BASE_FITNESS = 80.0
POSTLOAD_REJECTION_MAX_EXTRA_FITNESS = 20.0
SOLVED_FTL_REJECTION_BASE_FITNESS = 75.0
SOLVED_FTL_REJECTION_MAX_EXTRA_FITNESS = 20.0


@dataclass(frozen=True)
class GRTresnaGateRejection:
    """Uniform rejection payload for optimize/QD drivers."""

    reason: str
    fitness: float
    component_key: str
    metadata: dict[str, Any]
    notes: tuple[str, ...] = ()
    metrics: dict[str, Any] | None = None


@dataclass(frozen=True)
class GRTresnaPreEvolutionGateConfig:
    convergence_config: GRTresnaConvergenceConfig | None = None
    solved_ftl_enabled: bool = False
    solved_ftl_config: SolvedFtlGateConfig | None = None
    postload_enabled: bool = False
    postload_config: PostLoadGateConfig | None = None
    ftl_L: float | None = None


def reject_grtresna_convergence(
    convergence: Mapping[str, float] | None,
    *,
    config: GRTresnaConvergenceConfig | None,
) -> GRTresnaGateRejection | None:
    reason = convergence_rejection_reason(convergence, config=config)
    if reason is None:
        return None
    fitness = convergence_rejection_fitness(
        convergence,
        config=config,
        base=GRTRESNA_REJECTION_BASE_FITNESS,
        max_extra=GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
    )
    return GRTresnaGateRejection(
        reason=reason,
        fitness=fitness,
        component_key="grtresna_rejection",
        metadata={
            "grtresna_rejected": True,
            "grtresna_rejection_reason": reason,
            "grtresna_rejection_fitness": fitness,
        },
    )


def reject_solved_ftl(
    solved_ftl: Any,
    *,
    enabled: bool,
    config: SolvedFtlGateConfig | None,
) -> GRTresnaGateRejection | None:
    if not enabled or solved_ftl_has_signal(solved_ftl, config=config):
        return None
    fitness = solved_geometry_rejection_fitness(
        solved_ftl,
        config=config,
        base=SOLVED_FTL_REJECTION_BASE_FITNESS,
        max_extra=SOLVED_FTL_REJECTION_MAX_EXTRA_FITNESS,
    )
    solved_dict = (
        solved_geometry_ftl_to_dict(solved_ftl, config=config)
        if solved_ftl is not None
        else None
    )
    return GRTresnaGateRejection(
        reason="solved_ftl_rejected",
        fitness=fitness,
        component_key="solved_ftl_rejection",
        metadata={
            "solved_ftl_rejected": True,
            "solved_ftl_rejection_fitness": fitness,
            "solved_geometry_ftl": solved_dict,
        },
        notes=("solved geometry has no FTL/precursor signal",),
        metrics={"solved_geometry_ftl": solved_dict},
    )


def reject_postload_gate(
    gridinit_path: Path,
    *,
    episode: Episode,
    gte_overrides: Mapping[str, Any],
    cuda_devices: str | None,
    config: PostLoadGateConfig | None,
    dry_run: bool = False,
) -> GRTresnaGateRejection | None:
    gate = run_postload_gate(
        gridinit_path,
        out_dir=episode.path / "postload_gate",
        config=config,
        cuda_devices=cuda_devices,
        dry_run=dry_run,
        overrides=dict(gte_overrides),
    )
    update_metadata(
        episode,
        {
            "postload_gate": {
                "passed": gate.passed,
                "max_hamiltonian_l2": gate.max_hamiltonian_l2,
                "max_momentum_l2": gate.max_momentum_l2,
                "reason": gate.reason,
            }
        },
    )
    if gate.passed:
        return None

    fitness = POSTLOAD_REJECTION_BASE_FITNESS + POSTLOAD_REJECTION_MAX_EXTRA_FITNESS
    reason = gate.reason or "postload_constraints_exceeded"
    return GRTresnaGateRejection(
        reason=reason,
        fitness=fitness,
        component_key="postload_rejection",
        metadata={
            "postload_rejected": True,
            "postload_rejection_fitness": fitness,
            "postload_gate": {
                "passed": gate.passed,
                "max_hamiltonian_l2": gate.max_hamiltonian_l2,
                "max_momentum_l2": gate.max_momentum_l2,
                "reason": gate.reason,
            },
        },
        notes=(reason,),
    )


def apply_grtresna_pre_evolution_gates(
    *,
    episode: Episode,
    convergence: Mapping[str, float] | None,
    gridinit_path: Path,
    gte_overrides: Mapping[str, Any],
    cuda_devices: str | None,
    config: GRTresnaPreEvolutionGateConfig,
    dry_run: bool = False,
) -> GRTresnaGateRejection | None:
    """Run convergence, solved-FTL, and post-load gates in order."""
    rejection = reject_grtresna_convergence(
        convergence, config=config.convergence_config,
    )
    if rejection is not None:
        return rejection

    solved_ftl = compute_solved_geometry_ftl(gridinit_path, L=config.ftl_L)
    if solved_ftl is not None:
        update_metadata(
            episode,
            {
                "solved_geometry_ftl": solved_geometry_ftl_to_dict(
                    solved_ftl,
                    config=config.solved_ftl_config,
                )
            },
        )

    rejection = reject_solved_ftl(
        solved_ftl,
        enabled=config.solved_ftl_enabled,
        config=config.solved_ftl_config,
    )
    if rejection is not None:
        return rejection

    return reject_postload_gate(
        gridinit_path,
        episode=episode,
        gte_overrides=gte_overrides,
        cuda_devices=cuda_devices,
        config=config.postload_config,
        dry_run=dry_run,
    )
