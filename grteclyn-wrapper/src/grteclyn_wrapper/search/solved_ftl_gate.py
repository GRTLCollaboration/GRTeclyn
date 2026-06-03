"""Policy for the pre-GPU solved-geometry FTL gate.

The metric layer computes geometry diagnostics.  This module owns the search
policy that decides whether those diagnostics are useful enough to evolve.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from ..metrics.ftl_solved_geometry import SolvedGeometryFtl


@dataclass(frozen=True)
class SolvedFtlGateConfig:
    """Configurable thresholds for the solved-geometry FTL gate."""

    f_op_floor: float = 1.0e-4
    near_luminal_speed_floor: float = 0.99
    superluminal_speed_floor: float = 1.01
    superluminal_fraction_floor: float = 0.02
    max_physical_coord_speed: float = 8.0
    max_physical_f_op: float = 0.85
    rejection_speed_target: float = 1.01
    rejection_speed_width: float = 0.25
    rejection_f_op_target: float = 1.0e-4
    rejection_fraction_target: float = 0.02


DEFAULT_SOLVED_FTL_GATE_CONFIG = SolvedFtlGateConfig()


@dataclass(frozen=True)
class SolvedFtlGatePolicy:
    """Decision service for solved-geometry pre-evolution filtering."""

    config: SolvedFtlGateConfig = DEFAULT_SOLVED_FTL_GATE_CONFIG

    def is_degenerate(self, solved: SolvedGeometryFtl) -> bool:
        """True when the operational signal is a near-degenerate-cell artifact."""
        op = solved.operational
        if not math.isfinite(op.max_local_speed):
            return True
        if op.max_local_speed > self.config.max_physical_coord_speed:
            return True
        if op.f_op > self.config.max_physical_f_op:
            return True
        return False

    def has_signal(self, solved: SolvedGeometryFtl | None) -> bool:
        """True if the solved geometry shows a physical FTL signal or precursor."""
        if solved is None:
            return False
        if self.is_degenerate(solved):
            return False
        op = solved.operational
        if op.f_op > self.config.f_op_floor:
            return True
        if self.config.near_luminal_speed_floor <= op.max_local_speed < 1.0:
            return True
        if (
            op.max_local_speed >= self.config.superluminal_speed_floor
            and op.superluminal_fraction >= self.config.superluminal_fraction_floor
        ):
            return True
        if op.superluminal_fraction >= self.config.superluminal_fraction_floor:
            return True
        return False

    def rejection_fitness(
        self,
        solved: SolvedGeometryFtl | None,
        *,
        base: float = 75.0,
        max_extra: float = 20.0,
    ) -> float:
        """Graded penalty when solved geometry shows no FTL signal."""
        cfg = self.config
        if solved is None or self.is_degenerate(solved):
            return base + max_extra
        op = solved.operational
        f_op_deficit = max(
            0.0,
            (cfg.rejection_f_op_target - op.f_op)
            / max(cfg.rejection_f_op_target, 1.0e-12),
        )
        speed_deficit = max(
            0.0,
            (cfg.rejection_speed_target - op.max_local_speed)
            / max(cfg.rejection_speed_width, 1.0e-12),
        )
        frac_deficit = max(
            0.0,
            (cfg.rejection_fraction_target - op.superluminal_fraction)
            / max(cfg.rejection_fraction_target, 1.0e-12),
        )
        deficit = min(
            1.0,
            0.15 * f_op_deficit + 0.70 * speed_deficit + 0.15 * frac_deficit,
        )
        return base + max_extra * deficit

    def to_dict(self, solved: SolvedGeometryFtl) -> dict[str, float]:
        """JSON-serializable solved-FTL summary under this policy."""
        op = solved.operational
        mech = solved.mechanisms
        degenerate = self.is_degenerate(solved)
        return {
            "f_op": float(op.f_op),
            "f_op_physical": 0.0 if degenerate else float(op.f_op),
            "degenerate": float(degenerate),
            "t_min": float(op.t_min) if op.t_min is not None else float("nan"),
            "t_flat": float(op.t_flat),
            "max_local_speed": float(op.max_local_speed),
            "superluminal_fraction": float(op.superluminal_fraction),
            "path_offaxis": float(op.path_offaxis),
            "reachable": float(op.reachable),
            "shift_warp": mech.shift_warp,
            "throat_pinch": mech.throat_pinch,
            "portal_compression": mech.portal_compression,
            "mechanism_descriptor": mech.mechanism_descriptor,
            "integration_L": solved.integration_L,
        }


DEFAULT_SOLVED_FTL_GATE_POLICY = SolvedFtlGatePolicy()


def _policy_with_overrides(
    config: SolvedFtlGateConfig | None,
    *,
    f_op_floor: float | None = None,
    precursor_speed_floor: float | None = None,
    precursor_frac_floor: float | None = None,
) -> SolvedFtlGatePolicy:
    cfg = config or DEFAULT_SOLVED_FTL_GATE_CONFIG
    if f_op_floor is not None:
        cfg = replace(cfg, f_op_floor=f_op_floor)
    if precursor_speed_floor is not None:
        cfg = replace(cfg, near_luminal_speed_floor=precursor_speed_floor)
    if precursor_frac_floor is not None:
        cfg = replace(cfg, superluminal_fraction_floor=precursor_frac_floor)
    return SolvedFtlGatePolicy(cfg)


def solved_ftl_has_signal(
    solved: SolvedGeometryFtl | None,
    *,
    config: SolvedFtlGateConfig | None = None,
    f_op_floor: float | None = None,
    precursor_speed_floor: float | None = None,
    precursor_frac_floor: float | None = None,
) -> bool:
    """Backward-compatible functional wrapper around `SolvedFtlGatePolicy`."""
    return _policy_with_overrides(
        config,
        f_op_floor=f_op_floor,
        precursor_speed_floor=precursor_speed_floor,
        precursor_frac_floor=precursor_frac_floor,
    ).has_signal(solved)


def solved_geometry_rejection_fitness(
    solved: SolvedGeometryFtl | None,
    *,
    config: SolvedFtlGateConfig | None = None,
    base: float = 75.0,
    max_extra: float = 20.0,
) -> float:
    """Backward-compatible functional wrapper around `SolvedFtlGatePolicy`."""
    return SolvedFtlGatePolicy(
        config or DEFAULT_SOLVED_FTL_GATE_CONFIG
    ).rejection_fitness(solved, base=base, max_extra=max_extra)


def solved_geometry_ftl_to_dict(
    solved: SolvedGeometryFtl,
    *,
    config: SolvedFtlGateConfig | None = None,
) -> dict[str, float]:
    """Backward-compatible functional wrapper around `SolvedFtlGatePolicy`."""
    return SolvedFtlGatePolicy(config or DEFAULT_SOLVED_FTL_GATE_CONFIG).to_dict(solved)
