"""Quantum energy inequality proxy metrics."""

from __future__ import annotations

from ..probes.physical import PhysicalMetrics
from ..types.diagnostics import QeiMetrics


def read_qei_metrics(
    physical: PhysicalMetrics | None,
    *,
    trajectory_violation: float | None = None,
) -> QeiMetrics | None:
    spatial = getattr(physical, "qei_spatial_proxy", None) if physical else None
    if spatial is None and trajectory_violation is None:
        return None
    violation = trajectory_violation if trajectory_violation is not None else spatial
    s_qei = 1.0 / (1.0 + max(0.0, violation or 0.0)) if violation is not None else None
    return QeiMetrics(
        spatial_proxy=spatial,
        trajectory_violation=trajectory_violation,
        s_qei=s_qei,
    )
