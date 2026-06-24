"""Policy for the post-solve GRTresna Hamiltonian/momentum convergence gate."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

GRTRESNA_REJECTION_BASE_FITNESS = 100.0
GRTRESNA_REJECTION_MAX_EXTRA_FITNESS = 250.0


@dataclass(frozen=True)
class GRTresnaConvergenceConfig:
    """Configurable Ham/Mom residual thresholds for usable initial data."""

    max_ham_pct: float = 5.0
    max_mom_pct: float = 5.0


DEFAULT_GRTRESNA_CONVERGENCE_CONFIG = GRTresnaConvergenceConfig()


def convergence_rejection_reason(
    convergence: Mapping[str, float] | None,
    *,
    config: GRTresnaConvergenceConfig | None = None,
) -> str | None:
    """Return a rejection reason for unusable GRTresna initial data."""
    cfg = config or DEFAULT_GRTRESNA_CONVERGENCE_CONFIG
    if not convergence:
        return "missing GRTresna convergence diagnostics"
    try:
        ham = float(convergence["ham_pct"])
        mom = float(convergence["mom_pct"])
    except (KeyError, TypeError, ValueError):
        return "invalid GRTresna convergence diagnostics"
    # Static (time-symmetric) data has Pi=0 everywhere, so the momentum
    # constraint is trivially satisfied.  GRTresna may report NaN for Mom
    # in this case because it skips the momentum solve.  Treat NaN Mom as
    # 0% when Ham is finite (the constraint *is* satisfied analytically).
    if not math.isfinite(mom) and math.isfinite(ham):
        mom = 0.0
    if not math.isfinite(ham) or not math.isfinite(mom):
        return f"nonfinite GRTresna convergence: Ham={ham}, Mom={mom}"
    if ham > cfg.max_ham_pct or mom > cfg.max_mom_pct:
        return (
            f"GRTresna convergence too poor: Ham={ham:.6g}% "
            f"Mom={mom:.6g}% thresholds=({cfg.max_ham_pct:g}%, {cfg.max_mom_pct:g}%)"
        )
    return None


def convergence_rejection_fitness(
    convergence: Mapping[str, float] | None,
    *,
    config: GRTresnaConvergenceConfig | None = None,
    base: float = GRTRESNA_REJECTION_BASE_FITNESS,
    max_extra: float = GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
) -> float:
    """Graded minimization penalty so search can rank bad GRTresna candidates."""
    cfg = config or DEFAULT_GRTRESNA_CONVERGENCE_CONFIG
    if not convergence:
        return base + max_extra
    try:
        ham = float(convergence["ham_pct"])
        mom = float(convergence["mom_pct"])
    except (KeyError, TypeError, ValueError):
        return base + max_extra
    # Same NaN-Mom leniency as the rejection-reason gate above.
    if not math.isfinite(mom) and math.isfinite(ham):
        mom = 0.0
    if not math.isfinite(ham) or not math.isfinite(mom):
        return base + max_extra
    residual_penalty = (
        max(0.0, ham - cfg.max_ham_pct) + max(0.0, mom - cfg.max_mom_pct)
    )
    return base + min(max_extra, residual_penalty)
