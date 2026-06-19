"""Splash-aware early termination predicates."""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class SplashEarlyTermDecision:
    should_stop: bool
    reason: str | None = None


def evaluate_splash_early_term(
    *,
    t: float,
    rho: float | None,
    lapse: float | None,
    activity: float | None,
    peak_rho_so_far: float,
    previous_activity: float | None = None,
) -> SplashEarlyTermDecision:
    if rho is not None and not math.isfinite(rho):
        return SplashEarlyTermDecision(True, "numerical_failure")
    if lapse is not None and not math.isfinite(lapse):
        return SplashEarlyTermDecision(True, "numerical_failure")
    if activity is not None and not math.isfinite(activity):
        return SplashEarlyTermDecision(True, "numerical_failure")

    if lapse is not None and lapse < 0.005:
        return SplashEarlyTermDecision(True, "collapse_complete")

    if t > 6.0 and rho is not None and math.isfinite(rho) and peak_rho_so_far > 0.0:
        if rho < 0.1 * peak_rho_so_far:
            declining = (
                previous_activity is not None
                and activity is not None
                and math.isfinite(previous_activity)
                and math.isfinite(activity)
                and activity < previous_activity
            )
            if declining or (activity is not None and activity < 0.5 * previous_activity if previous_activity else True):
                return SplashEarlyTermDecision(True, "dispersion_complete")

    return SplashEarlyTermDecision(False)
