from __future__ import annotations

from ..types import CollapseMetrics

# A trapped surface is interior, so the minimum-expansion radius must sit well
# inside the domain.  If it lands beyond this fraction of the domain half-width
# the apparent-horizon proxy is miscentered (corner origin) or reading boundary
# noise, and its trapped verdict is rejected rather than silently vetoing the
# candidate.
HORIZON_OFFCENTER_FRACTION: float = 0.5
# Corroborated trapped surfaces that appear only in the trailing portion of a
# run (after the FTL measurement window) should not veto the candidate with the
# full binary penalty — the early-time theta+ blips are uncorroborated artifacts.
LATE_COLLAPSE_TIME_FRACTION: float = 0.75


def horizon_penalty_from_collapse(
    collapse: CollapseMetrics,
    *,
    domain_half_width: float | None,
) -> tuple[float, list[str]]:
    """Score-side trapped-surface proxy with corroboration and off-center guards."""
    notes: list[str] = []
    horizon = collapse.max_horizon_radius or 0.0
    r_at_min_theta = collapse.r_at_min_theta_plus
    horizon_offcenter = (
        horizon > 0.0
        and domain_half_width is not None
        and domain_half_width > 0.0
        and r_at_min_theta is not None
        and r_at_min_theta > HORIZON_OFFCENTER_FRACTION * domain_half_width
    )
    if horizon_offcenter:
        notes.append(
            "horizon proxy located at r="
            f"{r_at_min_theta:.1f} > {HORIZON_OFFCENTER_FRACTION:g}*half-width "
            f"({domain_half_width:.1f}); rejected as miscentered/boundary "
            "artifact (no interior trapped surface), horizon penalty suppressed"
        )
        return 0.0, notes

    if not collapse.corroborated_trapped:
        if horizon > 0.0:
            notes.append(
                "theta+<=0 without lapse-collapse corroboration at the same "
                "timestep; trapped-surface proxy suppressed"
            )
        return 0.0, notes

    final_time = collapse.final_time or 0.0
    first_time = collapse.first_corroborated_time or 0.0
    if (
        final_time > 0.0
        and first_time > LATE_COLLAPSE_TIME_FRACTION * final_time
    ):
        notes.append(
            "corroborated trapped surface only in the trailing "
            f"{100.0 * (1.0 - LATE_COLLAPSE_TIME_FRACTION):.0f}% of the run "
            f"(first at t={first_time:.2f} / {final_time:.2f}); "
            "late collapse penalty suppressed"
        )
        return 0.0, notes

    return -min(horizon, 1.0), notes
