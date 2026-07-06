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
# Kept as a hard suppression threshold: a collapse this late is likely a
# numerical end-of-run artifact, not a physical event.
LATE_COLLAPSE_TIME_FRACTION: float = 0.90


def horizon_penalty_from_collapse(
    collapse: CollapseMetrics,
    *,
    domain_half_width: float | None,
    target_stop_time: float | None = None,
) -> tuple[float, list[str]]:
    """Score-side trapped-surface proxy with corroboration and off-center guards.

    The penalty is **graded by collapse timing**: a configuration that survives
    longer before collapsing receives a smaller penalty than one that collapses
    immediately.  This gives the optimizer a smooth gradient to *delay* collapse
    rather than a binary veto that floors every collapsing candidate to the same
    score.

    Grading formula (when corroborated and not suppressed):
        penalty = -(1 - t_collapse / T)
    where ``t_collapse`` is ``first_corroborated_time`` and ``T`` is the
    ``target_stop_time`` (or ``final_time`` as fallback).  A collapse at t=0
    earns the full -1.0; a collapse at the end of the run earns ~0.0.
    """
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

    # Hard suppression for very late collapse (likely numerical end-of-run
    # artifact rather than physical event).
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

    # Grade the penalty by how early the collapse occurs.  A longer survival
    # time before collapse earns a smaller penalty, giving the optimizer a
    # smooth gradient to delay collapse rather than a binary floor.
    reference_time = target_stop_time if target_stop_time and target_stop_time > 0.0 else final_time
    if reference_time > 0.0 and first_time >= 0.0:
        survival_fraction = min(first_time / reference_time, 1.0)
        graded_penalty = -(1.0 - survival_fraction)
        notes.append(
            f"graded horizon penalty: collapse at t={first_time:.2f} / "
            f"{reference_time:.2f} (survived {survival_fraction:.1%}); "
            f"penalty={graded_penalty:.3f}"
        )
        return graded_penalty, notes

    # Fallback: no timing information available, use binary penalty.
    return -min(horizon, 1.0), notes
