from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from .episode_metrics import STATIONARY_BETA_EPS, EpisodeMetrics


DEFAULT_WEIGHTS: dict[str, float] = {
    # t=0 FTL signals are now only faint "looked promising at t=0" hints -- they
    # are downweighted hard because a t=0 shortcut that does not survive the
    # evolution is a gauge/initial-data artifact (see operational_ftl, which
    # rewards the *evolved* advantage instead).
    "ftl_shortcut": 1.0,
    "expansion_asymmetry": 0.75,
    "nonflat_geometry": 1.0,
    "comoving_stability": 2.5,
    "survival": 1.5,
    "constraint_health": 2.0,
    "lapse_health": 1.0,
    "horizon_penalty": 1.5,
    "energy_condition": 2.0,
    "stability": 0.5,
    "instability_penalty": 8.0,
    "nontrivial_geometry": 0.25,
    "initial_constraint_quality": 0.5,
    # Proposed extensions (Sec. "Proposed Extensions").
    "constraint_growth": 2.0,
    "anec_condition": 1.5,
    "tidal_comfort": 1.0,
    # Mechanism-agnostic operational FTL, measured on the EVOLVED spacetime
    # (Dijkstra arrival-time advantage that survives the dynamics) -- this is the
    # primary goal and dominates the reward budget.  A t=0-only shortcut scores
    # nothing here.  Curvature activity rewards genuinely non-trivial geometry
    # (and keeps flat space out of the running).
    "operational_ftl": 3.0,
    "operational_ftl_geodesic": 12.0,
    # Continuous FTL *precursor* (shaping gradient).  operational_ftl is a hard
    # gate -- it stays exactly 0 until a connected superluminal channel beats
    # the flat baseline end-to-end -- so an "almost there" geometry (light cones
    # tilted, a few locally superluminal cells, but no through-channel yet) gets
    # no signal and the optimizer cannot climb toward FTL.  The precursor rewards
    # the cone-tilting itself (max local coordinate light speed past 1 + the
    # superluminal-cell fraction), which rises smoothly *before* the gate fires.
    # Flat space has max_local_speed == 1 exactly, so the precursor is 0 in
    # flatness: it is a clean slope out of the flat-space basin toward FTL.
    "ftl_precursor": 3.0,
    # Direct frame-dragging / light-cone-tilt drive.  This gives the optimizer a
    # signal before the end-to-end FTL path exists, but unlike stability terms it
    # points at the mechanism we actually need.
    "shift_drive": 2.0,
    # Coupled path-progress shaping: rewards evolved geometries whose fastest
    # causal path is close to beating the flat baseline *and* carry both cone
    # opening (precursor) and frame-drag (shift) together -- not either alone.
    "channel_progress": 2.5,
    # Constraint-satisfying GRTresna geometry at t=0 (mechanism-agnostic Dijkstra).
    # Shaping gradient before GPU evolution; evolved operational_ftl remains primary.
    "operational_ftl_solved": 4.0,
    "curvature_activity": 0.5,
    # Exotic-matter PENALTY.  The objective is FTL *without* exotic matter, so
    # any negative-energy requirement is punished -- both the matter sourced at
    # t=0 by the constrained solve and the geometric exotic energy of the
    # evolved spacetime (effective NEC violation, T^eff = G/8pi).  The component
    # is negative; this weight sets how hard exoticity is penalized.
    "exotic_penalty": 8.0,
    "qei_penalty": 4.0,
    "boundary_penalty": 6.0,
    "transport_objective": 2.0,
}


# "Health/niceness" rewards that a trivial flat spacetime maximises by default.
# They are gated by a non-triviality factor (see ``score_episode``) so flat
# Minkowski cannot bank them and out-rank a genuinely exotic geometry.
HEALTH_COMPONENTS: frozenset[str] = frozenset({
    "survival",
    "constraint_health",
    "energy_condition",
    "initial_constraint_quality",
    "lapse_health",
    "stability",
    "comoving_stability",
    "constraint_growth",
    "anec_condition",
    "tidal_comfort",
})


@dataclass(frozen=True)
class Score:
    total: float
    components: dict[str, float]
    notes: list[str]


def _bounded_reward(value: float, scale: float) -> float:
    if not math.isfinite(value) or value < 0:
        return 0.0
    return 1.0 / (1.0 + value / scale)


def domain_half_width_from_overrides(
    overrides: Mapping[str, object] | None,
) -> float | None:
    """Half-width of the cubic domain (``L_full / 2``) for the horizon-finder
    sanity guard.  The domain is ``[0, L_full]`` with the geometry centered at
    ``L_full / 2``, so a genuine interior trapped surface must satisfy
    ``r << L_full / 2``."""
    if not overrides:
        return None
    try:
        l_full = float(overrides.get("L_full"))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return 0.5 * l_full if l_full > 0.0 else None


# A trapped surface is interior, so the minimum-expansion radius must sit well
# inside the domain.  If it lands beyond this fraction of the domain half-width
# the apparent-horizon proxy is miscentered (corner origin) or reading boundary
# noise, and its trapped verdict is rejected rather than silently vetoing the
# candidate.
HORIZON_OFFCENTER_FRACTION: float = 0.5


def score_episode(
    metrics: EpisodeMetrics,
    *,
    target_stop_time: float | None = None,
    weights: Mapping[str, float] | None = None,
    objective_mode: str = "weighted",
    domain_half_width: float | None = None,
) -> Score:
    w = dict(DEFAULT_WEIGHTS)
    if weights:
        w.update(weights)

    components: dict[str, float] = {}
    notes: list[str] = []

    final_time = None
    if metrics.collapse and metrics.collapse.final_time is not None:
        final_time = metrics.collapse.final_time
    elif metrics.constraints and metrics.constraints.final_time is not None:
        final_time = metrics.constraints.final_time

    if target_stop_time and target_stop_time > 0 and final_time is not None:
        components["survival"] = min(final_time / target_stop_time, 1.0)
    elif final_time is not None:
        components["survival"] = 1.0
    else:
        components["survival"] = 0.0
        notes.append("no time-series diagnostics were found")

    if metrics.constraints:
        ham = metrics.constraints.max_hamiltonian_l2 or 0.0
        mom = metrics.constraints.max_momentum_l2 or 0.0
        components["constraint_health"] = (
            0.5 * _bounded_reward(ham, 1.0e-2)
            + 0.5 * _bounded_reward(mom, 1.0e-2)
        )

        int_neg = metrics.constraints.integral_negative_rho
        min_rho = metrics.constraints.min_rho_required
        if int_neg is not None and min_rho is not None:
            if min_rho >= 0.0 and int_neg <= 1.0e-12:
                components["energy_condition"] = 1.0
            else:
                components["energy_condition"] = _bounded_reward(
                    int_neg if int_neg > 0 else abs(min_rho), 1.0e-2
                )
        else:
            components["energy_condition"] = 0.0
            notes.append("energy density diagnostics not available")

        if metrics.constraints.final_time is not None and final_time is not None:
            first_ham = metrics.constraints.max_hamiltonian_l2 or 0.0
            components["initial_constraint_quality"] = _bounded_reward(
                first_ham, 1.0e-3
            )
        else:
            components["initial_constraint_quality"] = 0.0
    else:
        components["constraint_health"] = 0.0
        components["energy_condition"] = 0.0
        components["initial_constraint_quality"] = 0.0
        notes.append("constraint_norms.dat missing")

    if metrics.collapse:
        lapse = metrics.collapse.min_lapse or 0.0
        horizon = metrics.collapse.max_horizon_radius or 0.0
        k_activity = metrics.collapse.max_abs_k or 0.0
        scalar_activity = max(
            metrics.collapse.scalar_phi_range or 0.0,
            metrics.collapse.scalar_pi_range or 0.0,
        )

        components["lapse_health"] = min(max(lapse / 1.0e-3, 0.0), 1.0)

        r_at_min_theta = metrics.collapse.r_at_min_theta_plus
        horizon_offcenter = (
            horizon > 0.0
            and domain_half_width is not None
            and domain_half_width > 0.0
            and r_at_min_theta is not None
            and r_at_min_theta > HORIZON_OFFCENTER_FRACTION * domain_half_width
        )
        if horizon_offcenter:
            components["horizon_penalty"] = 0.0
            notes.append(
                "horizon proxy located at r="
                f"{r_at_min_theta:.1f} > {HORIZON_OFFCENTER_FRACTION:g}*half-width "
                f"({domain_half_width:.1f}); rejected as miscentered/boundary "
                "artifact (no interior trapped surface), horizon penalty suppressed"
            )
        else:
            components["horizon_penalty"] = -min(horizon, 1.0)
        components["nontrivial_geometry"] = min(
            math.log1p(k_activity + scalar_activity), 1.0
        )
    else:
        components["lapse_health"] = 0.0
        components["horizon_penalty"] = 0.0
        components["nontrivial_geometry"] = 0.0
        notes.append("collapse_diagnostics.dat missing")

    if metrics.stability and metrics.stability.violation is not None:
        components["stability"] = _bounded_reward(metrics.stability.violation, 1.0)
        components["instability_penalty"] = -(1.0 - components["stability"])
        if components["stability"] < 0.25:
            notes.append("geometry changes rapidly over the evolution window (Eulerian)")
    else:
        components["stability"] = 0.0
        components["instability_penalty"] = 0.0
        notes.append("stability diagnostics not available")

    if metrics.comoving and metrics.comoving.stationary:
        components["comoving_stability"] = components["stability"]
        notes.append("co-moving stability uses Eulerian fallback (stationary geometry)")
    elif metrics.comoving and metrics.comoving.score is not None:
        components["comoving_stability"] = metrics.comoving.score
    else:
        components["comoving_stability"] = 0.0
        notes.append("co-moving stability diagnostics not available")

    if metrics.growth is not None and metrics.growth.s_growth is not None:
        components["constraint_growth"] = metrics.growth.s_growth
        if metrics.growth.lambda_effective is not None and metrics.growth.lambda_effective > 0.5:
            notes.append(
                "constraint/collapse series grow exponentially "
                f"(lambda={metrics.growth.lambda_effective:.3f}); slow-collapse penalized"
            )
    else:
        components["constraint_growth"] = 0.0

    if metrics.physical is not None:
        if metrics.physical.s_anec is not None:
            components["anec_condition"] = metrics.physical.s_anec
        else:
            components["anec_condition"] = 0.0
        if metrics.physical.s_tidal is not None:
            components["tidal_comfort"] = metrics.physical.s_tidal
        else:
            components["tidal_comfort"] = 0.0
        for note in metrics.physical.notes:
            notes.append(note)
    else:
        components["anec_condition"] = 0.0
        components["tidal_comfort"] = 0.0

    if metrics.ftl is not None:
        components["ftl_shortcut"] = metrics.ftl.f_log
        components["expansion_asymmetry"] = metrics.ftl.f_asymmetry
        components["nonflat_geometry"] = metrics.ftl.s_nonflat
        if metrics.ftl.f_shortcut <= 0.0 and metrics.ftl.s_nonflat < 0.05:
            notes.append("no FTL shortcut detected in t=0 profile")
        for note in metrics.ftl.notes:
            notes.append(note)
    else:
        components["ftl_shortcut"] = 0.0
        components["expansion_asymmetry"] = 0.0
        components["nonflat_geometry"] = 0.0
        notes.append("FTL profile metrics not available")

    # Mechanism-agnostic operational FTL.  CRITICAL: we reward the EVOLVED
    # arrival-time advantage (F_op^ev), not the t=0 slice.  A large t=0 shortcut
    # that the dynamics relax away (a gauge / initial-data artifact, max c -> 1)
    # is worthless; only a channel that *survives* the evolution counts.  The
    # t=0 value is kept solely as a diagnostic and as a faint "looked promising"
    # hint (small ftl_shortcut weight), never as the main FTL reward.
    f_op_t0 = (
        metrics.general_ftl.f_op
        if metrics.general_ftl is not None
        else 0.0
    )
    f_op_ev = (
        metrics.general_ftl_evolved.f_op
        if metrics.general_ftl_evolved is not None
        else 0.0
    )
    f_op_solved = (
        metrics.general_ftl_solved.f_op
        if metrics.general_ftl_solved is not None
        else 0.0
    )
    SOLVED_FTL_SCALE = 3.0e-3
    components["operational_ftl_solved"] = (
        min(math.log1p(f_op_solved / SOLVED_FTL_SCALE), 1.0)
        if math.isfinite(f_op_solved) and f_op_solved > 0
        else 0.0
    )
    if metrics.mechanism_descriptor is not None:
        components["mechanism_descriptor"] = float(
            min(max(metrics.mechanism_descriptor, 0.0), 1.0)
        )
    # A tiny coordinate-time win is useful as a shaping signal, but it should not
    # dominate the search.  The previous log scale saturated at F_op~0.01 and let
    # quasi-stationary "warp lens" artifacts score as solved FTL.  Reserve the
    # hard operational reward for a materially stronger evolved shortcut.
    OP_FTL_FLOOR = 3.0e-2
    OP_FTL_TARGET = 1.0e-1
    if math.isfinite(f_op_ev) and f_op_ev > OP_FTL_FLOOR:
        components["operational_ftl"] = min(
            (f_op_ev - OP_FTL_FLOOR) / (OP_FTL_TARGET - OP_FTL_FLOOR),
            1.0,
        )
    else:
        components["operational_ftl"] = 0.0

    GEO_FTL_FLOOR = 1.0e-3
    GEO_FTL_TARGET = 5.0e-2
    geo_report = metrics.geodesic_ftl
    f_geo = geo_report.f_geo if geo_report is not None else 0.0
    # A gauge-invariant shortcut is only trustworthy when the null-ray
    # integration stayed on the constraint surface (``h_quality_ok``) AND the
    # whole ray bundle reached the detector.  A high Hamiltonian-constraint
    # drift or a partial bundle means ``f_geo`` is integration noise / a caustic
    # rather than a certified shortcut, so it must never earn the dominant FTL
    # reward (this is what let a stationary warp-lens artifact top the ranking).
    geo_trustworthy = bool(
        geo_report is not None
        and geo_report.h_quality_ok
        and geo_report.n_rays > 0
        and geo_report.n_reached == geo_report.n_rays
    )
    if geo_trustworthy and math.isfinite(f_geo) and f_geo > GEO_FTL_FLOOR:
        components["operational_ftl_geodesic"] = min(
            (f_geo - GEO_FTL_FLOOR) / (GEO_FTL_TARGET - GEO_FTL_FLOOR),
            1.0,
        )
    else:
        components["operational_ftl_geodesic"] = 0.0
    if geo_report is not None:
        if f_geo > GEO_FTL_FLOOR and not geo_trustworthy:
            notes.append(
                "geodesic shortcut rejected as unreliable "
                f"(h_quality_ok={geo_report.h_quality_ok}, "
                f"rays {geo_report.n_reached}/{geo_report.n_rays}, "
                f"max H={geo_report.max_h_drift:.2e}, f_geo={f_geo:.3e})"
            )
        elif f_geo <= 0.0 and f_op_ev > OP_FTL_FLOOR:
            notes.append(
                "coordinate FTL channel is a gauge artifact "
                f"(F_op^ev={f_op_ev:.3e}, f_geo={f_geo:.3e})"
            )
        elif geo_trustworthy and f_geo > GEO_FTL_FLOOR:
            notes.append(
                f"gauge-invariant null-geodesic shortcut confirmed (f_geo={f_geo:.3e})"
            )
    # Sustained evolved FTL.  Preferred signal: the worst-case operational
    # shortcut across the last few retained plotfiles (needs consumer_keep_last
    # >= 2), normalized like operational_ftl.  This rewards a channel that holds
    # over the window and structurally rejects a one-frame gauge spike (whose
    # window minimum collapses to ~0).  Falls back to the old ratio diagnostic
    # (evolved/t=0) when no multi-frame window survived.
    if metrics.ftl_persistence is not None and metrics.ftl_persistence.f_op_min is not None:
        f_op_sustained = float(metrics.ftl_persistence.f_op_min)
        if math.isfinite(f_op_sustained) and f_op_sustained > OP_FTL_FLOOR:
            components["ftl_persistence"] = min(
                (f_op_sustained - OP_FTL_FLOOR) / (OP_FTL_TARGET - OP_FTL_FLOOR),
                1.0,
            )
        else:
            components["ftl_persistence"] = 0.0
    else:
        components["ftl_persistence"] = (
            float(min(max(f_op_ev / f_op_t0, 0.0), 1.0)) if f_op_t0 > 1.0e-9 else 0.0
        )
    stationary_geometry = bool(metrics.comoving and metrics.comoving.stationary)
    # The stationary-artifact verdict (and the gating of frozen-coordinate
    # shaping rewards) is deferred until ftl_precursor / channel_progress /
    # shift_drive have been computed below; see the "Stationary-artifact gate".
    if metrics.general_ftl_evolved is not None:
        c_ev = metrics.general_ftl_evolved.max_local_speed
        if c_ev > 1.0 and components["operational_ftl"] > 0.0:
            notes.append(
                "evolved geometry sustains a strong superluminal channel "
                f"(max c = {c_ev:.3f}, F_op^ev = {f_op_ev:.3e}); operational FTL"
            )
        elif c_ev > 1.0 and f_op_ev > 0.0:
            notes.append(
                "weak evolved coordinate shortcut treated as channel progress, "
                f"not solved FTL (max c = {c_ev:.3f}, F_op^ev = {f_op_ev:.3e})"
            )
        elif f_op_t0 > 0.0:
            notes.append(
                f"t=0 shortcut (F_op0={f_op_t0:.3e}, max c={metrics.general_ftl.max_local_speed:.3f}) "
                f"did not persist (F_op^ev={f_op_ev:.3e}); likely gauge/initial-data artifact"
            )
    elif metrics.general_ftl is not None and metrics.general_ftl.f_op <= 0.0:
        notes.append("no operational FTL shortcut on the t=0 slice")

    # Coordinate-invariant curvature activity: a bounded reward so flat space
    # scores ~0 while a structured warp/wormhole throat earns credit.  Keep the
    # scale deliberately broad; the previous log1p(x)->cap at 1.0 saturated the
    # whole GRTresna shell population and erased the precursor gradient.
    if metrics.curvature is not None and metrics.curvature.max_l2_ricci_scalar is not None:
        curvature = max(0.0, metrics.curvature.max_l2_ricci_scalar)
        CURVATURE_ACTIVITY_SCALE = 5.0
        components["curvature_activity"] = curvature / (curvature + CURVATURE_ACTIVITY_SCALE)
    else:
        components["curvature_activity"] = 0.0

    # FTL PRECURSOR -- the continuous shaping gradient that operational_ftl
    # lacks.  Two smooth signals that rise *before* a connected channel forms:
    #   * speed_term: how the fastest local coordinate light speed climbs toward
    #     and past the c=1 threshold (cones tilting superluminal), and
    #   * frac_term: what fraction of the slice is locally superluminal (the
    #     nascent channel growing toward an end-to-end shortcut).
    # We prefer the EVOLVED slice (the tilt must survive the dynamics) but fall
    # back to a half-credit t=0 reading so there is still a gradient before any
    # plotfile is available.
    #
    # Sub-luminal approach: the solved throats sit just *below* c=1 (they slow
    # light rather than speed it), so the old "only reward c>1" precursor was
    # dead-flat zero across the whole population -- leaving the search no FTL
    # direction and collapsing it onto stability.  Below the threshold we add a
    # structure-gated ramp toward c=1, multiplied by curvature activity so flat
    # space (c==1 trivially, zero curvature) earns nothing, and constructed so
    # the reward stays continuous and monotonic as a structured geometry crosses
    # into the genuinely superluminal (c>1) regime.
    PRECURSOR_SPEED_SCALE = 0.05  # (c - 1) reward saturation scale
    PRECURSOR_FRAC_SCALE = 0.05   # superluminal area-fraction reference
    PRECURSOR_SUBLUMINAL_FLOOR = 0.9   # below this c the geometry is too slow to count
    PRECURSOR_SUBLUMINAL_WEIGHT = 0.5  # max sub-luminal credit reached at the c=1 threshold

    structure = components["curvature_activity"]

    def _precursor(report) -> float:
        if report is None:
            return 0.0
        c = float(getattr(report, "max_local_speed", 0.0) or 0.0)
        frac = float(getattr(report, "superluminal_fraction", 0.0) or 0.0)
        sub_credit = PRECURSOR_SUBLUMINAL_WEIGHT * structure
        if c >= 1.0:
            over = c - 1.0
            superluminal = min(math.log1p(over / PRECURSOR_SPEED_SCALE), 1.0)
            speed_term = sub_credit + (1.0 - sub_credit) * superluminal
        elif c > PRECURSOR_SUBLUMINAL_FLOOR:
            ramp = (c - PRECURSOR_SUBLUMINAL_FLOOR) / (1.0 - PRECURSOR_SUBLUMINAL_FLOOR)
            speed_term = sub_credit * ramp
        else:
            speed_term = 0.0
        frac_term = min(frac / PRECURSOR_FRAC_SCALE, 1.0) if frac > 0.0 else 0.0
        return 0.7 * speed_term + 0.3 * frac_term

    precursor_ev = _precursor(metrics.general_ftl_evolved)
    precursor_t0 = _precursor(metrics.general_ftl)
    precursor_solved = _precursor(metrics.general_ftl_solved)
    components["ftl_precursor"] = max(
        precursor_ev, 0.5 * precursor_t0, 0.5 * precursor_solved
    )

    def _shift_drive(report) -> float:
        if report is None:
            return 0.0
        max_shift = float(getattr(report, "max_shift", 0.0) or 0.0)
        if not math.isfinite(max_shift) or max_shift <= 0.0:
            return 0.0
        SHIFT_DRIVE_SCALE = 0.25
        return min(math.log1p(max_shift / SHIFT_DRIVE_SCALE), 1.0)

    components["shift_drive"] = max(
        _shift_drive(metrics.general_ftl_evolved),
        0.5 * _shift_drive(metrics.general_ftl),
        0.5 * _shift_drive(metrics.general_ftl_solved),
    )

    # CHANNEL PROGRESS -- coupled shaping gradient that pushes the search away
    # from isolated precursor-only (eval 128) or shift-only (eval 57) basins.
    # Path closeness uses relative excess (t_min - t_flat) / t_flat so flat
    # space (t_min == t_flat) earns nothing unless precursor/shift are nonzero.
    CHANNEL_PATH_EXCESS_SCALE = 0.12

    def _path_closeness(report) -> float:
        if report is None or not getattr(report, "reachable", False):
            return 0.0
        t_min = getattr(report, "t_min", None)
        t_flat = float(getattr(report, "t_flat", 0.0) or 0.0)
        if t_min is None or not math.isfinite(t_min) or t_flat <= 0.0:
            return 0.0
        if t_min <= t_flat:
            return 1.0
        excess = (t_min - t_flat) / t_flat
        return max(0.0, 1.0 - excess / CHANNEL_PATH_EXCESS_SCALE)

    def _channel_progress(report) -> float:
        path = _path_closeness(report)
        if path <= 0.0:
            return 0.0
        prec = _precursor(report)
        shift = _shift_drive(report)
        mechanism = math.sqrt(max(prec, 0.0) * max(shift, 0.0))
        return path * mechanism

    channel_ev = _channel_progress(metrics.general_ftl_evolved)
    channel_t0 = _channel_progress(metrics.general_ftl)
    channel_solved = _channel_progress(metrics.general_ftl_solved)
    components["channel_progress"] = max(
        channel_ev, 0.5 * channel_t0, 0.5 * channel_solved
    )
    if (
        components["channel_progress"] > 0.05
        and components["operational_ftl"] <= 0.0
    ):
        notes.append(
            "coupled channel progress (path + precursor + shift): "
            f"{components['channel_progress']:.3f}"
        )

    if components["ftl_precursor"] > 0.05 and components["operational_ftl"] <= 0.0:
        notes.append(
            "FTL precursor active (cones climbing toward the c=1 threshold): "
            f"{components['ftl_precursor']:.3f}"
        )

    # Stationary-artifact gate.  A geometry with no mean shift in the bubble
    # (``comoving.stationary``) has no frame-dragging mechanism, so any
    # superluminal coordinate speed is a static "warp-lens" / gauge artifact,
    # not a propagating channel.  Unless such a geometry shows a *dynamical* FTL
    # signal we trust -- an evolved coordinate shortcut, a sustained channel, or
    # a reliability-gated gauge-invariant geodesic shortcut -- its FTL "signals"
    # (t=0 solved Dijkstra, cone-tilt precursor, channel progress, shift drive)
    # are frozen coordinate features rather than progress toward FTL.  Zeroing
    # those shaping rewards (instead of relying on a beatable additive penalty)
    # is what stops a static artifact from climbing the ranking; genuine
    # shift-driven candidates have ``beta_mean != 0`` and are never flagged
    # stationary, so their gradient is untouched.
    trustworthy_dynamical_ftl = (
        components["operational_ftl"] > 0.0
        or components["operational_ftl_geodesic"] > 0.0
        or components["ftl_persistence"] > 0.0
    )
    stationary_artifact = stationary_geometry and not trustworthy_dynamical_ftl
    if stationary_artifact:
        for key in (
            "operational_ftl_solved",
            "ftl_precursor",
            "channel_progress",
            "shift_drive",
        ):
            components[key] = 0.0
        notes.append(
            "stationary zero-shift geometry with no trustworthy dynamical FTL: "
            "coordinate-artifact shaping rewards gated out"
        )
    # Graded stationary penalty.  A *flat* -1.0 collapses every zero-net-shift
    # geometry to the same score, so the map has no slope and cannot tell
    # "almost propulsive" from "perfectly static" -- the whole stationary basin
    # is a flat floor the search cannot climb out of.  Instead we scale the
    # penalty by how much *coherent net axial shift* the bubble carries
    # (|beta_mean|, the directed frame-drag that an actual warp needs),
    # normalized to the stationarity threshold.  A truly static lens
    # (beta_mean -> 0) still earns the full -1.0, so the eval_083 artifact fix
    # is preserved, but a geometry developing directed shift sees the penalty
    # smoothly relax toward 0 as it approaches the threshold -- a continuous
    # downhill gradient out of the basin toward genuine shift-driven FTL.
    beta_mean_mag = (
        abs(metrics.comoving.beta_mean)
        if metrics.comoving is not None and metrics.comoving.beta_mean is not None
        else 0.0
    )
    shift_coherence = (
        min(beta_mean_mag / STATIONARY_BETA_EPS, 1.0)
        if STATIONARY_BETA_EPS > 0.0
        else 0.0
    )
    components["stationary_artifact_penalty"] = (
        -(1.0 - shift_coherence) if stationary_artifact else 0.0
    )

    # Exotic-matter penalty.  The goal is FTL *without* exotic matter, so any
    # negative-energy requirement is penalized.  Two independent probes feed it:
    #   1. matter sector at t=0: the constrained solve back-solves the matter
    #      from the proposed geometry; a non-zero integral of negative rho (or a
    #      negative peak min_rho_required) means exotic matter is needed up front.
    #   2. evolved geometry: the effective stress-energy T^eff = G/8pi of the
    #      evolved spacetime violating the NEC (geometric exotic energy that the
    #      matter-sector EC cannot see -- e.g. a shift-driven warp).
    #
    # Each probe is log-scaled against its own reference (the magnitudes differ
    # by orders of magnitude: matter integral ~O(1), geometric NEC ~O(1e-4)) so
    # the penalty is a *smooth gradient* across the observed range rather than a
    # saturated constant.  Steepest near zero, so "almost no exotic" is clearly
    # distinguished from "a little exotic"; heavy exotic saturates at -1.
    def _graded(x: float, scale: float, ceiling: float = 1.0) -> float:
        if not math.isfinite(x) or x <= 0.0:
            return 0.0
        return min(math.log1p(x / scale), ceiling)

    MATTER_INT_SCALE = 0.5     # integral of negative rho (volume measure)
    MATTER_PEAK_SCALE = 5.0e-2  # peak negative energy density
    GEO_NEC_SCALE = 5.0e-4      # effective-NEC violation of evolved geometry
    # The persistent-FTL winners all live in the int_neg ~ 0.5-2 band, where a
    # ceiling of 1.0 saturates and can no longer tell "a bit exotic" from "very
    # exotic".  Raise the matter ceiling so the penalty keeps a usable gradient
    # right through that band -- this is what lets the search minimise exotic
    # content *among* the geometries that actually sustain a channel.
    MATTER_CEILING = 1.6

    matter_exotic = 0.0
    geo_exotic = 0.0
    if metrics.constraints is not None:
        int_neg = metrics.constraints.integral_negative_rho
        min_rho = metrics.constraints.min_rho_required
        if int_neg is not None and int_neg > 0.0:
            matter_exotic = max(
                matter_exotic, _graded(int_neg, MATTER_INT_SCALE, MATTER_CEILING)
            )
        if min_rho is not None and min_rho < 0.0:
            matter_exotic = max(
                matter_exotic, _graded(-min_rho, MATTER_PEAK_SCALE, MATTER_CEILING)
            )
    if metrics.effective_ec is not None:
        eff = metrics.effective_ec
        eff_candidates = [
            v for v in (eff.nec_slack_min, eff.nec_min, eff.rho_eulerian_min)
            if v is not None
        ]
        eff_worst = min(eff_candidates) if eff_candidates else 0.0
        geo_exotic = _graded(-eff_worst, GEO_NEC_SCALE)

    exotic_term = max(matter_exotic, geo_exotic)
    components["exotic_penalty"] = -exotic_term
    if exotic_term > 0.05:
        notes.append(
            "exotic matter required "
            f"(matter={matter_exotic:.2f}, geometric={geo_exotic:.2f} of full penalty)"
        )

    qei_violation = 0.0
    if metrics.qei is not None:
        if metrics.qei.trajectory_violation is not None and metrics.qei.trajectory_violation > 0:
            qei_violation = max(qei_violation, metrics.qei.trajectory_violation)
        if metrics.qei.spatial_proxy is not None and metrics.qei.spatial_proxy > 0:
            qei_violation = max(qei_violation, metrics.qei.spatial_proxy)
    elif metrics.physical is not None and metrics.physical.qei_spatial_proxy is not None:
        qei_violation = max(0.0, metrics.physical.qei_spatial_proxy)
    components["qei_penalty"] = -min(math.log1p(qei_violation / 1.0e-2), 1.0) if qei_violation > 0 else 0.0

    if metrics.boundary_flux is not None and metrics.boundary_flux.reflection_contaminated:
        components["boundary_penalty"] = -1.0
        notes.append("late-time inward scalar flux at boundary (reflection contamination)")
    else:
        components["boundary_penalty"] = 0.0

    if metrics.transport is not None and metrics.transport.score is not None:
        components["transport_objective"] = metrics.transport.score
    else:
        components["transport_objective"] = 0.0

    # Non-triviality gate: a trivial flat spacetime aces every "health" reward
    # (survival, stability, clean constraints, no exotic energy) while scoring
    # zero on all exoticity terms, which previously let it out-rank genuine warp
    # / wormhole geometries (the flat-space attractor).  We gate the health
    # rewards by how non-trivial the geometry actually is, so they only count
    # once there is real structure to be healthy about.  The exoticity/FTL terms
    # are never gated -- they supply the gradient out of flatness.
    # Note: exotic energy is deliberately NOT a non-triviality signal anymore --
    # we do not want negative-energy content to "earn" the health rewards.  FTL
    # and curvature structure supply the gradient out of flatness instead.
    nontriviality = max(
        components.get("nonflat_geometry", 0.0),
        components.get("expansion_asymmetry", 0.0),
        components.get("nontrivial_geometry", 0.0),
        components.get("curvature_activity", 0.0),
        components.get("shift_drive", 0.0),
        components.get("channel_progress", 0.0),
        components.get("operational_ftl", 0.0),
        components.get("operational_ftl_geodesic", 0.0),
        components.get("ftl_persistence", 0.0),
        components.get("ftl_precursor", 0.0),
        components.get("operational_ftl_solved", 0.0),
        components.get("ftl_shortcut", 0.0),
    )
    nontriviality = float(min(max(nontriviality, 0.0), 1.0))
    components["nontriviality_gate"] = nontriviality
    if nontriviality < 0.05:
        notes.append(
            "near-trivial geometry: health rewards gated out (flat-space attractor guard)"
        )

    if objective_mode == "ftl_first":
        # Lexicographic-style scalarization: path-level FTL dominates; precursor
        # and shift are shaping only.  Trapped-surface proxies (horizon_penalty
        # < 0) veto high local-speed artifacts that are not traversable channels.
        health_gate = components.get("nontriviality_gate", 0.0)
        horizon = components.get("horizon_penalty", 0.0)
        trapped_surface = horizon < -0.05
        total = (
            1500.0 * components.get("operational_ftl_geodesic", 0.0)
            + 400.0 * components.get("operational_ftl", 0.0)
            + 300.0 * components.get("ftl_persistence", 0.0)
            + 350.0 * components.get("channel_progress", 0.0)
            + 180.0 * components.get("operational_ftl_solved", 0.0)
            + 60.0 * components.get("ftl_precursor", 0.0)
            + 40.0 * components.get("shift_drive", 0.0)
            + 50.0 * components.get("ftl_shortcut", 0.0)
            + 5.0 * components.get("nontrivial_geometry", 0.0)
            + health_gate * (
                20.0 * components.get("survival", 0.0)
                + 3.0 * components.get("constraint_health", 0.0)
                + 2.0 * components.get("stability", 0.0)
                + 3.0 * components.get("instability_penalty", 0.0)
                + 2.0 * components.get("comoving_stability", 0.0)
                + 1.0 * components.get("energy_condition", 0.0)
            )
            + 1.0 * components.get("exotic_penalty", 0.0)
            # Moderate stationary penalty: the FTL shaping rewards are already
            # zeroed for a stationary artifact (see the stationary-artifact gate)
            # and the geodesic shortcut is reliability-gated, so a static lens
            # can no longer masquerade as FTL.  The penalty therefore only needs
            # to keep a shift-free geometry ranked below a genuine shift-driven
            # one -- a catastrophic weight (it used to be 250) instead floored the
            # entire population negative and erased the QD gradient.  At this
            # weight a healthy static structure still scores mildly positive
            # while any real shift (graded by |beta|) lifts it further.
            + 8.0 * components.get("stationary_artifact_penalty", 0.0)
            + 500.0 * horizon
        )
        if trapped_surface:
            notes.append(
                f"trapped-surface proxy active (horizon_penalty={horizon:.3f}); "
                "local precursor/shift alone cannot rank this candidate highly"
            )
        notes.append(
            "objective_mode=ftl_first: channel/operational FTL dominate precursor/shift"
        )
    else:
        total = 0.0
        for key, value in components.items():
            if key == "nontriviality_gate":
                continue
            gate = nontriviality if key in HEALTH_COMPONENTS else 1.0
            total += w.get(key, 0.0) * value * gate
    return Score(total=total, components=components, notes=notes)
