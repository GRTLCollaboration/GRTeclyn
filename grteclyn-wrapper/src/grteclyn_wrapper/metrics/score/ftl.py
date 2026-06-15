from __future__ import annotations

import math

from .types import ScoringContext

SOLVED_FTL_SCALE = 3.0e-3
# Peak must clear c by >SOLVED_PEAK_FLOOR to earn any solved-FTL credit,
# full credit by SOLVED_PEAK_TARGET.
SOLVED_PEAK_FLOOR = 0.03
SOLVED_PEAK_TARGET = 0.20
# A superluminal region covering more than this fraction of the slice reads
# as a global coordinate/lapse offset, not a localized channel.
SOLVED_LOCALITY_CEILING = 0.5
OP_FTL_FLOOR = 3.0e-2
OP_FTL_TARGET = 1.0e-1
GEO_FTL_FLOOR = 1.0e-3
# Full marks require a genuinely dramatic gauge-invariant shortcut (~20% null
# arrival-time advantage over flat space).  A marginal few-percent shortcut --
# real and reliability-gated, but physically modest in a near-flat geometry --
# therefore earns only a small fraction here, so the scalar score tracks the
# *magnitude* of the shortcut instead of saturating the moment the floor is
# crossed.
GEO_FTL_TARGET = 2.0e-1
PRECURSOR_SPEED_SCALE = 0.15  # (c - 1) reward saturation scale
PRECURSOR_FRAC_SCALE = 0.15   # superluminal area-fraction reference
PRECURSOR_SUBLUMINAL_FLOOR = 0.9   # below this c the geometry is too slow to count
PRECURSOR_SUBLUMINAL_WEIGHT = 0.5  # max sub-luminal credit reached at the c=1 threshold
CHANNEL_PATH_EXCESS_SCALE = 0.12


def _geo_magnitude(value: float) -> float:
    return min(
        (value - GEO_FTL_FLOOR) / (GEO_FTL_TARGET - GEO_FTL_FLOOR),
        1.0,
    )


def _precursor(report, *, structure: float) -> float:
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


def _shift_drive(report) -> float:
    if report is None:
        return 0.0
    max_shift = float(getattr(report, "max_shift", 0.0) or 0.0)
    if not math.isfinite(max_shift) or max_shift <= 0.0:
        return 0.0
    shift_drive_scale = 0.25
    return min(math.log1p(max_shift / shift_drive_scale), 1.0)


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


def _channel_progress(report, *, structure: float) -> float:
    path = _path_closeness(report)
    if path <= 0.0:
        return 0.0
    prec = _precursor(report, structure=structure)
    shift = _shift_drive(report)
    mechanism = math.sqrt(max(prec, 0.0) * max(shift, 0.0))
    return path * mechanism


def compute_ftl_components(ctx: ScoringContext) -> None:
    metrics = ctx.metrics
    components = ctx.components
    notes = ctx.notes
    structural_persistence = ctx.structural_persistence

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
    solved_ftl_raw = (
        min(math.log1p(f_op_solved / SOLVED_FTL_SCALE), 1.0)
        if math.isfinite(f_op_solved) and f_op_solved > 0
        else 0.0
    )
    # Localization + peak-margin gate.  A near-uniform, marginally-superluminal
    # coordinate field (every cell at ~1.01, max c barely > 1) trivially fills
    # superluminal_fraction -> 1.0 and saturates this reward, yet it is a global
    # gauge/lapse offset, not a traversable warp channel.  A genuine channel is
    # *localized* (small superluminal fraction) and has a real peak above c, so
    # scale the reward by how far the peak clears c AND by how localized the
    # superluminal region is.  A sub-floor peak or a fraction >= the ceiling
    # collapses the reward to ~0; a tight high-speed bubble keeps it.
    solved_report = metrics.general_ftl_solved
    if solved_ftl_raw > 0.0 and solved_report is not None:
        peak = float(getattr(solved_report, "max_local_speed", 0.0) or 0.0)
        frac = float(getattr(solved_report, "superluminal_fraction", 0.0) or 0.0)
        peak_margin = max(
            0.0,
            min((peak - 1.0 - SOLVED_PEAK_FLOOR)
                / (SOLVED_PEAK_TARGET - SOLVED_PEAK_FLOOR), 1.0),
        )
        locality = max(
            0.0,
            min((SOLVED_LOCALITY_CEILING - frac) / SOLVED_LOCALITY_CEILING, 1.0),
        )
        solved_gate = peak_margin * locality
        components["operational_ftl_solved"] = solved_ftl_raw * solved_gate
        if solved_gate < 1.0 and solved_ftl_raw > 0.05:
            notes.append(
                "operational_ftl_solved down-gated as delocalized coordinate "
                f"offset (max c={peak:.3f}, superluminal_fraction={frac:.2f}, "
                f"gate={solved_gate:.2f})"
            )
    else:
        components["operational_ftl_solved"] = solved_ftl_raw
    if metrics.mechanism_descriptor is not None:
        components["mechanism_descriptor"] = float(
            min(max(metrics.mechanism_descriptor, 0.0), 1.0)
        )
    # A tiny coordinate-time win is useful as a shaping signal, but it should not
    # dominate the search.  The previous log scale saturated at F_op~0.01 and let
    # quasi-stationary "warp lens" artifacts score as solved FTL.  Reserve the
    # hard operational reward for a materially stronger evolved shortcut.
    if math.isfinite(f_op_ev) and f_op_ev > OP_FTL_FLOOR:
        components["operational_ftl"] = min(
            (f_op_ev - OP_FTL_FLOOR) / (OP_FTL_TARGET - OP_FTL_FLOOR),
            1.0,
        )
    else:
        components["operational_ftl"] = 0.0

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
    ctx.geo_trustworthy = geo_trustworthy
    ctx.f_geo = f_geo

    # Time-resolved headline.  The gauge-invariant shortcut peaks mid-run and
    # diffuses, so scoring only the final frame (``geo_report``) is half-blind:
    # it systematically under-credits a real-but-transient channel (f_geo>0 at
    # t=16 but ~0 by t=30) and would over-credit a late-blooming gauge artifact.
    # When the per-frame FTL stream is available we instead reward the MEAN over
    # the whole run of a per-frame gauge-invariant magnitude (each frame gated on
    # its own trustworthiness).  This is the avg/sum/divide design: a channel
    # that is FTL for most of the evolution scores high, a one-frame
    # Alcubierre-like transient averages down toward zero, and a diffused final
    # frame no longer zeroes a genuinely sustained shortcut.  The end-state
    # ``structural_persistence`` multiplier is kept so a fragmenting structure is
    # still discounted, consistent with the shaping rewards below.
    ts = metrics.ftl_timeseries
    geo_timeavg: float | None = None
    if ts is not None and ts.n_frames >= 2:
        per_frame = [
            _geo_magnitude(fg)
            if (trust and math.isfinite(fg) and fg > GEO_FTL_FLOOR)
            else 0.0
            for fg, trust in zip(ts.f_geo, ts.geo_trustworthy)
        ]
        geo_timeavg = sum(per_frame) / len(per_frame)
        components["ftl_geo_timeavg"] = geo_timeavg
        components["ftl_geo_peak"] = _geo_magnitude(ts.f_geo_peak) if ts.f_geo_peak > GEO_FTL_FLOOR else 0.0
        components["ftl_lifetime"] = ts.ftl_lifetime_fraction

    if geo_timeavg is not None:
        components["operational_ftl_geodesic"] = geo_timeavg * structural_persistence
        if geo_timeavg > 0.0:
            peak_at = (
                f" at t={ts.t_at_f_geo_peak:.2f}"
                if ts.t_at_f_geo_peak is not None
                else ""
            )
            notes.append(
                f"time-averaged gauge-invariant FTL over {ts.n_frames} frames: "
                f"mean_magnitude={geo_timeavg:.3e}, peak f_geo={ts.f_geo_peak:.3e}{peak_at}, "
                f"FTL-lifetime={ts.ftl_lifetime_fraction:.0%}"
            )
    elif geo_trustworthy and math.isfinite(f_geo) and f_geo > GEO_FTL_FLOOR:
        # Fallback: no per-frame stream -> single final-frame magnitude, scaled
        # by structural_persistence (defaults to 1.0 when matter slice missing).
        components["operational_ftl_geodesic"] = _geo_magnitude(f_geo) * structural_persistence
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
    # An evolved coordinate light-speed shortcut that the *trustworthy*
    # gauge-invariant geodesic probe contradicts (reliable ray bundle, yet
    # f_geo <= floor) is a coordinate/gauge artifact, not operational FTL.  The
    # gauge-invariant probe is the arbiter, so the coordinate signal must not
    # claim the primary operational-FTL reward -- otherwise a cone-tilt artifact
    # banks ~400 points and ranks beside a genuine gauge-invariant shortcut.
    if (
        geo_trustworthy
        and (not math.isfinite(f_geo) or f_geo <= GEO_FTL_FLOOR)
        and components["operational_ftl"] > 0.0
    ):
        notes.append(
            "operational_ftl zeroed: trustworthy geodesic probe found no "
            f"gauge-invariant shortcut (f_geo={f_geo:.3e}); coordinate channel "
            "is a gauge artifact"
        )
        components["operational_ftl"] = 0.0
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
        curvature_activity_scale = 5.0
        components["curvature_activity"] = curvature / (curvature + curvature_activity_scale)
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
    # Saturation scales.  These were previously so tight (0.05) that the
    # precursor became a near-binary 1.0 flag: any ``max_local_speed`` past
    # ~1.05 maxed the speed term and any superluminal area past ~5% maxed the
    # fraction term, so genuinely different geometries (c = 1.05 vs 1.47, 5% vs
    # 30% area) all scored ~1.0 and the reward carried no gradient.  Widening
    # them to span the observed evolved range (max_local_speed ~1.0-1.47,
    # superluminal area now noise-free after the margin fix) restores a smooth
    # slope the optimizer can climb instead of a cliff it instantly tops out on.
    structure = components["curvature_activity"]

    precursor_ev = _precursor(metrics.general_ftl_evolved, structure=structure)
    precursor_t0 = _precursor(metrics.general_ftl, structure=structure)
    precursor_solved = _precursor(metrics.general_ftl_solved, structure=structure)
    components["ftl_precursor"] = max(
        precursor_ev, 0.5 * precursor_t0, 0.5 * precursor_solved
    )

    components["shift_drive"] = max(
        _shift_drive(metrics.general_ftl_evolved),
        0.5 * _shift_drive(metrics.general_ftl),
        0.5 * _shift_drive(metrics.general_ftl_solved),
    )

    # CHANNEL PROGRESS -- coupled shaping gradient that pushes the search away
    # from isolated precursor-only (eval 128) or shift-only (eval 57) basins.
    # Path closeness uses relative excess (t_min - t_flat) / t_flat so flat
    # space (t_min == t_flat) earns nothing unless precursor/shift are nonzero.
    channel_ev = _channel_progress(metrics.general_ftl_evolved, structure=structure)
    channel_t0 = _channel_progress(metrics.general_ftl, structure=structure)
    channel_solved = _channel_progress(metrics.general_ftl_solved, structure=structure)
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
    ctx.stationary_artifact = stationary_artifact
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

    # Persistence gate on the FTL-shaping rewards.  ftl_precursor /
    # channel_progress / shift_drive credit a geometry for tilted light cones and
    # frame-drag, but a configuration that dissipates or fragments by the stop
    # time has not earned that "promising precursor" credit -- its cones tilt in
    # a structure that no longer exists (e.g. a coherent bubble that shatters
    # into turbulent lobes still shows local cone-tilt yet is not a propagating
    # channel).  Scaling these shaping rewards by the fraction of peak matter
    # energy density retained (structural_persistence) means only a structure
    # that actually holds together banks them, so a broken end-state can no
    # longer out-rank a coherent survivor.  Persistence defaults to 1.0 when the
    # matter-density series is unavailable, leaving the rewards untouched.
    persistence_gate = components.get("structural_persistence", 1.0)
    for key in ("ftl_precursor", "channel_progress", "shift_drive"):
        components[key] *= persistence_gate

    # Geodesic contradiction gate.  When a trustworthy null-ray probe ran and
    # found no gauge-invariant shortcut, coordinate shaping rewards are gauge
    # artifacts and must not inflate the archive (v12 eval 197 scored 130 from
    # shaping alone with f_geo=0).
    if (
        geo_trustworthy
        and (not math.isfinite(f_geo) or f_geo <= GEO_FTL_FLOOR)
    ):
        shaped_zeroed = False
        for key in (
            "operational_ftl_solved",
            "ftl_precursor",
            "channel_progress",
            "shift_drive",
        ):
            if components.get(key, 0.0) > 0.0:
                shaped_zeroed = True
            components[key] = 0.0
        if shaped_zeroed:
            notes.append(
                "FTL shaping zeroed: trustworthy geodesic probe found no "
                f"gauge-invariant shortcut (f_geo={f_geo:.3e}); coordinate "
                "precursor/channel rewards are artifacts"
            )
