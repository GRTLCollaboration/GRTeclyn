"""4D evolving null-geodesic FTL probe.

Traces null rays through a time-interpolated metric stack so the shortcut
measures end-to-end pulse transit through the *evolving* geometry, not a frozen
Cauchy slice.

Emission protocol
-----------------
* ``t_emit = times[0]`` (first available plotfile time).
* Rays launch at ``x^μ = (t_emit, x_start, y_mid, z_i)`` toward ``+x``.
* The metric is sampled at the ray's coordinate time ``x^0`` (linearly
  interpolated between plotfile slices).
* ``f_geo = max(0, (t_flat - t_travel) / t_flat)`` with
  ``t_travel = |x^0_arrival - t_emit|``.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import NDArray

from .geodesic import (
    H_REL_TOL,
    GeodesicFtlReport,
    NullRayResult,
    _geodesic_report_at_resolution,
    geodesic_report_from_metric_g,
    _hamiltonian_rhs,
    _null_relative_drift,
    compute_geodesic_ftl_from_plotfile,
    future_null_cov,
    null_hamiltonian,
    project_null,
    ray_is_captured,
    rays_complete,
)
from .metric_field import (
    EvolvingMetricField,
    MetricField,
    StaticMetricField,
    evolving_field_from_analytic_stack,
    evolving_field_from_plotfiles,
)
from .metric_stack_cache import (
    evolving_field_from_metric_stack_cache,
    list_slice_files,
    slice_time,
    subsample_slice_files,
)
from .evolving_geodesic_options import (
    EvolvingGeodesicOptions,
    HQ_OPTIONS,
    SEARCH_OPTIONS,
    evolving_geodesic_options_from_env,
    geo_directions_from_env,
)


_AXIS_LABELS = ("x", "y", "z")


@dataclass(frozen=True)
class EvolvingGeodesicFtlReport:
    """End-to-end evolving null-geodesic shortcut."""

    f_geo: float
    f_geo_frozen_peak: float | None
    t_emit: float
    t_arrival: float | None
    t_flat: float
    n_rays: int
    n_reached: int
    max_h_drift: float
    h_quality_ok: bool
    max_h_rel_drift: float = 0.0
    notes: tuple[str, ...] = ()
    # Continuous-emission sweep: (t_emit, f_geo, n_reached) per launch time.
    # Empty for a single-launch trace.  ``f_geo``/``t_emit`` above report the peak.
    emit_sweep: tuple[tuple[float, float, int], ...] = ()
    # Rays that fell into a puncture throat / horizon (excluded from f_geo).
    n_captured: int = 0


def _spatial_extent(
    field: MetricField, axis: int
) -> tuple[float, float]:
    """Return (min, max) coordinate along spatial axis ``axis`` (0=x, 1=y, 2=z)."""
    origin = field.origin[:3]
    spacing = field.spatial_spacing
    shape = field.spatial_shape
    lo = float(origin[axis])
    hi = float(origin[axis] + (shape[axis] - 1) * spacing[axis])
    return lo, hi


def integrate_null_ray_on_field(
    field: MetricField,
    *,
    x_start: float,
    x_end: float,
    y0: float,
    z0: float,
    t0: float = 0.0,
    axis: int = 0,
    max_steps: int = 50_000,
    ds_init: float = 0.05,
    h_tol: float = 1.0e-6,
    h_rel_abort: float | None = None,
    detect_capture: bool = True,
    allow_frozen_tail: bool = False,
) -> NullRayResult:
    """Trace one null ray through a possibly time-dependent metric field.

    ``allow_frozen_tail=True`` restores the pre-2026-07-28 behaviour of letting
    a ray complete through the (time-clamped) final slice after the stack ends.
    Only legitimate for genuinely stationary stacks (tests, analytic controls).
    """
    prop_idx = axis + 1
    pos = [0.0, 0.0, 0.0]
    pos[axis] = x_start
    transverse = [i for i in range(3) if i != axis]
    pos[transverse[0]] = y0
    pos[transverse[1]] = z0
    x = np.array([t0, pos[0], pos[1], pos[2]], dtype=float)
    t_flat = abs(x_end - x_start)
    n_hat = np.zeros(3, dtype=float)
    n_hat[axis] = 1.0

    g_pt, ginv_pt, _ = field.sample(x)
    k = future_null_cov(g_pt, n_hat)

    max_h = 0.0
    max_h_rel = 0.0
    ds = ds_init
    s_min, s_max = _spatial_extent(field, axis)

    # A ray still in flight past the last stored slice would silently complete
    # through a FROZEN final metric (EvolvingMetricField clamps its time
    # bracket).  On the candidate-146 RM stack a t_emit=12 launch "arrived" at
    # t=32.75 against a stack ending at t=30 -- 2.75 units of pretend geometry.
    # Such rays are integration failures, not measurements.
    t_stack_end = (
        float(field.times[-1])
        if isinstance(field, EvolvingMetricField) and not allow_frozen_tail
        else None
    )

    x_prev = x.copy()
    for _ in range(max_steps):
        if x[prop_idx] >= x_end:
            # Interpolate the crossing back onto the detector plane: the RK
            # step overshoots by up to one step, and reading the clock at the
            # overshot position biased t_coord late (f_geo low) by up to
            # ~ds ~ 0.05.
            frac = 1.0
            denom = float(x[prop_idx] - x_prev[prop_idx])
            if denom > 0.0 and x_prev[prop_idx] < x_end:
                frac = (x_end - float(x_prev[prop_idx])) / denom
            t_cross = float(x_prev[0]) + frac * float(x[0] - x_prev[0])
            if t_stack_end is not None and t_cross > t_stack_end + 1.0e-9:
                return NullRayResult(
                    reached=False,
                    t_coord=None,
                    t_flat=t_flat,
                    max_h_drift=max_h,
                    max_h_rel=max_h_rel,
                    notes=(
                        f"ray outlived metric stack (arrival t={t_cross:.2f} > "
                        f"last slice t={t_stack_end:.2f}; frozen-tail geometry)",
                    ),
                )
            t_coord = abs(t_cross - t0)
            return NullRayResult(
                reached=True,
                t_coord=t_coord,
                t_flat=t_flat,
                max_h_drift=max_h,
                max_h_rel=max_h_rel,
            )
        if t_stack_end is not None and float(x[0]) > t_stack_end + 1.0e-9:
            return NullRayResult(
                reached=False,
                t_coord=None,
                t_flat=t_flat,
                max_h_drift=max_h,
                max_h_rel=max_h_rel,
                notes=(
                    f"ray outlived metric stack (t={float(x[0]):.2f} > last "
                    f"slice t={t_stack_end:.2f}; frozen-tail geometry)",
                ),
            )

        g_pt, ginv_pt, dg_pt = field.sample(x)
        if detect_capture and ray_is_captured(g_pt, ginv_pt):
            return NullRayResult(
                reached=False,
                t_coord=None,
                t_flat=t_flat,
                max_h_drift=max_h,
                max_h_rel=max_h_rel,
                notes=("ray captured by puncture/horizon",),
                captured=True,
            )
        h = abs(null_hamiltonian(ginv_pt, k))
        max_h = max(max_h, h)
        max_h_rel = max(max_h_rel, _null_relative_drift(ginv_pt, k))
        if h_rel_abort is not None and max_h_rel > h_rel_abort:
            return NullRayResult(
                reached=False,
                t_coord=None,
                t_flat=t_flat,
                max_h_drift=max_h,
                max_h_rel=max_h_rel,
                notes=(f"h_rel abort ({max_h_rel:.2e}>{h_rel_abort:.2e})",),
            )
        if h > h_tol:
            k = project_null(g_pt, ginv_pt, k, dx_ref=(ginv_pt @ k)[1:])

        def rhs(xp: NDArray[np.float64], kp: NDArray[np.float64]) -> tuple[NDArray, NDArray]:
            _gp, ginvp, dgp = field.sample(xp)
            return _hamiltonian_rhs(ginvp, dgp, kp)

        k1x, k1k = rhs(x, k)
        k2x, k2k = rhs(x + 0.5 * ds * k1x, k + 0.5 * ds * k1k)
        k3x, k3k = rhs(x + 0.5 * ds * k2x, k + 0.5 * ds * k2k)
        k4x, k4k = rhs(x + ds * k3x, k + ds * k3k)
        x_prev = x
        x = x + (ds / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        k = k + (ds / 6.0) * (k1k + 2 * k2k + 2 * k3k + k4k)

        if x[prop_idx] < s_min or x[prop_idx] > s_max:
            return NullRayResult(
                reached=False,
                t_coord=None,
                t_flat=t_flat,
                max_h_drift=max_h,
                max_h_rel=max_h_rel,
                notes=("ray left grid",),
            )

    return NullRayResult(
        reached=False,
        t_coord=None,
        t_flat=t_flat,
        max_h_drift=max_h,
        max_h_rel=max_h_rel,
        notes=("max_steps exceeded",),
    )


def _ray_fan_geometry(
    field: MetricField, *, axis: int, n_rays: int
) -> tuple[float, float, float, float, float, NDArray[np.float64]]:
    """Return propagation endpoints, transverse centers, and fan offsets.

    ``integrate_null_ray_on_field`` maps the caller's ``y0`` onto the lower
    transverse axis (``transverse[0]``) and ``z0`` onto the higher one
    (``transverse[1]``).  ``transverse`` is built in ascending order, so the
    fan always sits on ``transverse[1]`` -- i.e. the offset is always applied
    to ``z0``.
    """
    shape = field.spatial_shape
    spacing = field.spatial_spacing
    origin = field.origin[:3]
    transverse = [i for i in range(3) if i != axis]
    fix_spatial, fan_spatial = transverse  # ascending: fan is the higher axis (z0 slot)

    prop_start = origin[axis] + 0.05 * (shape[axis] - 1) * spacing[axis]
    prop_end = origin[axis] + 0.95 * (shape[axis] - 1) * spacing[axis]
    t_flat = prop_end - prop_start

    y0_center = origin[fix_spatial] + 0.5 * (shape[fix_spatial] - 1) * spacing[fix_spatial]
    z0_center = origin[fan_spatial] + 0.5 * (shape[fan_spatial] - 1) * spacing[fan_spatial]
    fan_offsets = (
        np.linspace(-0.1, 0.1, max(1, n_rays)) * (shape[fan_spatial] - 1) * spacing[fan_spatial]
    )

    return prop_start, prop_end, y0_center, z0_center, t_flat, fan_offsets


def compute_evolving_geodesic_ftl(
    field: MetricField,
    *,
    t_emit: float | None = None,
    n_rays: int = 5,
    max_steps: int = 50_000,
    ds_init: float = 0.05,
    h_tol: float = 1.0e-6,
    h_rel_abort: float | None = None,
    frozen_peak: float | None = None,
    axis: int = 0,
    allow_frozen_tail: bool = False,
) -> EvolvingGeodesicFtlReport:
    """Run a fan of evolving null rays and return end-to-end ``f_geo``."""
    if isinstance(field, EvolvingMetricField):
        t_emit_val = float(field.times[0]) if t_emit is None else float(t_emit)
    else:
        t_emit_val = 0.0 if t_emit is None else float(t_emit)

    prop_start, prop_end, y0_center, z0_center, t_flat, fan_offsets = _ray_fan_geometry(
        field, axis=axis, n_rays=n_rays
    )

    results: list[NullRayResult] = []
    for offset in fan_offsets:
        y0, z0 = y0_center, z0_center + float(offset)
        results.append(
            integrate_null_ray_on_field(
                field,
                x_start=prop_start,
                x_end=prop_end,
                y0=y0,
                z0=z0,
                t0=t_emit_val,
                axis=axis,
                max_steps=max_steps,
                ds_init=ds_init,
                h_tol=h_tol,
                h_rel_abort=h_rel_abort,
                allow_frozen_tail=allow_frozen_tail,
            )
        )

    reached = [r for r in results if r.reached and r.t_coord is not None]
    n_captured = sum(1 for r in results if r.captured)
    if not reached:
        note = f"no rays reached detector (axis={_AXIS_LABELS[axis]})"
        if n_captured:
            note += f"; {n_captured}/{len(results)} captured by puncture/horizon"
        return EvolvingGeodesicFtlReport(
            f_geo=0.0,
            f_geo_frozen_peak=frozen_peak,
            t_emit=t_emit_val,
            t_arrival=None,
            t_flat=t_flat,
            n_rays=len(results),
            n_reached=0,
            max_h_drift=max((r.max_h_drift for r in results), default=0.0),
            h_quality_ok=False,
            max_h_rel_drift=max((r.max_h_rel for r in results), default=0.0),
            notes=(note,),
            n_captured=n_captured,
        )

    t_min = min(r.t_coord for r in reached if r.t_coord is not None)
    max_h = max(r.max_h_drift for r in results)
    max_h_rel = max(r.max_h_rel for r in results)
    h_ok = max_h_rel <= H_REL_TOL
    f_geo = max(0.0, (t_flat - t_min) / t_flat) if t_flat > 0 else 0.0
    t_arrival = t_emit_val + t_min

    notes: list[str] = []
    if not h_ok:
        notes.append(
            f"null constraint drift high (rel H={max_h_rel:.2e}, abs H={max_h:.2e})"
        )
    if n_captured:
        notes.append(
            f"{n_captured}/{len(results)} rays captured by puncture/horizon "
            "(excluded from f_geo)"
        )
    notes.append(f"evolving end-to-end trace t_emit={t_emit_val:.3f}")
    notes.append(f"probe_axis={_AXIS_LABELS[axis]}")

    return EvolvingGeodesicFtlReport(
        f_geo=f_geo,
        f_geo_frozen_peak=frozen_peak,
        t_emit=t_emit_val,
        t_arrival=t_arrival,
        t_flat=t_flat,
        n_rays=len(results),
        n_reached=len(reached),
        max_h_drift=max_h,
        h_quality_ok=h_ok,
        max_h_rel_drift=max_h_rel,
        notes=tuple(notes),
        n_captured=n_captured,
    )


def evolving_report_trustworthy(report: EvolvingGeodesicFtlReport) -> bool:
    """True when the 4D trace is certified: h-quality held and the bundle closed.

    The single trust bar for an evolving report.  Consumers that only have the
    scalar fields (``EvolvingGeodesicMetrics``, a JSON payload) should call
    :func:`rays_complete` with the same three counts rather than re-deriving.
    """
    return report.h_quality_ok and rays_complete(
        report.n_rays, report.n_reached, report.n_captured
    )


def _report_probe_score(report: EvolvingGeodesicFtlReport) -> float:
    if evolving_report_trustworthy(report):
        return report.f_geo
    return -1.0


def compute_evolving_geodesic_ftl_best_direction(
    field: MetricField,
    *,
    directions: Sequence[str],
    t_emit: float | None = None,
    n_rays: int = 5,
    max_steps: int = 50_000,
    ds_init: float = 0.05,
    h_tol: float = 1.0e-6,
    h_rel_abort: float | None = None,
    frozen_peak: float | None = None,
    allow_frozen_tail: bool = False,
) -> EvolvingGeodesicFtlReport:
    """Scan principal axes and return the report with the largest trustworthy ``f_geo``."""
    axis_map = {label: idx for idx, label in enumerate(_AXIS_LABELS)}
    axes = [axis_map[d] for d in directions if d in axis_map]
    if not axes:
        axes = [0]

    reports = [
        compute_evolving_geodesic_ftl(
            field,
            t_emit=t_emit,
            n_rays=n_rays,
            max_steps=max_steps,
            ds_init=ds_init,
            h_tol=h_tol,
            h_rel_abort=h_rel_abort,
            frozen_peak=frozen_peak if len(axes) == 1 else None,
            axis=axis,
            allow_frozen_tail=allow_frozen_tail,
        )
        for axis in axes
    ]
    best_axis_label = _AXIS_LABELS[axes[0]]
    best = reports[0]
    best_score = _report_probe_score(best)
    for axis, rep in zip(axes, reports):
        score = _report_probe_score(rep)
        if score > best_score:
            best_score = score
            best = rep
            best_axis_label = _AXIS_LABELS[axis]
    extra = tuple(n for n in best.notes if not n.startswith("best_direction="))
    return EvolvingGeodesicFtlReport(
        f_geo=best.f_geo,
        f_geo_frozen_peak=frozen_peak if frozen_peak is not None else best.f_geo_frozen_peak,
        t_emit=best.t_emit,
        t_arrival=best.t_arrival,
        t_flat=best.t_flat,
        n_rays=best.n_rays,
        n_reached=best.n_reached,
        max_h_drift=best.max_h_drift,
        h_quality_ok=best.h_quality_ok,
        max_h_rel_drift=best.max_h_rel_drift,
        notes=extra + (f"best_direction={best_axis_label}",),
        n_captured=best.n_captured,
    )


def _emission_times(
    field: MetricField, *, interval: float, max_emissions: int
) -> list[float] | None:
    """Launch times for a continuous-emission sweep, or ``None`` if disabled.

    Fires at ``times[0], +interval, +2*interval, ...`` up to ``max_emissions``
    launches, capped at the last available slice time (a ray launched later than
    that samples only the frozen final slice).
    """
    if interval <= 0.0 or max_emissions <= 1:
        return None
    if not isinstance(field, EvolvingMetricField):
        return None
    t0 = float(field.times[0])
    t_last = float(field.times[-1])
    times = [t0 + i * interval for i in range(max_emissions)]
    times = [t for t in times if t <= t_last + 1.0e-9]
    return times if len(times) > 1 else None


def _geodesic_emit_min_time() -> float | None:
    """Minimum eligible t_emit for the peak-f_geo selection.

    Preference order:
      1. ``GEODESIC_EMIT_MIN_TIME`` (explicit; use with an always-on pump)
      2. ``RL_PUMP_STOP_TIME`` (legacy igniter filter; preserves published numbers)

    Returns None when neither is set, so every launch remains eligible.
    """
    for key in ("GEODESIC_EMIT_MIN_TIME", "RL_PUMP_STOP_TIME"):
        raw = os.environ.get(key, "").strip()
        if not raw:
            continue
        try:
            stop = float(raw)
        except (TypeError, ValueError):
            continue
        if stop >= 0.0:
            return stop
    return None


def _pump_stop_time_for_geo() -> float | None:
    """Backward-compatible alias for :func:`_geodesic_emit_min_time`."""
    return _geodesic_emit_min_time()


def _eligible_emit_reports(
    reports: list[tuple[float, EvolvingGeodesicFtlReport]],
) -> list[tuple[float, EvolvingGeodesicFtlReport]]:
    """Keep launches with t_emit >= emit-min when a floor is configured.

    If every launch is below the floor (too short a run), fall back to all
    reports so the probe still returns a number; the peak_note flags the filter.
    """
    stop = _geodesic_emit_min_time()
    if stop is None:
        return reports
    kept = [(te, rep) for te, rep in reports if te + 1.0e-12 >= stop]
    return kept if kept else reports


def compute_evolving_geodesic_ftl_emission_sweep(
    field: MetricField,
    *,
    directions: Sequence[str],
    emit_times: Sequence[float],
    n_rays: int = 5,
    max_steps: int = 50_000,
    ds_init: float = 0.05,
    h_tol: float = 1.0e-6,
    h_rel_abort: float | None = None,
    frozen_peak: float | None = None,
) -> EvolvingGeodesicFtlReport:
    """Fire ray fans at successive launch times and report the peak-over-time.

    Implements the "surf the transient wave" probe: instead of a single launch at
    ``times[0]`` (before the warp has formed), it maps ``f_geo`` vs launch time so
    the strongest surfable window is captured.  The headline ``f_geo``/``t_emit``
    are the best trustworthy launch; the full map is in ``emit_sweep`` and notes.
    """
    reports: list[tuple[float, EvolvingGeodesicFtlReport]] = []
    for i, te in enumerate(emit_times):
        rep = compute_evolving_geodesic_ftl_best_direction(
            field,
            directions=directions,
            t_emit=float(te),
            n_rays=n_rays,
            max_steps=max_steps,
            ds_init=ds_init,
            h_tol=h_tol,
            h_rel_abort=h_rel_abort,
            # frozen_peak is launch-time independent; attach once.
            frozen_peak=frozen_peak if i == 0 else None,
        )
        reports.append((float(te), rep))

    sweep = tuple((te, float(rep.f_geo), int(rep.n_reached)) for te, rep in reports)
    # Optional emit floor (GEODESIC_EMIT_MIN_TIME, else RL_PUMP_STOP_TIME).
    eligible = _eligible_emit_reports(reports)
    best_te, best = max(eligible, key=lambda tr: _report_probe_score(tr[1]))
    sweep_note = "emit_sweep: " + ", ".join(
        f"t={te:.2f}->f={rep.f_geo:.3f}(n{rep.n_reached})" for te, rep in reports
    )
    peak_note = f"peak f_geo={best.f_geo:.3f} at t_emit={best_te:.2f} over {len(reports)} launches"
    emit_floor = _geodesic_emit_min_time()
    if emit_floor is not None and len(eligible) < len(reports):
        peak_note += (
            f" (emit floor: kept {len(eligible)}/{len(reports)} launches "
            f"with t_emit>={emit_floor:g})"
        )

    return EvolvingGeodesicFtlReport(
        f_geo=best.f_geo,
        f_geo_frozen_peak=frozen_peak if frozen_peak is not None else best.f_geo_frozen_peak,
        t_emit=best_te,
        t_arrival=best.t_arrival,
        t_flat=best.t_flat,
        n_rays=best.n_rays,
        n_reached=best.n_reached,
        max_h_drift=best.max_h_drift,
        h_quality_ok=best.h_quality_ok,
        max_h_rel_drift=best.max_h_rel_drift,
        notes=best.notes + (sweep_note, peak_note),
        emit_sweep=sweep,
        n_captured=best.n_captured,
    )


def _frozen_peak_from_g_slices(
    slices: Iterable[NDArray[np.float64]],
    origin: NDArray[np.float64],
    spacing: tuple[float, float, float],
    *,
    n_rays: int,
    h_tol: float,
) -> float | None:
    peak = 0.0
    for g in slices:
        rep = geodesic_report_from_metric_g(
            g, origin, spacing, n_rays=n_rays, h_tol=h_tol
        )
        if rep.h_quality_ok and rep.f_geo > peak:
            peak = rep.f_geo
    return peak if peak > 0.0 else None


def _frozen_peak_from_plotfiles(
    plotfiles: Sequence[str | Path],
    *,
    n: int,
    half_width: float | None,
    n_rays: int,
    h_tol: float,
) -> float | None:
    peak = 0.0
    for path in plotfiles:
        rep = _geodesic_report_at_resolution(
            path, n=n, half_width=half_width, n_rays=n_rays, h_tol=h_tol
        )
        if rep is not None and rep.h_quality_ok and rep.f_geo > peak:
            peak = rep.f_geo
    return peak if peak > 0.0 else None


def compute_evolving_geodesic_ftl_from_metric_stack_cache(
    cache_dir: Path,
    *,
    options: EvolvingGeodesicOptions | None = None,
    n_rays: int | None = None,
    h_tol: float = 1.0e-6,
    max_time: float | None = None,
    frozen_peak_override: float | None = None,
) -> EvolvingGeodesicFtlReport | None:
    """End-to-end evolving trace from cached per-plotfile metric slices.

    ``frozen_peak_override`` supplies a frozen-slice peak measured elsewhere
    (e.g. the consumer's per-plotfile ``f_geo`` column) when
    ``compute_frozen_peak`` is off; recomputing it from the cache costs two
    full-grid inversions plus gradients per slice, which at 257^3 dominates the
    entire pass while only reproducing a number of lower provenance.
    """
    opts = options or evolving_geodesic_options_from_env()
    ray_count = opts.n_rays if n_rays is None else n_rays
    field = evolving_field_from_metric_stack_cache(
        cache_dir,
        slice_stride=opts.slice_stride,
        max_slices=opts.max_slices,
        max_time=max_time,
    )
    if field is None:
        return None
    frozen_peak = frozen_peak_override
    if opts.compute_frozen_peak:
        all_files = list_slice_files(cache_dir)
        if max_time is not None:
            all_files = [p for p in all_files if slice_time(p) <= max_time + 1.0e-9]
        files = subsample_slice_files(
            all_files,
            stride=opts.slice_stride,
            max_slices=opts.max_slices,
        )
        spacing = field.spatial_spacing
        # Stream one slice at a time: preloading the whole list held a second
        # full-precision copy of the stack (~47 GB at 257^3) for the entire pass.
        slice_iter = (
            np.asarray(np.load(path)["g"], dtype=np.float64) for path in files
        )
        frozen_peak = _frozen_peak_from_g_slices(
            slice_iter,
            field.origin,
            spacing,
            n_rays=ray_count,
            h_tol=h_tol,
        )
    emit_times = _emission_times(
        field, interval=opts.emit_interval, max_emissions=opts.max_emissions
    )
    if emit_times is not None:
        return compute_evolving_geodesic_ftl_emission_sweep(
            field,
            directions=opts.directions,
            emit_times=emit_times,
            n_rays=ray_count,
            max_steps=opts.max_steps,
            ds_init=opts.ds_init,
            h_tol=h_tol,
            h_rel_abort=opts.h_rel_abort,
            frozen_peak=frozen_peak,
        )
    return compute_evolving_geodesic_ftl_best_direction(
        field,
        directions=opts.directions,
        n_rays=ray_count,
        max_steps=opts.max_steps,
        ds_init=opts.ds_init,
        h_tol=h_tol,
        h_rel_abort=opts.h_rel_abort,
        frozen_peak=frozen_peak,
    )


def compute_evolving_geodesic_ftl_from_plotfiles(
    plotfiles: Sequence[str | Path],
    *,
    n_space: int = 65,
    half_width: float | None = None,
    n_rays: int = 5,
    h_tol: float = 1.0e-6,
) -> EvolvingGeodesicFtlReport | None:
    """Episode entry point: evolving trace from ordered plotfiles."""
    paths = [str(p) for p in plotfiles]
    if len(paths) < 3:
        return None
    try:
        field = evolving_field_from_plotfiles(
            paths, n_space=n_space, half_width=half_width
        )
    except Exception:
        return None
    frozen_peak = _frozen_peak_from_plotfiles(
        paths, n=n_space, half_width=half_width, n_rays=n_rays, h_tol=h_tol
    )
    opts = evolving_geodesic_options_from_env()
    emit_times = _emission_times(
        field, interval=opts.emit_interval, max_emissions=opts.max_emissions
    )
    if emit_times is not None:
        return compute_evolving_geodesic_ftl_emission_sweep(
            field,
            directions=opts.directions,
            emit_times=emit_times,
            n_rays=n_rays,
            max_steps=opts.max_steps,
            ds_init=opts.ds_init,
            h_tol=h_tol,
            h_rel_abort=opts.h_rel_abort,
            frozen_peak=frozen_peak,
        )
    return compute_evolving_geodesic_ftl_best_direction(
        field,
        directions=opts.directions,
        n_rays=n_rays,
        h_tol=h_tol,
        frozen_peak=frozen_peak,
    )


def compute_evolving_geodesic_ftl_from_analytic(
    g: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    n_rays: int = 5,
    h_tol: float = 1.0e-6,
    directions: Sequence[str] | None = None,
) -> EvolvingGeodesicFtlReport:
    """Analytic-stack entry point (Alcubierre validation).

    allow_frozen_tail: analytic metrics are defined for all t; the stored stack
    is a sampling convenience, so completing a flight through the clamped final
    slice is acceptable here (and only here -- production cache paths stay
    strict).
    """
    field = evolving_field_from_analytic_stack(g, spacing)
    dirs = tuple(directions) if directions is not None else geo_directions_from_env()
    if len(dirs) == 1:
        axis = _AXIS_LABELS.index(dirs[0]) if dirs[0] in _AXIS_LABELS else 0
        return compute_evolving_geodesic_ftl(
            field, n_rays=n_rays, h_tol=h_tol, axis=axis, allow_frozen_tail=True
        )
    return compute_evolving_geodesic_ftl_best_direction(
        field, directions=dirs, n_rays=n_rays, h_tol=h_tol, allow_frozen_tail=True
    )


def _json_safe(value: object) -> object:
    """Coerce numpy scalars so ``json.dumps`` accepts probe reports."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def write_evolving_geodesic_json(path: Path, report: EvolvingGeodesicFtlReport) -> None:
    """Persist evolving geodesic report for scoring / inspection."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _json_safe(asdict(report))
    payload["notes"] = list(report.notes)
    payload["emit_sweep"] = [list(row) for row in report.emit_sweep]
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def read_evolving_geodesic_json(path: Path) -> EvolvingGeodesicFtlReport | None:
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return EvolvingGeodesicFtlReport(
            f_geo=float(data["f_geo"]),
            f_geo_frozen_peak=(
                float(data["f_geo_frozen_peak"])
                if data.get("f_geo_frozen_peak") is not None
                else None
            ),
            t_emit=float(data["t_emit"]),
            t_arrival=(
                float(data["t_arrival"]) if data.get("t_arrival") is not None else None
            ),
            t_flat=float(data["t_flat"]),
            n_rays=int(data["n_rays"]),
            n_reached=int(data["n_reached"]),
            max_h_drift=float(data["max_h_drift"]),
            h_quality_ok=bool(data["h_quality_ok"]),
            max_h_rel_drift=float(data.get("max_h_rel_drift", 0.0)),
            notes=tuple(data.get("notes", ())),
            emit_sweep=tuple(
                (float(r[0]), float(r[1]), int(r[2]))
                for r in data.get("emit_sweep", ())
                if len(r) >= 3
            ),
            n_captured=int(data.get("n_captured", 0)),
        )
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def patch_ftl_timeseries_evolving(
    path: Path,
    *,
    f_geo_evol: float,
    f_geo_evol_ok: bool,
) -> None:
    """Set evolving columns on the last data row of ``ftl_timeseries.dat``."""
    if not path.is_file():
        return
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return
    header = lines[0]
    data_idx = None
    for i in range(len(lines) - 1, 0, -1):
        line = lines[i].strip()
        if line and not line.startswith("#"):
            data_idx = i
            break
    if data_idx is None:
        return
    parts = lines[data_idx].split()
    if len(parts) < 12:
        return
    if len(parts) >= 14:
        parts[12] = f"{f_geo_evol:.16e}"
        parts[13] = "1" if f_geo_evol_ok else "0"
    else:
        parts.extend([f"{f_geo_evol:.16e}", "1" if f_geo_evol_ok else "0"])
    lines[data_idx] = "  ".join(parts)
    if "f_geo_evol" not in header:
        lines[0] = header.rstrip() + "  f_geo_evol  f_geo_evol_ok"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
