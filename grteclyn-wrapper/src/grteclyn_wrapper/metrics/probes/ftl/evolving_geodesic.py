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
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from .geodesic import (
    H_REL_TOL,
    GeodesicFtlReport,
    NullRayResult,
    _geodesic_report_at_resolution,
    _hamiltonian_rhs,
    _null_relative_drift,
    compute_geodesic_ftl_from_plotfile,
    future_null_cov,
    null_hamiltonian,
    project_null,
)
from .metric_field import (
    EvolvingMetricField,
    MetricField,
    StaticMetricField,
    evolving_field_from_analytic_stack,
    evolving_field_from_plotfiles,
)


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


def integrate_null_ray_on_field(
    field: MetricField,
    *,
    x_start: float,
    x_end: float,
    y0: float,
    z0: float,
    t0: float = 0.0,
    max_steps: int = 50_000,
    ds_init: float = 0.05,
    h_tol: float = 1.0e-6,
) -> NullRayResult:
    """Trace one null ray through a possibly time-dependent metric field."""
    x = np.array([t0, x_start, y0, z0], dtype=float)
    t_flat = abs(x_end - x_start)
    n_hat = np.array([1.0, 0.0, 0.0], dtype=float)

    g_pt, ginv_pt, _ = field.sample(x)
    k = future_null_cov(g_pt, n_hat)

    max_h = 0.0
    max_h_rel = 0.0
    ds = ds_init

    for _ in range(max_steps):
        if x[1] >= x_end:
            t_coord = abs(float(x[0] - t0))
            return NullRayResult(
                reached=True,
                t_coord=t_coord,
                t_flat=t_flat,
                max_h_drift=max_h,
                max_h_rel=max_h_rel,
            )

        g_pt, ginv_pt, dg_pt = field.sample(x)
        h = abs(null_hamiltonian(ginv_pt, k))
        max_h = max(max_h, h)
        max_h_rel = max(max_h_rel, _null_relative_drift(ginv_pt, k))
        if h > h_tol:
            k = project_null(g_pt, ginv_pt, k, dx_ref=(ginv_pt @ k)[1:])

        def rhs(xp: NDArray[np.float64], kp: NDArray[np.float64]) -> tuple[NDArray, NDArray]:
            _gp, ginvp, dgp = field.sample(xp)
            return _hamiltonian_rhs(ginvp, dgp, kp)

        k1x, k1k = rhs(x, k)
        k2x, k2k = rhs(x + 0.5 * ds * k1x, k + 0.5 * ds * k1k)
        k3x, k3k = rhs(x + 0.5 * ds * k2x, k + 0.5 * ds * k2k)
        k4x, k4k = rhs(x + ds * k3x, k + ds * k3k)
        x = x + (ds / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        k = k + (ds / 6.0) * (k1k + 2 * k2k + 2 * k3k + k4k)

        if x[1] < field.x_spatial_min() or x[1] > field.x_spatial_max():
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


def _ray_fan_geometry(field: MetricField, *, n_rays: int) -> tuple[float, float, float, float, float]:
    shape = field.spatial_shape
    cy = field.origin[1] + 0.5 * (shape[1] - 1) * field.spatial_spacing[1]
    cz = field.origin[2] + 0.5 * (shape[2] - 1) * field.spatial_spacing[2]
    x_start = field.origin[0] + 0.05 * (shape[0] - 1) * field.spatial_spacing[0]
    x_end = field.origin[0] + 0.95 * (shape[0] - 1) * field.spatial_spacing[0]
    t_flat = x_end - x_start
    return x_start, x_end, cy, cz, t_flat


def compute_evolving_geodesic_ftl(
    field: MetricField,
    *,
    t_emit: float | None = None,
    n_rays: int = 5,
    h_tol: float = 1.0e-6,
    frozen_peak: float | None = None,
) -> EvolvingGeodesicFtlReport:
    """Run a fan of evolving null rays and return end-to-end ``f_geo``."""
    if isinstance(field, EvolvingMetricField):
        t_emit_val = float(field.times[0]) if t_emit is None else float(t_emit)
    else:
        t_emit_val = 0.0 if t_emit is None else float(t_emit)

    x_start, x_end, cy, cz, t_flat = _ray_fan_geometry(field, n_rays=n_rays)
    shape = field.spatial_shape
    offsets = np.linspace(-0.1, 0.1, max(1, n_rays)) * (shape[2] - 1) * field.spatial_spacing[2]

    results: list[NullRayResult] = []
    for dz in offsets:
        results.append(
            integrate_null_ray_on_field(
                field,
                x_start=x_start,
                x_end=x_end,
                y0=cy,
                z0=cz + float(dz),
                t0=t_emit_val,
                h_tol=h_tol,
            )
        )

    reached = [r for r in results if r.reached and r.t_coord is not None]
    if not reached:
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
            notes=("no rays reached detector",),
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
    notes.append(f"evolving end-to-end trace t_emit={t_emit_val:.3f}")

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
    )


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
    return compute_evolving_geodesic_ftl(
        field, n_rays=n_rays, h_tol=h_tol, frozen_peak=frozen_peak
    )


def compute_evolving_geodesic_ftl_from_analytic(
    g: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    n_rays: int = 5,
    h_tol: float = 1.0e-6,
) -> EvolvingGeodesicFtlReport:
    """Analytic-stack entry point (Alcubierre validation)."""
    field = evolving_field_from_analytic_stack(g, spacing)
    return compute_evolving_geodesic_ftl(field, n_rays=n_rays, h_tol=h_tol)


def write_evolving_geodesic_json(path: Path, report: EvolvingGeodesicFtlReport) -> None:
    """Persist evolving geodesic report for scoring / inspection."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["notes"] = list(report.notes)
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
