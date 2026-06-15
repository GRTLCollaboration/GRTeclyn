"""Time-aware metric sampling for null geodesic integration.

``StaticMetricField`` wraps a single Cauchy slice (frozen snapshot).
``EvolvingMetricField`` linearly interpolates a stack of slices in simulation time.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, Sequence

import numpy as np
from numpy.typing import NDArray

from .geodesic_interp import (
    clamp_index,
    inverse_metric_4d,
    partial_inverse_metric,
    trilinear,
)


class MetricField(Protocol):
    """Sample ``g``, ``g^{-1}``, and ``∂_μ g^{ab}`` at a 4-position ``x^μ``."""

    origin: NDArray[np.float64]
    spatial_shape: tuple[int, int, int]
    spatial_spacing: tuple[float, float, float]

    def sample(
        self, x: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]: ...

    def x_spatial_min(self) -> float: ...

    def x_spatial_max(self) -> float: ...


@dataclass(frozen=True)
class StaticMetricField:
    """Single-slice metric ``g_{μν}(x,y,z)`` (time-independent)."""

    g: NDArray[np.float64]
    origin: NDArray[np.float64]
    spatial_spacing: tuple[float, float, float]
    _ginv: NDArray[np.float64] | None = None
    _dg_inv: NDArray[np.float64] | None = None

    @property
    def spatial_shape(self) -> tuple[int, int, int]:
        return self.g.shape[:3]

    def __post_init__(self) -> None:
        object.__setattr__(self, "_ginv", inverse_metric_4d(self.g))
        object.__setattr__(
            self, "_dg_inv", partial_inverse_metric(self.g, self.spatial_spacing)
        )

    def sample(
        self, x: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        assert self._ginv is not None and self._dg_inv is not None
        i0, frac = clamp_index(x, self.origin, self.spatial_spacing, self.g.shape[:3])
        g_pt = trilinear(self.g, i0, frac)
        ginv_pt = trilinear(self._ginv, i0, frac)
        dg_pt = trilinear(self._dg_inv, i0, frac)
        return g_pt, ginv_pt, dg_pt

    def x_spatial_min(self) -> float:
        return float(self.origin[0])

    def x_spatial_max(self) -> float:
        n0 = self.g.shape[0] - 1
        return float(self.origin[0] + n0 * self.spatial_spacing[0])


@dataclass(frozen=True)
class EvolvingMetricField:
    """Stack ``g_{μν}(t,x,y,z)`` with linear interpolation in simulation time."""

    g_stack: NDArray[np.float64]
    times: NDArray[np.float64]
    origin: NDArray[np.float64]
    spacing: tuple[float, float, float, float]

    @property
    def spatial_shape(self) -> tuple[int, int, int]:
        return self.g_stack.shape[1:4]

    @property
    def spatial_spacing(self) -> tuple[float, float, float]:
        return self.spacing[1], self.spacing[2], self.spacing[3]

    def _time_bracket(self, t: float) -> tuple[int, int, float]:
        times = self.times
        if t <= float(times[0]):
            return 0, 0, 0.0
        if t >= float(times[-1]):
            last = len(times) - 1
            return last, last, 0.0
        idx = int(np.searchsorted(times, t, side="right") - 1)
        idx = min(idx, len(times) - 2)
        t0, t1 = float(times[idx]), float(times[idx + 1])
        alpha = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
        return idx, idx + 1, alpha

    def _sample_g_slice(self, slice_idx: int, x: NDArray[np.float64]) -> NDArray[np.float64]:
        g_slice = self.g_stack[slice_idx]
        i0, frac = clamp_index(x, self.origin, self.spatial_spacing, g_slice.shape[:3])
        return trilinear(g_slice, i0, frac)

    def sample_g(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        i0, i1, alpha = self._time_bracket(float(x[0]))
        g0 = self._sample_g_slice(i0, x)
        if i0 == i1 or alpha <= 0.0:
            return g0
        g1 = self._sample_g_slice(i1, x)
        return (1.0 - alpha) * g0 + alpha * g1

    def sample(
        self, x: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        eps = np.array(
            [
                max(self.spacing[0] * 0.5, 1.0e-6),
                max(self.spatial_spacing[0] * 0.5, 1.0e-6),
                max(self.spatial_spacing[1] * 0.5, 1.0e-6),
                max(self.spatial_spacing[2] * 0.5, 1.0e-6),
            ],
            dtype=float,
        )
        g_pt = self.sample_g(x)
        ginv = inverse_metric_4d(g_pt)
        dg = np.zeros((4, 4, 4), dtype=float)
        for mu in range(4):
            xp = x.astype(float, copy=True)
            xm = x.astype(float, copy=True)
            xp[mu] += eps[mu]
            xm[mu] -= eps[mu]
            ginv_p = inverse_metric_4d(self.sample_g(xp))
            ginv_m = inverse_metric_4d(self.sample_g(xm))
            dg[..., mu] = (ginv_p - ginv_m) / (2.0 * eps[mu])
        return g_pt, ginv, dg

    def x_spatial_min(self) -> float:
        return float(self.origin[0])

    def x_spatial_max(self) -> float:
        n0 = self.g_stack.shape[1] - 1
        return float(self.origin[0] + n0 * self.spatial_spacing[0])


def evolving_field_from_plotfiles(
    plotfiles: Sequence[str | Path],
    *,
    n_space: int = 65,
    half_width: float | None = None,
) -> EvolvingMetricField:
    """Build an evolving field from time-ordered GRTeclyn plotfiles."""
    paths = [str(p) for p in plotfiles]
    if len(paths) < 3:
        raise ValueError("need >= 3 consecutive plotfiles for evolving geodesic trace")

    import yt  # local import: heavy optional dependency

    slices: list[NDArray[np.float64]] = []
    times: list[float] = []
    origin: NDArray[np.float64] | None = None
    spacing_xyz: tuple[float, float, float] | None = None

    from .geodesic import build_metric_3d_from_plotfile

    for path in paths:
        ds = yt.load(path)
        times.append(float(ds.current_time))
        g, orig, sp = build_metric_3d_from_plotfile(
            path, n=n_space, half_width=half_width
        )
        slices.append(g)
        if origin is None:
            origin = orig
            spacing_xyz = sp

    assert origin is not None and spacing_xyz is not None
    times_arr = np.asarray(times, dtype=float)
    dt = float(np.mean(np.diff(times_arr))) if len(times_arr) > 1 else 1.0
    if not np.isfinite(dt) or dt <= 0.0:
        dt = 1.0
    g_stack = np.stack(slices, axis=0)
    return EvolvingMetricField(
        g_stack=g_stack,
        times=times_arr,
        origin=origin,
        spacing=(dt, spacing_xyz[0], spacing_xyz[1], spacing_xyz[2]),
    )


def evolving_field_from_analytic_stack(
    g: NDArray[np.float64],
    spacing: Sequence[float],
    *,
    origin: NDArray[np.float64] | None = None,
) -> EvolvingMetricField:
    """Wrap an analytic ``(Nt,Nx,Ny,Nz,4,4)`` grid (e.g. Alcubierre)."""
    dt, dx, dy, dz = (float(spacing[0]), float(spacing[1]), float(spacing[2]), float(spacing[3]))
    n_time = g.shape[0]
    times = (np.arange(n_time) - n_time // 2) * dt
    if origin is None:
        half = 0.5 * (g.shape[1] - 1) * dx
        origin = np.array([-half, -half, -half], dtype=float)
    return EvolvingMetricField(
        g_stack=g,
        times=times.astype(float),
        origin=origin.astype(float),
        spacing=(dt, dx, dy, dz),
    )
