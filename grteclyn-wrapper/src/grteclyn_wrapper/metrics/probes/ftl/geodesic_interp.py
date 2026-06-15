"""Shared metric interpolation helpers for geodesic probes."""

from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from ..warpfactory import _d_dx


def inverse_metric_4d(g: NDArray[np.float64]) -> NDArray[np.float64]:
    """Invert a batch of 4x4 metrics.  Shape ``(..., 4, 4)`` -> same."""
    return np.linalg.inv(g)


def partial_inverse_metric(
    g: NDArray[np.float64],
    spacing: Sequence[float],
) -> NDArray[np.float64]:
    """``partial_mu g^{ab}`` for a static spatial grid (``∂_t = 0``)."""
    ginv = inverse_metric_4d(g)
    out = np.zeros(g.shape[:-2] + (4, 4, 4), dtype=float)
    for mu in range(1, 4):
        out[..., :, :, mu] = _d_dx(ginv, float(spacing[mu - 1]), axis=mu - 1)
    return out


def clamp_index(
    x: NDArray[np.float64],
    origin: NDArray[np.float64],
    spacing: Sequence[float],
    shape: Sequence[int],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Continuous grid coordinates and fractional cell position."""
    rel = (x[1:] - origin[:3]) / np.asarray(spacing, dtype=float)
    n = np.array(shape, dtype=float)
    rel_clamped = np.clip(rel, 0.0, n - 1.001)
    i0 = np.floor(rel_clamped).astype(int)
    frac = rel_clamped - i0
    return i0, frac


def trilinear(
    field: NDArray[np.float64],
    i0: NDArray[np.int64],
    frac: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Trilinear interpolation on a 3D field with trailing dimensions."""
    i0 = np.clip(i0, 0, np.array(field.shape[:3]) - 2)
    fx, fy, fz = frac[0], frac[1], frac[2]
    c000 = field[i0[0], i0[1], i0[2]]
    c100 = field[i0[0] + 1, i0[1], i0[2]]
    c010 = field[i0[0], i0[1] + 1, i0[2]]
    c110 = field[i0[0] + 1, i0[1] + 1, i0[2]]
    c001 = field[i0[0], i0[1], i0[2] + 1]
    c101 = field[i0[0] + 1, i0[1], i0[2] + 1]
    c011 = field[i0[0], i0[1] + 1, i0[2] + 1]
    c111 = field[i0[0] + 1, i0[1] + 1, i0[2] + 1]
    c00 = c000 * (1 - fx) + c100 * fx
    c10 = c010 * (1 - fx) + c110 * fx
    c01 = c001 * (1 - fx) + c101 * fx
    c11 = c011 * (1 - fx) + c111 * fx
    c0 = c00 * (1 - fy) + c10 * fy
    c1 = c01 * (1 - fy) + c11 * fy
    return c0 * (1 - fz) + c1 * fz
