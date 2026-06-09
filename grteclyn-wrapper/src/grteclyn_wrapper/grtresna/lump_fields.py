"""Analytic GRTresna scalar-lump profiles (shared by motif fit and gridinit export).

Mirrors ``MatterParams.hpp`` in GRTresna so painted per-lump ``phi``/``Pi``
fields match the independent-field stress tensor used in the constraint solve.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

EXOTIC_AMP_SCALE = 0.25
MAX_INDEPENDENT_LUMPS = 5


def shift_lump_centers_for_gridinit(
    lumps: Sequence[Mapping[str, Any]],
    *,
    grid_center: Sequence[float],
) -> list[dict[str, Any]]:
    """Translate shell-ansatz lump centres into gridinit coordinates.

    GRTresna shell lumps are specified as offsets from the throat.  The
    ``.gridinit`` header shifts the file origin so that the throat maps to
    ``grid_center`` in GRTeclyn coordinates, so per-lump scalar channels must be
    painted around that same center.
    """
    cx, cy, cz = (float(grid_center[0]), float(grid_center[1]), float(grid_center[2]))
    shifted: list[dict[str, Any]] = []
    for lump in lumps:
        center = lump.get("center", (0.0, 0.0, 0.0))
        new_lump = dict(lump)
        new_lump["center"] = (
            float(center[0]) + cx,
            float(center[1]) + cy,
            float(center[2]) + cz,
        )
        shifted.append(new_lump)
    return shifted


def effective_amp(lump: Mapping[str, Any]) -> float:
    amp = float(lump.get("amp", 0.0))
    if amp == 0.0:
        return 0.0
    return EXOTIC_AMP_SCALE * amp if int(lump.get("exotic", 0)) else amp


def angular_factor(mode: int, dx: float, dy: float, width: float) -> float:
    if mode == 1:
        return dx / width
    if mode == 2:
        return (dx * dx - dy * dy) / (width * width)
    return 1.0


def lump_phi_at(
    lump: Mapping[str, Any],
    point: Sequence[float],
) -> float:
    amp = float(lump.get("amp", 0.0))
    if amp == 0.0:
        return 0.0
    width = float(lump.get("width", 5.0))
    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=float)
    loc = np.asarray(point, dtype=float)
    dx = loc - center
    r2 = float(np.dot(dx, dx))
    env = math.exp(-r2 / (2.0 * width * width))
    mode = int(lump.get("mode", 0))
    return effective_amp(lump) * angular_factor(mode, dx[0], dx[1], width) * env


def lump_grad_at(
    lump: Mapping[str, Any],
    point: Sequence[float],
) -> np.ndarray:
    amp = float(lump.get("amp", 0.0))
    width = float(lump.get("width", 5.0))
    if amp == 0.0:
        return np.zeros(3, dtype=float)
    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=float)
    loc = np.asarray(point, dtype=float)
    eps = 1.0e-3 * width
    grad = np.zeros(3, dtype=float)
    for axis in range(3):
        lp = loc.copy()
        lm = loc.copy()
        lp[axis] += eps
        lm[axis] -= eps
        grad[axis] = (lump_phi_at(lump, lp) - lump_phi_at(lump, lm)) / (2.0 * eps)
    return grad


def lump_pi_at(
    lump: Mapping[str, Any],
    point: Sequence[float],
) -> float:
    amp = float(lump.get("amp", 0.0))
    velocity = tuple(float(v) for v in lump.get("velocity", (0.0, 0.0, 0.0)))
    omega = float(lump.get("omega", 0.0))
    if amp == 0.0:
        return 0.0
    if not any(abs(v) > 0.0 for v in velocity) and abs(omega) < 1.0e-12:
        return 0.0

    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=float)
    loc = np.asarray(point, dtype=float)
    grad = lump_grad_at(lump, point)
    boost = -float(np.dot(velocity, grad))
    dx = loc - center
    rot = -omega * (dx[0] * grad[1] - dx[1] * grad[0])
    return boost + rot


def lump_sign(lump: Mapping[str, Any]) -> int:
    return -1 if int(lump.get("exotic", 0)) else 1


def _lump_phi_grid(
    lump: Mapping[str, Any],
    px: NDArray,
    py: NDArray,
    pz: NDArray,
) -> NDArray:
    amp = float(lump.get("amp", 0.0))
    if amp == 0.0:
        return np.zeros_like(px, dtype=np.float64)
    width = float(lump.get("width", 5.0))
    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=np.float64)
    dx = px - center[0]
    dy = py - center[1]
    dz = pz - center[2]
    r2 = dx * dx + dy * dy + dz * dz
    env = np.exp(-r2 / (2.0 * width * width))
    mode = int(lump.get("mode", 0))
    angular = angular_factor(mode, dx, dy, width)
    return effective_amp(lump) * angular * env


def _lump_pi_grid(
    lump: Mapping[str, Any],
    px: NDArray,
    py: NDArray,
    pz: NDArray,
) -> NDArray:
    amp = float(lump.get("amp", 0.0))
    velocity = tuple(float(v) for v in lump.get("velocity", (0.0, 0.0, 0.0)))
    omega = float(lump.get("omega", 0.0))
    if amp == 0.0:
        return np.zeros_like(px, dtype=np.float64)
    if not any(abs(v) > 0.0 for v in velocity) and abs(omega) < 1.0e-12:
        return np.zeros_like(px, dtype=np.float64)

    width = float(lump.get("width", 5.0))
    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=np.float64)
    eps = 1.0e-3 * width

    phi_xp = _lump_phi_grid(lump, px + eps, py, pz)
    phi_xm = _lump_phi_grid(lump, px - eps, py, pz)
    phi_yp = _lump_phi_grid(lump, px, py + eps, pz)
    phi_ym = _lump_phi_grid(lump, px, py - eps, pz)
    phi_zp = _lump_phi_grid(lump, px, py, pz + eps)
    phi_zm = _lump_phi_grid(lump, px, py, pz - eps)
    grad_x = (phi_xp - phi_xm) / (2.0 * eps)
    grad_y = (phi_yp - phi_ym) / (2.0 * eps)
    grad_z = (phi_zp - phi_zm) / (2.0 * eps)

    boost = -(velocity[0] * grad_x + velocity[1] * grad_y + velocity[2] * grad_z)
    dx = px - center[0]
    dy = py - center[1]
    rot = -omega * (dx * grad_y - dy * grad_x)
    return boost + rot


def _coordinate_grids(
    nx: int,
    ny: int,
    nz: int,
    dx_xyz: NDArray,
    origin: NDArray,
) -> tuple[NDArray, NDArray, NDArray]:
    px = origin[0] + (np.arange(nx, dtype=np.float64) + 0.5) * float(dx_xyz[0])
    py = origin[1] + (np.arange(ny, dtype=np.float64) + 0.5) * float(dx_xyz[1])
    pz = origin[2] + (np.arange(nz, dtype=np.float64) + 0.5) * float(dx_xyz[2])
    pz_grid, py_grid, px_grid = np.meshgrid(pz, py, px, indexing="ij")
    return px_grid, py_grid, pz_grid


def _lump_momentum_density_grid(
    lump: Mapping[str, Any],
    px: NDArray,
    py: NDArray,
    pz: NDArray,
) -> tuple[NDArray, NDArray, NDArray]:
    """Scalar momentum density ``S_i = sign * Pi * d_i phi`` of one lump.

    This is the matter momentum flux that sources a coordinate shift: a moving
    or rotating lump (non-zero ``velocity`` / ``omega`` -> non-zero ``Pi``)
    carries a directed ``S_i`` that a warp-style shift should follow.
    """
    zero = np.zeros_like(px, dtype=np.float64)
    pi = _lump_pi_grid(lump, px, py, pz)
    if not np.any(pi):
        return zero, zero.copy(), zero.copy()
    width = float(lump.get("width", 5.0))
    eps = 1.0e-3 * width
    grad_x = (_lump_phi_grid(lump, px + eps, py, pz) - _lump_phi_grid(lump, px - eps, py, pz)) / (2.0 * eps)
    grad_y = (_lump_phi_grid(lump, px, py + eps, pz) - _lump_phi_grid(lump, px, py - eps, pz)) / (2.0 * eps)
    grad_z = (_lump_phi_grid(lump, px, py, pz + eps) - _lump_phi_grid(lump, px, py, pz - eps)) / (2.0 * eps)
    sign = float(lump_sign(lump))
    return sign * pi * grad_x, sign * pi * grad_y, sign * pi * grad_z


def paint_shift_seed_on_grid(
    data: NDArray,
    comp_names: list[str],
    dx_xyz: NDArray,
    origin: NDArray,
    lumps: Sequence[Mapping[str, Any]],
    *,
    shift_seed: float,
    max_lumps: int = MAX_INDEPENDENT_LUMPS,
) -> NDArray:
    """Seed an initial shift ``beta^i`` aligned with the matter momentum density.

    GRTresna writes a zero initial shift (a valid gauge choice), so the warp /
    channel mechanism that needs ``beta != 0`` is structurally unreachable at
    t=0.  The shift is a free, constraint-independent kinematic field, so we may
    seed it directly into the gridinit: ``beta^i = shift_seed * S_i / max|S|``,
    giving a peak shift magnitude of ``|shift_seed|`` pointed along (sign of
    ``shift_seed``) the matter momentum flux.  Returns ``data`` unchanged when
    ``shift_seed == 0``, there is no momentum flux, or the gridinit lacks shift
    channels.
    """
    if not lumps or not math.isfinite(shift_seed) or shift_seed == 0.0:
        return data
    name_to_idx = {name: i for i, name in enumerate(comp_names)}
    shift_cols = [name_to_idx.get(n) for n in ("shift1", "shift2", "shift3")]
    if any(col is None for col in shift_cols):
        return data

    nz, ny, nx, _ = data.shape
    px_grid, py_grid, pz_grid = _coordinate_grids(nx, ny, nz, dx_xyz, origin)
    s_x = np.zeros((nz, ny, nx), dtype=np.float64)
    s_y = np.zeros((nz, ny, nx), dtype=np.float64)
    s_z = np.zeros((nz, ny, nx), dtype=np.float64)
    for lump in list(lumps[:max_lumps]):
        lx, ly, lz = _lump_momentum_density_grid(lump, px_grid, py_grid, pz_grid)
        s_x += lx
        s_y += ly
        s_z += lz

    peak = float(np.sqrt(s_x * s_x + s_y * s_y + s_z * s_z).max())
    if peak <= 0.0:
        return data
    scale = float(shift_seed) / peak

    data = np.array(data, dtype=np.float64, copy=True)
    data[:, :, :, shift_cols[0]] = scale * s_x
    data[:, :, :, shift_cols[1]] = scale * s_y
    data[:, :, :, shift_cols[2]] = scale * s_z
    return data


def paint_lump_fields_on_grid(
    data: NDArray,
    comp_names: list[str],
    dx_xyz: NDArray,
    origin: NDArray,
    lumps: Sequence[Mapping[str, Any]],
    *,
    max_lumps: int = MAX_INDEPENDENT_LUMPS,
) -> tuple[NDArray, list[str]]:
    """Append per-lump ``phi_lumpK`` / ``Pi_lumpK`` channels to a uniform grid."""
    if not lumps:
        return data, comp_names

    nz, ny, nx, _ = data.shape
    selected = list(lumps[:max_lumps])
    new_names = list(comp_names)
    lump_arrays: list[tuple[NDArray, NDArray]] = []
    px_grid, py_grid, pz_grid = _coordinate_grids(nx, ny, nz, dx_xyz, origin)

    for idx, lump in enumerate(selected):
        phi_plane = _lump_phi_grid(lump, px_grid, py_grid, pz_grid)
        pi_plane = _lump_pi_grid(lump, px_grid, py_grid, pz_grid)
        lump_arrays.append((phi_plane, pi_plane))
        new_names.extend([f"phi_lump{idx}", f"Pi_lump{idx}"])

    extra = np.zeros((nz, ny, nx, len(lump_arrays) * 2), dtype=np.float64)
    for idx, (phi_plane, pi_plane) in enumerate(lump_arrays):
        extra[:, :, :, 2 * idx] = phi_plane
        extra[:, :, :, 2 * idx + 1] = pi_plane

    return np.concatenate([data, extra], axis=3), new_names
