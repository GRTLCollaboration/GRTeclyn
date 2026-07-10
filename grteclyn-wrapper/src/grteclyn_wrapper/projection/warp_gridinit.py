"""Analytic Alcubierre warp-bubble gridinit writer and post-solve painter.

Phase A of the Alcubierre validation baseline (see ``GeometryFirst.md``).

A0 — ``write_alcubierre_gridinit`` produces a *ground-truth* analytic
Alcubierre t=0 slice on a uniform grid (flat ``gamma_ij``, ``chi=1``,
``lapse=1``, painted axial shift ``beta^x = -v f(r_s)``, shift-consistent
``A_ij``, zero matter).  This is the calibration anchor for the
preservation checker and the FTL probe before any reconstruction is
attempted.

A1 — ``paint_alcubierre_warp_on_gridinit`` reads a *solved* gridinit and
overwrites the shift (and optionally ``A_ij``) with the analytic Alcubierre
values, leaving the solved ``chi`` / matter intact.  This is legitimate
because the shift is pure gauge: the constraints fix ``(gamma_ij, K_ij,
matter)`` and never ``beta^i``.  The real quantity to match is ``A_ij``
(which the solver sources from the matter momentum density ``S_i``).

Both functions reuse:
  - ``alcubierre_shape_function`` from ``initial_data/motif.py``
  - ``read_gridinit`` / ``write_gridinit`` / ``GRTECLYN_STATE_VARS`` from
    ``grtresna/io/gridinit.py``
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..grtresna.io.gridinit import (
    GRTECLYN_STATE_VARS,
    read_gridinit,
    write_gridinit,
)
from ..initial_data.motif import alcubierre_shape_function

# Default grid for the analytic anchor (matches the documented test config).
DEFAULT_ANALYTIC_N: int = 32
DEFAULT_ANALYTIC_L: float = 16.0


def _state_var_index(comp_names: list[str], name: str) -> int:
    """Return the component index for *name* or raise."""
    try:
        return comp_names.index(name)
    except ValueError as exc:  # pragma: no cover - defensive
        raise KeyError(f"gridinit is missing required component '{name}'") from exc


def _coordinate_grids(
    nx: int, ny: int, nz: int,
    dx_xyz: NDArray, origin: NDArray,
) -> tuple[NDArray, NDArray, NDArray]:
    """Physical coordinate arrays (px, py, pz) of shape (nz, ny, nx).

    Matches the convention in ``grtresna/fields/lump.py``: cell centres at
    ``origin + (i + 0.5) * dx`` and ``meshgrid(pz, py, px, indexing="ij")``.
    """
    px = float(origin[0]) + (np.arange(nx, dtype=np.float64) + 0.5) * float(dx_xyz[0])
    py = float(origin[1]) + (np.arange(ny, dtype=np.float64) + 0.5) * float(dx_xyz[1])
    pz = float(origin[2]) + (np.arange(nz, dtype=np.float64) + 0.5) * float(dx_xyz[2])
    pz_grid, py_grid, px_grid = np.meshgrid(pz, py, px, indexing="ij")
    return px_grid, py_grid, pz_grid


def alcubierre_analytic_fields(
    shape_xyz: tuple[int, int, int],
    dx_xyz: NDArray,
    origin: NDArray,
    *,
    velocity: float,
    bubble_radius: float,
    sigma: float,
    center: tuple[float, float, float],
) -> dict[str, NDArray]:
    """Evaluate the exact Alcubierre t=0 slice on a uniform grid.

    Returns a dict of arrays keyed by ``GRTECLYN_STATE_VARS`` names:
      ``chi = 1``, ``h11=h22=h33 = 1`` (off-diagonal ``h_ij = 0``),
      ``lapse = 1``, ``K = 0``,
      ``shift1 = -v * f(r_s)``, ``shift2 = shift3 = 0``,
      ``A_ij = (1/2)(d_i beta_j + d_j beta_i) - (1/3) delta_ij d_k beta_k``
    (with ``alpha=1`` and a flat conformal metric, ``K=0`` implies
    ``A_ij = K_ij``; the ADM relation ``K_ij = (D_i beta_j +
    D_j beta_i) / (2 alpha)`` reduces to the symmetric gradient above).
    Derivatives of ``beta`` are computed with ``np.gradient`` on the
    painted shift array (robust; no closed form needed at ``r_s = 0``).

    All matter fields (``phi``, ``Pi``, ``teo_*``) and the gauge auxiliaries
    (``Theta``, ``Gamma_i``, ``B_i``) are zero.
    """
    nx, ny, nz = shape_xyz
    px, py, pz = _coordinate_grids(nx, ny, nz, dx_xyz, origin)
    cx, cy, cz = float(center[0]), float(center[1]), float(center[2])

    # r_s = distance from the bubble centre
    dxs = px - cx
    dys = py - cy
    dzs = pz - cz
    r_s = np.sqrt(dxs * dxs + dys * dys + dzs * dzs)

    # Shift: beta^x = -v f(r_s), beta^y = beta^z = 0
    f = alcubierre_shape_function(r_s, bubble_radius=bubble_radius, sigma=sigma)
    beta_x = -float(velocity) * f
    beta_y = np.zeros_like(beta_x)
    beta_z = np.zeros_like(beta_x)

    # A_ij from the symmetric gradient of beta (K=0 => A_ij = K_ij).
    # np.gradient on an (nz, ny, nx) array returns (d/dz, d/dy, d/dx).
    d_beta_x_dz, d_beta_x_dy, d_beta_x_dx = np.gradient(
        beta_x, float(dx_xyz[2]), float(dx_xyz[1]), float(dx_xyz[0]),
    )
    # beta_y = beta_z = 0 => their gradients vanish
    div_beta = d_beta_x_dx  # + d_beta_y_dy + d_beta_z_dz  (those are zero)

    A11 = d_beta_x_dx - (1.0 / 3.0) * div_beta
    A22 = -(1.0 / 3.0) * div_beta
    A33 = -(1.0 / 3.0) * div_beta
    A12 = 0.5 * d_beta_x_dy
    A13 = 0.5 * d_beta_x_dz
    A23 = 0.0

    fields: dict[str, NDArray] = {}
    for name in GRTECLYN_STATE_VARS:
        fields[name] = np.zeros((nz, ny, nx), dtype=np.float64)
    # Flat spatial metric + lapse + chi
    fields["chi"] = np.ones((nz, ny, nx), dtype=np.float64)
    fields["h11"] = np.ones((nz, ny, nx), dtype=np.float64)
    fields["h22"] = np.ones((nz, ny, nx), dtype=np.float64)
    fields["h33"] = np.ones((nz, ny, nx), dtype=np.float64)
    fields["lapse"] = np.ones((nz, ny, nx), dtype=np.float64)
    # Shift
    fields["shift1"] = beta_x
    fields["shift2"] = beta_y
    fields["shift3"] = beta_z
    # Extrinsic curvature (K=0, A_ij from the shift gradient)
    fields["K"] = np.zeros((nz, ny, nx), dtype=np.float64)
    fields["A11"] = A11
    fields["A12"] = A12
    fields["A13"] = A13
    fields["A22"] = A22
    fields["A23"] = np.zeros((nz, ny, nx), dtype=np.float64)
    fields["A33"] = A33
    return fields


def write_alcubierre_gridinit(
    out_path: str | Path,
    *,
    n: int = DEFAULT_ANALYTIC_N,
    L: float = DEFAULT_ANALYTIC_L,
    velocity: float,
    bubble_radius: float,
    sigma: float,
) -> Path:
    """Write a full analytic Alcubierre gridinit (zero matter).

    The grid is ``n^3`` cells covering the box ``[0, 2L]^3`` with the bubble
    centred at ``(L, L, L)`` and ``origin = (0, 0, 0)`` — matching the
    convention used by ``grtresna/io/conversion.py`` for full-box solves
    (``target_center = (L/2, L/2, L/2)`` in physics coords, which maps to
    the gridinit centre).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    nx = ny = nz = n
    dx = (2.0 * float(L)) / float(n)
    dx_xyz = np.array([dx, dx, dx], dtype=np.float64)
    origin = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    center = (float(L), float(L), float(L))

    fields = alcubierre_analytic_fields(
        (nx, ny, nz), dx_xyz, origin,
        velocity=velocity, bubble_radius=bubble_radius, sigma=sigma,
        center=center,
    )
    # Stack into (nz, ny, nx, n_comp) in GRTECLYN_STATE_VARS order
    data = np.stack([fields[name] for name in GRTECLYN_STATE_VARS], axis=-1)
    return write_gridinit(data, list(GRTECLYN_STATE_VARS), dx_xyz, origin, out_path)


def paint_alcubierre_warp_on_gridinit(
    gridinit_path: str | Path,
    *,
    velocity: float,
    bubble_radius: float,
    sigma: float,
    paint_aij: bool = True,
) -> Path:
    """Overwrite the shift (and optionally ``A_ij``) on a solved gridinit.

    Reads the gridinit produced by GRTresna (which has ``shift = 0`` and a
    solver-produced ``A_ij``), paints the analytic Alcubierre shift
    ``beta^x = -v f(r_s)`` centred on the grid centre, and optionally
    overwrites ``A_ij`` with the shift-consistent analytic values.  The
    solved ``chi`` / matter / lapse are left intact.  Writes back in place
    and returns the path.

    With ``paint_aij=True`` the Hamiltonian constraint acquires an error
    because ``chi`` was solved with the *solver's* ``A_ij``, not the painted
    one — this trade-off is expected and should be logged by the caller.
    """
    gridinit_path = Path(gridinit_path)
    grid = read_gridinit(gridinit_path)
    nz, ny, nx, _ = grid.data.shape
    center = (
        float(grid.origin[0]) + 0.5 * nx * float(grid.dx_xyz[0]),
        float(grid.origin[1]) + 0.5 * ny * float(grid.dx_xyz[1]),
        float(grid.origin[2]) + 0.5 * nz * float(grid.dx_xyz[2]),
    )
    fields = alcubierre_analytic_fields(
        (nx, ny, nz), grid.dx_xyz, grid.origin,
        velocity=velocity, bubble_radius=bubble_radius, sigma=sigma,
        center=center,
    )
    data = grid.data.copy()
    for name in ("shift1", "shift2", "shift3"):
        idx = _state_var_index(grid.comp_names, name)
        data[:, :, :, idx] = fields[name]
    if paint_aij:
        for name in ("A11", "A12", "A13", "A22", "A23", "A33"):
            try:
                idx = _state_var_index(grid.comp_names, name)
            except KeyError:
                continue
            data[:, :, :, idx] = fields[name]
    return write_gridinit(data, list(grid.comp_names), grid.dx_xyz, grid.origin, gridinit_path)


def alcubierre_analytic_Si(
    *,
    velocity: float,
    bubble_radius: float,
    sigma: float,
    L: float,
    n: int = 64,
) -> tuple[NDArray, NDArray, NDArray, NDArray, NDArray]:
    """Analytic Eulerian momentum density ``S_i`` of an Alcubierre bubble.

    Returns ``(x_axis, z_axis, S_x, S_y, S_z)`` on an xz-slice (y=0) with
    ``x_axis`` and ``z_axis`` of length ``n``.  The momentum density is
    derived from the exact 4-metric via ``T^{0i}`` (equivalently
    ``S_i = -n_mu T^{mu i}`` with ``n^mu = (1, -beta)/alpha``); for the
    Alcubierre slice (``alpha=1``, flat ``gamma``) this reduces to
    ``S_i = -beta^j K_{ij}`` with ``K_{ij}`` the symmetric gradient of
    ``beta``.  The dominant component is ``S_x``; ``S_y`` and ``S_z`` are
    non-zero only through the off-diagonal ``K_{ij}`` terms.

    This is the target the matter fitter (A2) and the ``si_l2`` mismatch
    term compare against.
    """
    x_axis = np.linspace(-L, L, n)
    z_axis = np.linspace(-L, L, n)
    X, Z = np.meshgrid(x_axis, z_axis, indexing="ij")
    r_s = np.sqrt(X * X + Z * Z)

    f = alcubierre_shape_function(r_s, bubble_radius=bubble_radius, sigma=sigma)
    beta_x = -float(velocity) * f

    # df/dr_s via finite differences on a 1D radial grid (robust at r=0)
    r_1d = np.linspace(max(L / n, 1.0e-6), L * math.sqrt(2.0), n)
    f_1d = alcubierre_shape_function(r_1d, bubble_radius=bubble_radius, sigma=sigma)
    dr = r_1d[1] - r_1d[0]
    df_dr_1d = np.gradient(f_1d, dr, edge_order=2)
    df_dr = np.interp(r_s.ravel(), r_1d, df_dr_1d).reshape(n, n)

    r_safe = np.clip(r_s, 1.0e-10, None)
    # K_ij = (1/2)(d_i beta_j + d_j beta_i); beta = (beta_x, 0, 0)
    # On the xz-slice (y=0): d_x beta_x = -v df/dr * x/r, d_z beta_x = -v df/dr * z/r
    d_x_beta_x = -float(velocity) * df_dr * X / r_safe
    d_z_beta_x = -float(velocity) * df_dr * Z / r_safe
    K_xx = d_x_beta_x
    K_xz = 0.5 * d_z_beta_x
    K_zx = K_xz

    # S_i = -beta^j K_{ij} (alpha=1, gamma=delta => n_mu=(1,0,0,0), S_i = -n_mu T^{mu i}
    # with T^{0i} = -beta^j K_{ij}/(8 pi) ... but we want the *normalisation-free*
    # shape; the si_l2 term rescales anyway.  Use S_i = -beta^j K_{ij}.
    S_x = -beta_x * K_xx  # beta^z K_zx = 0 (beta_z=0)
    S_y = np.zeros((n, n), dtype=np.float64)  # beta_y = 0
    S_z = -beta_x * K_zx
    return x_axis, z_axis, S_x, S_y, S_z
