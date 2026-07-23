"""Teo rotating traversable wormhole initial data (Teo 1998, gr-qc/9803098).

Physics
-------
The Teo geometry is the stationary, axisymmetric generalization of the
Morris-Thorne wormhole:

    ds^2 = -N^2 dt^2 + (1 - b/r)^-1 dr^2
         + r^2 K^2 [dtheta^2 + sin^2(theta) (dphi - omega dt)^2]

with the worked example  N = K = 1 + (4 a cos(theta))^2 / r,  b = b0,
omega = 2 a / r^3.  This module samples that geometry on a uniform Cartesian
grid, performs the 3+1 (ADM) decomposition, and converts it to the CCZ4 state
GRTeclyn evolves.

Teo specifies a geometry and *infers* its stress-energy.  The effective source
``T_ab = G_ab / 8pi`` is therefore computed numerically from the 4-metric, not
from a fundamental matter field.  Callers that only want geometry-first controls
can request a zeroed source.

This module is intentionally pure (no argv/IO side effects).  The CLI launcher
lives in ``scripts/make_teo_wormhole_gridinit.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .adm_stationary import (
    ConstraintResiduals,
    StationaryGrid,
    christoffel,
    conformal_connection,
    constraint_residuals_from_adm,
    einstein_source,
    extrinsic_curvature,
    pack_ccz4_grid,
    ricci_tensor_3d,
    spatial_gradient,
)

TINY = 1.0e-12
INV_8PI = 1.0 / (8.0 * np.pi)

SourceModel = Literal["einstein", "zero"]

# Backward-compatible aliases used by older call sites / tests.
_grad = spatial_gradient
_christoffel = christoffel
_ricci_tensor_3d = ricci_tensor_3d


@dataclass
class TeoWormholeConfig:
    """Configuration for one Teo wormhole grid.

    The active grid (``nx, ny, nz`` over ``Lx, Ly, Lz``) must match the GRTeclyn
    box that loads the resulting ``.gridinit``.  ``origin`` is written into the
    file header so the loader maps absolute coordinates ``(i+0.5)*dx`` onto the
    same throat location used by the simulation ``center``.
    """

    nx: int = 64
    ny: int = 64
    nz: int = 64
    Lx: float = 64.0
    Ly: float = 64.0
    Lz: float = 64.0
    center: tuple[float, float, float] = (32.0, 32.0, 0.0)
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    throat_radius: float = 0.5
    spin: float = 0.05
    source: SourceModel = "einstein"

    @property
    def dx(self) -> tuple[float, float, float]:
        return (self.Lx / self.nx, self.Ly / self.ny, self.Lz / self.nz)


@dataclass
class TeoGrid:
    """Generated uniform grid plus summary diagnostics."""

    data: NDArray[np.float64]  # (nz, ny, nx, n_comp)
    comp_names: list[str]
    dx_xyz: NDArray[np.float64]  # (3,)
    origin: NDArray[np.float64]  # (3,)
    metrics: dict[str, float]


def build_adm(cfg: TeoWormholeConfig) -> tuple[NDArray, NDArray, NDArray]:
    """Sample the Teo metric and return (alpha, beta^i, gamma_ij).

    Uses the isotropic Schwarzschild-bridge radial map so the a=0 limit is the
    zero-tidal-force (Schwarzschild-like) wormhole of the Teo example.
    """
    dx = cfg.dx
    x = (np.arange(cfg.nx, dtype=np.float64) + 0.5) * dx[0]
    y = (np.arange(cfg.ny, dtype=np.float64) + 0.5) * dx[1]
    z = (np.arange(cfg.nz, dtype=np.float64) + 0.5) * dx[2]
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")

    rel_x = xx - cfg.center[0]
    rel_y = yy - cfg.center[1]
    rel_z = zz - cfg.center[2]
    rho = np.sqrt(rel_x * rel_x + rel_y * rel_y + rel_z * rel_z)
    rho_safe = np.maximum(rho, TINY)

    b0 = cfg.throat_radius
    psi = 1.0 + b0 / (4.0 * rho_safe)
    areal_r = rho_safe * psi * psi
    r_hat = np.maximum(areal_r / b0, TINY)
    cos_theta = rel_z / rho_safe

    teo_factor = 1.0 + (4.0 * cfg.spin * cos_theta) ** 2 / r_hat
    alpha = teo_factor
    K_potential = teo_factor
    omega = 2.0 * cfg.spin / (b0 * r_hat**3)

    radial_metric = psi**4
    angular_metric = areal_r * areal_r * K_potential * K_potential
    tangential_factor = angular_metric / np.maximum(rho_safe * rho_safe, TINY)

    n = np.stack([rel_x / rho_safe, rel_y / rho_safe, rel_z / rho_safe], axis=-1)
    ident = np.eye(3)
    gamma = np.zeros(rho.shape + (3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            gamma[..., i, j] = radial_metric * n[..., i] * n[..., j] + (
                tangential_factor * (ident[i, j] - n[..., i] * n[..., j])
            )

    # The angular decomposition is ill-conditioned on the puncture cell; use the
    # radial value isotropically there and rely on the chi/lapse floors.
    near_puncture = rho < 0.5 * min(dx)
    for i in range(3):
        for j in range(3):
            gamma[..., i, j] = np.where(
                near_puncture, radial_metric * ident[i, j], gamma[..., i, j]
            )

    beta = np.zeros(rho.shape + (3,), dtype=np.float64)
    beta[..., 0] = omega * rel_y
    beta[..., 1] = -omega * rel_x
    return alpha, beta, gamma


def constraint_residuals(
    cfg: TeoWormholeConfig, *, trim_cells: int = 3, puncture_cells: float = 2.0
) -> ConstraintResiduals:
    """Evaluate ADM constraint residuals for the Teo data."""
    alpha, beta, gamma = build_adm(cfg)
    return constraint_residuals_from_adm(
        alpha,
        beta,
        gamma,
        cfg.dx,
        center=cfg.center,
        trim_cells=trim_cells,
        puncture_cells=puncture_cells,
        origin=cfg.origin,
    )


def build_grid(cfg: TeoWormholeConfig) -> TeoGrid:
    """Assemble the full CCZ4 + effective-source state for a Teo wormhole."""
    alpha, beta, gamma = build_adm(cfg)
    packed: StationaryGrid = pack_ccz4_grid(
        alpha,
        beta,
        gamma,
        cfg.dx,
        cfg.origin,
        include_einstein_source=(cfg.source == "einstein"),
    )
    return TeoGrid(
        data=packed.data,
        comp_names=packed.comp_names,
        dx_xyz=packed.dx_xyz,
        origin=packed.origin,
        metrics=packed.metrics,
    )


__all__ = [
    "ConstraintResiduals",
    "SourceModel",
    "TINY",
    "INV_8PI",
    "TeoGrid",
    "TeoWormholeConfig",
    "build_adm",
    "build_grid",
    "christoffel",
    "conformal_connection",
    "constraint_residuals",
    "einstein_source",
    "extrinsic_curvature",
    "ricci_tensor_3d",
    "spatial_gradient",
    "_christoffel",
    "_grad",
    "_ricci_tensor_3d",
]
