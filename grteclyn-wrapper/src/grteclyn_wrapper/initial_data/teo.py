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

from ..grtresna.io import GRTECLYN_STATE_VARS

TINY = 1.0e-12
INV_8PI = 1.0 / (8.0 * np.pi)

SourceModel = Literal["einstein", "zero"]


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


@dataclass
class ConstraintResiduals:
    """ADM Hamiltonian / momentum constraint residuals of the generated data.

    These are a *self-consistency* check: the effective source is T_ab = G_ab/8pi
    from the same metric, so the constraints hold by construction up to
    finite-difference truncation error. Small, resolution-decreasing residuals
    indicate the K_ij / source derivatives are mutually consistent; large ones
    flag a bug or under-resolution.
    """

    ham_l2: float
    ham_max: float
    mom_l2: float
    mom_max: float
    valid_cells: int
    valid_fraction: float

    def as_dict(self) -> dict[str, float]:
        return {
            "ham_l2": self.ham_l2,
            "ham_max": self.ham_max,
            "mom_l2": self.mom_l2,
            "mom_max": self.mom_max,
            "valid_cells": float(self.valid_cells),
            "valid_fraction": self.valid_fraction,
        }


def _grad(field: NDArray, dx: tuple[float, float, float]) -> list[NDArray]:
    """Return [d/dx, d/dy, d/dz] for an array indexed (z, y, x, ...)."""
    edge_order = 2 if min(field.shape[:3]) > 2 else 1
    d_dz, d_dy, d_dx = np.gradient(field, dx[2], dx[1], dx[0], edge_order=edge_order)
    return [d_dx, d_dy, d_dz]


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


def _christoffel(metric: NDArray, metric_inv: NDArray, dx, dim: int) -> NDArray:
    """Levi-Civita connection for a (dim x dim) metric on a spatial grid.

    For 4-metrics, time derivatives vanish (stationarity), so only spatial
    derivatives contribute.
    """
    dmetric: list = []
    for a in range(dim):
        row = []
        for b in range(dim):
            row.append(_grad(metric[..., a, b], dx))
        dmetric.append(row)

    def d(idx_lo: int, a: int, b: int) -> NDArray | float:
        # spatial coordinate index idx_lo: 0->x,1->y,2->z for 3D; for 4D the
        # coordinate 0 is time (stationary -> 0) and 1..3 map to x,y,z.
        if dim == 4:
            if idx_lo == 0:
                return 0.0
            return dmetric[a][b][idx_lo - 1]
        return dmetric[a][b][idx_lo]

    chris = np.zeros(metric.shape[:-2] + (dim, dim, dim), dtype=np.float64)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                acc = 0.0
                for e in range(dim):
                    acc = acc + metric_inv[..., a, e] * (
                        d(c, e, b) + d(b, e, c) - d(e, b, c)
                    )
                chris[..., a, b, c] = 0.5 * acc
    return chris


def extrinsic_curvature(alpha, beta, gamma, dx) -> tuple[NDArray, NDArray]:
    """Stationary K_ij = (D_i beta_j + D_j beta_i) / (2 alpha); return (K, A_phys)."""
    gamma_inv = np.linalg.inv(gamma)
    beta_low = np.einsum("...ij,...j->...i", gamma, beta)
    dbeta = [_grad(beta_low[..., i], dx) for i in range(3)]
    chris = _christoffel(gamma, gamma_inv, dx, dim=3)

    Kij = np.zeros_like(gamma)
    for i in range(3):
        for j in range(3):
            cov = dbeta[j][i] + dbeta[i][j]
            for k in range(3):
                cov = cov - 2.0 * chris[..., k, i, j] * beta_low[..., k]
            Kij[..., i, j] = 0.5 * cov / alpha

    K = np.einsum("...ij,...ij->...", gamma_inv, Kij)
    A_phys = Kij - (K[..., None, None] / 3.0) * gamma
    return K, A_phys


def conformal_connection(h: NDArray, h_inv: NDArray, dx) -> NDArray:
    """Conformal connection Gamma^i = h^{jk} Gamma^i_{jk} from the unit-det metric."""
    chris = _christoffel(h, h_inv, dx, dim=3)
    return np.einsum("...jk,...ijk->...i", h_inv, chris)


def einstein_source(alpha, beta, gamma, dx) -> tuple[NDArray, NDArray, NDArray]:
    """Effective (rho, j_i, S_ij) = (G_ab / 8pi) projected into the 3+1 frame."""
    beta_low = np.einsum("...ij,...j->...i", gamma, beta)
    beta_sq = np.einsum("...i,...i->...", beta_low, beta)

    g4 = np.zeros(gamma.shape[:-2] + (4, 4), dtype=np.float64)
    g4[..., 0, 0] = beta_sq - alpha * alpha
    for i in range(3):
        g4[..., 0, i + 1] = g4[..., i + 1, 0] = beta_low[..., i]
        for j in range(3):
            g4[..., i + 1, j + 1] = gamma[..., i, j]
    g4_inv = np.linalg.inv(g4)

    chris = _christoffel(g4, g4_inv, dx, dim=4)

    def dchris(idx_lo: int, a: int, b: int, c: int) -> NDArray | float:
        if idx_lo == 0:
            return 0.0
        return _grad(chris[..., a, b, c], dx)[idx_lo - 1]

    trace = np.zeros(g4.shape[:-2] + (4,), dtype=np.float64)
    for a in range(4):
        for c in range(4):
            trace[..., a] += chris[..., c, a, c]

    def dtrace(idx_lo: int, a: int) -> NDArray | float:
        if idx_lo == 0:
            return 0.0
        return _grad(trace[..., a], dx)[idx_lo - 1]

    ricci = np.zeros_like(g4)
    for a in range(4):
        for b in range(4):
            term1 = 0.0
            for c in range(4):
                term1 = term1 + dchris(c, c, a, b)
            term2 = dtrace(b, a)
            term3 = 0.0
            term4 = 0.0
            for c in range(4):
                for e in range(4):
                    term3 = term3 + chris[..., c, a, b] * chris[..., e, c, e]
                    term4 = term4 + chris[..., e, a, c] * chris[..., c, b, e]
            ricci[..., a, b] = term1 - term2 + term3 - term4

    scalar = np.einsum("...ab,...ab->...", g4_inv, ricci)
    einstein = ricci - 0.5 * g4 * scalar[..., None, None]

    n_up = np.zeros(beta.shape[:-1] + (4,), dtype=np.float64)
    n_up[..., 0] = 1.0 / alpha
    n_up[..., 1:] = -beta / alpha[..., None]

    rho = np.einsum("...a,...ab,...b->...", n_up, einstein, n_up) * INV_8PI
    j = np.zeros(beta.shape, dtype=np.float64)
    for i in range(3):
        j[..., i] = -np.einsum("...a,...a->...", n_up, einstein[..., i + 1, :]) * INV_8PI
    S = einstein[..., 1:, 1:] * INV_8PI
    return rho, j, S


def _ricci_tensor_3d(gamma: NDArray, gamma_inv: NDArray, dx) -> NDArray:
    """Spatial Ricci tensor R_ij from the 3-metric (second derivatives)."""
    chris = _christoffel(gamma, gamma_inv, dx, dim=3)
    ricci = np.zeros_like(gamma)
    for i in range(3):
        for j in range(3):
            t1 = 0.0  # d_k Gamma^k_ij
            t2 = 0.0  # d_j Gamma^k_ik
            for k in range(3):
                t1 = t1 + _grad(chris[..., k, i, j], dx)[k]
                t2 = t2 + _grad(chris[..., k, i, k], dx)[j]
            t3 = 0.0  # Gamma^k_ij Gamma^m_km
            t4 = 0.0  # Gamma^m_ik Gamma^k_jm
            for k in range(3):
                for m in range(3):
                    t3 = t3 + chris[..., k, i, j] * chris[..., m, k, m]
                    t4 = t4 + chris[..., m, i, k] * chris[..., k, j, m]
            ricci[..., i, j] = t1 - t2 + t3 - t4
    return ricci


def constraint_residuals(
    cfg: TeoWormholeConfig, *, trim_cells: int = 3, puncture_cells: float = 2.0
) -> ConstraintResiduals:
    """Evaluate ADM constraint residuals for the Teo data.

    Hamiltonian:  H = R + K^2 - K_ij K^ij - 16 pi rho
    Momentum:     M_i = D_j (K^j_i - delta^j_i K) - 8 pi j_i

    The boundary ``trim_cells`` layers (where finite differences are one-sided)
    and a ``puncture_cells``-radius ball around the throat are excluded before
    reducing to L2 (RMS) and max norms.
    """
    alpha, beta, gamma = build_adm(cfg)
    dx = cfg.dx
    gamma_inv = np.linalg.inv(gamma)

    K_trace, A_phys = extrinsic_curvature(alpha, beta, gamma, dx)
    Kij = A_phys + (K_trace[..., None, None] / 3.0) * gamma
    rho, j_low, _S = einstein_source(alpha, beta, gamma, dx)

    # Hamiltonian constraint.
    ricci = _ricci_tensor_3d(gamma, gamma_inv, dx)
    R_scalar = np.einsum("...ij,...ij->...", gamma_inv, ricci)
    K_up = np.einsum("...ia,...jb,...ab->...ij", gamma_inv, gamma_inv, Kij)
    KijKij = np.einsum("...ij,...ij->...", Kij, K_up)
    ham = R_scalar + K_trace * K_trace - KijKij - 16.0 * np.pi * rho

    # Momentum constraint M_i = D_j K^j_i - D_i K - 8 pi j_i.
    chris = _christoffel(gamma, gamma_inv, dx, dim=3)
    K_mixed = np.einsum("...jk,...ki->...ji", gamma_inv, Kij)  # K^j_i
    mom = np.zeros(beta.shape, dtype=np.float64)
    for i in range(3):
        div = 0.0
        for jj in range(3):
            div = div + _grad(K_mixed[..., jj, i], dx)[jj]
            for ll in range(3):
                div = div + chris[..., jj, jj, ll] * K_mixed[..., ll, i]
                div = div - chris[..., ll, jj, i] * K_mixed[..., jj, ll]
        di_K = _grad(K_trace, dx)[i]
        mom[..., i] = div - di_K - 8.0 * np.pi * j_low[..., i]

    # Build the valid-region mask: trim boundaries and the puncture ball.
    nz, ny, nx = alpha.shape
    mask = np.zeros(alpha.shape, dtype=bool)
    t = trim_cells
    mask[t : nz - t, t : ny - t, t : nx - t] = True

    x = (np.arange(nx) + 0.5) * dx[0]
    y = (np.arange(ny) + 0.5) * dx[1]
    z = (np.arange(nz) + 0.5) * dx[2]
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    radius = np.sqrt(
        (xx - cfg.center[0]) ** 2
        + (yy - cfg.center[1]) ** 2
        + (zz - cfg.center[2]) ** 2
    )
    mask &= radius > puncture_cells * min(dx)

    valid = int(np.count_nonzero(mask))
    if valid == 0:
        raise ValueError("Constraint mask excluded every cell; grid too small")

    ham_v = ham[mask]
    mom2 = np.einsum("...i,...i->...", mom, mom)[mask]
    return ConstraintResiduals(
        ham_l2=float(np.sqrt(np.mean(ham_v * ham_v))),
        ham_max=float(np.max(np.abs(ham_v))),
        mom_l2=float(np.sqrt(np.mean(mom2))),
        mom_max=float(np.sqrt(np.max(mom2))),
        valid_cells=valid,
        valid_fraction=valid / alpha.size,
    )


def build_grid(cfg: TeoWormholeConfig) -> TeoGrid:
    """Assemble the full CCZ4 + effective-source state for a Teo wormhole."""
    alpha, beta, gamma = build_adm(cfg)
    dx = cfg.dx
    K, A_phys = extrinsic_curvature(alpha, beta, gamma, dx)

    det_gamma = np.linalg.det(gamma)
    chi = np.power(np.maximum(det_gamma, TINY), -1.0 / 3.0)
    h = chi[..., None, None] * gamma
    h_inv = np.linalg.inv(h)
    A = chi[..., None, None] * A_phys
    Gamma = conformal_connection(h, h_inv, dx)

    if cfg.source == "einstein":
        rho, j, S = einstein_source(alpha, beta, gamma, dx)
    else:
        rho = np.zeros_like(alpha)
        j = np.zeros(beta.shape, dtype=np.float64)
        S = np.zeros(gamma.shape, dtype=np.float64)

    names = list(GRTECLYN_STATE_VARS)
    idx = {name: i for i, name in enumerate(names)}
    data = np.zeros((cfg.nz, cfg.ny, cfg.nx, len(names)), dtype=np.float64)

    data[..., idx["chi"]] = chi
    data[..., idx["h11"]] = h[..., 0, 0]
    data[..., idx["h12"]] = h[..., 0, 1]
    data[..., idx["h13"]] = h[..., 0, 2]
    data[..., idx["h22"]] = h[..., 1, 1]
    data[..., idx["h23"]] = h[..., 1, 2]
    data[..., idx["h33"]] = h[..., 2, 2]
    data[..., idx["K"]] = K
    data[..., idx["A11"]] = A[..., 0, 0]
    data[..., idx["A12"]] = A[..., 0, 1]
    data[..., idx["A13"]] = A[..., 0, 2]
    data[..., idx["A22"]] = A[..., 1, 1]
    data[..., idx["A23"]] = A[..., 1, 2]
    data[..., idx["A33"]] = A[..., 2, 2]
    data[..., idx["Gamma1"]] = Gamma[..., 0]
    data[..., idx["Gamma2"]] = Gamma[..., 1]
    data[..., idx["Gamma3"]] = Gamma[..., 2]
    data[..., idx["lapse"]] = alpha
    data[..., idx["shift1"]] = beta[..., 0]
    data[..., idx["shift2"]] = beta[..., 1]
    data[..., idx["shift3"]] = beta[..., 2]
    data[..., idx["teo_rho"]] = rho
    data[..., idx["teo_j1"]] = j[..., 0]
    data[..., idx["teo_j2"]] = j[..., 1]
    data[..., idx["teo_j3"]] = j[..., 2]
    data[..., idx["teo_S11"]] = S[..., 0, 0]
    data[..., idx["teo_S12"]] = S[..., 0, 1]
    data[..., idx["teo_S13"]] = S[..., 0, 2]
    data[..., idx["teo_S22"]] = S[..., 1, 1]
    data[..., idx["teo_S23"]] = S[..., 1, 2]
    data[..., idx["teo_S33"]] = S[..., 2, 2]

    if not np.all(np.isfinite(data)):
        raise ValueError("Generated Teo grid contains non-finite values")
    if np.min(data[..., idx["chi"]]) <= 0.0:
        raise ValueError("Generated Teo grid contains non-positive chi")
    if np.min(data[..., idx["lapse"]]) <= 0.0:
        raise ValueError("Generated Teo grid contains non-positive lapse")

    metrics = {
        "min_chi": float(np.min(chi)),
        "max_chi": float(np.max(chi)),
        "min_lapse": float(np.min(alpha)),
        "max_lapse": float(np.max(alpha)),
        "max_abs_shift": float(np.max(np.abs(beta))),
        "max_abs_K": float(np.max(np.abs(K))),
        "min_teo_rho": float(np.min(rho)),
        "max_teo_rho": float(np.max(rho)),
    }
    return TeoGrid(
        data=data,
        comp_names=names,
        dx_xyz=np.array(dx, dtype=np.float64),
        origin=np.array(cfg.origin, dtype=np.float64),
        metrics=metrics,
    )
