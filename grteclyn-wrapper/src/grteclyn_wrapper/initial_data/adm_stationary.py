"""Reusable stationary ADM / CCZ4 helpers for geometry-first pipelines.

Extracted from the Teo wormhole generator so arbitrary stationary metrics
(alpha, beta^i, gamma_ij) can share the same extrinsic curvature, Einstein
source, constraint residual, and CCZ4 packing logic.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ..grtresna.io import GRTECLYN_STATE_VARS

TINY = 1.0e-12
INV_8PI = 1.0 / (8.0 * np.pi)


@dataclass
class ConstraintResiduals:
    """ADM Hamiltonian / momentum constraint residuals."""

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


@dataclass
class StationaryGrid:
    """Uniform CCZ4 + effective-source state on a Cartesian grid."""

    data: NDArray[np.float64]  # (nz, ny, nx, n_comp)
    comp_names: list[str]
    dx_xyz: NDArray[np.float64]
    origin: NDArray[np.float64]
    metrics: dict[str, float]
    alpha: NDArray[np.float64]
    beta: NDArray[np.float64]
    gamma: NDArray[np.float64]
    rho: NDArray[np.float64]
    j: NDArray[np.float64]
    S: NDArray[np.float64]
    K: NDArray[np.float64]
    A: NDArray[np.float64]


def spatial_gradient(field: NDArray, dx: tuple[float, float, float]) -> list[NDArray]:
    """Return [d/dx, d/dy, d/dz] for an array indexed (z, y, x, ...)."""
    edge_order = 2 if min(field.shape[:3]) > 2 else 1
    d_dz, d_dy, d_dx = np.gradient(field, dx[2], dx[1], dx[0], edge_order=edge_order)
    return [d_dx, d_dy, d_dz]


def christoffel(metric: NDArray, metric_inv: NDArray, dx, dim: int) -> NDArray:
    """Levi-Civita connection for a (dim x dim) metric on a spatial grid.

    For 4-metrics, time derivatives vanish (stationarity), so only spatial
    derivatives contribute.
    """
    dmetric: list = []
    for a in range(dim):
        row = []
        for b in range(dim):
            row.append(spatial_gradient(metric[..., a, b], dx))
        dmetric.append(row)

    def d(idx_lo: int, a: int, b: int) -> NDArray | float:
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


def extrinsic_curvature(
    alpha: NDArray, beta: NDArray, gamma: NDArray, dx
) -> tuple[NDArray, NDArray]:
    """Stationary K_ij = (D_i beta_j + D_j beta_i) / (2 alpha); return (K, A_phys)."""
    gamma_inv = np.linalg.inv(gamma)
    beta_low = np.einsum("...ij,...j->...i", gamma, beta)
    dbeta = [spatial_gradient(beta_low[..., i], dx) for i in range(3)]
    chris = christoffel(gamma, gamma_inv, dx, dim=3)

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
    chris = christoffel(h, h_inv, dx, dim=3)
    return np.einsum("...jk,...ijk->...i", h_inv, chris)


def einstein_source(
    alpha: NDArray, beta: NDArray, gamma: NDArray, dx
) -> tuple[NDArray, NDArray, NDArray]:
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

    chris = christoffel(g4, g4_inv, dx, dim=4)

    def dchris(idx_lo: int, a: int, b: int, c: int) -> NDArray | float:
        if idx_lo == 0:
            return 0.0
        return spatial_gradient(chris[..., a, b, c], dx)[idx_lo - 1]

    trace = np.zeros(g4.shape[:-2] + (4,), dtype=np.float64)
    for a in range(4):
        for c in range(4):
            trace[..., a] += chris[..., c, a, c]

    def dtrace(idx_lo: int, a: int) -> NDArray | float:
        if idx_lo == 0:
            return 0.0
        return spatial_gradient(trace[..., a], dx)[idx_lo - 1]

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


def ricci_tensor_3d(gamma: NDArray, gamma_inv: NDArray, dx) -> NDArray:
    """Spatial Ricci tensor R_ij from the 3-metric."""
    chris = christoffel(gamma, gamma_inv, dx, dim=3)
    ricci = np.zeros_like(gamma)
    for i in range(3):
        for j in range(3):
            t1 = 0.0
            t2 = 0.0
            for k in range(3):
                t1 = t1 + spatial_gradient(chris[..., k, i, j], dx)[k]
                t2 = t2 + spatial_gradient(chris[..., k, i, k], dx)[j]
            t3 = 0.0
            t4 = 0.0
            for k in range(3):
                for m in range(3):
                    t3 = t3 + chris[..., k, i, j] * chris[..., m, k, m]
                    t4 = t4 + chris[..., m, i, k] * chris[..., k, j, m]
            ricci[..., i, j] = t1 - t2 + t3 - t4
    return ricci


def constraint_residuals_from_adm(
    alpha: NDArray,
    beta: NDArray,
    gamma: NDArray,
    dx: tuple[float, float, float],
    *,
    center: tuple[float, float, float] | None = None,
    trim_cells: int = 3,
    puncture_cells: float = 0.0,
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> ConstraintResiduals:
    """Evaluate ADM constraint residuals for stationary ADM data."""
    gamma_inv = np.linalg.inv(gamma)
    K_trace, A_phys = extrinsic_curvature(alpha, beta, gamma, dx)
    Kij = A_phys + (K_trace[..., None, None] / 3.0) * gamma
    rho, j_low, _S = einstein_source(alpha, beta, gamma, dx)

    ricci = ricci_tensor_3d(gamma, gamma_inv, dx)
    R_scalar = np.einsum("...ij,...ij->...", gamma_inv, ricci)
    K_up = np.einsum("...ia,...jb,...ab->...ij", gamma_inv, gamma_inv, Kij)
    KijKij = np.einsum("...ij,...ij->...", Kij, K_up)
    ham = R_scalar + K_trace * K_trace - KijKij - 16.0 * np.pi * rho

    chris = christoffel(gamma, gamma_inv, dx, dim=3)
    K_mixed = np.einsum("...jk,...ki->...ji", gamma_inv, Kij)
    mom = np.zeros(beta.shape, dtype=np.float64)
    for i in range(3):
        div = 0.0
        for jj in range(3):
            div = div + spatial_gradient(K_mixed[..., jj, i], dx)[jj]
            for ll in range(3):
                div = div + chris[..., jj, jj, ll] * K_mixed[..., ll, i]
                div = div - chris[..., ll, jj, i] * K_mixed[..., jj, ll]
        di_K = spatial_gradient(K_trace, dx)[i]
        mom[..., i] = div - di_K - 8.0 * np.pi * j_low[..., i]

    nz, ny, nx = alpha.shape
    mask = np.zeros(alpha.shape, dtype=bool)
    t = trim_cells
    mask[t : nz - t, t : ny - t, t : nx - t] = True

    if puncture_cells > 0.0 and center is not None:
        x = origin[0] + (np.arange(nx) + 0.5) * dx[0]
        y = origin[1] + (np.arange(ny) + 0.5) * dx[1]
        z = origin[2] + (np.arange(nz) + 0.5) * dx[2]
        zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
        radius = np.sqrt(
            (xx - center[0]) ** 2 + (yy - center[1]) ** 2 + (zz - center[2]) ** 2
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


def pack_ccz4_grid(
    alpha: NDArray,
    beta: NDArray,
    gamma: NDArray,
    dx: tuple[float, float, float],
    origin: tuple[float, float, float] | NDArray,
    *,
    include_einstein_source: bool = True,
    kij_extra: NDArray | None = None,
) -> StationaryGrid:
    """Assemble full CCZ4 + optional effective-source state from ADM fields.

    ``kij_extra``, when provided, is added to the stationary extrinsic
    curvature before conformal packing.  Frozen 4-metric probes still depend
    only on ``(alpha, beta, gamma)``; the extra ``K_ij`` is for Stage-2
    handoff diversity.
    """
    K_trace_stat, A_phys = extrinsic_curvature(alpha, beta, gamma, dx)
    gamma_inv = np.linalg.inv(gamma)
    Kij = A_phys + (K_trace_stat[..., None, None] / 3.0) * gamma
    if kij_extra is not None:
        Kij = Kij + np.asarray(kij_extra, dtype=np.float64)
        K = np.einsum("...ij,...ij->...", gamma_inv, Kij)
        A_phys = Kij - (K[..., None, None] / 3.0) * gamma
    else:
        K = K_trace_stat

    det_gamma = np.linalg.det(gamma)
    chi = np.power(np.maximum(det_gamma, TINY), -1.0 / 3.0)
    h = chi[..., None, None] * gamma
    h_inv = np.linalg.inv(h)
    A = chi[..., None, None] * A_phys
    Gamma = conformal_connection(h, h_inv, dx)

    if include_einstein_source:
        rho, j, S = einstein_source(alpha, beta, gamma, dx)
    else:
        rho = np.zeros_like(alpha)
        j = np.zeros(beta.shape, dtype=np.float64)
        S = np.zeros(gamma.shape, dtype=np.float64)

    names = list(GRTECLYN_STATE_VARS)
    idx = {name: i for i, name in enumerate(names)}
    nz, ny, nx = alpha.shape
    data = np.zeros((nz, ny, nx, len(names)), dtype=np.float64)

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
        raise ValueError("Generated stationary grid contains non-finite values")
    if np.min(data[..., idx["chi"]]) <= 0.0:
        raise ValueError("Generated stationary grid contains non-positive chi")
    if np.min(data[..., idx["lapse"]]) <= 0.0:
        raise ValueError("Generated stationary grid contains non-positive lapse")

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
    return StationaryGrid(
        data=data,
        comp_names=names,
        dx_xyz=np.asarray(dx, dtype=np.float64),
        origin=np.asarray(origin, dtype=np.float64),
        metrics=metrics,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        rho=rho,
        j=j,
        S=S,
        K=K,
        A=A,
    )
