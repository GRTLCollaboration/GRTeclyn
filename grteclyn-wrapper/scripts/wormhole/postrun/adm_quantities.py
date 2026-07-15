#!/usr/bin/env python3
"""ADM mass and angular momentum from a GRTeclyn .gridinit (t=0).

Surface integrals on a coordinate sphere of radius R (default 24):

    M_ADM = (1/16π) ∮ (∂_j γ_ij − ∂_i γ_jj) n^i dA

    J_z   = (1/8π)  ∮ ε_{zjk} x^j K^k_l n^l dA

Uses the CCZ4/BSSN conformal variables stored in the gridinit:
    γ_ij = χ^{-1} h̃_ij,   K_ij = χ^{-1} (Ã_ij + (1/3) K h̃_ij).

Usage:
    python adm_quantities.py path/to/initial_data.gridinit [--radius 24]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Allow running without installing the wrapper package.
_WRAPPER_SRC = Path(__file__).resolve().parents[3] / "src"
if str(_WRAPPER_SRC) not in sys.path:
    sys.path.insert(0, str(_WRAPPER_SRC))

from grteclyn_wrapper.grtresna.io import read_gridinit  # noqa: E402


def _comp_index(names: list[str], name: str) -> int:
    try:
        return names.index(name)
    except ValueError as exc:
        raise KeyError(f"gridinit missing component {name!r}; have {names}") from exc


def _finite_diff_grad(arr: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Central differences; arr shaped (nz, ny, nx). Returns (dz, dy, dx)."""
    gz = np.gradient(arr, dx, axis=0, edge_order=1)
    gy = np.gradient(arr, dx, axis=1, edge_order=1)
    gx = np.gradient(arr, dx, axis=2, edge_order=1)
    return gz, gy, gx


def adm_from_gridinit(
    gridinit_path: str | Path,
    radius: float = 24.0,
    n_theta: int = 64,
    n_phi: int = 128,
    center: tuple[float, float, float] | None = None,
) -> dict[str, float]:
    g = read_gridinit(gridinit_path)
    names = list(g.comp_names)
    data = g.data  # (nz, ny, nx, ncomp)
    dx = float(np.asarray(g.dx_xyz).flat[0])
    origin = np.asarray(g.origin, dtype=float)
    nz, ny, nx, _ = data.shape

    # Cell-center coordinates.
    xs = origin[0] + (np.arange(nx) + 0.5) * dx
    ys = origin[1] + (np.arange(ny) + 0.5) * dx
    zs = origin[2] + (np.arange(nz) + 0.5) * dx
    if center is None:
        center = (float(xs.mean()), float(ys.mean()), float(zs.mean()))
    cx, cy, cz = center

    chi = np.clip(data[:, :, :, _comp_index(names, "chi")], 1e-12, None)
    h = {
        (1, 1): data[:, :, :, _comp_index(names, "h11")],
        (1, 2): data[:, :, :, _comp_index(names, "h12")],
        (1, 3): data[:, :, :, _comp_index(names, "h13")],
        (2, 2): data[:, :, :, _comp_index(names, "h22")],
        (2, 3): data[:, :, :, _comp_index(names, "h23")],
        (3, 3): data[:, :, :, _comp_index(names, "h33")],
    }
    A = {
        (1, 1): data[:, :, :, _comp_index(names, "A11")],
        (1, 2): data[:, :, :, _comp_index(names, "A12")],
        (1, 3): data[:, :, :, _comp_index(names, "A13")],
        (2, 2): data[:, :, :, _comp_index(names, "A22")],
        (2, 3): data[:, :, :, _comp_index(names, "A23")],
        (3, 3): data[:, :, :, _comp_index(names, "A33")],
    }
    Ktr = data[:, :, :, _comp_index(names, "K")]

    # Physical metric γ_ij = χ^{-1} h̃_ij.
    inv_chi = 1.0 / chi
    gamma = {ij: inv_chi * h[ij] for ij in h}

    # Gradients of γ_ij (Cartesian): g_ij,k
    dgamma: dict[tuple[int, int, int], np.ndarray] = {}
    for (i, j), arr in gamma.items():
        gz, gy, gx = _finite_diff_grad(arr, dx)
        dgamma[(i, j, 1)] = gx
        dgamma[(i, j, 2)] = gy
        dgamma[(i, j, 3)] = gz

    def gcomp(i: int, j: int) -> np.ndarray:
        a, b = (i, j) if i <= j else (j, i)
        return gamma[(a, b)]

    def dg(i: int, j: int, k: int) -> np.ndarray:
        a, b = (i, j) if i <= j else (j, i)
        return dgamma[(a, b, k)]

    # Extrinsic curvature K_ij = χ^{-1} (Ã_ij + K/3 h̃_ij).
    Kij = {}
    for ij in A:
        Kij[ij] = inv_chi * (A[ij] + (Ktr / 3.0) * h[ij])

    def Kcomp(i: int, j: int) -> np.ndarray:
        a, b = (i, j) if i <= j else (j, i)
        return Kij[(a, b)]

    # Sphere samples (trilinear via round-to-nearest for robustness / speed).
    theta = np.linspace(0.0, np.pi, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    dtheta = np.pi / (n_theta - 1)
    dphi = 2.0 * np.pi / n_phi
    w = np.sin(TH) * dtheta * dphi  # solid-angle weight

    nx_hat = np.sin(TH) * np.cos(PH)
    ny_hat = np.sin(TH) * np.sin(PH)
    nz_hat = np.cos(TH)
    px = cx + radius * nx_hat
    py = cy + radius * ny_hat
    pz = cz + radius * nz_hat

    def sample(arr: np.ndarray) -> np.ndarray:
        ix = np.clip(np.rint((px - origin[0]) / dx - 0.5).astype(int), 0, nx - 1)
        iy = np.clip(np.rint((py - origin[1]) / dx - 0.5).astype(int), 0, ny - 1)
        iz = np.clip(np.rint((pz - origin[2]) / dx - 0.5).astype(int), 0, nz - 1)
        return arr[iz, iy, ix]

    # ADM mass integrand: (∂_j γ_ij − ∂_i γ_jj) n^i
    # i summed with n^i; j summed.
    integrand_M = np.zeros_like(TH)
    n_hat = {1: nx_hat, 2: ny_hat, 3: nz_hat}
    for i in (1, 2, 3):
        term = np.zeros_like(TH)
        for j in (1, 2, 3):
            term = term + sample(dg(i, j, j)) - sample(dg(j, j, i))
        integrand_M = integrand_M + term * n_hat[i]

    M_adm = float(np.sum(integrand_M * w) * radius**2 / (16.0 * np.pi))

    # J_z = (1/8π) ∮ ε_zjk x^j K^k_l n^l dA
    # ε_zxy=1, ε_zyx=-1 → x K^y_l n^l − y K^x_l n^l
    # rebuild coordinate arrays matching data layout (z, y, x):
    Zg, Yg, Xg = np.meshgrid(zs, ys, xs, indexing="ij")

    x_s = sample(Xg) - cx
    y_s = sample(Yg) - cy

    Ky_n = (
        sample(Kcomp(2, 1)) * nx_hat
        + sample(Kcomp(2, 2)) * ny_hat
        + sample(Kcomp(2, 3)) * nz_hat
    )
    Kx_n = (
        sample(Kcomp(1, 1)) * nx_hat
        + sample(Kcomp(1, 2)) * ny_hat
        + sample(Kcomp(1, 3)) * nz_hat
    )
    integrand_Jz = x_s * Ky_n - y_s * Kx_n
    J_adm = float(np.sum(integrand_Jz * w) * radius**2 / (8.0 * np.pi))

    # Sanity: conformal-factor form for nearly conformally flat data.
    # ψ = χ^{-1/4}; M_ψ = −(1/2π) ∮ ∂_r ψ r^2 dΩ
    psi = chi ** (-0.25)
    dpsi_z, dpsi_y, dpsi_x = _finite_diff_grad(psi, dx)
    dpsi_r = (
        sample(dpsi_x) * nx_hat + sample(dpsi_y) * ny_hat + sample(dpsi_z) * nz_hat
    )
    M_psi = float(-np.sum(dpsi_r * w) * radius**2 / (2.0 * np.pi))

    return {
        "M_ADM": M_adm,
        "M_psi": M_psi,
        "J_ADM_z": J_adm,
        "radius": float(radius),
        "center": (cx, cy, cz),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("gridinit", type=Path)
    ap.add_argument("--radius", type=float, default=24.0)
    ap.add_argument("--n-theta", type=int, default=64)
    ap.add_argument("--n-phi", type=int, default=128)
    ap.add_argument(
        "--center",
        type=float,
        nargs=3,
        default=None,
        metavar=("X", "Y", "Z"),
    )
    ap.add_argument(
        "--radii",
        type=float,
        nargs="+",
        default=None,
        help="If set, evaluate M/J at several radii (convergence check).",
    )
    args = ap.parse_args()

    radii = args.radii or [args.radius]
    print(f"gridinit: {args.gridinit}")
    for R in radii:
        out = adm_from_gridinit(
            args.gridinit,
            radius=R,
            n_theta=args.n_theta,
            n_phi=args.n_phi,
            center=tuple(args.center) if args.center else None,
        )
        print(
            f"  R={R:g}:  M_ADM={out['M_ADM']:.6g}  M_psi={out['M_psi']:.6g}  "
            f"J_ADM_z={out['J_ADM_z']:.6g}  center={out['center']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
