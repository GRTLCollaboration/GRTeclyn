"""Boundary reflection auditing via scalar-field energy flux."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class BoundaryFluxMetrics:
    final_time: float | None
    net_outward_flux: float | None
    late_inward_flux: float | None
    reflection_contaminated: bool
    psi4_boundary_amp: float | None = None


def extract_scalar_boundary_flux(
    plotfile: str | Path,
    *,
    shell_fraction: float = 0.92,
) -> tuple[float, float] | None:
    """Integrate signed radial scalar energy flux on a near-boundary sphere.

    Returns ``(net_outward_flux, psi4_boundary_amp)`` where outward flux is
    positive.  Uses ``T^t_r ~ -Pi * partial_r phi / alpha`` for a minimally
    coupled scalar.
    """
    import yt  # local import

    ds = yt.load(str(plotfile))
    center = np.array([float(c.to_value()) for c in ds.domain_center])
    dleft = np.array([float(c.to_value()) for c in ds.domain_left_edge])
    dright = np.array([float(c.to_value()) for c in ds.domain_right_edge])
    half = 0.5 * (dright - dleft)
    r_shell = shell_fraction * float(np.min(half))

    n_theta, n_phi = 24, 48
    theta = np.linspace(0.0, math.pi, n_theta)
    phi = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    x = center[0] + r_shell * np.sin(TH) * np.cos(PH)
    y = center[1] + r_shell * np.sin(TH) * np.sin(PH)
    z = center[2] + r_shell * np.cos(TH)

    pts = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)

    def sample(name: str) -> np.ndarray:
        try:
            return np.asarray(ds.sample(name, pts)["boxlib", name], dtype=float)
        except Exception:  # noqa: BLE001
            return np.asarray(ds.sample(name, pts)[name], dtype=float)

    try:
        phi_f = sample("phi").reshape(x.shape)
        pi_f = sample("Pi").reshape(x.shape)
        alpha = np.clip(sample("lapse").reshape(x.shape), 1.0e-10, None)
        shift1 = sample("shift1").reshape(x.shape)
        shift2 = sample("shift2").reshape(x.shape)
        shift3 = sample("shift3").reshape(x.shape)
    except Exception:
        return None

    dx = float(ds.index.grids[0].dds[0].to_value())
    dr = max(dx, 1.0e-6)
    phi_rp = np.asarray(
        ds.sample("phi", pts + np.array([dr, 0.0, 0.0]))["boxlib", "phi"],
        dtype=float,
    ).reshape(x.shape)
    phi_rm = np.asarray(
        ds.sample("phi", pts - np.array([dr, 0.0, 0.0]))["boxlib", "phi"],
        dtype=float,
    ).reshape(x.shape)
    dphi_dr = (phi_rp - phi_rm) / (2.0 * dr)

    nx = (x - center[0]) / np.maximum(r_shell, 1.0e-12)
    ny = (y - center[1]) / np.maximum(r_shell, 1.0e-12)
    nz = (z - center[2]) / np.maximum(r_shell, 1.0e-12)

    beta_dot_n = shift1 * nx + shift2 * ny + shift3 * nz
    flux = (-pi_f * dphi_dr / alpha + beta_dot_n * pi_f * pi_f / (alpha * alpha))

    dtheta = theta[1] - theta[0] if n_theta > 1 else math.pi
    dphi = 2.0 * math.pi / n_phi
    weights = (r_shell ** 2) * np.sin(TH) * dtheta * dphi
    net = float(np.sum(flux * weights))

    psi4_amp = None
    try:
        w_re = sample("Weyl4_Re")
        w_im = sample("Weyl4_Im")
        psi4_amp = float(np.sqrt(np.mean(w_re ** 2 + w_im ** 2)))
    except Exception:
        psi4_amp = None

    return net, psi4_amp


def read_boundary_flux_metrics(path: Path) -> BoundaryFluxMetrics | None:
    rows: list[list[float]] = []
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                rows.append([float(v) for v in stripped.split()])
            except ValueError:
                continue
    if not rows:
        return None
    rows.sort(key=lambda r: r[0])
    final = rows[-1]
    late = rows[-max(1, len(rows) // 5):]
    late_inward = min((r[1] for r in late if len(r) >= 2), default=0.0)
    net = final[1] if len(final) >= 2 else None
    psi4 = final[2] if len(final) >= 3 else None
    contaminated = late_inward < -1.0e-8
    return BoundaryFluxMetrics(
        final_time=final[0],
        net_outward_flux=net,
        late_inward_flux=late_inward,
        reflection_contaminated=contaminated,
        psi4_boundary_amp=psi4,
    )
