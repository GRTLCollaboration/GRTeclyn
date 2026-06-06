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

    def _coerce_scalar_samples(values) -> np.ndarray:
        arr = np.asarray(values)
        if arr.dtype != object and arr.ndim == 1:
            return arr.astype(float, copy=False)
        out = np.empty(len(values), dtype=float)
        for jj, v in enumerate(values):
            vv = np.asarray(v, dtype=float).reshape(-1)
            out[jj] = vv[0] if vv.size else np.nan
        return out

    def sample_at(names: list[str], points: np.ndarray) -> dict[str, np.ndarray] | None:
        """Robustly sample fields at points using the yt point API.

        ``ds.sample`` does not exist; ``find_field_values_at_points`` is fast but
        can choke on AMR points, so fall back to per-point ``ds.point``.
        """
        try:
            vals = ds.find_field_values_at_points(
                [("boxlib", n) for n in names], points
            )
            return {n: _coerce_scalar_samples(vals[i]) for i, n in enumerate(names)}
        except Exception:  # noqa: BLE001
            out = {n: np.full(len(points), np.nan, dtype=float) for n in names}
            for ii in range(len(points)):
                pt = ds.point([points[ii, 0], points[ii, 1], points[ii, 2]])
                for n in names:
                    v = np.asarray(pt[("boxlib", n)], dtype=float).reshape(-1)
                    out[n][ii] = float(v[0]) if v.size else np.nan
            return out

    fields = sample_at(["phi", "Pi", "lapse", "shift1", "shift2", "shift3"], pts)
    if fields is None:
        return None
    phi_f = fields["phi"].reshape(x.shape)
    pi_f = fields["Pi"].reshape(x.shape)
    alpha = np.clip(fields["lapse"].reshape(x.shape), 1.0e-10, None)
    shift1 = fields["shift1"].reshape(x.shape)
    shift2 = fields["shift2"].reshape(x.shape)
    shift3 = fields["shift3"].reshape(x.shape)
    if not np.any(np.isfinite(phi_f)):
        return None

    dx = float(ds.index.grids[0].dds[0].to_value())
    dr = max(dx, 1.0e-6)
    phi_rp = sample_at(["phi"], pts + np.array([dr, 0.0, 0.0]))["phi"].reshape(x.shape)
    phi_rm = sample_at(["phi"], pts - np.array([dr, 0.0, 0.0]))["phi"].reshape(x.shape)
    dphi_dr = (phi_rp - phi_rm) / (2.0 * dr)

    nx = (x - center[0]) / np.maximum(r_shell, 1.0e-12)
    ny = (y - center[1]) / np.maximum(r_shell, 1.0e-12)
    nz = (z - center[2]) / np.maximum(r_shell, 1.0e-12)

    beta_dot_n = shift1 * nx + shift2 * ny + shift3 * nz
    flux = (-pi_f * dphi_dr / alpha + beta_dot_n * pi_f * pi_f / (alpha * alpha))

    dtheta = theta[1] - theta[0] if n_theta > 1 else math.pi
    dphi = 2.0 * math.pi / n_phi
    weights = (r_shell ** 2) * np.sin(TH) * dtheta * dphi
    net = float(np.nansum(flux * weights))

    psi4_amp = None
    try:
        weyl = sample_at(["Weyl4_Re", "Weyl4_Im"], pts)
        w_re = weyl["Weyl4_Re"]
        w_im = weyl["Weyl4_Im"]
        finite = np.isfinite(w_re) & np.isfinite(w_im)
        if np.any(finite):
            psi4_amp = float(
                np.sqrt(np.mean(w_re[finite] ** 2 + w_im[finite] ** 2))
            )
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
