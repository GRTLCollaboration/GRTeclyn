from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

from ..sphere import spin_weighted_sph_harm_s2_l2_m0, spin_weighted_sph_harm_s2_l2


def _coerce_scalar_samples(values) -> np.ndarray:
    arr = np.asarray(values)
    if arr.dtype != object and arr.ndim == 1:
        return arr.astype(float, copy=False)

    out = np.empty(len(values), dtype=float)
    for jj, v in enumerate(values):
        vv = np.asarray(v, dtype=float).reshape(-1)
        out[jj] = vv[0] if vv.size else np.nan
    return out


def _extract_mode_amps_l2m0(
    ds,
    radii: Sequence[float],
    n_points: int,
    center: Sequence[float],
) -> List[complex]:
    """
    Extract (l=2,m=0) mode amplitudes for multiple radii in ONE batched yt query.

    This is much faster than querying radii one-by-one because AMR indexing and
    disk reads are amortized.
    """
    if not radii:
        return []

    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    center = np.asarray(center, dtype=float)

    # Precompute unit-sphere grid (theta,phi) and weights (independent of radius).
    theta = np.linspace(0.0, np.pi, n_points)
    phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
    sinT = np.sin(THETA)
    X1 = sinT * np.cos(PHI)
    Y1 = sinT * np.sin(PHI)
    Z1 = np.cos(THETA)
    dtheta = np.pi / (n_points - 1)
    dphi = 2.0 * np.pi / n_points
    W = sinT * dtheta * dphi
    Wf = W.ravel()
    Y20 = spin_weighted_sph_harm_s2_l2_m0(THETA)
    shape = X1.shape

    # Build one big list of points for all radii.
    starts: List[int] = []
    ends: List[int] = []
    sx_all: List[np.ndarray] = []
    sy_all: List[np.ndarray] = []
    sz_all: List[np.ndarray] = []
    cursor = 0
    for r in radii:
        starts.append(cursor)
        n = shape[0] * shape[1]
        cursor += n
        ends.append(cursor)
        sx_all.append((float(r) * X1).ravel() + center[0])
        sy_all.append((float(r) * Y1).ravel() + center[1])
        sz_all.append((float(r) * Z1).ravel() + center[2])

    sx = np.concatenate(sx_all)
    sy = np.concatenate(sy_all)
    sz = np.concatenate(sz_all)

    in_domain = (
        (sx >= left[0])
        & (sx <= right[0])
        & (sy >= left[1])
        & (sy <= right[1])
        & (sz >= left[2])
        & (sz <= right[2])
    )
    idxs = np.where(in_domain)[0]
    if idxs.size < 4:
        raise RuntimeError("Too few sphere points inside domain (all radii).")

    weyl = np.full(sx.shape, np.nan + 1j * np.nan, dtype=np.complex128)
    pts = np.column_stack((sx[idxs], sy[idxs], sz[idxs]))

    try:
        vals = ds.find_field_values_at_points(
            [("boxlib", "Weyl4_Re"), ("boxlib", "Weyl4_Im")], pts
        )
        re_vals = _coerce_scalar_samples(vals[0])
        im_vals = _coerce_scalar_samples(vals[1])
        weyl[idxs] = re_vals + 1j * im_vals
    except Exception:
        # Fallback: per-point sampling (slow, but should still work)
        for ii, i in enumerate(idxs):
            pt = ds.point([pts[ii, 0], pts[ii, 1], pts[ii, 2]])
            re = np.asarray(pt[("boxlib", "Weyl4_Re")])
            im = np.asarray(pt[("boxlib", "Weyl4_Im")])
            if re.size and im.size:
                weyl[i] = float(re.flat[0]) + 1j * float(im.flat[0])

    amps: List[complex] = []
    for r, s, e in zip(radii, starts, ends):
        weyl_r = weyl[s:e]
        mask = np.isfinite(weyl_r.real) & np.isfinite(weyl_r.imag)
        if np.sum(mask) < 4:
            raise RuntimeError(f"Too few valid Weyl4 samples at R={r}.")

        omega_tot = np.sum(Wf[mask])
        if omega_tot <= 0.0:
            raise RuntimeError(f"Invalid integration weights at R={r} (omega_tot <= 0).")

        W_renorm = (Wf / omega_tot * 4.0 * np.pi).reshape(shape)
        psi4 = np.where(mask, weyl_r, 0.0).reshape(shape)
        amp = np.sum(psi4 * np.conj(Y20) * W_renorm) * float(r)
        amps.append(complex(amp))

    return amps


class DirectionalPsi4Metrics:
    """Container for directional Psi4 decomposition per radius."""

    def __init__(self) -> None:
        self.modes: Dict[int, List[complex]] = {m: [] for m in (-2, -1, 0, 1, 2)}
        self.p_total: List[float] = []
        self.p_z_beam: List[float] = []
        self.beam_ratio: List[float] = []

    def append_radius(self, modes: Dict[int, complex]) -> None:
        p_total = sum(abs(c) ** 2 for c in modes.values())
        p_z_beam = abs(modes[-2]) ** 2 + abs(modes[2]) ** 2
        self.p_total.append(float(p_total))
        self.p_z_beam.append(float(p_z_beam))
        self.beam_ratio.append(float(p_z_beam / p_total) if p_total > 0.0 else 0.0)
        for m, c in modes.items():
            self.modes[m].append(complex(c))


def _extract_mode_amps_l2_all(
    ds,
    radii: Sequence[float],
    n_points: int,
    center: Sequence[float],
) -> Tuple[List[complex], DirectionalPsi4Metrics]:
    """Extract all five (l=2,m) s=-2 modes and directional power metrics.

    Returns the legacy (l=2,m=0) amplitudes for backward compatibility plus a
    ``DirectionalPsi4Metrics`` object holding the full mode decomposition and the
    z-axis beaming proxy.
    """
    if not radii:
        return [], DirectionalPsi4Metrics()

    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    center = np.asarray(center, dtype=float)

    # Precompute unit-sphere grid and weights.
    theta = np.linspace(0.0, np.pi, n_points)
    phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
    sinT = np.sin(THETA)
    X1 = sinT * np.cos(PHI)
    Y1 = sinT * np.sin(PHI)
    Z1 = np.cos(THETA)
    dtheta = np.pi / (n_points - 1)
    dphi = 2.0 * np.pi / n_points
    W = sinT * dtheta * dphi
    Wf = W.ravel()
    Ylm = spin_weighted_sph_harm_s2_l2(THETA, PHI)
    shape = X1.shape

    starts: List[int] = []
    ends: List[int] = []
    sx_all: List[np.ndarray] = []
    sy_all: List[np.ndarray] = []
    sz_all: List[np.ndarray] = []
    cursor = 0
    for r in radii:
        starts.append(cursor)
        n = shape[0] * shape[1]
        cursor += n
        ends.append(cursor)
        sx_all.append((float(r) * X1).ravel() + center[0])
        sy_all.append((float(r) * Y1).ravel() + center[1])
        sz_all.append((float(r) * Z1).ravel() + center[2])

    sx = np.concatenate(sx_all)
    sy = np.concatenate(sy_all)
    sz = np.concatenate(sz_all)

    in_domain = (
        (sx >= left[0])
        & (sx <= right[0])
        & (sy >= left[1])
        & (sy <= right[1])
        & (sz >= left[2])
        & (sz <= right[2])
    )
    idxs = np.where(in_domain)[0]
    if idxs.size < 4:
        raise RuntimeError("Too few sphere points inside domain (all radii).")

    weyl = np.full(sx.shape, np.nan + 1j * np.nan, dtype=np.complex128)
    pts = np.column_stack((sx[idxs], sy[idxs], sz[idxs]))

    try:
        vals = ds.find_field_values_at_points(
            [("boxlib", "Weyl4_Re"), ("boxlib", "Weyl4_Im")], pts
        )
        re_vals = _coerce_scalar_samples(vals[0])
        im_vals = _coerce_scalar_samples(vals[1])
        weyl[idxs] = re_vals + 1j * im_vals
    except Exception:
        for ii, i in enumerate(idxs):
            pt = ds.point([pts[ii, 0], pts[ii, 1], pts[ii, 2]])
            re = np.asarray(pt[("boxlib", "Weyl4_Re")])
            im = np.asarray(pt[("boxlib", "Weyl4_Im")])
            if re.size and im.size:
                weyl[i] = float(re.flat[0]) + 1j * float(im.flat[0])

    l2m0_amps: List[complex] = []
    directional = DirectionalPsi4Metrics()
    for r, s, e in zip(radii, starts, ends):
        weyl_r = weyl[s:e]
        mask = np.isfinite(weyl_r.real) & np.isfinite(weyl_r.imag)
        if np.sum(mask) < 4:
            raise RuntimeError(f"Too few valid Weyl4 samples at R={r}.")

        omega_tot = np.sum(Wf[mask])
        if omega_tot <= 0.0:
            raise RuntimeError(f"Invalid integration weights at R={r} (omega_tot <= 0).")

        W_renorm = (Wf / omega_tot * 4.0 * np.pi).reshape(shape)
        psi4 = np.where(mask, weyl_r, 0.0).reshape(shape)

        modes: Dict[int, complex] = {}
        for m, ylm in Ylm.items():
            amp = np.sum(psi4 * np.conj(ylm) * W_renorm) * float(r)
            modes[m] = complex(amp)
        directional.append_radius(modes)
        l2m0_amps.append(modes[0])

    return l2m0_amps, directional
