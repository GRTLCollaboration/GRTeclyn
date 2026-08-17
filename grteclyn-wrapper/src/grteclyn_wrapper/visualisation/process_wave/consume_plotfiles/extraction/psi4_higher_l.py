"""Higher-multipole Psi4 extraction (l >= 3).

Why this is a separate module rather than an extension of ``psi4.py``
---------------------------------------------------------------------
``psi4.py`` writes ``psi4_mode_l2m0.dat``, ``psi4_mode_l2_all.dat`` and
``psi4_directional.dat``, whose column layouts are consumed by the article
figure code.  Those contracts must not widen.  This module therefore owns its
own stream (``psi4_mode_higher_l.dat``) and is **off by default**, enabled per
campaign with ``GRTECLYN_PSI4_HIGHER_L=1``.

Physics motivation (Debug.md item C)
------------------------------------
The runaway is an x-axis phenomenon while the mode decomposition is taken about
z, so the momentum-flux beaming along the runaway axis shows up in the
interference between odd and even l.  With l = 2 alone that interference is
invisible.  Note that psi4 has spin weight -2, so **l = 1 does not exist** --
there is no "dipole mode" of psi4 to extract, whatever the naive
negative-mass-binary expectation suggests.

The projection convention matches ``psi4.py`` exactly:

    A_lm = r * sum( psi4 * conj(sYlm) * W_renorm )

so an l = 2 column produced here is directly comparable to the legacy stream --
which is what ``selftest_against_l2`` checks.
"""

from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np

# Deliberately the *general* harmonic, not ``spin_weighted_sph_harm_s2``.
# The latter routes l = 2 to the legacy closed forms, whose m = +-1 modes carry
# the opposite overall sign.  That is harmless within a single l (it cancels in
# |amplitude|), but this stream exists precisely to measure odd-l/even-l
# interference, where a relative sign between l = 2 and l = 3 is the signal.
from ..sphere import spin_weighted_sph_harm

# l values this module knows how to project onto.  sphere.py implements the
# s=-2 harmonics for l=2,3,4 only; asking for anything else is a hard error
# rather than a silently skipped multipole.
SUPPORTED_ELLS = (2, 3, 4)
DEFAULT_ELLS = (3, 4)


def _coerce_scalar_samples(values) -> np.ndarray:
    arr = np.asarray(values)
    if arr.dtype != object and arr.ndim == 1:
        return arr.astype(float, copy=False)

    out = np.empty(len(values), dtype=float)
    for jj, v in enumerate(values):
        vv = np.asarray(v, dtype=float).reshape(-1)
        out[jj] = vv[0] if vv.size else np.nan
    return out


def parse_ells(raw: Sequence[int] | str | None) -> tuple[int, ...]:
    """Normalise an ell specification into a validated, ordered tuple."""
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        return DEFAULT_ELLS
    if isinstance(raw, str):
        items = [tok for tok in raw.replace(",", " ").split() if tok]
    else:
        items = list(raw)
    ells: List[int] = []
    for tok in items:
        val = int(tok)
        if val not in SUPPORTED_ELLS:
            raise ValueError(
                f"unsupported l={val} for s=-2 harmonics; supported: {SUPPORTED_ELLS}"
            )
        if val not in ells:
            ells.append(val)
    return tuple(sorted(ells))


def extract_higher_l_modes(
    ds,
    radii: Sequence[float],
    n_points: int,
    center: Sequence[float],
    ells: Sequence[int] = DEFAULT_ELLS,
) -> Dict[int, List[Dict[int, complex]]]:
    """Project Psi4 onto s=-2 harmonics for each requested l.

    Returns ``{l: [ {m: amplitude} for each radius ]}``.

    The sphere is sampled once and reused for every l, so requesting l = 3 and
    l = 4 together costs one field query, not two.
    """
    ells = tuple(ells)
    if not radii or not ells:
        return {l: [] for l in ells}

    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    center = np.asarray(center, dtype=float)

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
    shape = X1.shape

    # One harmonic table per l, computed once for the whole plotfile.
    harmonics: Dict[int, Dict[int, np.ndarray]] = {
        l: {m: spin_weighted_sph_harm(-2, l, m, THETA, PHI) for m in range(-l, l + 1)}
        for l in ells
    }

    starts: List[int] = []
    ends: List[int] = []
    sx_all: List[np.ndarray] = []
    sy_all: List[np.ndarray] = []
    sz_all: List[np.ndarray] = []
    cursor = 0
    for r in radii:
        starts.append(cursor)
        cursor += shape[0] * shape[1]
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
        weyl[idxs] = _coerce_scalar_samples(vals[0]) + 1j * _coerce_scalar_samples(vals[1])
    except Exception:
        for ii, i in enumerate(idxs):
            pt = ds.point([pts[ii, 0], pts[ii, 1], pts[ii, 2]])
            re = np.asarray(pt[("boxlib", "Weyl4_Re")])
            im = np.asarray(pt[("boxlib", "Weyl4_Im")])
            if re.size and im.size:
                weyl[i] = float(re.flat[0]) + 1j * float(im.flat[0])

    out: Dict[int, List[Dict[int, complex]]] = {l: [] for l in ells}
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

        for l in ells:
            modes: Dict[int, complex] = {}
            for m, ylm in harmonics[l].items():
                modes[m] = complex(np.sum(psi4 * np.conj(ylm) * W_renorm) * float(r))
            out[l].append(modes)

    return out


def higher_l_header(ells: Sequence[int], radii: Sequence[float]) -> str:
    """Column header for ``psi4_mode_higher_l.dat``.

    Layout: time, then for each radius, for each l, for each m in -l..+l, a
    Re/Im pair.  Column names carry the radius so the file is self-describing.
    """
    cols = ["time"]
    for r in radii:
        for l in ells:
            for m in range(-l, l + 1):
                cols.append(f"R{float(r):g}_l{l}_m{m}_re")
                cols.append(f"R{float(r):g}_l{l}_m{m}_im")
    return "# " + "  ".join(cols)


def higher_l_line(
    t: float,
    modes_by_l: Dict[int, List[Dict[int, complex]]],
    ells: Sequence[int],
    n_radii: int,
) -> str:
    """Format one time row matching :func:`higher_l_header`."""
    parts = [f"{t:.16e}"]
    for ir in range(n_radii):
        for l in ells:
            per_radius = modes_by_l.get(l) or []
            modes = per_radius[ir] if ir < len(per_radius) else {}
            for m in range(-l, l + 1):
                c = modes.get(m, 0j)
                parts.append(f"{c.real:.16e}")
                parts.append(f"{c.imag:.16e}")
    return "  ".join(parts)
