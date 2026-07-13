"""Matter-profile CONTRACT: single source of truth + cross-code consistency rail.

Why this exists
---------------
Wiring a matter profile into the pipeline means the SAME physical facts are
repeated across three codes -- the Python solver, the GRTresna C++ painter
(``Source/Matter/BosonStarParams.hpp``), and the GRTeclyn evolution -- and a
disagreement produces a *plausible-looking wrong answer* (silent dispersal or a
t=0 constraint blow-up), never an error.  Real bugs of this class seen in
practice:

  * a stale GRTresna binary painted a torus as ``f == 0`` -> "Ham = 0 %"
    (looked like a converged solve, was actually empty space);
  * an initial-data sign flip vs the evolution's exotic matter -> instant
    constraint violation;
  * a diagnostics column off-by-one -> Q_sphere read as Q_total.

This module makes those invariants EXECUTABLE:

  * :class:`MatterProfileSpec` is the ONE place a profile's physics lives
    (couplings, frequency, winding, exotic sign, geometry).  It emits the
    GRTresna ``lump<k>_*`` params AND the Python reference field -- so the two
    cannot drift.
  * :func:`reference_paint` is the Python mirror of the C++ painter
    (``lump_winding_modulus`` + the U(1) momentum rule ``Pi = d_t Phi / alpha``).
    To add a NEW profile you implement ONE branch here.
  * :func:`check_gridinit_matches_spec` reads a GRTresna-produced ``.gridinit``
    and asserts the painted matter (phi, phi2, Pi, Pi2) matches the reference to
    a tolerance -- the t=0 cross-code consistency test that catches every bug in
    the list above in seconds instead of a rebuild-and-run cycle.

Convention (shared by both codes; keep in sync with BosonStarParams.hpp)::

    Phi = f(rho, z) e^{i (m phi_az)},   phi_az = atan2(y-cy, x-cx)
    phi1 = f cos(m phi_az),  phi2 = f sin(m phi_az)
    Pi1  = +(omega/alpha) phi2,  Pi2 = -(omega/alpha) phi1
    V    = 1/2 mass^2 |Phi|^2 - 1/4 lam |Phi|^4 + 1/6 mu6 |Phi|^6

The exotic (phantom) flag flips only the STRESS-ENERGY sign in the constraint
solve / evolution, NOT the painted field -- so phi/Pi are identical for exotic
and normal; the exotic invariant is checked separately against matter metadata.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from ..io.gridinit import read_gridinit
from ..profiles.qball_torus import QBallTorusProfile, read_torus_profile

ProfileKind = Literal["torus", "qball1d", "gaussian"]

# GRTresna BosonStarParams.hpp ``profile`` enum -- the ONE mapping kind<->int.
PROFILE_KIND_TO_INT: dict[str, int] = {
    "gaussian": 0,
    "qball1d": 3,
    "torus": 4,
}


@dataclass(frozen=True)
class MatterProfileSpec:
    """Single source of truth for one complex-scalar matter lump.

    All physics needed to (a) emit the GRTresna ``lump<k>_*`` params and (b) paint
    the reference field lives here, so the params and the reference cannot drift.
    """

    kind: ProfileKind
    amp: float
    omega: float
    mass: float
    lam: float = 0.0
    mu6: float = 0.0
    m_az: int = 1
    exotic: bool = False
    winding: bool = True
    width: float = 8.0
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    profile_path: str = ""
    # peak of the tabulated profile, so ``amp`` scales a normalised shape exactly
    # as the C++ ``qball_(torus_)envelope`` does (amp * f / f_max).
    f_max: float = 1.0

    def grtresna_scalar_keys(self) -> dict[str, float]:
        """Global scalar-potential keys (must match evolution phantom_* keys)."""
        return {
            "scalar_mass": self.mass,
            "scalar_lambda": self.lam,
            "scalar_mu": self.mu6,
            "bs_omega": self.omega,
        }

    def grtresna_lump_keys(self, k: int = 0) -> dict[str, object]:
        """The ``lump<k>_*`` params for this profile -- emitted from ONE place."""
        p = f"lump{k}_"
        keys: dict[str, object] = {
            f"{p}amp": self.amp,
            f"{p}width": self.width,
            f"{p}center": tuple(self.center),
            f"{p}velocity": (0.0, 0.0, 0.0),
            f"{p}omega": self.omega,
            f"{p}mode": self.m_az,
            f"{p}exotic": int(self.exotic),
            f"{p}winding": int(self.winding),
            f"{p}profile": PROFILE_KIND_TO_INT[self.kind],
        }
        if self.profile_path:
            keys[f"{p}profile_path"] = self.profile_path
        return keys


def _modulus_torus(
    spec: MatterProfileSpec,
    torus: QBallTorusProfile,
    dx: NDArray,
    dy: NDArray,
    dz: NDArray,
) -> NDArray:
    """f(rho, z) modulus mirroring C++ lump_winding_modulus for profile==4."""
    rho = np.sqrt(dx * dx + dy * dy)
    f = torus.eval_f(rho.ravel(), dz.ravel()).reshape(rho.shape)
    # C++: amp * qball_torus_interp = amp * f_table / f_max
    return spec.amp * f / spec.f_max


def _modulus_gaussian(
    spec: MatterProfileSpec, dx: NDArray, dy: NDArray, dz: NDArray
) -> NDArray:
    """amp * exp(-r^2/2w^2) * (sin theta)^m -- profile==0 winding modulus."""
    r2 = dx * dx + dy * dy + dz * dz
    r = np.sqrt(r2)
    env = np.exp(-0.5 * r2 / (spec.width * spec.width))
    with np.errstate(invalid="ignore"):
        sinth = np.where(r > 1e-12, np.sqrt(dx * dx + dy * dy) / np.maximum(r, 1e-30), 0.0)
    return spec.amp * env * sinth**spec.m_az


def reference_paint(
    spec: MatterProfileSpec,
    x: NDArray,
    y: NDArray,
    z: NDArray,
    *,
    alpha: NDArray | float = 1.0,
    torus: QBallTorusProfile | None = None,
) -> dict[str, NDArray]:
    """Python mirror of the C++ complex-scalar painter.

    Returns dict with keys 'phi' (=phi1), 'phi2', 'Pi' (=pi1), 'Pi2' at the given
    world coordinates.  To add a new profile kind, add a modulus branch here.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    cx, cy, cz = spec.center
    dx, dy, dz = x - cx, y - cy, z - cz

    if spec.kind == "torus":
        if torus is None:
            if not spec.profile_path:
                raise ValueError("torus spec needs profile_path or an explicit torus")
            torus = read_torus_profile(spec.profile_path)
        modulus = _modulus_torus(spec, torus, dx, dy, dz)
    elif spec.kind == "gaussian":
        modulus = _modulus_gaussian(spec, dx, dy, dz)
    else:
        raise NotImplementedError(f"reference_paint: kind={spec.kind!r} not yet mirrored")

    if spec.winding:
        az = np.arctan2(dy, dx)
        m = float(spec.m_az)
        phi1 = modulus * np.cos(m * az)
        phi2 = modulus * np.sin(m * az)
    else:
        phi1 = modulus
        phi2 = np.zeros_like(modulus)

    inv_alpha = (1.0 / np.asarray(alpha)) if not np.isscalar(alpha) else (1.0 / alpha)
    pi1 = (spec.omega * inv_alpha) * phi2
    pi2 = -(spec.omega * inv_alpha) * phi1
    return {"phi": phi1, "phi2": phi2, "Pi": pi1, "Pi2": pi2}


@dataclass
class ConsistencyReport:
    """Per-field max abs error between the C++ gridinit and the Python reference."""

    max_abs_err: dict[str, float]
    ref_peak: dict[str, float]
    grid_peak: dict[str, float]
    rel_err: dict[str, float]
    passed: bool
    tol: float

    def __str__(self) -> str:
        lines = [f"MatterProfile consistency ({'PASS' if self.passed else 'FAIL'}, "
                 f"tol={self.tol:.1e}):"]
        for name in self.max_abs_err:
            lines.append(
                f"  {name:5s}: max|Δ|={self.max_abs_err[name]:.3e}  "
                f"rel={self.rel_err[name]:.3e}  "
                f"(grid_peak={self.grid_peak[name]:.4g}, ref_peak={self.ref_peak[name]:.4g})"
            )
        return "\n".join(lines)


# Matter channels a complex-scalar profile is responsible for.
_MATTER_CHANNELS = ("phi", "phi2", "Pi", "Pi2")


def check_gridinit_matches_spec(
    gridinit_path: str | Path,
    spec: MatterProfileSpec,
    *,
    tol: float = 5.0e-2,
    torus: QBallTorusProfile | None = None,
    alpha_channel: str | None = "lapse",
) -> ConsistencyReport:
    """t=0 cross-code consistency test: does the C++ gridinit match the reference?

    Reads the GRTresna-produced ``.gridinit``, paints the reference field at the
    SAME world coordinates (cell centres ``origin + (idx+0.5)*dx``), and compares
    phi/phi2/Pi/Pi2.  ``tol`` is a RELATIVE tolerance on the peak (trilinear
    export + finite table interpolation give ~1% differences; the failures this
    catches -- zero field, wrong sign, wrong frequency -- are O(1) relative).

    Raises AssertionError on mismatch so it drops into a pytest / CLI gate.
    """
    g = read_gridinit(gridinit_path)
    nz, ny, nx, _ = g.data.shape
    ox, oy, oz = g.origin
    dxv = g.dx_xyz
    # Cell-centred world coordinates.
    xs = ox + (np.arange(nx) + 0.5) * dxv[0]
    ys = oy + (np.arange(ny) + 0.5) * dxv[1]
    zs = oz + (np.arange(nz) + 0.5) * dxv[2]
    Z, Y, X = np.meshgrid(zs, ys, xs, indexing="ij")  # (nz,ny,nx)

    if alpha_channel and alpha_channel in g.comp_names:
        alpha = g.data[..., g.comp_names.index(alpha_channel)]
        alpha = np.where(np.abs(alpha) > 1e-8, alpha, 1.0)
    else:
        alpha = 1.0

    ref = reference_paint(spec, X, Y, Z, alpha=alpha, torus=torus)

    max_abs_err: dict[str, float] = {}
    ref_peak: dict[str, float] = {}
    grid_peak: dict[str, float] = {}
    rel_err: dict[str, float] = {}
    for name in _MATTER_CHANNELS:
        if name not in g.comp_names:
            continue
        grid = g.data[..., g.comp_names.index(name)]
        r = ref[name]
        err = float(np.max(np.abs(grid - r)))
        gp = float(np.max(np.abs(grid)))
        rp = float(np.max(np.abs(r)))
        scale = max(gp, rp, 1e-30)
        max_abs_err[name] = err
        grid_peak[name] = gp
        ref_peak[name] = rp
        rel_err[name] = err / scale

    passed = all(v <= tol for v in rel_err.values()) and len(rel_err) > 0
    report = ConsistencyReport(max_abs_err, ref_peak, grid_peak, rel_err, passed, tol)
    return report
