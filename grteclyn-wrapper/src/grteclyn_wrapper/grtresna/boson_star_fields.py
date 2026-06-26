"""Paint boson-star complex scalar fields onto a gridinit array."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray

from .boson_star_profile import BosonStarProfile
from .lump_fields import _coordinate_grids, angular_factor

# GRTeclyn RadialRecipe state layout (must match StateVariables.hpp).
_NUM_CCZ4_VARS = 25
# Canonical complex field Phi+ (slots 25-28) and phantom field Phi- (29-32).
_C_PHI_P = _NUM_CCZ4_VARS + 0
_C_PI_P = _NUM_CCZ4_VARS + 1
_C_PHI2_P = _NUM_CCZ4_VARS + 2
_C_PI2_P = _NUM_CCZ4_VARS + 3
_C_PHI_M = _NUM_CCZ4_VARS + 4
_C_PI_M = _NUM_CCZ4_VARS + 5
_C_PHI2_M = _NUM_CCZ4_VARS + 6
_C_PI2_M = _NUM_CCZ4_VARS + 7
_BICOMPLEX_NCOMP = _NUM_CCZ4_VARS + 8

# GRTresna HDF5 names -> GRTeclyn RadialRecipe names
_COMPLEX_SCALAR_RENAMES = {
    "phi_re": "phi",
    "Pi_re": "Pi",
    "phi_im": "phi2",
    "Pi_im": "Pi2",
}


def rename_complex_scalar_components(comp_names: list[str]) -> list[str]:
    return [_COMPLEX_SCALAR_RENAMES.get(name, name) for name in comp_names]


def apply_post_solve_lapse_correction(
    data: NDArray, comp_names: Sequence[str]
) -> NDArray:
    """Divide the painted Pi2 by the solved lapse.

    GRTresna paints the imaginary momentum as ``Pi_im = -omega*phi0`` assuming
    ``alpha = 1``.  A stationary boson star needs ``Pi_im = -omega*phi0/alpha``,
    so once the constraint solve has fixed the lapse the evolved initial data
    must rescale ``Pi2 /= alpha``.  Component names must already be renamed to
    the GRTeclyn convention (``Pi2``, ``lapse``).
    """
    idx = {name: i for i, name in enumerate(comp_names)}
    if "Pi2" not in idx or "lapse" not in idx:
        return data
    out = np.array(data, dtype=np.float64, copy=True)
    lapse = out[:, :, :, idx["lapse"]]
    with np.errstate(divide="ignore", invalid="ignore"):
        out[:, :, :, idx["Pi2"]] = out[:, :, :, idx["Pi2"]] / np.maximum(
            lapse, 1.0e-12
        )
    return out


def paint_boson_star_fields_on_grid(
    data: NDArray,
    comp_names: list[str],
    dx: NDArray,
    origin: NDArray,
    profile: BosonStarProfile,
    *,
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
    apply_post_solve_pi2: bool = True,
) -> tuple[NDArray, list[str]]:
    """Paint phi, phi2, Pi, Pi2 from a radial BS profile onto the grid."""
    names = rename_complex_scalar_components(list(comp_names))
    name_to_idx = {name: i for i, name in enumerate(names)}

    required = ("phi", "Pi", "phi2", "Pi2")
    for key in required:
        if key not in name_to_idx:
            raise ValueError(f"Missing component {key!r} in gridinit names: {names}")

    nz, ny, nx, _ = data.shape
    px_grid, py_grid, pz_grid = _coordinate_grids(nx, ny, nz, dx, origin)
    cx, cy, cz = center
    rr = np.sqrt(
        (px_grid - cx) ** 2 + (py_grid - cy) ** 2 + (pz_grid - cz) ** 2
    )

    phi1 = np.asarray(profile.eval_phi0(rr), dtype=np.float64)
    phi2 = np.zeros_like(phi1)
    pi1 = np.zeros_like(phi1)

    if apply_post_solve_pi2 and "lapse" in name_to_idx:
        lapse = data[:, :, :, name_to_idx["lapse"]]
        with np.errstate(divide="ignore", invalid="ignore"):
            pi2 = -profile.omega * phi1 / np.maximum(lapse, 1.0e-12)
    else:
        pi2 = -profile.omega * phi1

    painted = np.array(data, dtype=np.float64, copy=True)
    painted[:, :, :, name_to_idx["phi"]] = phi1
    painted[:, :, :, name_to_idx["phi2"]] = phi2
    painted[:, :, :, name_to_idx["Pi"]] = pi1
    painted[:, :, :, name_to_idx["Pi2"]] = pi2

    return painted, names


def _raw_lump_phi_grid(
    lump: Mapping[str, Any], px: NDArray, py: NDArray, pz: NDArray,
) -> NDArray:
    """phi of one lump using the RAW amplitude (no exotic 0.25 scaling).

    Matches GRTresna's ``BosonStarParams::lump_phi1`` so the painted Phi+/Phi-
    fields agree with the per-lump-signed constraint solve.  The exotic sign is
    carried by the phantom FIELD (Phi-), not by an amplitude rescale.
    """
    amp = float(lump.get("amp", 0.0))
    if amp == 0.0:
        return np.zeros_like(px, dtype=np.float64)
    width = float(lump.get("width", 5.0))
    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=np.float64)
    dx = px - center[0]
    dy = py - center[1]
    dz = pz - center[2]
    r2 = dx * dx + dy * dy + dz * dz
    profile = int(lump.get("profile", 0))
    if profile == 1:
        r = np.sqrt(r2)
        soft = 0.25 * width
        env = 0.5 * (1.0 - np.tanh((r - width) / soft))
    elif profile == 2:
        # Bound boson lump sech(r/width) (PROFILE_SECH_BOUND) -- correct
        # exponential tail; must match GRTresna lump_envelope profile==2.
        r = np.sqrt(r2)
        env = 1.0 / np.cosh(r / width)
    else:
        env = np.exp(-r2 / (2.0 * width * width))
    return amp * angular_factor(int(lump.get("mode", 0)), dx, dy, width) * env


def _boosted_lump_fields(
    lump: Mapping[str, Any],
    px: NDArray, py: NDArray, pz: NDArray,
    omega: float,
    inv_lapse: NDArray,
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """Compute (phi_re, phi_im, Pi_re, Pi_im) for a Lorentz-boosted Q-ball.

    For a Q-ball at rest the field is ``Φ(x,t) = φ₀(r) exp(+iωt)`` (matching
    the GRTeclyn convention ``Pi2 = -ω φ₀/α``).  Under a Lorentz boost with
    velocity *v* in the lab frame, at t=0:

        Φ(x, 0) = φ₀(ρ) exp(−i θ)
        ∂_t Φ|₀ = [φ₀'(ρ) D + iωγ φ₀(ρ)] exp(−i θ)

    where ρ = √(δx_⊥² + γ² δx_∥²) is the Lorentz-contracted radius,
    θ = ωγ (v · δx) is the de Broglie phase, D = −γ²|v| δx_∥ / ρ is the
    advection factor, and γ = 1/√(1−|v|²).

    Returns zero fields when lump amplitude is zero.  Falls back to the
    standard at-rest initialization when velocity is negligible (|v| < 1e-10).
    """
    amp = float(lump.get("amp", 0.0))
    zero = np.zeros_like(px, dtype=np.float64)
    if amp == 0.0:
        return zero, zero.copy(), zero.copy(), zero.copy()

    width = float(lump.get("width", 5.0))
    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=np.float64)
    velocity = np.asarray(lump.get("velocity", (0.0, 0.0, 0.0)), dtype=np.float64)

    # Displacement from lump center.
    ddx = px - center[0]
    ddy = py - center[1]
    ddz = pz - center[2]

    v_mag = float(np.sqrt(np.dot(velocity, velocity)))

    if v_mag < 1.0e-10:
        # At rest: standard initialization.
        r = np.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
        profile_type = int(lump.get("profile", 0))
        if profile_type == 2:
            phi0 = amp / np.cosh(r / width)
        elif profile_type == 1:
            soft = 0.25 * width
            phi0 = amp * 0.5 * (1.0 - np.tanh((r - width) / soft))
        else:
            phi0 = amp * np.exp(-r * r / (2.0 * width * width))
        return phi0, zero.copy(), zero.copy(), -omega * phi0 * inv_lapse

    # --- Lorentz boost ---
    gamma = 1.0 / np.sqrt(1.0 - v_mag * v_mag)
    v_hat = velocity / v_mag

    # Parallel displacement (scalar projection onto velocity direction).
    dx_par = ddx * v_hat[0] + ddy * v_hat[1] + ddz * v_hat[2]
    # Perpendicular squared: |δx|² - δx_∥²
    r2 = ddx * ddx + ddy * ddy + ddz * ddz
    dx_perp2 = np.maximum(r2 - dx_par * dx_par, 0.0)

    # Lorentz-contracted radial distance.
    rho = np.sqrt(dx_perp2 + gamma * gamma * dx_par * dx_par)

    # Profile and its radial derivative at contracted radius.
    profile_type = int(lump.get("profile", 0))
    if profile_type == 2:
        phi0_rho = amp / np.cosh(rho / width)
        dphi0_rho = -phi0_rho * np.tanh(rho / width) / width
    elif profile_type == 1:
        soft = 0.25 * width
        arg = (rho - width) / soft
        phi0_rho = amp * 0.5 * (1.0 - np.tanh(arg))
        dphi0_rho = -amp * 0.5 / (soft * np.cosh(arg) ** 2)
    else:
        phi0_rho = amp * np.exp(-rho * rho / (2.0 * width * width))
        dphi0_rho = -phi0_rho * rho / (width * width)

    # Phase: θ = ω γ (v · δx).
    v_dot_dx = velocity[0] * ddx + velocity[1] * ddy + velocity[2] * ddz
    theta = omega * gamma * v_dot_dx

    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    # Boosted field components at t=0: Φ = φ₀(ρ) exp(−iθ).
    phi_re = phi0_rho * cos_theta
    phi_im = -phi0_rho * sin_theta

    # Time derivative for Pi = −∂_t φ / α.
    # D = ∂ρ/∂t|_{t=0} = −γ² |v| δx_∥ / ρ
    rho_safe = np.maximum(rho, 1.0e-30)
    D = -gamma * gamma * v_mag * dx_par / rho_safe

    # ∂_t Φ|₀ = [φ₀'(ρ)·D + iωγ φ₀(ρ)] · exp(−iθ)
    # Real: φ₀'D cosθ + ωγ φ₀ sinθ
    # Imag: ωγ φ₀ cosθ − φ₀'D sinθ
    omega_gamma = omega * gamma
    dt_phi_re = dphi0_rho * D * cos_theta + omega_gamma * phi0_rho * sin_theta
    dt_phi_im = omega_gamma * phi0_rho * cos_theta - dphi0_rho * D * sin_theta

    Pi_re = -dt_phi_re * inv_lapse
    Pi_im = -dt_phi_im * inv_lapse

    return phi_re, phi_im, Pi_re, Pi_im


def paint_bicomplex_fields_on_grid(
    data: NDArray,
    comp_names: list[str],
    dx: NDArray,
    origin: NDArray,
    lumps: Sequence[Mapping[str, Any]],
    *,
    omega: float,
) -> tuple[NDArray, list[str]]:
    """Paint TWO complex fields for the ``grtresna_bicomplex_scalar`` model.

    Canonical lumps (``exotic`` == 0) are superposed into the CANONICAL field
    Phi+ (state slots 25-28); exotic lumps into the PHANTOM field Phi- (slots
    29-32).  Each field carries the full Q-ball state:

    - **At rest** (velocity=0): ``phi_im = 0``, ``Pi_im = -omega*phi_re/lapse``
    - **Boosted** (velocity≠0): full Lorentz-boosted complex field with
      contracted profile, de Broglie phase oscillation, and advection momentum.

    The metric (slots 0-24, incl. the solved lapse) is taken from GRTresna
    unchanged.  The momentum constraint is not re-solved after the boost; at
    Q-ball amplitudes (~0.15) the momentum source is O(0.3%·v) which is
    sub-dominant to the solve tolerance.
    """
    nz, ny, nx, ncomp = data.shape
    # Ensure room for the full canonical+phantom matter block.
    if ncomp < _BICOMPLEX_NCOMP:
        pad = np.zeros((nz, ny, nx, _BICOMPLEX_NCOMP - ncomp), dtype=np.float64)
        data = np.concatenate([data, pad], axis=3)
        ncomp = _BICOMPLEX_NCOMP

    out = np.array(data, dtype=np.float64, copy=True)

    name_to_idx = {name: i for i, name in enumerate(comp_names)}
    lapse_idx = name_to_idx.get("lapse", None)
    lapse = (
        out[:, :, :, lapse_idx]
        if lapse_idx is not None
        else np.ones((nz, ny, nx), dtype=np.float64)
    )
    inv_lapse = 1.0 / np.maximum(lapse, 1.0e-12)

    px_grid, py_grid, pz_grid = _coordinate_grids(nx, ny, nz, dx, origin)

    # Accumulate per-group (canonical / phantom) via boosted field computation.
    phi_re_p = np.zeros((nz, ny, nx), dtype=np.float64)
    phi_im_p = np.zeros((nz, ny, nx), dtype=np.float64)
    Pi_re_p = np.zeros((nz, ny, nx), dtype=np.float64)
    Pi_im_p = np.zeros((nz, ny, nx), dtype=np.float64)

    phi_re_m = np.zeros((nz, ny, nx), dtype=np.float64)
    phi_im_m = np.zeros((nz, ny, nx), dtype=np.float64)
    Pi_re_m = np.zeros((nz, ny, nx), dtype=np.float64)
    Pi_im_m = np.zeros((nz, ny, nx), dtype=np.float64)

    for lump in lumps:
        is_exotic = bool(int(lump.get("exotic", 0)))
        fr, fi, pr, pi = _boosted_lump_fields(
            lump, px_grid, py_grid, pz_grid, omega, inv_lapse,
        )
        if is_exotic:
            phi_re_m += fr
            phi_im_m += fi
            Pi_re_m += pr
            Pi_im_m += pi
        else:
            phi_re_p += fr
            phi_im_p += fi
            Pi_re_p += pr
            Pi_im_p += pi

    # Zero out the whole matter block first, then write the two fields.
    out[:, :, :, _NUM_CCZ4_VARS:_BICOMPLEX_NCOMP] = 0.0
    # Canonical Phi+.
    out[:, :, :, _C_PHI_P] = phi_re_p
    out[:, :, :, _C_PHI2_P] = phi_im_p
    out[:, :, :, _C_PI_P] = Pi_re_p
    out[:, :, :, _C_PI2_P] = Pi_im_p
    # Phantom Phi-.
    out[:, :, :, _C_PHI_M] = phi_re_m
    out[:, :, :, _C_PHI2_M] = phi_im_m
    out[:, :, :, _C_PI_M] = Pi_re_m
    out[:, :, :, _C_PI2_M] = Pi_im_m

    names = list(comp_names)
    matter_names = [
        "phi", "Pi", "phi2", "Pi2",
        "phi_lump1", "Pi_lump1", "phi_lump2", "Pi_lump2",
    ]
    for slot, nm in zip(range(_NUM_CCZ4_VARS, _BICOMPLEX_NCOMP), matter_names):
        if slot < len(names):
            names[slot] = nm
        else:
            names.append(nm)
    return out, names
