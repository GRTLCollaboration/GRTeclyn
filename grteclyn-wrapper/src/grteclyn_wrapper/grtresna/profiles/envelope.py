"""Shared lump profile envelopes (mirrors GRTresna MatterParams.hpp)."""

from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from .boson_star import (
    PROFILE_ODE_BOUND,
    PROFILE_SECH_BOUND,
    PROFILE_SELFGRAV_BOUND,
)

EXOTIC_AMP_SCALE = 0.25


def angular_factor(mode: int, dx: float | NDArray, dy: float | NDArray, width: float) -> float | NDArray:
    if mode == 1:
        return dx / width
    if mode == 2:
        return (dx * dx - dy * dy) / (width * width)
    return 1.0 if isinstance(dx, float) else np.ones_like(dx, dtype=np.float64)


def effective_amp(lump: Mapping[str, Any], *, raw: bool = False) -> float:
    amp = float(lump.get("amp", 0.0))
    if amp == 0.0:
        return 0.0
    if raw:
        return amp
    return EXOTIC_AMP_SCALE * amp if int(lump.get("exotic", 0)) else amp


def envelope_scalar(profile: int, r2: float, width: float) -> float:
    """Radial envelope at squared radius r2 (scalar path)."""
    if profile == 1:
        r = math.sqrt(r2)
        soft = 0.25 * width
        return 0.5 * (1.0 - math.tanh((r - width) / soft))
    if profile == PROFILE_SECH_BOUND:
        r = math.sqrt(r2)
        return 1.0 / math.cosh(r / width)
    return math.exp(-r2 / (2.0 * width * width))


def envelope_grid(profile: int, r2: NDArray, width: float) -> NDArray:
    """Radial envelope on a grid (vectorized)."""
    if profile == 1:
        r = np.sqrt(r2)
        soft = 0.25 * width
        return 0.5 * (1.0 - np.tanh((r - width) / soft))
    if profile == PROFILE_SECH_BOUND:
        r = np.sqrt(r2)
        return 1.0 / np.cosh(r / width)
    return np.exp(-r2 / (2.0 * width * width))


def phi0_at_radius(r: NDArray, lump: Mapping[str, Any], *, raw_amp: bool) -> NDArray:
    """Radial profile phi0(r) including ODE-bound dispatch."""
    amp = float(lump.get("amp", 0.0))
    if amp == 0.0:
        return np.zeros_like(r, dtype=np.float64)
    width = float(lump.get("width", 5.0))
    profile_type = int(lump.get("profile", 0))
    scale = amp if raw_amp else effective_amp(lump)
    if profile_type == PROFILE_ODE_BOUND:
        from .qball_ode import profile_for_lump

        radial = profile_for_lump(
            dict(lump),
            mass=float(lump.get("qball_mass", 1.0)),
            lam=float(lump.get("qball_lam", 0.0)),
            mu=float(lump.get("qball_mu", 0.0)),
            omega=float(lump.get("qball_omega", 0.0)),
        )
        return np.asarray(radial.eval_phi0(r), dtype=np.float64)
    if profile_type == PROFILE_SELFGRAV_BOUND:
        # Gravitationally dressed star: amp == the star's own phi_c by
        # construction, so no amp rescale applies.  Falling through to the
        # Gaussian here would paint evolution matter that contradicts the
        # constraint solve's tabulated star.  Dispatch must mirror
        # _write_selfgrav_profile so both sides hit the same cache entry.
        q_mass = float(lump.get("qball_mass", 1.0))
        q_lam = float(lump.get("qball_lam", 0.0))
        q_mu = float(lump.get("qball_mu", 0.0))
        q_omega = float(lump.get("qball_omega", 0.0))
        if q_omega > 0.0 and q_lam > 0.0 and q_mu > 0.0:
            from .boson_star_ode import cached_selfgrav_at_omega

            gravity_sign = -1.0 if int(lump.get("exotic", 0)) else 1.0
            star = cached_selfgrav_at_omega(q_mass, q_lam, q_mu, q_omega, gravity_sign)
        else:
            from .boson_star_ode import cached_selfgrav_profile

            star = cached_selfgrav_profile(q_mass, q_lam, q_mu, amp)
        return np.asarray(star.eval_phi0(r), dtype=np.float64)
    if profile_type == PROFILE_SECH_BOUND:
        return scale / np.cosh(r / width)
    if profile_type == 1:
        soft = 0.25 * width
        return scale * 0.5 * (1.0 - np.tanh((r - width) / soft))
    return scale * np.exp(-r * r / (2.0 * width * width))
