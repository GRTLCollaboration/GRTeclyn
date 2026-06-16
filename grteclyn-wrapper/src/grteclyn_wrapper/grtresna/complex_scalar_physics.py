"""Reference physics for canonical complex scalar matter (boson-star path).

Formulas match GRTeclyn ScalarField 1/2 convention extended to two real
components with a shared potential V(|Phi|^2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ComplexScalarState:
    """Pointwise matter + metric data for EM tensor evaluation."""

    phi1: float
    phi2: float
    pi1: float
    pi2: float
    chi: float
    h_uu: np.ndarray  # shape (3, 3) inverse metric
    h_ll: np.ndarray  # shape (3, 3) covariant metric
    dphi1: np.ndarray  # shape (3,)
    dphi2: np.ndarray  # shape (3,)


@dataclass(frozen=True)
class EMTensor:
    rho: float
    j: np.ndarray  # shape (3,)
    S: np.ndarray  # shape (3, 3)
    trS: float


def phi_modulus_sq(phi1: float, phi2: float) -> float:
    return phi1 * phi1 + phi2 * phi2


def potential_and_derivatives(
    phi1: float,
    phi2: float,
    mass: float,
    lam: float,
) -> tuple[float, float, float]:
    """Return V, dV/dphi1, dV/dphi2 for V = 1/2 m^2 |Phi|^2 - lambda/4 |Phi|^4."""
    mod2 = phi_modulus_sq(phi1, phi2)
    mphi = mass * math.sqrt(mod2) if mod2 > 0.0 else 0.0
    v = 0.5 * mphi * mphi - 0.25 * lam * mod2 * mod2
    if mod2 <= 0.0:
        return v, 0.0, 0.0
    dv_dmod2 = 0.5 * mass * mass - 0.5 * lam * mod2
    dphi1 = 2.0 * phi1 * dv_dmod2
    dphi2 = 2.0 * phi2 * dv_dmod2
    return v, dphi1, dphi2


def _component_emtensor(
    pi: float,
    dphi: np.ndarray,
    chi: float,
    h_uu: np.ndarray,
    h_ll: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Single real-component kinetic block (ScalarField convention)."""
    grad_sq = float(dphi @ (chi * (h_uu @ dphi)))
    vt = -pi * pi + grad_sq
    rho_k = pi * pi + 0.5 * vt
    j = -pi * dphi
    s = np.zeros((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            s[i, j] = (
                -0.5 * h_ll[i, j] * vt / chi
                + dphi[i] * dphi[j]
            )
    return rho_k, j, s


def compute_emtensor(
    state: ComplexScalarState,
    mass: float = 0.0,
    lam: float = 0.0,
) -> EMTensor:
    """Sum two ScalarField blocks + one shared V(|Phi|^2)."""
    v, _, _ = potential_and_derivatives(state.phi1, state.phi2, mass, lam)

    rho1, j1, s1 = _component_emtensor(
        state.pi1, state.dphi1, state.chi, state.h_uu, state.h_ll
    )
    rho2, j2, s2 = _component_emtensor(
        state.pi2, state.dphi2, state.chi, state.h_uu, state.h_ll
    )

    rho = rho1 + rho2 + v
    j = j1 + j2
    s = s1 + s2
    h_ll = state.h_ll
    for i in range(3):
        s[i, i] -= h_ll[i, i] * v / state.chi

    tr_s = state.chi * float(np.trace(state.h_uu @ s)) - 3.0 * v
    return EMTensor(rho=rho, j=j, S=s, trS=tr_s)


def conserved_charge_density(phi1: float, phi2: float, pi1: float, pi2: float) -> float:
    return phi1 * pi2 - phi2 * pi1


def integrate_charge(
    phi1: np.ndarray,
    phi2: np.ndarray,
    pi1: np.ndarray,
    pi2: np.ndarray,
    sqrt_det: np.ndarray,
) -> float:
    dens = phi1 * pi2 - phi2 * pi1
    return float(np.sum(dens * sqrt_det))


def phase_velocity_pi2(phi0: float, omega: float, lapse: float) -> float:
    """Pi2 for Phi ~ exp(-i omega t) with d_t phi = alpha Pi (wormhole convention)."""
    if lapse == 0.0:
        raise ValueError("lapse must be non-zero")
    return -omega * phi0 / lapse


def flat_minkowski_h() -> tuple[np.ndarray, np.ndarray]:
    h_uu = np.eye(3, dtype=float)
    h_ll = np.eye(3, dtype=float)
    return h_uu, h_ll


def random_test_state(rng: np.random.Generator) -> ComplexScalarState:
    h_uu, h_ll = flat_minkowski_h()
    return ComplexScalarState(
        phi1=float(rng.uniform(-0.5, 0.5)),
        phi2=float(rng.uniform(-0.5, 0.5)),
        pi1=float(rng.uniform(-0.3, 0.3)),
        pi2=float(rng.uniform(-0.3, 0.3)),
        chi=1.0,
        h_uu=h_uu,
        h_ll=h_ll,
        dphi1=rng.uniform(-0.2, 0.2, size=3),
        dphi2=rng.uniform(-0.2, 0.2, size=3),
    )
