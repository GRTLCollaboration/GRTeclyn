"""Unit tests for the shared-sum lambda*phi^4 potential convention."""

from __future__ import annotations

import math


def _potential(s: float, mass: float, lam: float) -> float:
    mphi = mass * s
    return 0.5 * mphi * mphi - 0.25 * lam * (s ** 4)


def _dpotential_ds(s: float, mass: float, lam: float) -> float:
    mphi = mass * s
    return mass * mphi - lam * (s ** 3)


def test_lambda_zero_matches_mass_only() -> None:
    s, mass = 0.12, 0.6
    assert _potential(s, mass, 0.0) == 0.5 * (mass * s) ** 2
    assert _dpotential_ds(s, mass, 0.0) == mass * mass * s


def test_quartic_lowers_potential_at_large_amplitude() -> None:
    mass, lam = 0.6, 0.05
    s = 2.0
    mass_only = 0.5 * (mass * s) ** 2
    assert _potential(s, mass, lam) < mass_only
    assert _dpotential_ds(0.0, mass, lam) == 0.0


def test_grtresna_params_include_scalar_lambda_default() -> None:
    from pathlib import Path

    from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, write_grtresna_params

    cfg = GRTresnaConfig(scalar_lambda=0.0)
    path = Path("/tmp/test_scalar_lambda_params.txt")
    write_grtresna_params(cfg, path)
    text = path.read_text()
    assert "scalar_lambda = 0.0" in text

    cfg_lam = GRTresnaConfig(scalar_lambda=0.05)
    write_grtresna_params(cfg_lam, path)
    assert "scalar_lambda = 0.05" in path.read_text()
