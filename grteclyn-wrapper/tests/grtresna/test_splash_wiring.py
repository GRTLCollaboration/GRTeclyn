"""Tests for splash-only GRTeclyn evolution overrides."""

from __future__ import annotations

from grteclyn_wrapper.grtresna.matter_wiring import plot_vars_for_independent_scalars
from grteclyn_wrapper.grtresna.splash_wiring import splash_evolution_overrides


def test_default_scalar_plot_vars_exclude_ham_abs() -> None:
    names = plot_vars_for_independent_scalars(2)
    assert "Ham_abs_terms" not in names


def test_splash_evolution_overrides_add_ham_abs() -> None:
    overrides = splash_evolution_overrides()
    derive = overrides["amr.derive_plot_vars"]
    assert "Ham_abs_terms" in derive
    assert "Mom_abs_terms" in derive
    assert "rho_req" in derive
