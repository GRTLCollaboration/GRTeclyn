"""Tests for splash-only GRTeclyn evolution overrides."""

from __future__ import annotations

from grteclyn_wrapper.grtresna.matter.wiring import plot_vars_for_complex_scalar
from grteclyn_wrapper.grtresna.matter.splash import (
    append_gw_proxy_plot_vars,
    apply_splash_overrides,
    splash_evolution_overrides,
)


def test_default_scalar_plot_vars_exclude_ham_abs() -> None:
    names = plot_vars_for_complex_scalar()
    assert "Ham_abs_terms" not in names


def test_splash_evolution_overrides_add_ham_abs() -> None:
    overrides = splash_evolution_overrides()
    derive = overrides["amr.derive_plot_vars"]
    assert "Ham_abs_terms" in derive
    assert "Mom_abs_terms" in derive
    assert "rho_req" in derive
    assert "Weyl4" not in derive


def test_splash_plot_vars_include_aij_for_gw_proxy() -> None:
    overrides = splash_evolution_overrides()
    plot_vars = overrides["amr.plot_vars"]
    assert "A11" in plot_vars
    assert "A12" in plot_vars
    assert "A22" in plot_vars


def test_append_gw_proxy_plot_vars_deduplicates() -> None:
    base = ("chi", "K", "A11")
    merged = append_gw_proxy_plot_vars(base)
    assert merged.count("A11") == 1
    assert "A12" in merged


def test_apply_splash_overrides_merges_existing_plot_vars() -> None:
    merged = apply_splash_overrides({"amr.plot_vars": plot_vars_for_complex_scalar()})
    assert "phi" in merged["amr.plot_vars"]
    assert "A22" in merged["amr.plot_vars"]
