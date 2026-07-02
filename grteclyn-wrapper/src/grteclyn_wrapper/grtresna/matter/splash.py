"""Splash-only GRTeclyn evolution overrides."""

from __future__ import annotations

from typing import Any

from .wiring import plot_vars_for_complex_scalar

_GW_PROXY_PLOT_VARS = ("A11", "A12", "A22")


def append_gw_proxy_plot_vars(plot_vars: Any) -> tuple[str, ...]:
    """Append native A_ij components needed for GW strain proxy extraction."""
    if isinstance(plot_vars, str):
        names = tuple(plot_vars.split())
    elif plot_vars is None:
        names = ()
    else:
        names = tuple(plot_vars)
    out = list(names)
    for name in _GW_PROXY_PLOT_VARS:
        if name not in out:
            out.append(name)
    return tuple(out)


def splash_evolution_overrides() -> dict[str, Any]:
    """Constraint + GW-proxy plot vars for splash campaigns only."""
    return {
        "amr.derive_plot_vars": "Ham_abs_terms Mom_abs_terms rho_req",
        "amr.plot_vars": append_gw_proxy_plot_vars(plot_vars_for_complex_scalar()),
    }


def apply_splash_overrides(overrides: dict[str, Any]) -> dict[str, Any]:
    """Merge splash overrides into episode GRTeclyn params."""
    merged = dict(overrides)
    splash = splash_evolution_overrides()
    merged.update(splash)
    if "amr.plot_vars" in merged:
        merged["amr.plot_vars"] = append_gw_proxy_plot_vars(merged["amr.plot_vars"])
    return merged
