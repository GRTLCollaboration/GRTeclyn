"""Splash-only GRTeclyn evolution overrides."""

from __future__ import annotations

from typing import Any


def splash_evolution_overrides() -> dict[str, Any]:
    """Constraint field plot vars for origin quality checks (splash campaigns only)."""
    return {
        "amr.derive_plot_vars": "Ham_abs_terms Mom_abs_terms rho_req",
    }
