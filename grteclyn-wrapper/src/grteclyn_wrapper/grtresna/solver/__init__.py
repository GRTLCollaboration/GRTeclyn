"""GRTresna elliptic solver orchestration."""

from .config import GRTresnaConfig, apply_exotic_safe_solver, config_has_exotic_lump
from .convergence import parse_convergence
from .params import _lump_lines, write_grtresna_params
from .runner import solve

__all__ = [
    "GRTresnaConfig",
    "apply_exotic_safe_solver",
    "config_has_exotic_lump",
    "parse_convergence",
    "solve",
    "write_grtresna_params",
    "_lump_lines",
]
