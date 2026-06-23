"""Closed-loop RL control for GRTeclyn FTL campaigns."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from grteclyn_wrapper.rl.env import SpacetimeFtlEnv

__all__ = ["SpacetimeFtlEnv"]


def __getattr__(name: str):
    if name == "SpacetimeFtlEnv":
        from grteclyn_wrapper.rl.env import SpacetimeFtlEnv

        return SpacetimeFtlEnv
    raise AttributeError(name)
