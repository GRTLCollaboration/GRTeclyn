"""Profiles for 4D evolving null-geodesic integration (search vs HQ verify)."""

from __future__ import annotations

import os
from dataclasses import dataclass, replace


@dataclass(frozen=True)
class EvolvingGeodesicOptions:
    """Knobs for end-of-run 4D null-ray integration."""

    compute_frozen_peak: bool = True
    slice_stride: int = 1
    max_slices: int | None = None
    n_rays: int = 5
    max_steps: int = 50_000
    h_rel_abort: float | None = None
    ds_init: float = 0.05
    directions: tuple[str, ...] = ("x",)


SEARCH_OPTIONS = EvolvingGeodesicOptions(
    compute_frozen_peak=False,
    slice_stride=2,
    max_slices=40,
    n_rays=3,
    max_steps=15_000,
    h_rel_abort=0.5,
)

HQ_OPTIONS = EvolvingGeodesicOptions()


def geo_directions_from_env() -> tuple[str, ...]:
    """Principal axes to scan for the end-of-run 4D trace (``GRTECLYN_GEO_DIRECTIONS``)."""
    raw = os.environ.get("GRTECLYN_GEO_DIRECTIONS", "x").strip()
    dirs = tuple(tok for tok in raw.split() if tok in {"x", "y", "z"})
    return dirs or ("x",)


def evolving_geodesic_options_from_env() -> EvolvingGeodesicOptions:
    """Resolve integration profile from ``GRTECLYN_EVOLVING_GEODESIC_MODE``."""
    mode = os.environ.get("GRTECLYN_EVOLVING_GEODESIC_MODE", "search").strip().lower()
    directions = geo_directions_from_env()
    base = HQ_OPTIONS if mode in {"hq", "full", "verify", "promote"} else SEARCH_OPTIONS
    return replace(base, directions=directions)


def metric_stack_n_space_from_env(*, default: int = 65) -> int:
    """Spatial resolution for cached metric slices (``GRTECLYN_METRIC_STACK_N_SPACE``)."""
    raw = os.environ.get("GRTECLYN_METRIC_STACK_N_SPACE", "").strip()
    if not raw:
        mode = os.environ.get("GRTECLYN_EVOLVING_GEODESIC_MODE", "search").strip().lower()
        return 33 if mode == "search" else default
    return max(17, int(raw))
