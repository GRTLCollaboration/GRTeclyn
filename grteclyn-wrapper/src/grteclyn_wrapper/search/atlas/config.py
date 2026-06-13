"""Atlas sampling configuration and path helpers."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from ...core.config import ExampleConfig


DEFAULT_RANGES: dict[str, tuple[float, float]] = {
    "wormhole_phi_perturbation_amplitude": (-0.04, 0.04),
    "wormhole_support_strength": (0.2, 1.0),
    "wormhole_phi_perturbation_width": (0.25, 1.0),
    "wormhole_phi_monopole_amplitude": (-0.02, 0.02),
}

DEFAULT_RECIPE_RANGES: dict[str, tuple[float, float]] = {
    "recipe_chi_coeff_0": (-0.5, 0.1),
    "recipe_chi_coeff_1": (-0.2, 0.2),
    "recipe_phi_coeff_0": (-0.1, 0.1),
    "recipe_basis_width": (0.5, 2.0),
}

DEFAULT_LOW_RES_OVERRIDES: dict[str, Any] = {
    "N_full": 32,
    "max_level": 0,
    "stop_time": 0.04,
    "plot_interval": 1000,
    "checkpoint_interval": 1000,
}


def atlas_ranges_for_example(example: ExampleConfig) -> dict[str, tuple[float, float]]:
    if example.name == "RadialRecipe":
        return DEFAULT_RECIPE_RANGES
    return DEFAULT_RANGES


@dataclass(frozen=True)
class AtlasThresholds:
    constraint_blowup: float = 1.0e-1
    lapse_collapse: float = 1.0e-3
    horizon_radius: float = 1.0e-8
    trivial_geometry: float = 1.0e-3


@dataclass(frozen=True)
class AtlasPaths:
    root: Path
    jsonl: Path
    csv: Path
    summary: Path


def make_atlas_dir(runs_dir: Path, name: str | None = None) -> AtlasPaths:
    if name is None:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"atlas_{timestamp}"
    root = (runs_dir / name).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=False)
    return AtlasPaths(
        root=root,
        jsonl=root / "atlas.jsonl",
        csv=root / "atlas.csv",
        summary=root / "summary.json",
    )


def sample_overrides(
    rng,
    *,
    base: Mapping[str, Any] | None = None,
    ranges: Mapping[str, tuple[float, float]] | None = None,
) -> dict[str, Any]:
    overrides = dict(DEFAULT_LOW_RES_OVERRIDES)
    if base:
        overrides.update(base)

    for key, (lo, hi) in (ranges or DEFAULT_RANGES).items():
        overrides.setdefault(key, rng.uniform(lo, hi))
    return overrides
