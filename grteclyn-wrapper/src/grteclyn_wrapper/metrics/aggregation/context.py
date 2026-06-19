"""Resolved episode paths and inputs for metric collection."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class EpisodeContext:
    episode_dir: Path
    ftl_L: float | None
    collapse_path: Path
    constraint_path: Path
    areal_path: Path
    energy_conditions_path: Path
    curvature_path: Path
    boundary_flux_path: Path
    boundary_flux_fallback_path: Path
    gridinit_path: Path
    central_timeseries_path: Path
    central_radial_profile_path: Path
    ftl_timeseries_path: Path


def build_episode_context(episode_dir: Path, *, ftl_L: float | None = None) -> EpisodeContext:
    data_dir = episode_dir / "data"
    small_data_dir = episode_dir / "small_data"

    collapse_path = data_dir / "collapse_diagnostics.dat"
    constraint_path = data_dir / "constraint_norms.dat"
    areal_path = small_data_dir / "areal_radius.dat"
    if not collapse_path.exists():
        collapse_path = episode_dir / "collapse_diagnostics.dat"
    if not constraint_path.exists():
        constraint_path = episode_dir / "constraint_norms.dat"
    if not areal_path.exists():
        areal_path = episode_dir / "areal_radius.dat"

    energy_conditions_path = data_dir / "energy_conditions.dat"
    if not energy_conditions_path.exists():
        energy_conditions_path = episode_dir / "energy_conditions.dat"
    curvature_path = data_dir / "curvature_invariants.dat"
    if not curvature_path.exists():
        curvature_path = episode_dir / "curvature_invariants.dat"

    return EpisodeContext(
        episode_dir=episode_dir,
        ftl_L=ftl_L,
        collapse_path=collapse_path,
        constraint_path=constraint_path,
        areal_path=areal_path,
        energy_conditions_path=energy_conditions_path,
        curvature_path=curvature_path,
        boundary_flux_path=data_dir / "boundary_flux.dat",
        boundary_flux_fallback_path=episode_dir / "boundary_flux.dat",
        gridinit_path=episode_dir / "initial_data.gridinit",
        central_timeseries_path=small_data_dir / "central_timeseries.dat",
        central_radial_profile_path=small_data_dir / "central_radial_profile.dat",
        ftl_timeseries_path=small_data_dir / "ftl_timeseries.dat",
    )
