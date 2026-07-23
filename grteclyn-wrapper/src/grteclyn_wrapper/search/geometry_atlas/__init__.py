"""Pure-geometry MAP-Elites atlas (Stage-1 inverse-design scout).

Searches smooth stationary asymptotically flat 4-metrics, scores frozen
null-geodesic ``f_geo`` and stationary free-fall ``f_ff``, and archives
shortcut strength against exotic-energy cost.  Matter fitting (GRTresna) and
GRTeclyn evolution are later promotion stages.
"""

from __future__ import annotations

from .calibrate import calibrate_atlas_probe
from .config import GeometryAtlasConfig
from .driver import run_geometry_atlas
from .genome import (
    GeometryGenome,
    GeometryGenomeConfig,
    mutate_genome,
    sample_genome,
    zero_genome,
)
from .refine import run_geometry_cmaes, seed_alcubierre_genome
from .render import RenderConfig, render_and_write, render_genome
from .score import GeometryAtlasEvaluation, evaluate_genome

__all__ = [
    "GeometryAtlasConfig",
    "GeometryAtlasEvaluation",
    "GeometryGenome",
    "GeometryGenomeConfig",
    "RenderConfig",
    "calibrate_atlas_probe",
    "evaluate_genome",
    "mutate_genome",
    "render_and_write",
    "render_genome",
    "run_geometry_atlas",
    "run_geometry_cmaes",
    "sample_genome",
    "seed_alcubierre_genome",
    "zero_genome",
]
