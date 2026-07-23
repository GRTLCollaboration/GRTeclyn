"""Configuration for the pure-geometry MAP-Elites atlas."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

from .genome import GeometryGenomeConfig
from .render import RenderConfig


@dataclass
class GeometryAtlasConfig:
    """Campaign knobs for Stage-1 stationary geometry search."""

    runs_dir: Path = Path("runs/geometry_atlas")
    name: str | None = None
    target_evals: int = 32
    bins: int = 8
    seed: int = 7
    batch_size: int = 4
    n_rays: int = 3
    compute_ff: bool = True
    resume: bool = False
    # Fraction of each batch drawn as fresh random samples (rest = mutate elites).
    random_fraction: float = 0.35
    # Localise emitter/detector about the genome support (recommended).
    localise_probe: bool = True
    # Soft cost: score -= exotic_penalty * integral_negative_rho.
    exotic_penalty: float = 0.0
    # Soft reward: score += exotic_bonus * integral_negative_rho (trust exotic).
    exotic_bonus: float = 0.0
    # Hard reject if integral_negative_rho > exotic_ban (0 disables).
    exotic_ban: float = 0.0
    # JSON genome paths injected into the archive at seeding (e.g. a CMA elite).
    warm_start_genomes: tuple[str, ...] = ()
    genome: GeometryGenomeConfig = field(default_factory=GeometryGenomeConfig)
    render: RenderConfig = field(default_factory=lambda: RenderConfig(n=24, L=64.0))

    def to_dict(self) -> dict:
        return {
            "runs_dir": str(self.runs_dir),
            "name": self.name,
            "target_evals": self.target_evals,
            "bins": self.bins,
            "seed": self.seed,
            "batch_size": self.batch_size,
            "n_rays": self.n_rays,
            "compute_ff": self.compute_ff,
            "resume": self.resume,
            "random_fraction": self.random_fraction,
            "localise_probe": self.localise_probe,
            "exotic_penalty": self.exotic_penalty,
            "exotic_bonus": self.exotic_bonus,
            "exotic_ban": self.exotic_ban,
            "warm_start_genomes": list(self.warm_start_genomes),
            "genome": asdict(self.genome),
            "render": {"n": self.render.n, "L": self.render.L},
        }
