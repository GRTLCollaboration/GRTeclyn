"""Domain policy shared by GRTresna solves and GRTeclyn evolution params."""

from __future__ import annotations

from dataclasses import dataclass

from .solver import GRTresnaConfig


@dataclass(frozen=True)
class GRTresnaDomainConfig:
    """Configurable domain mapping for GRTresna -> GRTeclyn `.gridinit` runs."""

    full_z: bool = False
    l_full: float = 64.0
    n_full: int = 64
    grtresna_l: float = 128.0
    grtresna_nx: int = 64
    grtresna_ny: int = 64
    grtresna_nz: int | None = None
    gridinit_nx: int = 64
    gridinit_ny: int = 64
    gridinit_nz: int = 64
    center: tuple[float, float, float] | None = None
    lo_boundary: tuple[int, int, int] | None = None
    hi_boundary: tuple[int, int, int] = (0, 0, 0)
    evolution_lo_boundary: tuple[int, int, int] | None = None

    @property
    def target_center(self) -> tuple[float, float, float]:
        if self.center is not None:
            return self.center
        half = 0.5 * float(self.l_full)
        return (half, half, half if self.full_z else 0.0)

    @property
    def solve_n(self) -> tuple[int, int, int]:
        nz = self.grtresna_nz
        if nz is None:
            nz = self.grtresna_nx if self.full_z else max(1, self.grtresna_nx // 2)
        return (int(self.grtresna_nx), int(self.grtresna_ny), int(nz))

    @property
    def solve_lo_boundary(self) -> tuple[int, int, int]:
        if self.lo_boundary is not None:
            return self.lo_boundary
        return (0, 0, 0) if self.full_z else (0, 0, 1)

    @property
    def evolution_center_override(self) -> str:
        return " ".join(f"{value:g}" for value in self.target_center)

    @property
    def evolution_lo_boundary_override(self) -> str:
        boundary = self.evolution_lo_boundary
        if boundary is None:
            boundary = (1, 1, 1) if self.full_z else (1, 1, 2)
        return " ".join(str(value) for value in boundary)

    def apply_to_solver(self, cfg: GRTresnaConfig) -> GRTresnaConfig:
        cfg.N = self.solve_n
        cfg.L = float(self.grtresna_l)
        cfg.lo_boundary = self.solve_lo_boundary
        cfg.hi_boundary = self.hi_boundary
        cfg.gridinit_nx = int(self.gridinit_nx)
        cfg.gridinit_ny = int(self.gridinit_ny)
        cfg.gridinit_nz = int(self.gridinit_nz)
        cfg.target_center = self.target_center
        return cfg

    def evolution_overrides(self) -> dict[str, object]:
        return {
            "L_full": float(self.l_full),
            "N_full": int(self.n_full),
            "center": self.evolution_center_override,
            "lo_boundary": self.evolution_lo_boundary_override,
        }
