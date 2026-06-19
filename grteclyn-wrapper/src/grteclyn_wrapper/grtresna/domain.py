"""Domain policy shared by GRTresna solves and GRTeclyn evolution params."""

from __future__ import annotations

from dataclasses import dataclass

from .gridinit_export import GridinitExportSpec
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

    @property
    def solve_center(self) -> tuple[float, float, float]:
        """Matter centre in GRTresna native coordinates."""
        half = 0.5 * float(self.grtresna_l)
        if self.full_z:
            return (half, half, half)
        return (half, half, 0.0)

    def evolution_extent(self) -> tuple[float, float, float]:
        """GRTeclyn evolution box size per axis (code units)."""
        length = float(self.l_full)
        if self.full_z:
            return (length, length, length)
        return (length, length, 0.5 * length)

    def gridinit_export_spec(self) -> GridinitExportSpec:
        """Export window matched to the evolution box and ``N_full`` spacing."""
        if self.full_z:
            lx, ly, lz = self.evolution_extent()
            source_origin = tuple(
                self.solve_center[i] - 0.5 * extent
                for i, extent in enumerate((lx, ly, lz))
            )
        else:
            # Legacy half-z: export the full solve domain; origin alignment handles
            # the GRTeclyn reflective plane (see ``convert_chombo_to_gridinit``).
            span = float(self.grtresna_l)
            lx = ly = lz = span
            source_origin = (0.0, 0.0, 0.0)
        return GridinitExportSpec(
            nx=int(self.gridinit_nx),
            ny=int(self.gridinit_ny),
            nz=int(self.gridinit_nz),
            lx=lx,
            ly=ly,
            lz=lz,
            target_center=self.target_center,
            source_origin=source_origin,
        )

    def apply_to_solver(self, cfg: GRTresnaConfig) -> GRTresnaConfig:
        cfg.N = self.solve_n
        cfg.L = float(self.grtresna_l)
        cfg.lo_boundary = self.solve_lo_boundary
        cfg.hi_boundary = self.hi_boundary
        cfg.gridinit_nx = int(self.gridinit_nx)
        cfg.gridinit_ny = int(self.gridinit_ny)
        cfg.gridinit_nz = int(self.gridinit_nz)
        cfg.target_center = self.target_center
        cfg.gridinit_export = self.gridinit_export_spec()
        return cfg

    def evolution_overrides(self) -> dict[str, object]:
        return {
            "L_full": float(self.l_full),
            "N_full": int(self.n_full),
            "center": self.evolution_center_override,
            "lo_boundary": self.evolution_lo_boundary_override,
        }
