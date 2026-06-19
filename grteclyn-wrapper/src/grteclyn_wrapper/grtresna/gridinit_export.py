"""Gridinit export policy: align GRTresna solve data with GRTeclyn evolution spacing."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GridinitExportSpec:
    """Physical window resampled from a GRTresna solve into a ``.gridinit`` file.

    ``source_origin`` is the minimum corner of the export window in GRTresna
    native coordinates.  ``origin`` is the gridinit header origin so GRTeclyn's
    evolution ``target_center`` maps to the array centre when loaded.
    """

    nx: int
    ny: int
    nz: int
    lx: float
    ly: float
    lz: float
    target_center: tuple[float, float, float]
    source_origin: tuple[float, float, float] = (0.0, 0.0, 0.0)

    @property
    def dx_xyz(self) -> tuple[float, float, float]:
        return (self.lx / self.nx, self.ly / self.ny, self.lz / self.nz)

    @property
    def origin(self) -> tuple[float, float, float]:
        return (
            self.target_center[0] - 0.5 * self.lx,
            self.target_center[1] - 0.5 * self.ly,
            self.target_center[2] - 0.5 * self.lz,
        )

    def matches_evolution_spacing(
        self,
        *,
        l_full: float,
        n_full: int,
        rtol: float = 1.0e-9,
    ) -> bool:
        """True when per-axis gridinit spacing equals ``l_full / n_full``."""
        if n_full <= 0:
            return False
        evo_dx = float(l_full) / float(n_full)
        return all(abs(dx - evo_dx) <= rtol * max(evo_dx, 1.0e-15) for dx in self.dx_xyz)
