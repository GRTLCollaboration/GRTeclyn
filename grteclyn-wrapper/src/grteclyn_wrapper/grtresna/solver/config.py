"""GRTresna solver configuration and exotic-matter policy."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

from ...core.config import REPO_ROOT
from ..domain.export import GridinitExportSpec

_DEFAULT_GRTRESNA_ROOT = Path(
    os.environ.get("GRTRESNA_ROOT", str(REPO_ROOT.parent / "GRTresna"))
)
_DEFAULT_EXAMPLE = "ScalarFieldBH"


@dataclass
class GRTresnaConfig:
    """All knobs for a GRTresna solve."""

    grtresna_root: Path = _DEFAULT_GRTRESNA_ROOT
    example: str = _DEFAULT_EXAMPLE
    mpi_ranks: int = 8
    max_NL_iterations: int = 50
    nl_exit_tolerance: float = 1.0
    nl_stall_tolerance: float = 0.02
    write_diagnostics: bool = False
    timeout: int = 3600
    mpirun: str | None = None

    N: tuple[int, int, int] = (64, 64, 32)
    L: float = 128.0
    max_level: int = 3
    block_factor: int = 16
    max_grid_size: int = 16
    refine_threshold: float = 0.5
    regrid_radius: float = 10.0
    coefficient_average_type: str = "harmonic"

    bh1_bare_mass: float = 1.0
    bh1_spin: tuple[float, float, float] = (0.0, 0.0, 0.5)
    bh1_momentum: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh1_offset: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_bare_mass: float = 0.0
    bh2_spin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_momentum: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_offset: tuple[float, float, float] = (0.0, 0.0, 0.0)

    G_Newton: float = 1.0
    phi_0: float = 0.0
    dphi: float = 0.1
    dphi_length: float = 5.0
    pi_0: float = 0.0
    dpi: float = 0.1
    dpi_length: float = 5.0
    scalar_mass: float = 0.1
    scalar_lambda: float = 0.0
    scalar_mu: float = 0.0

    matter_model: str = ""

    bs_phi_c: float = 0.08
    bs_profile_width: float = 8.0
    bs_omega: float = 0.0
    scalar_sign: int = 1

    lump_amp: float = 0.0
    lump_width: float = 5.0
    lump_center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    lump_velocity: tuple[float, float, float] = (0.0, 0.0, 0.0)
    lump_omega: float = 0.0
    lump_mode: int = 0
    lump_exotic: int = 0

    lumps: list[dict] = field(default_factory=list)
    shift_seed: float = 0.0

    regularised_part_psi: float = 1.0
    sign_of_K: int = 1

    maximal_slicing: bool = False
    psi_relaxation: float = 1.0
    psi_floor: float = -1.0
    maximal_jacobian_cap: float = -1.0

    use_compact_Vi_ansatz: int = 1
    hi_boundary: tuple[int, int, int] = (0, 0, 0)
    lo_boundary: tuple[int, int, int] = (0, 0, 1)

    gridinit_nx: int = 64
    gridinit_ny: int = 64
    gridinit_nz: int = 64
    gridinit_workers: int = 0

    target_center: tuple[float, float, float] = (32.0, 32.0, 0.0)
    gridinit_export: GridinitExportSpec | None = None

    cleanup: bool = True
    cleanup_workdir: bool = False


def apply_exotic_safe_solver(cfg: GRTresnaConfig) -> GRTresnaConfig:
    """Switch ``cfg`` to the K=0 maximal-slicing York/Lichnerowicz path in place."""
    if not cfg.maximal_slicing:
        cfg.maximal_slicing = True
    if cfg.psi_relaxation == 1.0:
        cfg.psi_relaxation = 0.6
    if cfg.psi_floor <= 0.0:
        cfg.psi_floor = 0.1
    if cfg.maximal_jacobian_cap <= 0.0:
        cfg.maximal_jacobian_cap = 25.0
    if cfg.coefficient_average_type == "harmonic":
        cfg.coefficient_average_type = "arithmetic"
    return cfg


def config_has_exotic_lump(cfg: GRTresnaConfig) -> bool:
    """True when any configured matter is flagged exotic (negative energy)."""
    if cfg.matter_model == "grtresna_complex_scalar" and cfg.scalar_sign < 0:
        return True
    if cfg.lumps:
        return any(int(lump.get("exotic", 0)) for lump in cfg.lumps)
    return bool(cfg.lump_exotic)
