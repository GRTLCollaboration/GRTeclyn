"""Orchestrate GRTresna to produce constraint-satisfying initial data.

Given a parameter vector theta (BH masses, spins, scalar field config),
this module:
1. Writes a GRTresna params.txt
2. Runs GRTresna via MPI
3. Converts the Chombo HDF5 output to a .gridinit file for GRTeclyn

The resulting .gridinit can be loaded by GRTeclyn's
ExternalGridInitialData to start an evolution from a fully
constraint-satisfied initial state.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from ..core.config import REPO_ROOT
from .io import convert_chombo_to_gridinit

# Paths auto-detected or overridden by environment variables
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
    write_diagnostics: bool = False

    # Grid
    N: tuple[int, int, int] = (64, 64, 32)
    L: float = 128.0
    max_level: int = 3
    block_factor: int = 16
    max_grid_size: int = 16

    # BH params
    bh1_bare_mass: float = 1.0
    bh1_spin: tuple[float, float, float] = (0.0, 0.0, 0.5)
    bh1_momentum: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh1_offset: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_bare_mass: float = 0.0
    bh2_spin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_momentum: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bh2_offset: tuple[float, float, float] = (0.0, 0.0, 0.0)

    # Scalar field
    G_Newton: float = 1.0
    phi_0: float = 0.0
    dphi: float = 0.1
    dphi_length: float = 5.0
    pi_0: float = 0.0
    dpi: float = 0.1
    dpi_length: float = 5.0
    scalar_mass: float = 0.1

    # Conformal / K
    regularised_part_psi: float = 1.0
    sign_of_K: int = 1

    # Boundary conditions
    use_compact_Vi_ansatz: int = 1
    hi_boundary: tuple[int, int, int] = (0, 0, 0)
    lo_boundary: tuple[int, int, int] = (0, 0, 1)

    # Output target resolution for .gridinit
    gridinit_N: int = 64


def _fmt(v: Any) -> str:
    """Format a value for GRTresna params.txt."""
    if isinstance(v, (tuple, list)):
        return " ".join(str(x) for x in v)
    return str(v)


def write_grtresna_params(cfg: GRTresnaConfig, path: Path) -> None:
    """Write a GRTresna params.txt."""
    lines = [
        f"output_path = Outputs/",
        f"output_filename = InitialDataFinal.3d.hdf5",
        f"error_filename = Ham_and_Mom_errors",
        f"write_diagnostics = {1 if cfg.write_diagnostics else 0}",
        f"",
        f"N = {_fmt(cfg.N)}",
        f"L = {cfg.L}",
        f"max_level = {cfg.max_level}",
        f"block_factor = {cfg.block_factor}",
        f"max_grid_size = {cfg.max_grid_size}",
        f"regrid_radius = 10",
        f"",
        f"is_periodic = 0 0 0",
        f"use_compact_Vi_ansatz = {cfg.use_compact_Vi_ansatz}",
        f"hi_boundary = {_fmt(cfg.hi_boundary)}",
        f"lo_boundary = {_fmt(cfg.lo_boundary)}",
        f"",
        f"G_Newton = {cfg.G_Newton}",
        f"phi_0 = {cfg.phi_0}",
        f"dphi = {cfg.dphi}",
        f"dphi_length = {cfg.dphi_length}",
        f"pi_0 = {cfg.pi_0}",
        f"dpi = {cfg.dpi}",
        f"dpi_length = {cfg.dpi_length}",
        f"scalar_mass = {cfg.scalar_mass}",
        f"regularised_part_psi = {cfg.regularised_part_psi}",
        f"sign_of_K = {cfg.sign_of_K}",
        f"",
        f"bh1_bare_mass = {cfg.bh1_bare_mass}",
        f"bh1_spin = {_fmt(cfg.bh1_spin)}",
        f"bh1_momentum = {_fmt(cfg.bh1_momentum)}",
        f"bh1_offset = {_fmt(cfg.bh1_offset)}",
        f"bh2_bare_mass = {cfg.bh2_bare_mass}",
        f"bh2_spin = {_fmt(cfg.bh2_spin)}",
        f"bh2_momentum = {_fmt(cfg.bh2_momentum)}",
        f"bh2_offset = {_fmt(cfg.bh2_offset)}",
        f"",
        f"max_NL_iterations = {cfg.max_NL_iterations}",
        f"deactivate_zero_mode = 0",
    ]
    path.write_text("\n".join(lines) + "\n")


def _find_executable(cfg: GRTresnaConfig) -> Path:
    """Locate the compiled GRTresna executable."""
    example_dir = cfg.grtresna_root / "Examples" / cfg.example
    candidates = list(example_dir.glob("Main_*3d.*.ex"))
    if not candidates:
        raise FileNotFoundError(
            f"No GRTresna executable found in {example_dir}. "
            "Build it first with `make all`."
        )
    return candidates[0]


def solve(
    cfg: GRTresnaConfig,
    work_dir: Path | None = None,
    gridinit_path: Path | None = None,
) -> Path:
    """Run GRTresna and return the path to the .gridinit file.

    Parameters
    ----------
    cfg : solver configuration
    work_dir : directory to run in (created if None)
    gridinit_path : where to write the .gridinit (defaults to work_dir)

    Returns
    -------
    Path to the .gridinit file
    """
    if work_dir is None:
        work_dir = Path(tempfile.mkdtemp(prefix="grtresna_"))
    work_dir.mkdir(parents=True, exist_ok=True)

    exe = _find_executable(cfg)

    outputs_dir = work_dir / "Outputs"
    outputs_dir.mkdir(exist_ok=True)
    (work_dir / "pout").mkdir(exist_ok=True)

    params_path = work_dir / "params.txt"
    write_grtresna_params(cfg, params_path)

    oversubscribe = "--oversubscribe" if cfg.mpi_ranks > 1 else ""
    cmd = (
        f"mpirun {oversubscribe} -np {cfg.mpi_ranks} "
        f"{exe} {params_path}"
    ).split()

    result = subprocess.run(
        cmd,
        cwd=str(work_dir),
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"GRTresna failed (rc={result.returncode}):\n"
            f"stdout: {result.stdout[-2000:]}\n"
            f"stderr: {result.stderr[-2000:]}"
        )

    chombo_hdf5 = outputs_dir / "InitialDataFinal.3d.hdf5"
    if not chombo_hdf5.exists():
        raise FileNotFoundError(
            f"GRTresna did not produce {chombo_hdf5}. "
            f"Check {work_dir / 'Ham_and_Mom_errors.txt'}."
        )

    if gridinit_path is None:
        gridinit_path = work_dir / "initial_data.gridinit"

    convert_chombo_to_gridinit(
        chombo_hdf5,
        gridinit_path,
        N=cfg.gridinit_N,
        L=cfg.L,
    )
    return gridinit_path
