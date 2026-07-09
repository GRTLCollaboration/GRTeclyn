"""Quantitative geometry-mismatch fitness for the iterative projection loop.

Computes an L2 distance between target analytic profiles (from a GeometryMotif)
and the profiles actually realized by a GRTresna-solved .gridinit, plus exotic-
mass and constraint-convergence penalties.  Lower fitness is better (CMA-ES
convention).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from ..grtresna.io import read_gridinit
from ..grtresna.solver.convergence import parse_convergence
from ..initial_data.motif import GeometryMotif
from ..metrics.probes.ftl.analytic import _axis_profiles
from ..metrics.probes.ftl.solved import build_xz_slice_from_gridinit

# Default weight schedule (all positive; fitness = sum of weighted terms).
W_CHI: float = 1.0
W_BETA: float = 1.0
W_EXOTIC: float = 0.3
W_CONVERGENCE: float = 0.5

# Gate fitness returned when a solve fails entirely (analogous to
# cmaes_adapter.DEFAULT_GATE_FITNESS).
GATE_FITNESS: float = 100.0

# Profile comparison resolution (number of 1-D sample points along axis).
N_PROFILE_POINTS: int = 256


@dataclass(frozen=True)
class MismatchReport:
    """Per-eval fitness breakdown for the iteration loop."""

    fitness: float
    chi_l2: float
    beta_l2: float
    exotic_penalty: float
    convergence_penalty: float
    solve_failed: bool
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _target_profiles(
    motif: GeometryMotif,
    *,
    n_points: int = N_PROFILE_POINTS,
    L: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Extract 1-D chi(x) and beta_x(x) target profiles from motif overrides."""
    overrides = motif.overrides
    eff_L = L if L is not None else float(overrides.get("recipe_basis_radius_max", 8.0)) * 2.0
    x, _r, chi, _alpha, beta_x = _axis_profiles(overrides, L=eff_L, n_points=n_points)
    return x, chi, beta_x


def _solved_midline_profiles(
    gridinit_path: Path,
    *,
    n_points: int = N_PROFILE_POINTS,
    L: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]] | None:
    """Extract 1-D midline chi(x) and beta_x(x) from a solved .gridinit."""
    try:
        grid = read_gridinit(gridinit_path)
        alpha, beta, gamma, spacing = build_xz_slice_from_gridinit(grid, n=n_points, L=L)
    except (OSError, ValueError, KeyError):
        return None

    mid_z = alpha.shape[1] // 2
    # chi = 1/gamma_{xx} along the x-axis midline (z = center)
    chi_solved = 1.0 / np.clip(gamma[:, mid_z, 0, 0], 1.0e-10, None)
    beta_x_solved = beta[:, mid_z, 0]
    # Construct matching x-axis coordinates
    n = alpha.shape[0]
    half = 0.5 * n * spacing[0]
    x = np.linspace(-half, half, n)
    return x, chi_solved, beta_x_solved


def _l2_norm(a: NDArray[np.float64], b: NDArray[np.float64]) -> float:
    """Normalized L2 distance (RMS)."""
    diff = a - b
    return float(np.sqrt(np.mean(diff * diff)))


def _resample(
    x_src: NDArray[np.float64],
    y_src: NDArray[np.float64],
    x_dst: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Linearly interpolate y_src(x_src) onto x_dst grid."""
    return np.interp(x_dst, x_src, y_src)


def _exotic_mass_proxy(lumps: list[dict[str, Any]]) -> float:
    """Simple proxy for total exotic matter: sum of |amp| for exotic lumps."""
    total = 0.0
    for lump in lumps:
        if int(lump.get("exotic", 0)):
            total += abs(float(lump.get("amp", 0.0)))
    return total


def compute_mismatch(
    motif: GeometryMotif,
    gridinit_path: str | Path,
    *,
    lumps: list[dict[str, Any]] | None = None,
    grtresna_work_dir: str | Path | None = None,
    L: float | None = None,
    n_points: int = N_PROFILE_POINTS,
    w_chi: float = W_CHI,
    w_beta: float = W_BETA,
    w_exotic: float = W_EXOTIC,
    w_convergence: float = W_CONVERGENCE,
) -> MismatchReport:
    """Compute geometry-mismatch fitness between motif target and solved gridinit.

    Parameters
    ----------
    motif : GeometryMotif
        The geometry target (carries overrides for analytic profile generation).
    gridinit_path : path
        The solved .gridinit from GRTresna.
    lumps : list of dicts, optional
        Lump parameters (for exotic penalty). If None, exotic penalty is 0.
    grtresna_work_dir : path, optional
        GRTresna work directory (for convergence penalty via parse_convergence).
    L : float, optional
        Integration half-width override for profile extraction.
    """
    gridinit_path = Path(gridinit_path)
    notes: list[str] = []

    # Target profiles from motif
    x_target, chi_target, beta_target = _target_profiles(motif, n_points=n_points, L=L)

    # Solved profiles from gridinit
    solved = _solved_midline_profiles(gridinit_path, n_points=n_points, L=L)
    if solved is None:
        return MismatchReport(
            fitness=GATE_FITNESS,
            chi_l2=0.0,
            beta_l2=0.0,
            exotic_penalty=0.0,
            convergence_penalty=0.0,
            solve_failed=True,
            notes=("failed to read solved gridinit profiles",),
        )

    x_solved, chi_solved, beta_solved = solved

    # Resample solved profiles onto the target x-grid for L2 comparison
    chi_on_target = _resample(x_solved, chi_solved, x_target)
    beta_on_target = _resample(x_solved, beta_solved, x_target)

    chi_l2 = _l2_norm(chi_on_target, chi_target)
    beta_l2 = _l2_norm(beta_on_target, beta_target)

    # Exotic mass penalty
    exotic_penalty = 0.0
    if lumps is not None:
        exotic_penalty = _exotic_mass_proxy(lumps)

    # Convergence penalty from Ham/Mom residuals (parse_convergence returns
    # keys: iteration, ham_pct, mom_pct — percentages of the norm).
    convergence_penalty = 0.0
    if grtresna_work_dir is not None:
        conv = parse_convergence(Path(grtresna_work_dir))
        if conv is not None:
            ham = conv.get("ham_pct", 0.0)
            mom = conv.get("mom_pct", 0.0)
            # Penalty grows linearly above 0.1% threshold
            convergence_penalty = max(0.0, ham - 0.1) + max(0.0, mom - 0.1)
            if ham > 5.0 or mom > 5.0:
                notes.append(f"high constraint residual: Ham={ham:.2f}% Mom={mom:.2f}%")
        else:
            notes.append("no convergence data available")

    fitness = (
        w_chi * chi_l2
        + w_beta * beta_l2
        + w_exotic * exotic_penalty
        + w_convergence * convergence_penalty
    )

    return MismatchReport(
        fitness=fitness,
        chi_l2=chi_l2,
        beta_l2=beta_l2,
        exotic_penalty=exotic_penalty,
        convergence_penalty=convergence_penalty,
        solve_failed=False,
        notes=tuple(notes),
    )
