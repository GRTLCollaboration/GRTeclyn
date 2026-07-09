"""Quantitative geometry-mismatch fitness for the iterative projection loop.

Computes an L2 distance between target analytic profiles (from a GeometryMotif)
and the profiles actually realized by a GRTresna-solved .gridinit, plus exotic-
mass and constraint-convergence penalties.  Lower fitness is better (CMA-ES
convention).

Fitness uses **two-phase** shaping:
  - Phase 1 (infeasible): Ham > ``FEASIBILITY_THRESHOLD`` — fitness is dominated
    by ``tanh``-saturated convergence penalty so CMA-ES finds the convergence
    basin first without the geometry signal being drowned.
  - Phase 2 (feasible): Ham ≤ threshold — geometry mismatch (chi + beta L2) is
    the primary signal; convergence and exotic penalties are small corrections.

Phase 4 extensions:
  - **2D slice mismatch**: chi and beta are compared on the full xz-plane slice
    (y = domain centre), not just the 1D midline.  This captures angular
    structure that the radial profile misses (e.g. ring-distributed lumps).
  - **K_ij matching**: the traceless extrinsic curvature A_ij and trace K are
    included in the mismatch when available in the gridinit.
  - **Feasibility pre-check**: ``feasibility_precheck`` estimates whether the
    target geometry is likely to be feasible before running GRTresna, based on
    the required matter density magnitude.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from ..grtresna.io import read_gridinit
from ..grtresna.solver.convergence import parse_convergence
from ..initial_data.motif import GeometryMotif, alcubierre_shape_function, _alcubierre_eulerian_rho
from ..metrics.probes.ftl.analytic import _axis_profiles
from ..metrics.probes.ftl.solved import build_xz_slice_from_gridinit

# Default weight schedule (all positive; fitness = sum of weighted terms).
W_CHI: float = 1.0
W_BETA: float = 1.0
W_EXOTIC: float = 0.3
W_CONVERGENCE: float = 0.5

# Phase 4: weights for 2D and K_ij terms.
W_CHI_2D: float = 2.0       # 2D chi slice (captures angular structure)
W_BETA_2D: float = 1.5      # 2D beta slice
W_KIJ: float = 0.5          # Extrinsic curvature trace K
W_AIJ: float = 0.3          # Traceless extrinsic curvature A_ij

# Gate fitness returned when a solve fails entirely (analogous to
# cmaes_adapter.DEFAULT_GATE_FITNESS).
GATE_FITNESS: float = 100.0

# Profile comparison resolution (number of 1-D sample points along axis).
N_PROFILE_POINTS: int = 256

# 2D slice resolution (n × n grid for the xz-plane).
N_SLICE_POINTS: int = 64

# Two-phase feasibility shaping.
FEASIBILITY_THRESHOLD: float = 5.0   # Ham% below which we switch to geometry-mode
CONV_TANH_SCALE: float = 20.0        # tanh(ham_pct / scale) — saturates ~1.0 at ham≈60%
FEASIBILITY_WEIGHT: float = 50.0     # max contribution of the feasibility term

# Phase 4c: feasibility pre-check thresholds.
# Empirically, GRTresna converges well when |rho_req| < ~0.5 and struggles
# when |rho_req| > ~2.0.  These are heuristic, derived from test runs.
FEASIBILITY_RHO_SAFE: float = 0.5
FEASIBILITY_RHO_HARD: float = 2.0


@dataclass(frozen=True)
class MismatchReport:
    """Per-eval fitness breakdown for the iteration loop."""

    fitness: float
    chi_l2: float
    beta_l2: float
    exotic_penalty: float
    convergence_penalty: float
    solve_failed: bool
    # Phase 4: additional mismatch terms.
    chi_2d_l2: float = 0.0
    beta_2d_l2: float = 0.0
    kij_l2: float = 0.0
    aij_l2: float = 0.0
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class FeasibilityEstimate:
    """Phase 4c: pre-solve feasibility estimate from matter requirements."""

    feasible: bool               # True if likely to converge (Ham < threshold)
    rho_peak: float              # Peak |rho_req| from motif
    rho_integral: float          # Integrated |rho_req| (proxy for total matter)
    risk_level: str              # "safe", "marginal", or "hard"
    notes: tuple[str, ...] = ()


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


def _target_2d_slice(
    motif: GeometryMotif,
    *,
    n: int = N_SLICE_POINTS,
    L: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Evaluate target chi and beta on a 2D xz-slice through the transport axis.

    The RadialRecipe basis is spherically symmetric (depends on r = |x| only),
    so the 2D target is obtained by evaluating the basis on r = sqrt(x² + z²).

    Returns (x_axis, z_axis, chi_2d, beta_x_2d, beta_z_2d).
    For a spherically symmetric target, beta_z = 0 everywhere.
    """
    from ..initial_data.constrained_recipe import RecipeBasis

    overrides = motif.overrides
    eff_L = L if L is not None else float(overrides.get("recipe_basis_radius_max", 8.0)) * 2.0

    num_bases = int(overrides.get("recipe_num_bases", 4))
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=float(overrides.get("recipe_basis_width", 1.0)),
        basis_radius_max=float(overrides.get("recipe_basis_radius_max", 8.0)),
    )

    chi_coeffs = [float(overrides.get(f"recipe_chi_coeff_{n}", 0.0)) for n in range(num_bases)]
    beta_coeffs = [float(overrides.get(f"recipe_beta_coeff_{n}", 0.0)) for n in range(num_bases)]
    chi_asym = float(overrides.get("recipe_chi_asymptotic", 1.0))
    beta_asym = float(overrides.get("recipe_beta_asymptotic", 0.0))

    x_axis = np.linspace(-eff_L, eff_L, n)
    z_axis = np.linspace(-eff_L, eff_L, n)
    X, Z = np.meshgrid(x_axis, z_axis, indexing="ij")
    R = np.sqrt(X**2 + Z**2)

    chi_2d = basis.evaluate(R.ravel(), chi_asym, chi_coeffs).reshape(n, n)
    chi_2d = np.clip(chi_2d, 1.0e-10, None)

    # For spherically symmetric recipe, beta is radial.  Along the x-axis
    # (transport axis), beta_x = beta(r).  In 2D, the radial direction is
    # (x, z)/r, so beta_x = beta(r) * x/r and beta_z = beta(r) * z/r.
    beta_r = basis.evaluate(R.ravel(), beta_asym, beta_coeffs).reshape(n, n)
    # Avoid division by zero at r=0
    r_safe = np.clip(R, 1.0e-10, None)
    beta_x_2d = beta_r * X / r_safe
    beta_z_2d = beta_r * Z / r_safe

    return x_axis, z_axis, chi_2d, beta_x_2d, beta_z_2d


def _is_warp_motif(motif: GeometryMotif) -> bool:
    """Detect an Alcubierre warp-drive motif (toroidal template + flat chi)."""
    return (
        motif.momentum_target.template == "toroidal"
        and motif.episode_path.startswith("alcubierre:")
    )


def _warp_params_from_motif(motif: GeometryMotif) -> dict[str, float]:
    """Extract (velocity, bubble_radius, sigma) from an Alcubierre motif."""
    parts = motif.episode_path.split(":")
    params: dict[str, float] = {}
    for part in parts[1:]:
        if "=" in part:
            key, val = part.split("=", 1)
            try:
                params[key.strip()] = float(val)
            except ValueError:
                pass
    return params


def _target_2d_slice_warp(
    motif: GeometryMotif,
    *,
    n: int = N_SLICE_POINTS,
    L: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Evaluate Alcubierre target chi, beta, K, A_ij on a 2D xz-slice.

    For an Alcubierre bubble:
      - chi = 1 everywhere (flat spatial metric)
      - beta_x = -v * f(r_s), beta_z = 0 (axial shift only)
      - K (trace) = 0 (maximal slicing is NOT assumed; Alcubierre has K=0
        because the spatial metric is flat and the trace of K_ij for the
        standard Alcubierre metric is zero)
      - A_ij (traceless extrinsic curvature) is non-zero: it comes from the
        shift gradient.  We use |A_ij|^2/2 as a scalar proxy, computed from
        the analytic K_ij of the Alcubierre metric.

    Returns (x_axis, z_axis, chi_2d, beta_x_2d, beta_z_2d, K_2d, A_trace_2d).
    """
    params = _warp_params_from_motif(motif)
    velocity = params.get("v", 0.5)
    bubble_radius = params.get("R", 2.0)
    sigma = params.get("sigma", 2.0)

    overrides = motif.overrides
    eff_L = L if L is not None else float(overrides.get("recipe_basis_radius_max", 8.0)) * 2.0

    x_axis = np.linspace(-eff_L, eff_L, n)
    z_axis = np.linspace(-eff_L, eff_L, n)
    X, Z = np.meshgrid(x_axis, z_axis, indexing="ij")
    r_s = np.sqrt(X**2 + Z**2)

    # Flat chi
    chi_2d = np.ones((n, n), dtype=np.float64)

    # Axial shift: beta_x = -v * f(r_s), beta_z = 0
    f = alcubierre_shape_function(r_s, bubble_radius=bubble_radius, sigma=sigma)
    beta_x_2d = -velocity * f
    beta_z_2d = np.zeros((n, n), dtype=np.float64)

    # K (trace) = 0 for Alcubierre (flat spatial metric, maximal slice)
    K_2d = np.zeros((n, n), dtype=np.float64)

    # A_ij proxy from the shift: the extrinsic curvature for a shift-only
    # metric (alpha=1, gamma=delta) is K_ij = (1/2)(d_i beta_j + d_j beta_i).
    # For beta = (-v f(r_s), 0, 0), the dominant components are:
    #   K_xx = d_x beta_x = -v * f' * x/r_s
    #   K_xz = (1/2) d_z beta_x = -v * f' * z / (2 r_s)
    #   K_zz = 0 (beta_z = 0)
    # |A_ij|^2/2 ~ |K_ij|^2/2 (since K=0, A_ij = K_ij).
    r_safe = np.clip(r_s, 1.0e-10, None)
    # df/dr_s via finite differences on a 1D radial grid
    r_1d = np.linspace(max(eff_L / n, 1.0e-6), eff_L * math.sqrt(2.0), n)
    f_1d = alcubierre_shape_function(r_1d, bubble_radius=bubble_radius, sigma=sigma)
    dr = r_1d[1] - r_1d[0]
    df_dr_1d = np.gradient(f_1d, dr, edge_order=2)
    df_dr = np.interp(r_s.ravel(), r_1d, df_dr_1d).reshape(n, n)

    K_xx = -velocity * df_dr * X / r_safe
    K_xz = -velocity * df_dr * Z / (2.0 * r_safe)
    K_zz = np.zeros((n, n), dtype=np.float64)
    # |K_ij|^2 = K_xx^2 + 2*K_xz^2 + K_zz^2 (symmetric tensor)
    k_sq = K_xx**2 + 2.0 * K_xz**2 + K_zz**2
    A_trace_2d = np.sqrt(k_sq / 2.0)

    return x_axis, z_axis, chi_2d, beta_x_2d, beta_z_2d, K_2d, A_trace_2d


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


def _solved_2d_slice(
    gridinit_path: Path,
    *,
    n: int = N_SLICE_POINTS,
    L: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]] | None:
    """Extract 2D xz-slice chi, beta, K, A_ij from a solved .gridinit.

    Returns (x_axis, z_axis, chi_2d, beta_x_2d, beta_z_2d, K_2d, A_trace_2d)
    or None if the gridinit can't be read.  K and A_trace are 0.0 if those
    components are not present in the gridinit.
    """
    try:
        grid = read_gridinit(gridinit_path)
        alpha, beta, gamma, spacing = build_xz_slice_from_gridinit(grid, n=n, L=L)
    except (OSError, ValueError, KeyError):
        return None

    # chi = 1 / gamma_xx (conformal factor)
    chi_2d = 1.0 / np.clip(gamma[:, :, 0, 0], 1.0e-10, None)
    beta_x_2d = beta[:, :, 0]
    beta_z_2d = beta[:, :, 1]

    # Extract K and A_ij from the raw grid data
    data = grid.data
    names = list(grid.comp_names)
    nz, ny, nx, _ = data.shape
    dx, dy, dz = grid.dx_xyz
    ox, oy, oz = grid.origin

    cx = ox + 0.5 * nx * dx
    cy = oy + 0.5 * ny * dy
    cz = oz + 0.5 * nz * dz

    half_x = 0.49 * nx * dx
    half_z = 0.49 * nz * dz
    if L is not None:
        half_x = min(half_x, float(L))
        half_z = min(half_z, float(L))

    x_axis = np.linspace(cx - half_x, cx + half_x, n)
    z_axis = np.linspace(cz - half_z, cz + half_z, n)

    # Nearest-neighbour indices (same logic as build_xz_slice_from_gridinit)
    ii = np.clip(np.round((x_axis - ox) / dx - 0.5).astype(np.int64), 0, nx - 1)
    kk = np.clip(np.round((z_axis - oz) / dz - 0.5).astype(np.int64), 0, nz - 1)
    jj = int(np.clip(round((cy - oy) / dy - 0.5), 0, ny - 1))
    ki = kk[None, :].repeat(n, axis=0)
    xi = ii[:, None].repeat(n, axis=1)

    def _comp(name: str) -> NDArray[np.float64] | None:
        try:
            idx = names.index(name)
        except ValueError:
            return None
        return data[ki, jj, xi, idx]

    # K (trace of extrinsic curvature)
    K_2d = np.zeros((n, n), dtype=np.float64)
    K_arr = _comp("K")
    if K_arr is not None:
        K_2d = K_arr

    # A_ij traceless part — use |A_ij|^2 / 2 as a scalar proxy
    A_trace_2d = np.zeros((n, n), dtype=np.float64)
    a_comps = {}
    for name in ("A11", "A12", "A13", "A22", "A23", "A33"):
        arr = _comp(name)
        if arr is not None:
            a_comps[name] = arr
    if a_comps:
        a_sq = np.zeros((n, n), dtype=np.float64)
        for name, arr in a_comps.items():
            a_sq += arr * arr
        A_trace_2d = np.sqrt(a_sq / 2.0)  # |A_ij|^2 / 2

    return x_axis, z_axis, chi_2d, beta_x_2d, beta_z_2d, K_2d, A_trace_2d


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


def _resample_2d(
    x_src: NDArray[np.float64],
    z_src: NDArray[np.float64],
    field_src: NDArray[np.float64],
    x_dst: NDArray[np.float64],
    z_dst: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Bilinearly interpolate a 2D field from (x_src, z_src) to (x_dst, z_dst)."""
    from scipy.interpolate import RegularGridInterpolator

    interp = RegularGridInterpolator((x_src, z_src), field_src, bounds_error=False, fill_value=None)
    X, Z = np.meshgrid(x_dst, z_dst, indexing="ij")
    return interp(np.stack([X.ravel(), Z.ravel()], axis=-1)).reshape(len(x_dst), len(z_dst))


def _exotic_mass_proxy(lumps: list[dict[str, Any]]) -> float:
    """Simple proxy for total exotic matter: sum of |amp| for exotic lumps."""
    total = 0.0
    for lump in lumps:
        if int(lump.get("exotic", 0)):
            total += abs(float(lump.get("amp", 0.0)))
    return total


def feasibility_precheck(motif: GeometryMotif) -> FeasibilityEstimate:
    """Phase 4c: estimate whether the target geometry is likely feasible.

    Uses the peak |rho_req| from the motif's support regions and the
    ``min_rho_required`` from the constraint file to classify the target
    as safe, marginal, or hard — *before* running GRTresna.

    This is a heuristic: the actual convergence depends on the matter ansatz,
    grid resolution, and solver settings.  But it provides a quick filter to
    skip obviously infeasible targets and warn on marginal ones.
    """
    notes: list[str] = []

    # Peak |rho| from support regions
    rho_peak = 0.0
    for region in motif.support_regions:
        rho_peak = max(rho_peak, abs(region.peak_rho))

    # Also check min_rho_required from constraint file
    if motif.min_rho_required is not None:
        rho_peak = max(rho_peak, abs(motif.min_rho_required))

    # Integral proxy: sum of |peak_rho| * width^3 (very rough)
    rho_integral = 0.0
    for region in motif.support_regions:
        rho_integral += abs(region.peak_rho) * region.width**3

    if rho_peak < FEASIBILITY_RHO_SAFE:
        risk = "safe"
        feasible = True
    elif rho_peak < FEASIBILITY_RHO_HARD:
        risk = "marginal"
        feasible = True  # still try, but expect slower convergence
        notes.append(f"marginal rho_peak={rho_peak:.4f} — convergence may be slow")
    else:
        risk = "hard"
        feasible = False
        notes.append(f"hard rho_peak={rho_peak:.4f} — likely infeasible for Gaussian ansatz")

    if motif.exotic_needed:
        notes.append("exotic matter required — convergence basin is narrower")

    return FeasibilityEstimate(
        feasible=feasible,
        rho_peak=rho_peak,
        rho_integral=rho_integral,
        risk_level=risk,
        notes=tuple(notes),
    )


def compute_mismatch(
    motif: GeometryMotif,
    gridinit_path: str | Path,
    *,
    lumps: list[dict[str, Any]] | None = None,
    grtresna_work_dir: str | Path | None = None,
    L: float | None = None,
    n_points: int = N_PROFILE_POINTS,
    n_slice: int = N_SLICE_POINTS,
    w_chi: float = W_CHI,
    w_beta: float = W_BETA,
    w_exotic: float = W_EXOTIC,
    w_convergence: float = W_CONVERGENCE,
    w_chi_2d: float = W_CHI_2D,
    w_beta_2d: float = W_BETA_2D,
    w_kij: float = W_KIJ,
    w_aij: float = W_AIJ,
    use_2d: bool = True,
    use_kij: bool = True,
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
    n_points : int
        Number of 1D profile sample points.
    n_slice : int
        Number of points per axis for 2D slice comparison.
    use_2d : bool
        If True (default), include 2D slice mismatch in the fitness.
    use_kij : bool
        If True (default), include K_ij mismatch in the fitness.
    """
    gridinit_path = Path(gridinit_path)
    notes: list[str] = []

    # --- 1D midline profiles (always computed — cheap and used as fallback) ---
    x_target, chi_target, beta_target = _target_profiles(motif, n_points=n_points, L=L)

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

    # --- Phase 4a: 2D slice mismatch ---
    chi_2d_l2 = 0.0
    beta_2d_l2 = 0.0
    if use_2d:
        try:
            is_warp = _is_warp_motif(motif)
            if is_warp:
                warp_2d = _target_2d_slice_warp(motif, n=n_slice, L=L)
                x_t, z_t, chi_2d_target, beta_x_2d_target, beta_z_2d_target = warp_2d[:5]
                warp_K_target, warp_A_target = warp_2d[5], warp_2d[6]
            else:
                x_t, z_t, chi_2d_target, beta_x_2d_target, beta_z_2d_target = _target_2d_slice(
                    motif, n=n_slice, L=L,
                )
                warp_K_target = None
                warp_A_target = None
            solved_2d = _solved_2d_slice(gridinit_path, n=n_slice, L=L)
            if solved_2d is not None:
                x_s, z_s, chi_2d_solved, beta_x_2d_solved, beta_z_2d_solved, _, _ = solved_2d
                # Resample solved 2D onto target grid
                chi_2d_resampled = _resample_2d(x_s, z_s, chi_2d_solved, x_t, z_t)
                beta_x_2d_resampled = _resample_2d(x_s, z_s, beta_x_2d_solved, x_t, z_t)
                beta_z_2d_resampled = _resample_2d(x_s, z_s, beta_z_2d_solved, x_t, z_t)

                chi_2d_l2 = _l2_norm(chi_2d_resampled, chi_2d_target)
                # Combine beta_x and beta_z into a single beta 2D mismatch
                beta_x_2d_l2 = _l2_norm(beta_x_2d_resampled, beta_x_2d_target)
                beta_z_2d_l2 = _l2_norm(beta_z_2d_resampled, beta_z_2d_target)
                beta_2d_l2 = math.sqrt(beta_x_2d_l2**2 + beta_z_2d_l2**2)
        except Exception as exc:
            notes.append(f"2D slice mismatch failed: {exc}")

    # --- Phase 4b: K_ij matching ---
    kij_l2 = 0.0
    aij_l2 = 0.0
    if use_kij:
        try:
            is_warp = _is_warp_motif(motif)
            if is_warp and warp_K_target is not None:
                # Warp target: K and A_ij come from the Alcubierre metric
                # (non-zero A_ij from the shift gradient, K=0).
                target_K = warp_K_target
                target_A_trace = warp_A_target
            else:
                # Standard recipe target: K=0 (maximal slicing), A_ij=0
                # (conformally flat).
                target_K = 0.0
                target_A_trace = 0.0

            solved_2d = _solved_2d_slice(gridinit_path, n=n_slice, L=L)
            if solved_2d is not None:
                _, _, _, _, _, K_2d_solved, A_trace_2d_solved = solved_2d
                # Resample solved K and A onto the target grid if warp
                if is_warp and warp_K_target is not None:
                    K_2d_solved = _resample_2d(x_s, z_s, K_2d_solved, x_t, z_t)
                    A_trace_2d_solved = _resample_2d(x_s, z_s, A_trace_2d_solved, x_t, z_t)
                # L2 of K against target
                kij_l2 = float(np.sqrt(np.mean((K_2d_solved - target_K)**2)))
                # L2 of |A_ij| against target
                aij_l2 = float(np.sqrt(np.mean((A_trace_2d_solved - target_A_trace)**2)))
        except Exception as exc:
            notes.append(f"K_ij matching failed: {exc}")

    # Exotic mass penalty
    exotic_penalty = 0.0
    if lumps is not None:
        exotic_penalty = _exotic_mass_proxy(lumps)

    # Convergence penalty from Ham/Mom residuals (parse_convergence returns
    # keys: iteration, ham_pct, mom_pct — percentages of the norm).
    convergence_penalty = 0.0
    ham_pct = 0.0
    if grtresna_work_dir is not None:
        conv = parse_convergence(Path(grtresna_work_dir))
        if conv is not None:
            ham_pct = conv.get("ham_pct", 0.0)
            mom_pct = conv.get("mom_pct", 0.0)
            # tanh-saturated penalty: provides gradient toward feasibility
            # without drowning geometry signal once solve is "close".
            convergence_penalty = (
                math.tanh(max(0.0, ham_pct) / CONV_TANH_SCALE)
                + math.tanh(max(0.0, mom_pct) / CONV_TANH_SCALE)
            )
            if ham_pct > 5.0 or mom_pct > 5.0:
                notes.append(f"high constraint residual: Ham={ham_pct:.2f}% Mom={mom_pct:.2f}%")
        else:
            notes.append("no convergence data available")
            convergence_penalty = 1.0  # assume bad if no data

    # Two-phase fitness shaping:
    #  - Infeasible (Ham > threshold): feasibility dominates so CMA-ES finds the
    #    convergence basin first.  Geometry terms still contribute a small gradient.
    #  - Feasible (Ham ≤ threshold): geometry mismatch is the primary signal;
    #    convergence and exotic are small corrections.
    geometry_fitness = (
        w_chi * chi_l2
        + w_beta * beta_l2
        + w_chi_2d * chi_2d_l2
        + w_beta_2d * beta_2d_l2
        + w_kij * kij_l2
        + w_aij * aij_l2
        + w_exotic * exotic_penalty
    )
    feasibility_fitness = FEASIBILITY_WEIGHT * convergence_penalty

    if ham_pct > FEASIBILITY_THRESHOLD:
        # Phase 1: feasibility-first (geometry contributes 10% for gradient)
        fitness = feasibility_fitness + 0.1 * geometry_fitness
    else:
        # Phase 2: geometry-first (convergence is a small correction)
        fitness = geometry_fitness + w_convergence * convergence_penalty

    return MismatchReport(
        fitness=fitness,
        chi_l2=chi_l2,
        beta_l2=beta_l2,
        exotic_penalty=exotic_penalty,
        convergence_penalty=convergence_penalty,
        solve_failed=False,
        chi_2d_l2=chi_2d_l2,
        beta_2d_l2=beta_2d_l2,
        kij_l2=kij_l2,
        aij_l2=aij_l2,
        notes=tuple(notes),
    )


@dataclass(frozen=True)
class WarpCrossCheck:
    """Cross-check of reconstructed rho against the analytic Alcubierre T_{mu nu}.

    Compares the Eulerian energy density from the motif's rho_req (derived
    from K_ij/shift) against the exact ``T = G/8 pi`` from the warp-factory
    analytic metric, and reports the energy-condition diagnostics.
    """

    rho_l2: float               # L2 of (rho_motif - rho_warpfactory) on the slice
    rho_min_motif: float        # min rho from the motif reconstruction
    rho_min_warpfactory: float  # min rho from the exact T = G/8pi
    nec_violation_fraction: float
    wec_violation_fraction: float
    exotic_energy: float        # |integral of negative rho| from warpfactory
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def warp_factory_cross_check(
    motif: GeometryMotif,
    *,
    n: int = N_SLICE_POINTS,
    L: float | None = None,
) -> WarpCrossCheck:
    """Cross-check the motif's reconstructed rho against the exact Alcubierre T.

    Builds the analytic Alcubierre 4-metric via ``warpfactory.alcubierre_metric``,
    computes ``T_{mu nu} = G_{mu nu} / 8 pi`` and the Eulerian energy density,
    then compares it against the motif's rho_req (from ``_alcubierre_eulerian_rho``).
    Also reports the energy-condition violation fractions and exotic energy budget.
    """
    from ..metrics.probes.warpfactory import (
        alcubierre_metric,
        eulerian_energy_density,
        evaluate_four_metric,
        exotic_energy_budget,
    )

    params = _warp_params_from_motif(motif)
    velocity = params.get("v", 0.5)
    bubble_radius = params.get("R", 2.0)
    sigma = params.get("sigma", 2.0)

    eff_L = L if L is not None else float(
        motif.overrides.get("recipe_basis_radius_max", 8.0)
    ) * 2.0

    # --- Motif rho (analytic formula) ---
    _, _, _, rho_motif = _alcubierre_eulerian_rho(
        velocity=velocity, bubble_radius=bubble_radius, sigma=sigma, L=eff_L, n=n,
    )

    # --- Warp-factory rho (T = G/8pi from the 4-metric) ---
    # Use a grid that covers the same spatial extent; the warp-factory builds
    # a (t,x,y,z) grid.  We evaluate the central time slice.
    n_wf = max(n, 24)  # warp-factory needs enough points for FD stencils
    g, spacing = alcubierre_metric(
        velocity=velocity,
        bubble_radius=bubble_radius,
        sigma=sigma,
        half_width=eff_L,
        n_space=n_wf,
        dt=0.2,
    )
    rho_wf_all = eulerian_energy_density(g, spacing)
    nt = g.shape[0]
    tc = nt // 2
    crop = 4  # discard FD boundary cells
    rho_wf = rho_wf_all[tc, crop:n_wf - crop, crop:n_wf - crop, crop:n_wf - crop]

    # Energy conditions and exotic budget
    ec_report = evaluate_four_metric(g, spacing, crop=crop)
    budget = exotic_energy_budget(g, spacing, crop=crop)

    # Resample warp-factory rho onto the motif grid for L2 comparison.
    # The WF grid is cubic [crop:n_wf-crop]^3; we take the y=centre slice.
    wf_axis = np.linspace(-eff_L, eff_L, n_wf)[crop:n_wf - crop]
    wf_slice = rho_wf[:, rho_wf.shape[1] // 2, :]  # xz-slice at y=centre
    motif_axis = np.linspace(-eff_L, eff_L, n)
    wf_resampled = _resample_2d(wf_axis, wf_axis, wf_slice, motif_axis, motif_axis)

    rho_l2 = _l2_norm(rho_motif, wf_resampled)

    notes: list[str] = []
    if ec_report.nec_min < -1.0e-8:
        notes.append(f"NEC violated (min={ec_report.nec_min:.3e}) — exotic matter confirmed")
    if ec_report.wec_min < -1.0e-8:
        notes.append(f"WEC violated (min={ec_report.wec_min:.3e})")
    notes.append(f"exotic energy = {budget.exotic_energy:.4f} (geometric units)")

    return WarpCrossCheck(
        rho_l2=rho_l2,
        rho_min_motif=float(rho_motif.min()),
        rho_min_warpfactory=float(rho_wf.min()),
        nec_violation_fraction=ec_report.nec_violation_fraction,
        wec_violation_fraction=ec_report.wec_violation_fraction,
        exotic_energy=budget.exotic_energy,
        notes=tuple(notes),
    )
