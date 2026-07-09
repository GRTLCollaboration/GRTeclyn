"""CMA-ES iteration loop for geometry-first matter adjustment.

Given a GeometryMotif target, this module repeatedly adjusts scalar-lump
parameters (amp, width, center, velocity, omega) and invokes GRTresna to find
the matter configuration that best reproduces the target geometry while
satisfying the Hamiltonian and momentum constraints.

The loop is cheap (GRTresna-only, no GPU evolution) and fits between the
geometry-first scout and the existing QD/CMA-ES/HQ campaign stages.

Improvements over the initial implementation:
  - **Amplitude pre-conditioning**: before CMA-ES, a bisection on a global
    amplitude scale factor finds the largest *feasible* (Ham < 10%) amplitude,
    so CMA-ES starts inside the convergence basin instead of wasting evals
    finding it.
  - **Solver fallback ladder**: on crash or high residual, retry once with
    safer relaxation / more iterations / reduced amplitude.
  - **Tighter bounds**: lump centers are bounded to the motif support-region
    radius ± 2×width, not the full box.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ..grtresna.fit.motif import (
    FittedMatter,
    MAX_LUMP_AMP,
    MAX_LUMP_WIDTH,
    MAX_VELOCITY,
    MIN_LUMP_WIDTH,
    build_grtresna_config_from_fitted,
    fitted_matter_to_dict,
    write_fitted_matter_json,
)
from ..grtresna.solver import GRTresnaConfig, solve
from ..grtresna.solver.convergence import parse_convergence
from ..initial_data.motif import GeometryMotif, MomentumTarget, SupportRegion
from ..projection.motif_preservation import compare_motif_preservation
from .mismatch import GATE_FITNESS, MismatchReport, compute_mismatch, feasibility_precheck

logger = logging.getLogger(__name__)

# Per-lump parameter order (fixed lump count / exotic / mode flags).
_LUMP_PARAMS = ("amp", "width", "cx", "cy", "cz", "vx", "vy", "vz", "omega")
_PARAMS_PER_LUMP = len(_LUMP_PARAMS)

# Default bounds for each per-lump parameter (centers will be tightened per-motif).
# Note: amp lower bound is 0.0 so CMA-ES can "turn off" unwanted lumps,
# effectively selecting the lump count as part of the optimisation.
_LOWER = np.array([0.0, MIN_LUMP_WIDTH, -64.0, -64.0, -64.0, -MAX_VELOCITY, -MAX_VELOCITY, -MAX_VELOCITY, -0.5])
_UPPER = np.array([MAX_LUMP_AMP, MAX_LUMP_WIDTH, 64.0, 64.0, 64.0, MAX_VELOCITY, MAX_VELOCITY, MAX_VELOCITY, 0.5])

# Sparsity penalty: rewards using fewer active lumps (amp > threshold).
# This lets CMA-ES naturally select the lump count by zeroing out
# unwanted lumps instead of being forced to use all max_lumps.
# The weight must be large enough that turning off a lump (saving ~0.02
# in geometry mismatch) is worth the penalty reduction.
SPARSITY_AMP_THRESHOLD: float = 0.005   # amps below this count as "off"
SPARSITY_WEIGHT: float = 0.05           # per-active-lump penalty

# --- Phase 1b: amplitude pre-conditioning ------------------------------------
PRECOND_HAM_THRESHOLD: float = 10.0   # Ham% below which amplitude is "feasible"
PRECOND_MAX_STEPS: int = 5             # max bisection steps (×1, ×0.5, ×0.25, ...)
PRECOND_MIN_SCALE: float = 0.0625      # don't go below 1/16 of original amplitude

# --- Phase 1c: solver fallback ladder ----------------------------------------
FALLBACK_HAM_THRESHOLD: float = 50.0   # retry if Ham% above this
FALLBACK_AMPLITUDE_FACTOR: float = 0.7  # reduce amplitude on retry


@dataclass(frozen=True)
class IterationConfig:
    """Configuration for the iterative adjustment loop."""

    max_evals: int = 80
    popsize: int = 8
    sigma0: float = 0.2
    max_concurrent_grtresna: int = 6
    mpi_ranks: int = 8
    plateau_generations: int = 5
    plateau_tolerance: float = 1e-4
    retention_threshold: float = 0.5
    L: float | None = None
    gridinit_n: int = 64
    grtresna_L: float = 128.0
    seed: int | None = None
    # Phase 1b: amplitude pre-conditioning
    precondition: bool = True
    precond_ham_threshold: float = PRECOND_HAM_THRESHOLD
    precond_max_steps: int = PRECOND_MAX_STEPS
    # Phase 1c: solver fallback ladder
    fallback_on_crash: bool = True
    fallback_on_high_residual: bool = True
    # Variable lump count: sparsity penalty lets CMA-ES turn off lumps
    sparsity_weight: float = SPARSITY_WEIGHT


@dataclass
class IterationResult:
    """Output summary of the iteration loop."""

    best_fitness: float
    best_vector: list[float]
    best_fitted_matter: FittedMatter
    best_gridinit_path: Path | None
    best_preservation_score: float
    total_evals: int
    generations: int
    converged: bool
    fitness_history: list[float]
    notes: list[str]


def _vector_from_lumps(lumps: Sequence[Mapping[str, Any]]) -> np.ndarray:
    """Encode lump parameters into a flat vector (warm-start)."""
    parts = []
    for lump in lumps:
        parts.extend([
            float(lump.get("amp", 0.1)),
            float(lump.get("width", 3.0)),
            float(lump.get("center", (0.0, 0.0, 0.0))[0]),
            float(lump.get("center", (0.0, 0.0, 0.0))[1]),
            float(lump.get("center", (0.0, 0.0, 0.0))[2]),
            float(lump.get("velocity", (0.0, 0.0, 0.0))[0]),
            float(lump.get("velocity", (0.0, 0.0, 0.0))[1]),
            float(lump.get("velocity", (0.0, 0.0, 0.0))[2]),
            float(lump.get("omega", 0.0)),
        ])
    return np.array(parts, dtype=float)


def _lumps_from_vector(
    vector: np.ndarray,
    *,
    template_lumps: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Decode a flat vector back into lump dicts (preserving exotic/mode flags)."""
    n_lumps = len(template_lumps)
    lumps = []
    for i in range(n_lumps):
        offset = i * _PARAMS_PER_LUMP
        v = vector[offset: offset + _PARAMS_PER_LUMP]
        template = template_lumps[i]
        lumps.append({
            "amp": float(v[0]),
            "width": float(v[1]),
            "center": (float(v[2]), float(v[3]), float(v[4])),
            "velocity": (float(v[5]), float(v[6]), float(v[7])),
            "omega": float(v[8]),
            "mode": int(template.get("mode", 0)),
            "exotic": int(template.get("exotic", 0)),
        })
    return lumps


def _clip_vector(vector: np.ndarray, n_lumps: int, lower: np.ndarray | None = None, upper: np.ndarray | None = None) -> np.ndarray:
    """Clip vector to per-lump bounds."""
    lo = lower if lower is not None else np.tile(_LOWER, n_lumps)
    hi = upper if upper is not None else np.tile(_UPPER, n_lumps)
    return np.clip(vector, lo, hi)


# --- Phase 1d: tighter bounds from motif support regions ---------------------

def _compute_tight_bounds(
    motif: GeometryMotif,
    n_lumps: int,
    *,
    grtresna_L: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-lump bounds tightened around motif support regions.

    Centers are bounded to support-region center ± 2×width (or a fallback
    radius if no support regions exist).  Amplitude/width/velocity bounds
    stay at the global defaults.
    """
    lower = np.tile(_LOWER, n_lumps).copy()
    upper = np.tile(_UPPER, n_lumps).copy()

    # Determine center bounds from support regions
    if motif.support_regions:
        # Use the union of support region extents
        centers = []
        max_extent = 0.0
        for region in motif.support_regions:
            centers.append(np.array(region.center, dtype=float))
            extent = float(region.width) * 2.0
            max_extent = max(max_extent, extent)

        # Use the centroid of support regions as the center of the allowed box
        centroid = np.mean(centers, axis=0)
        half_box = max(max_extent, 8.0)  # at least ±8 in each direction
        # Don't exceed the grid box
        half_box = min(half_box, grtresna_L * 0.4)

        for i in range(n_lumps):
            offset = i * _PARAMS_PER_LUMP
            lower[offset + 2] = centroid[0] - half_box
            upper[offset + 2] = centroid[0] + half_box
            lower[offset + 3] = centroid[1] - half_box
            upper[offset + 3] = centroid[1] + half_box
            lower[offset + 4] = centroid[2] - half_box
            upper[offset + 4] = centroid[2] + half_box
    # else: keep default ±64 bounds

    return lower, upper


# --- Phase 1b: amplitude pre-conditioning ------------------------------------

def _scale_lump_amplitudes(lumps: list[dict[str, Any]], scale: float) -> list[dict[str, Any]]:
    """Return a copy of lumps with all amplitudes multiplied by *scale*."""
    return [
        {**lump, "amp": float(lump.get("amp", 0.1)) * scale}
        for lump in lumps
    ]


def _run_precondition_solve(
    fitted: FittedMatter,
    base_cfg: GRTresnaConfig,
    work_dir: Path,
    gridinit_path: Path,
) -> tuple[float, Path | None]:
    """Run a single GRTresna solve and return (ham_pct, gridinit_path).

    Returns (inf, None) if the solve crashes.
    """
    cfg = build_grtresna_config_from_fitted(fitted, base=base_cfg)
    grtresna_dir = work_dir / "grtresna"
    try:
        solve(cfg, work_dir=grtresna_dir, gridinit_path=gridinit_path)
    except (RuntimeError, FileNotFoundError, OSError) as exc:
        logger.debug("precondition solve failed: %s", exc)
        return float("inf"), None

    conv = parse_convergence(grtresna_dir)
    if conv is not None:
        return float(conv.get("ham_pct", 100.0)), gridinit_path
    return 100.0, gridinit_path


def amplitude_precondition(
    fitted: FittedMatter,
    base_cfg: GRTresnaConfig,
    work_dir: Path,
    *,
    ham_threshold: float = PRECOND_HAM_THRESHOLD,
    max_steps: int = PRECOND_MAX_STEPS,
) -> tuple[FittedMatter, float, list[str]]:
    """Bisection on a global amplitude scale to find the largest feasible config.

    Returns (preconditioned_fitted, best_scale, notes).  If the original
    amplitude is already feasible, returns it unchanged with scale=1.0.
    """
    notes: list[str] = []
    work_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: try original amplitude
    gridinit = work_dir / "initial_data.gridinit"
    ham, _ = _run_precondition_solve(fitted, base_cfg, work_dir / "step_1x", gridinit)
    logger.info("precondition: scale=1.0, Ham=%.2f%%", ham)

    if ham <= ham_threshold:
        notes.append(f"precondition: original amplitude feasible (Ham={ham:.2f}%)")
        return fitted, 1.0, notes

    # Step 2: bisection — try halving until feasible
    best_scale = 1.0
    best_fitted = fitted
    scale = 0.5
    for step in range(max_steps):
        if scale < PRECOND_MIN_SCALE:
            break

        scaled_lumps = _scale_lump_amplitudes(list(fitted.lumps), scale)
        scaled_fitted = FittedMatter(
            lumps=tuple(scaled_lumps),
            scalar_mass=fitted.scalar_mass,
            maximal_slicing=fitted.maximal_slicing,
            static_lens_only=fitted.static_lens_only,
            momentum_target=fitted.momentum_target,
            notes=fitted.notes,
        )

        step_dir = work_dir / f"step_{scale:.4f}"
        gridinit = step_dir / "initial_data.gridinit"
        ham, _ = _run_precondition_solve(scaled_fitted, base_cfg, step_dir, gridinit)
        logger.info("precondition: scale=%.4f, Ham=%.2f%%", scale, ham)

        if ham <= ham_threshold:
            best_scale = scale
            best_fitted = scaled_fitted
            notes.append(f"precondition: feasible at scale={scale:.4f} (Ham={ham:.2f}%)")
            break
        scale *= 0.5

    if best_scale == 1.0:
        notes.append(
            f"precondition: no feasible amplitude found in {max_steps} steps; "
            f"using original (best Ham={ham:.2f}%)"
        )

    return best_fitted, best_scale, notes


# --- Phase 1c: solver fallback ladder ----------------------------------------

def _solve_with_fallback(
    cfg: GRTresnaConfig,
    work_dir: Path,
    gridinit_path: Path,
    *,
    fallback_on_crash: bool = True,
    fallback_on_high_residual: bool = True,
) -> tuple[Path | None, str | None]:
    """Run GRTresna with a fallback ladder on crash or high residual.

    Returns (gridinit_path, error_message).  On success, gridinit_path is
    set and error_message is None.  On failure, gridinit_path is None and
    error_message describes the failure.
    """
    # First attempt
    try:
        solve(cfg, work_dir=work_dir, gridinit_path=gridinit_path)
    except (RuntimeError, FileNotFoundError, OSError) as exc:
        if not fallback_on_crash:
            return None, f"solve failed: {exc}"

        # Fallback: safer relaxation, more iterations, reduced amplitude
        logger.debug("solve crashed, trying fallback: %s", exc)
        fb_cfg = dataclasses.replace(cfg)
        fb_cfg.psi_relaxation = min(fb_cfg.psi_relaxation, 0.4)
        fb_cfg.max_NL_iterations = max(fb_cfg.max_NL_iterations, 150)
        # Reduce all lump amplitudes
        if fb_cfg.lumps:
            fb_cfg.lumps = [
                {**l, "amp": float(l.get("amp", 0.1)) * FALLBACK_AMPLITUDE_FACTOR}
                for l in fb_cfg.lumps
            ]
        fb_dir = work_dir.parent / (work_dir.name + "_fallback")
        try:
            solve(fb_cfg, work_dir=fb_dir, gridinit_path=gridinit_path)
            logger.info("fallback solve succeeded after crash")
            return gridinit_path, None
        except (RuntimeError, FileNotFoundError, OSError) as exc2:
            return None, f"solve failed (fallback also failed): {exc2}"

    # Check residual — if too high, retry with safer settings
    if fallback_on_high_residual:
        conv = parse_convergence(work_dir)
        if conv is not None:
            ham = float(conv.get("ham_pct", 0.0))
            if ham > FALLBACK_HAM_THRESHOLD:
                logger.debug("high residual Ham=%.1f%%, trying fallback", ham)
                fb_cfg = dataclasses.replace(cfg)
                fb_cfg.psi_relaxation = min(fb_cfg.psi_relaxation, 0.4)
                fb_cfg.max_NL_iterations = max(fb_cfg.max_NL_iterations, 150)
                if fb_cfg.lumps:
                    fb_cfg.lumps = [
                        {**l, "amp": float(l.get("amp", 0.1)) * FALLBACK_AMPLITUDE_FACTOR}
                        for l in fb_cfg.lumps
                    ]
                fb_dir = work_dir.parent / (work_dir.name + "_fallback")
                try:
                    solve(fb_cfg, work_dir=fb_dir, gridinit_path=gridinit_path)
                    fb_conv = parse_convergence(fb_dir)
                    if fb_conv is not None:
                        fb_ham = float(fb_conv.get("ham_pct", 0.0))
                        if fb_ham < ham:
                            logger.info(
                                "fallback improved Ham: %.1f%% → %.1f%%", ham, fb_ham
                            )
                            return gridinit_path, None
                except (RuntimeError, FileNotFoundError, OSError):
                    pass  # keep the original solve result

    return gridinit_path, None


def _sparsity_penalty(vector: np.ndarray, n_lumps: int) -> float:
    """Count active lumps (amp > threshold) and return a penalty.

    This lets CMA-ES naturally select the lump count: lumps with amp ≈ 0
    are effectively "off" and don't contribute to the matter distribution.
    The penalty rewards simpler configs (fewer active lumps) when the
    geometry mismatch is similar.
    """
    n_active = 0
    for i in range(n_lumps):
        offset = i * _PARAMS_PER_LUMP
        if vector[offset] > SPARSITY_AMP_THRESHOLD:
            n_active += 1
    return float(n_active)


def _evaluate_candidate(
    *,
    vector: np.ndarray,
    motif: GeometryMotif,
    template_lumps: Sequence[Mapping[str, Any]],
    base_fitted: FittedMatter,
    base_cfg: GRTresnaConfig,
    work_dir: Path,
    L: float | None,
    fallback_on_crash: bool = True,
    fallback_on_high_residual: bool = True,
    sparsity_weight: float = SPARSITY_WEIGHT,
) -> tuple[float, MismatchReport, Path | None]:
    """Single candidate evaluation: vector → lumps → solve → mismatch."""
    n_lumps = len(template_lumps)
    lumps = _lumps_from_vector(vector, template_lumps=template_lumps)

    # Build updated FittedMatter (preserving non-lump fields)
    fitted = FittedMatter(
        lumps=tuple(lumps),
        scalar_mass=base_fitted.scalar_mass,
        maximal_slicing=base_fitted.maximal_slicing,
        static_lens_only=base_fitted.static_lens_only,
        momentum_target=base_fitted.momentum_target,
        notes=base_fitted.notes,
    )

    cfg = build_grtresna_config_from_fitted(fitted, base=base_cfg)
    gridinit_path = work_dir / "initial_data.gridinit"
    grtresna_dir = work_dir / "grtresna"

    gridinit_path_result, error_msg = _solve_with_fallback(
        cfg,
        work_dir=grtresna_dir,
        gridinit_path=gridinit_path,
        fallback_on_crash=fallback_on_crash,
        fallback_on_high_residual=fallback_on_high_residual,
    )

    if gridinit_path_result is None:
        logger.debug("GRTresna solve failed for %s: %s", work_dir.name, error_msg)
        report = MismatchReport(
            fitness=GATE_FITNESS,
            chi_l2=0.0,
            beta_l2=0.0,
            exotic_penalty=0.0,
            convergence_penalty=0.0,
            solve_failed=True,
            notes=(error_msg or "solve failed",),
        )
        return GATE_FITNESS, report, None

    # Use the grtresna_dir for convergence data (fallback may have written there)
    conv_dir = grtresna_dir
    fb_dir = work_dir.parent / (grtresna_dir.name + "_fallback")
    if fb_dir.exists():
        conv_dir = fb_dir

    report = compute_mismatch(
        motif,
        gridinit_path,
        lumps=lumps,
        grtresna_work_dir=conv_dir,
        L=L,
    )

    # Add sparsity penalty so CMA-ES can select lump count by zeroing
    # unwanted lumps.  Only applies in the feasible phase (low convergence
    # penalty) — in the infeasible phase, the feasibility term dominates.
    if report.convergence_penalty < 0.5:
        sparsity = _sparsity_penalty(vector, n_lumps)
        fitness = report.fitness + sparsity_weight * sparsity
    else:
        fitness = report.fitness

    return fitness, report, gridinit_path


def run_iterate(
    motif: GeometryMotif,
    initial_fitted: FittedMatter,
    *,
    out_dir: Path,
    config: IterationConfig | None = None,
    base_grtresna_cfg: GRTresnaConfig | None = None,
) -> IterationResult:
    """Run the CMA-ES iteration loop to adjust matter for a geometry target.

    Parameters
    ----------
    motif : GeometryMotif
        Target geometry (carries overrides for profile generation).
    initial_fitted : FittedMatter
        Starting point from one-shot fit (warm-start for CMA-ES).
    out_dir : Path
        Output directory for logs, best gridinit, and fitted_matter JSONs.
    config : IterationConfig, optional
        Loop parameters (defaults are sensible for 1-3 lumps at n=64).
    base_grtresna_cfg : GRTresnaConfig, optional
        Base solver config (overridden by fitted-matter fields each eval).
    """
    try:
        import cma
    except ImportError as exc:
        raise ImportError(
            "Iteration loop requires the 'cma' package. Install: pip install cma"
        ) from exc

    cfg = config or IterationConfig()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    n_lumps = len(initial_fitted.lumps)
    if n_lumps == 0:
        raise ValueError("initial_fitted has no lumps to optimize")

    # Base GRTresna config
    if base_grtresna_cfg is None:
        n_grid = cfg.gridinit_n
        base_grtresna_cfg = GRTresnaConfig(
            mpi_ranks=cfg.mpi_ranks,
            N=(n_grid, n_grid, n_grid // 2),
            L=cfg.grtresna_L,
            gridinit_nx=n_grid,
            gridinit_ny=n_grid,
            gridinit_nz=n_grid,
            scalar_mass=initial_fitted.scalar_mass,
            dphi=0.0,
            dpi=0.0,
            bh1_bare_mass=0.0,
            bh1_spin=(0.0, 0.0, 0.0),
            # For exotic matter the default stall tolerance (2%) gives up too
            # early — Ham can plateau at ~96% then slowly creep down.  Give
            # the solver more iterations and a tighter stall threshold.
            max_NL_iterations=200,
            nl_stall_tolerance=0.005,
        )

    notes: list[str] = []

    # --- Phase 4c: feasibility pre-check -------------------------------------
    fe_check = feasibility_precheck(motif)
    if fe_check.notes:
        notes.extend(fe_check.notes)
    if not fe_check.feasible:
        notes.append(
            f"feasibility pre-check: target classified as '{fe_check.risk_level}' "
            f"(rho_peak={fe_check.rho_peak:.4f}) — proceeding but expect poor convergence"
        )
    logger.info(
        "feasibility pre-check: risk=%s, rho_peak=%.4f, feasible=%s",
        fe_check.risk_level, fe_check.rho_peak, fe_check.feasible,
    )

    # --- Phase 1b: amplitude pre-conditioning --------------------------------
    working_fitted = initial_fitted
    if cfg.precondition:
        precond_dir = out_dir / "precondition"
        working_fitted, best_scale, precond_notes = amplitude_precondition(
            initial_fitted,
            base_grtresna_cfg,
            precond_dir,
            ham_threshold=cfg.precond_ham_threshold,
            max_steps=cfg.precond_max_steps,
        )
        notes.extend(precond_notes)
        # Clean up precondition work dirs to save disk
        shutil.rmtree(precond_dir, ignore_errors=True)

    template_lumps = list(working_fitted.lumps)
    x0 = _vector_from_lumps(template_lumps)

    # --- Phase 1d: tighter bounds from motif support regions -----------------
    lower, upper = _compute_tight_bounds(motif, n_lumps, grtresna_L=cfg.grtresna_L)
    x0 = _clip_vector(x0, n_lumps, lower=lower, upper=upper)
    bounds = [lower.tolist(), upper.tolist()]

    opts: dict[str, Any] = {
        "popsize": cfg.popsize,
        "bounds": bounds,
        "maxfevals": cfg.max_evals,
        "verbose": -9,  # suppress CMA-ES stdout
        "CMA_stds": (cfg.sigma0 * (upper - lower)).tolist(),
    }
    if cfg.seed is not None:
        opts["seed"] = cfg.seed

    es = cma.CMAEvolutionStrategy(x0.tolist(), cfg.sigma0, opts)

    # Logging
    log_path = out_dir / "iteration_log.jsonl"
    log_handle = open(log_path, "w", encoding="utf-8")

    best_fitness = GATE_FITNESS
    best_vector = x0.copy()
    best_gridinit: Path | None = None
    best_dir = out_dir / "best"
    best_dir.mkdir(exist_ok=True)
    fitness_history: list[float] = []
    total_evals = 0
    generation = 0
    plateau_count = 0
    prev_best = GATE_FITNESS

    try:
        while not es.stop() and total_evals < cfg.max_evals:
            generation += 1
            solutions = es.ask()
            n_pop = len(solutions)
            clipped = [_clip_vector(np.array(s), n_lumps, lower=lower, upper=upper) for s in solutions]

            # Parallel GRTresna evaluations
            fitnesses = [GATE_FITNESS] * n_pop
            with ThreadPoolExecutor(max_workers=cfg.max_concurrent_grtresna) as pool:
                futures = {}
                for idx, vec in enumerate(clipped):
                    eval_id = total_evals + idx
                    eval_dir = out_dir / f"eval_{eval_id:06d}"
                    eval_dir.mkdir(exist_ok=True)
                    fut = pool.submit(
                        _evaluate_candidate,
                        vector=vec,
                        motif=motif,
                        template_lumps=template_lumps,
                        base_fitted=working_fitted,
                        base_cfg=dataclasses.replace(base_grtresna_cfg),
                        work_dir=eval_dir,
                        L=cfg.L,
                        fallback_on_crash=cfg.fallback_on_crash,
                        fallback_on_high_residual=cfg.fallback_on_high_residual,
                        sparsity_weight=cfg.sparsity_weight,
                    )
                    futures[fut] = idx

                for fut in as_completed(futures):
                    idx = futures[fut]
                    eval_id = total_evals + idx
                    try:
                        fitness, report, gridinit_path = fut.result()
                    except Exception as exc:
                        logger.warning("eval %d raised: %s", eval_id, exc)
                        fitness = GATE_FITNESS
                        report = MismatchReport(
                            fitness=GATE_FITNESS,
                            chi_l2=0.0,
                            beta_l2=0.0,
                            exotic_penalty=0.0,
                            convergence_penalty=0.0,
                            solve_failed=True,
                            notes=(f"exception: {exc}",),
                        )
                        gridinit_path = None

                    fitnesses[idx] = fitness

                    # Log this eval (report.to_dict() has raw fitness;
                    # overwrite with total fitness including sparsity)
                    log_entry = {
                        "eval_id": eval_id,
                        "generation": generation,
                        **report.to_dict(),
                        "fitness": fitness,  # total fitness (with sparsity)
                    }
                    log_handle.write(json.dumps(log_entry) + "\n")
                    log_handle.flush()

                    # Track best
                    if fitness < best_fitness:
                        best_fitness = fitness
                        best_vector = clipped[idx].copy()
                        # Copy best gridinit to stable location
                        if gridinit_path is not None and gridinit_path.exists():
                            dst = best_dir / "initial_data.gridinit"
                            shutil.copy2(gridinit_path, dst)
                            best_gridinit = dst

            total_evals += n_pop
            es.tell(solutions, fitnesses)

            fitness_history.append(best_fitness)
            logger.info(
                "gen %d: best_fitness=%.6f, evals=%d/%d",
                generation, best_fitness, total_evals, cfg.max_evals,
            )

            # Plateau detection
            if abs(prev_best - best_fitness) < cfg.plateau_tolerance:
                plateau_count += 1
            else:
                plateau_count = 0
            prev_best = best_fitness

            if plateau_count >= cfg.plateau_generations:
                notes.append(
                    f"plateau stop after {cfg.plateau_generations} generations "
                    f"with no improvement > {cfg.plateau_tolerance}"
                )
                break

            # Prune non-best eval dirs to save disk
            for idx in range(n_pop):
                eval_id = total_evals - n_pop + idx
                eval_dir = out_dir / f"eval_{eval_id:06d}"
                if eval_dir.exists() and eval_dir != best_dir:
                    gridinit_in_dir = eval_dir / "initial_data.gridinit"
                    if not gridinit_in_dir.exists() or fitnesses[idx] > best_fitness:
                        shutil.rmtree(eval_dir, ignore_errors=True)

    finally:
        log_handle.close()

    # Build best FittedMatter for output
    best_lumps = _lumps_from_vector(best_vector, template_lumps=template_lumps)
    best_fitted = FittedMatter(
        lumps=tuple(best_lumps),
        scalar_mass=working_fitted.scalar_mass,
        maximal_slicing=working_fitted.maximal_slicing,
        static_lens_only=working_fitted.static_lens_only,
        momentum_target=working_fitted.momentum_target,
        notes=working_fitted.notes + (f"iterate: {total_evals} evals, best_fitness={best_fitness:.6f}",),
    )
    write_fitted_matter_json(best_fitted, out_dir / "best_fitted_matter.json")

    # Preservation check on best gridinit
    preservation_score = 0.0
    if best_gridinit is not None and best_gridinit.exists():
        try:
            preservation = compare_motif_preservation(
                motif, best_gridinit, ftl_L=cfg.L,
            )
            preservation_score = preservation.retention_score
        except Exception as exc:
            notes.append(f"preservation check failed: {exc}")

    converged = preservation_score >= cfg.retention_threshold

    # Write summary
    summary = {
        "best_fitness": best_fitness,
        "best_preservation_score": preservation_score,
        "total_evals": total_evals,
        "generations": generation,
        "converged": converged,
        "fitness_history": fitness_history,
        "notes": notes,
    }
    (out_dir / "iterate_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    return IterationResult(
        best_fitness=best_fitness,
        best_vector=best_vector.tolist(),
        best_fitted_matter=best_fitted,
        best_gridinit_path=best_gridinit,
        best_preservation_score=preservation_score,
        total_evals=total_evals,
        generations=generation,
        converged=converged,
        fitness_history=fitness_history,
        notes=notes,
    )
