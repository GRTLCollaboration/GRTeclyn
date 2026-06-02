"""CMA-ES optimization driver for the RadialRecipe metric search.

Replaces random sampling with gradient-free optimization over the
Gaussian basis coefficients.  Each evaluation runs a full GRTeclyn
episode and returns the negative score (CMA-ES minimizes).

Supports multi-GPU parallel evaluation within each generation.
"""

from __future__ import annotations

import json
import math
import os
import random
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..core.config import (
    DEFAULT_RADIAL_RECIPE_TEMPLATE,
    ExecutableConfig,
    ExampleConfig,
    resolve_example,
)
from ..core.episode import Episode, create_episode, update_metadata, write_json
from ..core.params import write_params
from ..core.runner import run_episode
from ..initial_data.constrained_recipe import constrained_overrides
from ..initial_data.preflight import preflight_check
from ..metrics.episode_metrics import dataclass_to_dict, read_episode_metrics
from ..metrics.score import Score, score_episode
from .surrogate import RBFSurrogate, screen_candidates


ObjectiveMode = str

try:
    import cma
except ImportError:
    cma = None  # type: ignore[assignment]

try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]


@dataclass
class SearchDimension:
    """One axis in the optimizer's search space."""

    param_key: str
    lower: float
    upper: float
    initial: float | None = None

    @property
    def center(self) -> float:
        if self.initial is not None:
            return self.initial
        return 0.5 * (self.lower + self.upper)

    @property
    def range(self) -> float:
        return self.upper - self.lower


DEFAULT_SEARCH_SPACE: list[SearchDimension] = [
    SearchDimension("recipe_chi_coeff_0", -0.5, 0.1, -0.1),
    SearchDimension("recipe_chi_coeff_1", -0.3, 0.3, 0.0),
    SearchDimension("recipe_chi_coeff_2", -0.2, 0.2, 0.0),
    SearchDimension("recipe_chi_coeff_3", -0.2, 0.2, 0.0),
    SearchDimension("recipe_basis_width", 0.3, 3.0, 1.0),
    SearchDimension("recipe_alpha_coeff_0", -0.3, 0.3, 0.0),
    SearchDimension("recipe_alpha_coeff_1", -0.3, 0.3, 0.0),
    SearchDimension("recipe_beta_coeff_0", -0.5, 0.5, 0.0),
    SearchDimension("recipe_beta_coeff_1", -0.5, 0.5, 0.0),
]


# --- Non-spherical (angular) search ----------------------------------------
#
# The radial recipe above is spherically symmetric, so the optimizer can only
# ever discover breathing shells / radial channels.  To open up *directional*
# FTL geometries we add axisymmetric Legendre angular modes -- but ONLY on the
# lapse (alpha) and shift (beta) fields.
#
# This is a deliberate, constraint-preserving choice:
#   * alpha and beta are pure GAUGE: they do not appear in the t=0 Hamiltonian
#     or momentum constraints, so making them aspherical adds NO exotic matter
#     and does NOT break the 1D constrained (phi-from-chi) solve.
#   * Angular structure on the shift tilts/reshapes the local light cones (the
#     same mechanism Alcubierre uses), creating a *directional* operational-FTL
#     channel that the search can sculpt.
#   * Angular modes on chi or K, by contrast, WOULD violate the spherical 1D
#     constraint solve and surface as constraint error / forced exotic content,
#     so they are intentionally left out of this stage.
#
# Only the mode amplitudes are searched; the angular order (ell), radial center
# and width are fixed constants (ANGULAR_BASE_OVERRIDES) to keep the dimension
# count low.  ell=1 gives a fore/aft dipole (a genuine "direction"), ell=2 a
# quadrupolar pinch.
# Amplitude AND radial placement (center rc, width rw) of every angular mode
# are searched, so the optimizer can decide not just how strong each lobe is
# but *where along the radius* the directional deformation sits and how sharp
# it is -- a much richer directional-channel family than fixed-placement modes.
ANGULAR_SEARCH_SPACE: list[SearchDimension] = [
    # lapse mode 0 (ell=1 dipole)
    SearchDimension("recipe_lapse_mode_amp_0", -0.25, 0.25, 0.0),
    SearchDimension("recipe_lapse_mode_rc_0", 0.5, 6.0, 2.5),
    SearchDimension("recipe_lapse_mode_rw_0", 0.8, 4.0, 2.0),
    # lapse mode 1 (ell=2 quadrupole)
    SearchDimension("recipe_lapse_mode_amp_1", -0.25, 0.25, 0.0),
    SearchDimension("recipe_lapse_mode_rc_1", 0.5, 6.0, 2.5),
    SearchDimension("recipe_lapse_mode_rw_1", 0.8, 4.0, 2.0),
    # shift mode 0 (ell=1 dipole)
    SearchDimension("recipe_beta_mode_amp_0", -0.6, 0.6, 0.0),
    SearchDimension("recipe_beta_mode_rc_0", 0.5, 6.0, 2.5),
    SearchDimension("recipe_beta_mode_rw_0", 0.8, 4.0, 2.0),
    # shift mode 1 (ell=2 quadrupole)
    SearchDimension("recipe_beta_mode_amp_1", -0.6, 0.6, 0.0),
    SearchDimension("recipe_beta_mode_rc_1", 0.5, 6.0, 2.5),
    SearchDimension("recipe_beta_mode_rw_1", 0.8, 4.0, 2.0),
]

# Fixed (non-searched) parameters that activate the angular modes above: the
# mode counts and the angular order (ell).  Radial center/width are searched
# (see ANGULAR_SEARCH_SPACE) and therefore intentionally NOT pinned here.
ANGULAR_BASE_OVERRIDES: dict[str, Any] = {
    "recipe_num_lapse_angular_modes": 2,
    "recipe_lapse_mode_ell_0": 1,
    "recipe_lapse_mode_ell_1": 2,
    "recipe_num_beta_angular_modes": 2,
    "recipe_beta_mode_ell_0": 1,
    "recipe_beta_mode_ell_1": 2,
}


# --- GRTresna (constraint-satisfying, momentum-carrying matter) search -------
#
# When GRTresna is in the loop the initial data is produced by a full 3D
# elliptic constraint solve rather than the 1D radial recipe.  This unlocks a
# genuinely new FTL axis the recipe cannot represent: a scalar-field cloud that
# carries net linear and/or angular MOMENTUM, with the momentum constraint
# solved (S_i = -Pi d_i phi != 0 -> matter-momentum-driven frame dragging).
#
# These dimensions use the ``grtresna_`` prefix so the objective can route them
# to a GRTresnaConfig (which writes GRTresna's params.txt and runs the solver)
# instead of writing them into GRTeclyn's params.txt.  Everything downstream
# (GPU evolution + scoring) is unchanged: GRTresna just supplies the .gridinit.
#
# The matter source is a basis of K scalar "lumps"; each lump k contributes 11
# searched dimensions (grtresna_lump{k}_*).  A superposition paints an arbitrary
# energy/momentum distribution, which the elliptic solve maps to a rich family
# of constraint-satisfying geometries (chi wells + momentum-driven A_ij).
GRTRESNA_DEFAULT_NUM_LUMPS = 5


def grtresna_search_space(num_lumps: int = GRTRESNA_DEFAULT_NUM_LUMPS) -> list[SearchDimension]:
    """Build the K-lump GRTresna matter search space.

    Initial lump centres are staggered near the origin so the lumps start
    distinct without painting broad matter across the entire visual frame.
    """
    dims: list[SearchDimension] = []
    for k in range(num_lumps):
        # staggered initial x-centre, symmetric about 0
        cx0 = (k - (num_lumps - 1) / 2.0) * 2.5
        dims += [
            # The FTL channel here is shift-driven frame dragging: lump
            # velocity/omega -> conjugate momentum Pi -> momentum density S_i ->
            # A_ij -> evolved shift, and the cones only tilt superluminal once
            # the momentum density is strong.  The right lever for that is
            # VELOCITY/OMEGA, which we widen aggressively -- NOT amplitude.  The
            # elliptic Hamiltonian solve only converges for modest matter
            # amplitude (the prior search lived at amp~0.1-0.15 with clean
            # constraints; amp>~0.4 degrades to a ~99% Hamiltonian residual =
            # unphysical initial data), so amplitude is kept close to that proven
            # range while the momentum knobs are opened up.
            SearchDimension(f"grtresna_lump{k}_amp", 0.0, 0.35, 0.15),
            SearchDimension(f"grtresna_lump{k}_width", 1.5, 4.5, 3.0),
            SearchDimension(f"grtresna_lump{k}_center_x", -6.0, 6.0, cx0),
            SearchDimension(f"grtresna_lump{k}_center_y", -6.0, 6.0, 0.0),
            SearchDimension(f"grtresna_lump{k}_center_z", -3.0, 3.0, 0.0),
            SearchDimension(f"grtresna_lump{k}_velocity_x", -0.9, 0.9, 0.0),
            SearchDimension(f"grtresna_lump{k}_velocity_y", -0.9, 0.9, 0.0),
            SearchDimension(f"grtresna_lump{k}_velocity_z", -0.9, 0.9, 0.0),
            SearchDimension(f"grtresna_lump{k}_omega", -0.6, 0.6, 0.0),
            # azimuthal mode m (rounded to int): 0 axisymmetric, >=1 enables L_z
            SearchDimension(f"grtresna_lump{k}_mode", 0.0, 2.0, 1.0),
            # phantom/ghost flag (rounded to int): 0 canonical (+rho), 1 EXOTIC
            # (-rho).  The search can flip individual lumps exotic to source the
            # NEC violation a persistent FTL channel needs, while keeping others
            # canonical.  Starts canonical (0.0) so the optimizer departs from
            # known-good positive matter and discovers exotic mixes from there.
            SearchDimension(f"grtresna_lump{k}_exotic", 0.0, 1.0, 0.0),
        ]
    return dims


# Default GRTresna search space (used as a fallback when none is supplied).
GRTRESNA_SEARCH_SPACE: list[SearchDimension] = grtresna_search_space()

_LUMP_KEY_RE = re.compile(r"^grtresna_lump(\d+)_(\w+)$")

GRTRESNA_MAX_HAM_PCT = 5.0
GRTRESNA_MAX_MOM_PCT = 5.0
GRTRESNA_REJECTION_BASE_FITNESS = 100.0
GRTRESNA_REJECTION_MAX_EXTRA_FITNESS = 250.0


def build_search_space(
    nonspherical: bool = False,
    grtresna: bool = False,
    grtresna_lumps: int = GRTRESNA_DEFAULT_NUM_LUMPS,
) -> list[SearchDimension]:
    """Return the optimizer search space.

    When ``nonspherical`` is True the spherically-symmetric radial space is
    extended with axisymmetric Legendre angular modes on the lapse and shift
    (see ANGULAR_SEARCH_SPACE).  The caller is responsible for also merging
    ANGULAR_BASE_OVERRIDES into the simulation base overrides so the modes are
    actually activated in params.txt.

    When ``grtresna`` is True the GRTresna momentum-carrying-matter dimensions
    (a ``grtresna_lumps``-lump scalar basis) REPLACE the radial recipe space,
    because the initial data is then produced entirely by the GRTresna solve.
    """
    if grtresna:
        return grtresna_search_space(grtresna_lumps)
    space = list(DEFAULT_SEARCH_SPACE)
    if nonspherical:
        space += ANGULAR_SEARCH_SPACE
    return space


def build_grtresna_config(
    overrides: Mapping[str, Any], base: "GRTresnaConfig | None" = None
) -> "GRTresnaConfig":
    """Build a GRTresnaConfig from a candidate's ``grtresna_*`` overrides.

    ``base`` provides the fixed solver settings (grid, MPI ranks, BH content,
    cleanup); the searched ``grtresna_lump{k}_*`` keys are assembled into the
    ``lumps`` basis.  The legacy un-indexed ``grtresna_lump_*`` keys are also
    honoured (assembled as a single lump) for backward compatibility.
    """
    import dataclasses

    from ..grtresna.solver import GRTresnaConfig

    cfg = dataclasses.replace(base) if base is not None else GRTresnaConfig()

    def _enable_exotic_safe_solver() -> None:
        """Switch to the K=0 maximal-slicing York/Lichnerowicz path.

        Only used for candidates that contain at least one EXOTIC lump: the
        default CTTKHybrid K = sign*sqrt(24 pi G rho) ansatz produces NaN for
        rho < 0, so exotic configs need the maximal-slicing path (matter sourced
        elliptically, with under-relaxation + a psi-positivity floor for the
        indefinite operator). PURELY CANONICAL candidates keep the standard
        ansatz, which is markedly more robust for multi-lump / angular-mode /
        higher-amplitude positive matter (the maximal-slicing path can diverge
        to NaN on those). ``base`` values, if explicitly set, win.
        """
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

    # Group indexed lump keys by index.
    by_index: dict[int, dict[str, float]] = {}
    for key, val in overrides.items():
        m = _LUMP_KEY_RE.match(str(key))
        if m:
            by_index.setdefault(int(m.group(1)), {})[m.group(2)] = float(val)

    if by_index:
        lumps: list[dict] = []
        for k in sorted(by_index):
            f = by_index[k]
            lumps.append({
                "amp": f.get("amp", 0.0),
                "width": f.get("width", 5.0),
                "center": (f.get("center_x", 0.0), f.get("center_y", 0.0),
                           f.get("center_z", 0.0)),
                "velocity": (f.get("velocity_x", 0.0), f.get("velocity_y", 0.0),
                             f.get("velocity_z", 0.0)),
                "omega": f.get("omega", 0.0),
                "mode": int(round(f.get("mode", 0.0))),
                "exotic": int(round(f.get("exotic", 0.0))),
            })
        cfg.lumps = lumps
        if any(lump["exotic"] for lump in lumps):
            _enable_exotic_safe_solver()
        return cfg

    # Backward-compatible single (un-indexed) lump.
    def _get(key: str, default: float) -> float:
        return float(overrides.get(key, default))

    if any(str(k).startswith("grtresna_lump_") for k in overrides):
        cfg.lump_amp = _get("grtresna_lump_amp", cfg.lump_amp)
        cfg.lump_width = _get("grtresna_lump_width", cfg.lump_width)
        cfg.lump_velocity = (
            _get("grtresna_lump_velocity_x", cfg.lump_velocity[0]),
            _get("grtresna_lump_velocity_y", cfg.lump_velocity[1]),
            _get("grtresna_lump_velocity_z", cfg.lump_velocity[2]),
        )
        cfg.lump_omega = _get("grtresna_lump_omega", cfg.lump_omega)
        if "grtresna_lump_mode" in overrides:
            cfg.lump_mode = int(round(float(overrides["grtresna_lump_mode"])))
        if "grtresna_lump_exotic" in overrides:
            cfg.lump_exotic = int(round(float(overrides["grtresna_lump_exotic"])))
    if cfg.lump_exotic:
        _enable_exotic_safe_solver()
    return cfg


def parse_convergence_safe(work_dir: Path) -> dict[str, float] | None:
    """Best-effort read of GRTresna Ham/Mom convergence (never raises)."""
    try:
        from ..grtresna.solver import parse_convergence

        return parse_convergence(work_dir)
    except Exception:
        return None


def _grtresna_convergence_rejection_reason(
    convergence: Mapping[str, float] | None,
    *,
    max_ham_pct: float = GRTRESNA_MAX_HAM_PCT,
    max_mom_pct: float = GRTRESNA_MAX_MOM_PCT,
) -> str | None:
    """Return a rejection reason for unusable GRTresna initial data."""
    if not convergence:
        return "missing GRTresna convergence diagnostics"
    try:
        ham = float(convergence["ham_pct"])
        mom = float(convergence["mom_pct"])
    except (KeyError, TypeError, ValueError):
        return "invalid GRTresna convergence diagnostics"
    if not math.isfinite(ham) or not math.isfinite(mom):
        return f"nonfinite GRTresna convergence: Ham={ham}, Mom={mom}"
    if ham > max_ham_pct or mom > max_mom_pct:
        return (
            f"GRTresna convergence too poor: Ham={ham:.6g}% "
            f"Mom={mom:.6g}% thresholds=({max_ham_pct:g}%, {max_mom_pct:g}%)"
        )
    return None


def _grtresna_rejection_fitness(
    convergence: Mapping[str, float] | None,
    *,
    base: float = GRTRESNA_REJECTION_BASE_FITNESS,
    max_extra: float = GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
) -> float:
    """Graded minimization penalty so CMA can rank bad GRTresna candidates."""
    if not convergence:
        return base + max_extra
    try:
        ham = float(convergence["ham_pct"])
        mom = float(convergence["mom_pct"])
    except (KeyError, TypeError, ValueError):
        return base + max_extra
    if not math.isfinite(ham) or not math.isfinite(mom):
        return base + max_extra
    residual_penalty = max(0.0, ham - GRTRESNA_MAX_HAM_PCT) + max(0.0, mom - GRTRESNA_MAX_MOM_PCT)
    return base + min(max_extra, residual_penalty)


@dataclass(frozen=True)
class OptimizeResult:
    best_params: dict[str, float]
    best_score: float
    best_episode: str
    generations: int
    evaluations: int
    trajectory: list[dict[str, Any]]


def _vector_to_overrides(
    x: Sequence[float],
    dims: Sequence[SearchDimension],
    base: Mapping[str, Any],
) -> dict[str, Any]:
    overrides = dict(base)
    for xi, dim in zip(x, dims):
        clamped = max(dim.lower, min(dim.upper, xi))
        overrides[dim.param_key] = clamped
    return overrides


def _clip_vector(x: Sequence[float], dims: Sequence[SearchDimension]) -> list[float]:
    """Clamp an optimizer vector to the configured search bounds."""
    return [
        max(dim.lower, min(dim.upper, float(xi)))
        for xi, dim in zip(x, dims)
    ]


def _overrides_to_vector(
    overrides: Mapping[str, Any],
    dims: Sequence[SearchDimension],
) -> list[float]:
    """Map a trajectory overrides dict back onto the current search vector."""
    vals: list[float] = []
    for dim in dims:
        raw = overrides.get(dim.param_key, dim.center)
        try:
            value = float(raw)
        except (TypeError, ValueError):
            value = dim.center
        vals.append(max(dim.lower, min(dim.upper, value)))
    return vals


def _load_warm_start_vectors(
    trajectories: Sequence[Path],
    dims: Sequence[SearchDimension],
    top_k: int,
) -> list[list[float]]:
    """Load top-scoring vectors from one or more prior trajectory.jsonl files."""
    records: list[tuple[float, Mapping[str, Any]]] = []
    for path in trajectories:
        path = path.expanduser().resolve()
        if not path.exists():
            continue
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if rec.get("surrogate_predicted"):
                    continue
                score = rec.get("score")
                overrides = rec.get("overrides")
                if not isinstance(overrides, Mapping):
                    continue
                try:
                    score_f = float(score)
                except (TypeError, ValueError):
                    continue
                if math.isfinite(score_f):
                    records.append((score_f, overrides))
    records.sort(key=lambda item: item[0], reverse=True)
    return [_overrides_to_vector(overrides, dims) for _, overrides in records[:top_k]]


def _random_vector(
    dims: Sequence[SearchDimension],
    rng: random.Random,
) -> list[float]:
    return [rng.uniform(dim.lower, dim.upper) for dim in dims]


def _jitter_vector(
    x: Sequence[float],
    dims: Sequence[SearchDimension],
    rng: random.Random,
    scale: float,
) -> list[float]:
    jittered = [
        float(xi) + rng.gauss(0.0, scale * dim.range)
        for xi, dim in zip(x, dims)
    ]
    return _clip_vector(jittered, dims)


def _force_exotic_template(
    x: Sequence[float],
    dims: Sequence[SearchDimension],
    rng: random.Random,
    template_index: int,
) -> list[float]:
    """Impose a categorical exotic/multi-lump pattern on a continuous vector."""
    vec = _clip_vector(x, dims)
    index = {dim.param_key: i for i, dim in enumerate(dims)}
    lump_ids = sorted({
        int(m.group(1))
        for dim in dims
        if (m := _LUMP_KEY_RE.match(dim.param_key))
    })
    if not lump_ids:
        return vec

    n = len(lump_ids)
    if template_index % 4 == 0:
        exotic = {lump_ids[n // 2]}
    elif template_index % 4 == 1:
        exotic = {lump_ids[0], lump_ids[-1]}
    elif template_index % 4 == 2:
        exotic = {k for i, k in enumerate(lump_ids) if i % 2 == 0}
    else:
        exotic = set(lump_ids)

    radius = rng.uniform(2.0, 5.0)
    for pos, k in enumerate(lump_ids):
        angle = 2.0 * math.pi * pos / max(n, 1)
        tangential = 0.45 if k in exotic else 0.30
        assignments = {
            "amp": rng.uniform(0.12, 0.24),
            "width": rng.uniform(1.5, 3.5),
            "center_x": radius * math.cos(angle),
            "center_y": radius * math.sin(angle),
            "center_z": rng.uniform(-2.0, 2.0),
            "velocity_x": -tangential * math.sin(angle),
            "velocity_y": tangential * math.cos(angle),
            "velocity_z": rng.uniform(-0.15, 0.15),
            "omega": rng.uniform(-0.15, 0.15),
            "mode": 1.0 if pos % 2 == 0 else 2.0,
            "exotic": 1.0 if k in exotic else 0.0,
        }
        for suffix, value in assignments.items():
            key = f"grtresna_lump{k}_{suffix}"
            if key in index:
                i = index[key]
                dim = dims[i]
                vec[i] = max(dim.lower, min(dim.upper, value))
    return vec


def _assign_gpu(index: int, gpu_ids: Sequence[int]) -> str:
    """Round-robin GPU assignment for parallel evaluation."""
    return str(gpu_ids[index % len(gpu_ids)])


def _collect_training(
    trajectory: Sequence[Mapping[str, Any]],
    dims: Sequence[SearchDimension],
):
    """Build (X, y) arrays of evaluated candidates for surrogate fitting."""
    xs: list[list[float]] = []
    ys: list[float] = []
    for rec in trajectory:
        if rec.get("preflight_rejected") or rec.get("surrogate_predicted"):
            continue
        overrides = rec.get("overrides")
        score = rec.get("score")
        if not isinstance(overrides, dict) or score is None:
            continue
        try:
            xs.append([float(overrides[d.param_key]) for d in dims])
            ys.append(float(score))
        except (KeyError, TypeError, ValueError):
            continue
    if not xs:
        return None, None
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def _objective(
    x: Sequence[float],
    *,
    dims: Sequence[SearchDimension],
    base_overrides: Mapping[str, Any],
    opt_dir: Path,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    eval_counter: list[int],
    constrained: bool,
    phantom: bool,
    use_preflight: bool,
    cuda_devices: str | None,
    check_params: bool,
    dry_run: bool,
    trajectory: list[dict[str, Any]],
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    objective_mode: ObjectiveMode = "weighted",
    ftl_L: float | None = None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
    consumer_keep_last: int = 1,
    grtresna: bool = False,
    grtresna_base: "GRTresnaConfig | None" = None,
) -> float:
    """Evaluate one candidate.  Returns negative score (CMA-ES minimizes)."""
    eval_counter[0] += 1
    idx = eval_counter[0]

    overrides = _vector_to_overrides(x, dims, base_overrides)

    # In GRTresna mode the constraint solve replaces the 1D radial recipe, so
    # the recipe-specific constrained/preflight steps are skipped.
    if constrained and not grtresna:
        constrained_overrides(overrides, phantom=phantom)

    if use_preflight and not grtresna:
        pf = preflight_check(overrides, phantom=phantom)
        if not pf.passed:
            record = {
                "eval": idx,
                "score": 0.0,
                "preflight_rejected": True,
                "reason": pf.reason,
                "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
            }
            trajectory.append(record)
            return 100.0

    episode = create_episode(
        opt_dir,
        name=f"eval_{idx:06d}",
        metadata={
            "mode": "optimize",
            "example": example.name,
            "eval_index": idx,
            "overrides": overrides,
        },
    )

    # GRTeclyn params get everything except the grtresna_* search keys (which
    # drive the upstream solver, not GRTeclyn's params.txt).
    gte_overrides = {
        k: v for k, v in overrides.items() if not str(k).startswith("grtresna_")
    }

    if grtresna and not dry_run:
        from ..grtresna.solver import solve as grtresna_solve

        cfg = build_grtresna_config(overrides, grtresna_base)
        try:
            gridinit = grtresna_solve(
                cfg,
                work_dir=episode.path / "grtresna",
                gridinit_path=episode.path / "initial_data.gridinit",
            )
            gte_overrides["recipe_initial_data_file"] = gridinit
            convergence = parse_convergence_safe(episode.path / "grtresna")
            update_metadata(episode, {"grtresna_convergence": convergence})
            rejection_reason = _grtresna_convergence_rejection_reason(convergence)
            if rejection_reason is not None:
                fitness = _grtresna_rejection_fitness(convergence)
                update_metadata(episode, {
                    "grtresna_rejected": True,
                    "grtresna_rejection_reason": rejection_reason,
                    "grtresna_rejection_fitness": fitness,
                    "simulation_exit_code": None,
                })
                record = {
                    "eval": idx,
                    "episode": str(episode.path),
                    "score": -fitness,
                    "fitness": fitness,
                    "grtresna_rejected": True,
                    "reason": rejection_reason,
                    "components": {"grtresna_rejection": -fitness},
                    "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
                }
                trajectory.append(record)
                return fitness
        except Exception as exc:  # solver failure -> penalise, keep searching
            update_metadata(episode, {"grtresna_error": repr(exc)})
            record = {
                "eval": idx,
                "episode": str(episode.path),
                "score": 0.0,
                "grtresna_failed": True,
                "reason": repr(exc),
                "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
            }
            trajectory.append(record)
            return 100.0

    write_params(
        template, episode.params_path,
        episode_dir=episode.path, example=example, overrides=gte_overrides,
    )

    exit_code: int | None = None
    if dry_run:
        update_metadata(episode, {"dry_run": True})
    else:
        if executable is None:
            raise ValueError("executable is required unless dry_run=True")
        try:
            result = run_episode(
                episode, executable,
                check_params=check_params,
                cuda_devices=cuda_devices,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_delete=True,
                consumer_keep_last=consumer_keep_last,
            )
            exit_code = result.returncode
        except Exception as exc:
            exit_code = 1
            update_metadata(episode, {
                "simulation_error": repr(exc),
                "simulation_exit_code": exit_code,
            })

    metrics = read_episode_metrics(episode.path, ftl_L=ftl_L)
    score = score_episode(
        metrics,
        target_stop_time=target_stop_time,
        weights=score_weights,
        objective_mode=objective_mode,
    )

    write_json(episode.score_path, {
        "score": dataclass_to_dict(score),
        "metrics": dataclass_to_dict(metrics),
    })

    record = {
        "eval": idx,
        "episode": str(episode.path),
        "exit_code": exit_code,
        "score": score.total,
        "components": score.components,
        "overrides": {d.param_key: overrides.get(d.param_key) for d in dims},
    }
    trajectory.append(record)
    return -score.total


def run_optimize(
    *,
    runs_dir: Path,
    executable: ExecutableConfig | None = None,
    max_generations: int = 50,
    population_size: int | None = None,
    sigma0: float = 0.3,
    seed: int | None = None,
    base_overrides: Mapping[str, Any] | None = None,
    search_space: Sequence[SearchDimension] | None = None,
    template: Path | None = None,
    example: ExampleConfig | str = "RadialRecipe",
    name: str | None = None,
    dry_run: bool = False,
    constrained: bool = True,
    phantom: bool = True,
    use_preflight: bool = True,
    cuda_devices: str | None = None,
    gpu_ids: Sequence[int] | None = None,
    check_params: bool = True,
    score_weights: Mapping[str, float] | None = None,
    objective_mode: ObjectiveMode = "weighted",
    ftl_L: float | None = None,
    x0: Sequence[float] | None = None,
    consume_plotfiles: bool = True,
    consumer_radii: Sequence[float] = (4.0, 8.0),
    consumer_keep_last: int = 1,
    surrogate: bool = False,
    surrogate_keep_fraction: float = 0.5,
    surrogate_warmup: int | None = None,
    grtresna: bool = False,
    grtresna_config: "GRTresnaConfig | None" = None,
    warm_start_trajectories: Sequence[Path] = (),
    warm_start_top_k: int = 8,
    warm_start_jitter: float = 0.08,
    random_injection_fraction: float = 0.0,
    exotic_injection_fraction: float = 0.0,
) -> OptimizeResult:
    """Run multi-GPU CMA-ES optimization loop.

    Parameters
    ----------
    gpu_ids : list of int, optional
        Available GPU indices for parallel evaluation. If None, uses
        cuda_devices for sequential mode. If provided, each member of
        the CMA-ES population is assigned a GPU round-robin.
    """
    if cma is None:
        raise ImportError(
            "CMA-ES optimization requires the 'cma' package. "
            "Install it with: pip install cma"
        )
    if np is None:
        raise ImportError("numpy is required for optimization.")

    example_cfg = example if isinstance(example, ExampleConfig) else resolve_example(example)
    dims = list(search_space or (GRTRESNA_SEARCH_SPACE if grtresna else DEFAULT_SEARCH_SPACE))
    tpl = template or example_cfg.template
    base = dict(base_overrides or {})
    base.setdefault("N_full", 64)
    base.setdefault("max_level", 2)
    base.setdefault("stop_time", 2.0)
    base.setdefault("plot_interval", 10)
    base.setdefault("checkpoint_interval", -1)
    base.setdefault("dt_multiplier", 0.02)
    base.setdefault("regrid_threshold", 0.01)
    max_lvl = int(base["max_level"])
    if "regrid_interval" not in base and max_lvl > 0:
        intervals = [16] * min(max_lvl, 2) + [8] * max(0, max_lvl - 2)
        base["regrid_interval"] = intervals

    target_stop_time = float(base["stop_time"])

    if name is None:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        name = f"optimize_{timestamp}"
    opt_dir = (runs_dir / name).expanduser().resolve()
    opt_dir.mkdir(parents=True, exist_ok=False)

    warm_vectors = _load_warm_start_vectors(
        warm_start_trajectories, dims, max(0, warm_start_top_k)
    )

    if x0 is not None:
        initial = list(x0)
    elif warm_vectors:
        initial = list(warm_vectors[0])
    else:
        initial = [d.center for d in dims]

    bounds = [[d.lower for d in dims], [d.upper for d in dims]]

    effective_popsize = population_size
    if effective_popsize is None and gpu_ids is not None:
        effective_popsize = len(gpu_ids)

    opts: dict[str, Any] = {
        "maxiter": max_generations,
        "bounds": bounds,
        "CMA_stds": [sigma0 * d.range for d in dims],
        "verbose": -1,
    }
    if effective_popsize is not None:
        opts["popsize"] = effective_popsize
    if seed is not None:
        opts["seed"] = seed

    eval_counter = [0]
    trajectory: list[dict[str, Any]] = []

    write_json(opt_dir / "metadata.json", {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "example": example_cfg.name,
        "max_generations": max_generations,
        "population_size": effective_popsize,
        "sigma0": sigma0,
        "seed": seed,
        "constrained": constrained,
        "phantom": phantom,
        "use_preflight": use_preflight,
        "dry_run": dry_run,
        "gpu_ids": list(gpu_ids) if gpu_ids else None,
        "consume_plotfiles": consume_plotfiles,
        "grtresna": grtresna,
        "objective_mode": objective_mode,
        "warm_start_trajectories": [str(p) for p in warm_start_trajectories],
        "warm_start_top_k": warm_start_top_k,
        "warm_start_jitter": warm_start_jitter,
        "random_injection_fraction": random_injection_fraction,
        "exotic_injection_fraction": exotic_injection_fraction,
        "base_overrides": base,
        "search_space": [
            {"key": d.param_key, "lower": d.lower, "upper": d.upper, "initial": d.center}
            for d in dims
        ],
    })

    es = cma.CMAEvolutionStrategy(initial, sigma0, opts)
    rng = random.Random(seed)

    best_score = -math.inf
    best_params: dict[str, float] = {}
    best_episode = ""
    gen = 0

    print(f"[optimize] Starting CMA-ES: {len(dims)}D, popsize={es.popsize}, "
          f"max_gen={max_generations}, GPUs={gpu_ids or cuda_devices}")
    if warm_vectors:
        print(f"[optimize] loaded {len(warm_vectors)} warm-start vectors")

    def _evaluate_subset(subset: list) -> list[float]:
        if gpu_ids is not None and len(gpu_ids) > 1 and not dry_run:
            return _evaluate_generation_parallel(
                subset,
                dims=dims,
                base_overrides=base,
                opt_dir=opt_dir,
                example=example_cfg,
                template=tpl,
                executable=executable,
                eval_counter=eval_counter,
                constrained=constrained,
                phantom=phantom,
                use_preflight=use_preflight,
                gpu_ids=gpu_ids,
                check_params=check_params,
                trajectory=trajectory,
                target_stop_time=target_stop_time,
                score_weights=score_weights,
                objective_mode=objective_mode,
                ftl_L=ftl_L,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_keep_last=consumer_keep_last,
                grtresna=grtresna,
                grtresna_config=grtresna_config,
            )
        out: list[float] = []
        for sol in subset:
            out.append(_objective(
                sol,
                dims=dims,
                base_overrides=base,
                opt_dir=opt_dir,
                example=example_cfg,
                template=tpl,
                executable=executable,
                eval_counter=eval_counter,
                constrained=constrained,
                phantom=phantom,
                use_preflight=use_preflight,
                cuda_devices=cuda_devices,
                check_params=check_params,
                dry_run=dry_run,
                trajectory=trajectory,
                target_stop_time=target_stop_time,
                score_weights=score_weights,
                objective_mode=objective_mode,
                ftl_L=ftl_L,
                consume_plotfiles=consume_plotfiles,
                consumer_radii=consumer_radii,
                consumer_keep_last=consumer_keep_last,
                grtresna=grtresna,
                grtresna_base=grtresna_config,
            ))
        return out

    warmup = surrogate_warmup if surrogate_warmup is not None else 2 * es.popsize

    while not es.stop():
        gen += 1
        solutions = [list(sol) for sol in es.ask()]
        n_solutions = len(solutions)

        if gen == 1 and warm_vectors:
            n_warm = min(len(warm_vectors), n_solutions)
            for i in range(n_warm):
                source = warm_vectors[i]
                if i == 0:
                    solutions[i] = _clip_vector(source, dims)
                else:
                    solutions[i] = _jitter_vector(
                        source, dims, rng, warm_start_jitter
                    )

        n_random = max(0, min(
            n_solutions,
            int(round(n_solutions * max(0.0, random_injection_fraction))),
        ))
        n_exotic = max(0, min(
            n_solutions - n_random,
            int(round(n_solutions * max(0.0, exotic_injection_fraction))),
        ))
        inject_pos = n_solutions - n_random - n_exotic
        for j in range(n_exotic):
            base_vec = (
                _jitter_vector(warm_vectors[j % len(warm_vectors)], dims, rng, warm_start_jitter)
                if warm_vectors else _random_vector(dims, rng)
            )
            solutions[inject_pos + j] = _force_exotic_template(
                base_vec, dims, rng, template_index=gen + j
            )
        for j in range(n_random):
            solutions[n_solutions - n_random + j] = _random_vector(dims, rng)

        traj_before = len(trajectory)

        x_train, y_train = _collect_training(trajectory, dims)
        use_surrogate_now = (
            surrogate
            and x_train is not None
            and x_train.shape[0] >= warmup
        )

        if use_surrogate_now:
            model = RBFSurrogate(
                lower=np.asarray([d.lower for d in dims], dtype=float),
                upper=np.asarray([d.upper for d in dims], dtype=float),
            ).fit(x_train, y_train)
            sols_arr = np.asarray(solutions, dtype=float)
            eval_idx, predicted = screen_candidates(
                model, sols_arr,
                keep_fraction=surrogate_keep_fraction,
                min_eval=1,
            )
            subset = [solutions[i] for i in eval_idx]
            sub_fit = _evaluate_subset(subset)
            fitnesses = [-float(predicted[i]) for i in range(len(solutions))]
            for k, i in enumerate(eval_idx):
                fitnesses[i] = sub_fit[k]
            n_skipped = len(solutions) - len(eval_idx)
            for i in range(len(solutions)):
                if i not in set(eval_idx):
                    trajectory.append({
                        "eval": None,
                        "surrogate_predicted": True,
                        "score": float(predicted[i]),
                        "overrides": {d.param_key: float(solutions[i][j]) for j, d in enumerate(dims)},
                    })
            print(f"[optimize] gen {gen}: surrogate screened "
                  f"{len(eval_idx)}/{len(solutions)} evaluated on GPU "
                  f"({n_skipped} predicted)")
        else:
            fitnesses = _evaluate_subset(list(solutions))

        es.tell(solutions, fitnesses)

        evaluated_records = [
            rec for rec in trajectory[traj_before:]
            if not rec.get("surrogate_predicted")
        ]
        gen_scores = [rec.get("score", -math.inf) for rec in evaluated_records]
        gen_best = max(gen_scores) if gen_scores else -math.inf
        gen_mean = sum(gen_scores) / len(gen_scores) if gen_scores else 0.0

        for rec in evaluated_records:
            sc = rec.get("score", -math.inf)
            if sc > best_score:
                best_score = sc
                best_params = dict(rec.get("overrides", {}))
                best_episode = rec.get("episode", "")

        print(f"[optimize] gen {gen}/{max_generations}: "
              f"best={gen_best:.4f} mean={gen_mean:.4f} "
              f"all-time-best={best_score:.4f} evals={eval_counter[0]}")

        with (opt_dir / "trajectory.jsonl").open("w", encoding="utf-8") as fh:
            for rec in trajectory:
                fh.write(json.dumps(rec, sort_keys=True) + "\n")

    result = OptimizeResult(
        best_params=best_params,
        best_score=best_score,
        best_episode=best_episode,
        generations=gen,
        evaluations=eval_counter[0],
        trajectory=trajectory,
    )

    write_json(opt_dir / "result.json", {
        "best_params": result.best_params,
        "best_score": result.best_score,
        "best_episode": result.best_episode,
        "generations": result.generations,
        "evaluations": result.evaluations,
    })

    print(f"\n[optimize] Done. {gen} generations, {eval_counter[0]} evaluations.")
    print(f"[optimize] Best score: {best_score:.4f}")
    print(f"[optimize] Best episode: {best_episode}")
    print(f"[optimize] Results: {opt_dir}")

    return result


def _evaluate_generation_parallel(
    solutions: list,
    *,
    dims: Sequence[SearchDimension],
    base_overrides: Mapping[str, Any],
    opt_dir: Path,
    example: ExampleConfig,
    template: Path,
    executable: ExecutableConfig | None,
    eval_counter: list[int],
    constrained: bool,
    phantom: bool,
    use_preflight: bool,
    gpu_ids: Sequence[int],
    check_params: bool,
    trajectory: list[dict[str, Any]],
    target_stop_time: float | None,
    score_weights: Mapping[str, float] | None,
    objective_mode: ObjectiveMode,
    ftl_L: float | None,
    consume_plotfiles: bool,
    consumer_radii: Sequence[float],
    consumer_keep_last: int = 1,
    grtresna: bool = False,
    grtresna_config: "GRTresnaConfig | None" = None,
) -> list[float]:
    """Evaluate an entire CMA-ES generation in parallel across GPUs.

    Each solution is assigned a GPU round-robin, and all are launched
    concurrently via threads (GIL released during subprocess waits).
    """
    import threading

    fitnesses: list[float | None] = [None] * len(solutions)
    lock = threading.Lock()

    def _eval_one(idx_in_gen: int, sol) -> None:
        gpu = _assign_gpu(idx_in_gen, gpu_ids)
        f = _objective(
            sol,
            dims=dims,
            base_overrides=base_overrides,
            opt_dir=opt_dir,
            example=example,
            template=template,
            executable=executable,
            eval_counter=eval_counter,
            constrained=constrained,
            phantom=phantom,
            use_preflight=use_preflight,
            cuda_devices=gpu,
            check_params=check_params,
            dry_run=False,
            trajectory=trajectory,
            target_stop_time=target_stop_time,
            score_weights=score_weights,
            objective_mode=objective_mode,
            ftl_L=ftl_L,
            consume_plotfiles=consume_plotfiles,
            consumer_radii=consumer_radii,
            consumer_keep_last=consumer_keep_last,
            grtresna=grtresna,
            grtresna_base=grtresna_config,
        )
        with lock:
            fitnesses[idx_in_gen] = f

    threads = []
    for i, sol in enumerate(solutions):
        t = threading.Thread(target=_eval_one, args=(i, sol))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()

    return [f if f is not None else 100.0 for f in fitnesses]
