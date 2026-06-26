"""CMA-ES optimization driver for the RadialRecipe metric search.

Replaces random sampling with gradient-free optimization over the
Gaussian basis coefficients.  Each evaluation runs a full GRTeclyn
episode and returns the negative score (CMA-ES minimizes).

Supports multi-GPU parallel evaluation within each generation.
"""

from __future__ import annotations

from ..grtresna_convergence_gate import (
    DEFAULT_GRTRESNA_CONVERGENCE_CONFIG,
    GRTRESNA_REJECTION_BASE_FITNESS,
    GRTRESNA_REJECTION_MAX_EXTRA_FITNESS,
    GRTresnaConvergenceConfig,
    convergence_rejection_fitness as _grtresna_rejection_fitness,
    convergence_rejection_reason as _grtresna_convergence_rejection_reason,
)
from .candidates import (
    _clip_vector,
    _force_exotic_template,
    _jitter_vector,
    _load_warm_start_vectors,
    _overrides_to_vector,
    _random_vector,
    _vector_to_overrides,
)
from .config import build_grtresna_config, parse_convergence_safe
from .dimension import SearchDimension
from .driver import OptimizeResult, run_optimize
from .eval import (
    SOLVED_FTL_REJECTION_BASE_FITNESS,
    SOLVED_FTL_REJECTION_MAX_EXTRA_FITNESS,
    ObjectiveMode,
)
from .spaces import (
    ANGULAR_BASE_OVERRIDES,
    ANGULAR_SEARCH_SPACE,
    DEFAULT_SEARCH_SPACE,
    GRTRESNA_DEFAULT_NUM_LUMPS,
    GRTRESNA_SEARCH_SPACE,
    SH_DEFAULT_ELL_MAX,
    SH_PROFILE_CHOICES,
    SHELL_PROFILE_CHOICES,
    TRAJECTORY_BOSON_PROFILE_CHOICES,
    TRAJECTORY_DEFAULT_NUM_LUMPS,
    TRAJECTORY_PROFILE_CHOICES,
    build_search_space,
    grtresna_boson_splash_search_space,
    grtresna_boson_star_search_space,
    grtresna_ring_search_space,
    grtresna_search_space,
    grtresna_sh_search_space,
    grtresna_shell_search_space,
    grtresna_trajectory_boson_search_space,
    grtresna_trajectory_search_space,
)

GRTRESNA_MAX_HAM_PCT = DEFAULT_GRTRESNA_CONVERGENCE_CONFIG.max_ham_pct
GRTRESNA_MAX_MOM_PCT = DEFAULT_GRTRESNA_CONVERGENCE_CONFIG.max_mom_pct

__all__ = [
    "ANGULAR_BASE_OVERRIDES",
    "ANGULAR_SEARCH_SPACE",
    "DEFAULT_GRTRESNA_CONVERGENCE_CONFIG",
    "DEFAULT_SEARCH_SPACE",
    "GRTRESNA_DEFAULT_NUM_LUMPS",
    "GRTRESNA_MAX_HAM_PCT",
    "GRTRESNA_MAX_MOM_PCT",
    "GRTRESNA_REJECTION_BASE_FITNESS",
    "GRTRESNA_REJECTION_MAX_EXTRA_FITNESS",
    "GRTRESNA_SEARCH_SPACE",
    "GRTresnaConvergenceConfig",
    "ObjectiveMode",
    "OptimizeResult",
    "SH_DEFAULT_ELL_MAX",
    "SH_PROFILE_CHOICES",
    "SHELL_PROFILE_CHOICES",
    "SOLVED_FTL_REJECTION_BASE_FITNESS",
    "SOLVED_FTL_REJECTION_MAX_EXTRA_FITNESS",
    "SearchDimension",
    "TRAJECTORY_BOSON_PROFILE_CHOICES",
    "TRAJECTORY_DEFAULT_NUM_LUMPS",
    "TRAJECTORY_PROFILE_CHOICES",
    "_force_exotic_template",
    "_grtresna_convergence_rejection_reason",
    "_grtresna_rejection_fitness",
    "_load_warm_start_vectors",
    "build_grtresna_config",
    "build_search_space",
    "grtresna_boson_splash_search_space",
    "grtresna_boson_star_search_space",
    "grtresna_ring_search_space",
    "grtresna_search_space",
    "grtresna_sh_search_space",
    "grtresna_shell_search_space",
    "grtresna_trajectory_boson_search_space",
    "grtresna_trajectory_search_space",
    "parse_convergence_safe",
    "run_optimize",
]
