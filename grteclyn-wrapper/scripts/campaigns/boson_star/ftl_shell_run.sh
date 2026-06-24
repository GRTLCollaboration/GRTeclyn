#!/usr/bin/env bash
# Bosonic shell MAP-Elites with FTL scoring (for RL chassis / general_ftl elites).
#
# Same matter as splash (BosonStarBH superposed shell, grtresna_complex_scalar) but
# OBJECTIVE_MODE=general_ftl + 4D geodesic — NOT critical_collapse.
#
# Exotic matter: ON by default — searches grtresna_shell_exotic_fraction/phase
# (same idiom as real-scalar wormhole QD). Splash pins canonical-only; this path
# does not.
#
# Contrast:
#   splash/run.sh      → boson shell + critical_collapse (grlab), exotic off
#   boson_star/run.sh  → boson_star sector + ftl_first (7-D centered OR shell)
#   general_ftl/       → real scalar lumps + general_ftl (wormhole QD)
#
# Example (200 evals, frames, t=16, ~6 plot dumps):
#   cd grteclyn-wrapper
#   QD_NAME=boson_shell_ftl_rl_v1 QD_TARGET_EVALS=200 \
#     GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=1 \
#     bash scripts/campaigns/boson_star/ftl_shell_run.sh

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=shell
export GRTRESNA_SHELL_PROFILE="${GRTRESNA_SHELL_PROFILE:-compact}"
export GRTRESNA_MATTER_COUPLING="${GRTRESNA_MATTER_COUPLING:-canonical}"
export GRTRESNA_BOSON_ALLOW_EXOTIC="${GRTRESNA_BOSON_ALLOW_EXOTIC:-1}"
export GRTRESNA_BOSON_MATCHED_BOUNDS="${GRTRESNA_BOSON_MATCHED_BOUNDS:-0}"
export GRTRESNA_FULL_Z=1

export OBJECTIVE_MODE="${OBJECTIVE_MODE:-general_ftl}"
export DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-ftl_lifetime}"

# Static shell; do NOT pin scalar_sign (exotic wedge search active).
export PIN_DIMS="${PIN_DIMS:-grtresna_shell_static=1}"

export STOP_TIME="${STOP_TIME:-16.0}"
# dt=0.01 → 1600 steps; interval 320 → 6 plotfiles over [0, 16].
export PLOT_INTERVAL="${PLOT_INTERVAL:-320}"

export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-1}"
export GRTECLYN_SPLASH_EARLY_TERM=0

export ITERATIONS="${ITERATIONS:-80}"
export POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-3e-2}"
export USE_PIPELINE="${USE_PIPELINE:-1}"
export MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-5}"

exec bash "${SCRIPT_DIR}/../qd/run.sh"
