#!/usr/bin/env bash
# Real-scalar shell MAP-Elites with FTL scoring — boson comparison arm.
#
# Matched to boson_star/ftl_shell_run.sh (same objective, grid, t=16, plots,
# frames, postload gate, pipeline). Only matter sector differs:
#   scalar + grtresna_independent_scalars  vs  boson + grtresna_complex_scalar
#
# Example (200 evals, after boson_shell_ftl_rl_v1 finishes):
#   cd grteclyn-wrapper
#   QD_NAME=scalar_shell_ftl_rl_v1 QD_TARGET_EVALS=200 \
#     GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=1 \
#     bash scripts/campaigns/general_ftl/scalar_shell_ftl_run.sh

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=scalar
export GRTRESNA_ANSATZ=shell
export GRTRESNA_SHELL_PROFILE="${GRTRESNA_SHELL_PROFILE:-compact}"
export GRTRESNA_MATTER_COUPLING="${GRTRESNA_MATTER_COUPLING:-canonical}"
export GRTRESNA_FULL_Z=1

export OBJECTIVE_MODE="${OBJECTIVE_MODE:-general_ftl}"
export DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-ftl_lifetime}"

# Static shell; exotic wedge in default scalar shell search space (no scalar_sign pin).
export PIN_DIMS="${PIN_DIMS:-grtresna_shell_static=1}"

export STOP_TIME="${STOP_TIME:-16.0}"
export PLOT_INTERVAL="${PLOT_INTERVAL:-320}"

export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-1}"
export GRTECLYN_SPLASH_EARLY_TERM=0

export ITERATIONS="${ITERATIONS:-80}"
export POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-3e-2}"
export USE_PIPELINE="${USE_PIPELINE:-1}"
export MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-5}"

exec bash "${SCRIPT_DIR}/../qd/run.sh"
