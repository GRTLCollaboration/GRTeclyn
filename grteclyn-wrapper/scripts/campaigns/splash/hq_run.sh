#!/usr/bin/env bash
# HQ promotion for spacetime_splash QD elites (critical_collapse, no 4D geodesic).
#
# Example (top-1 eval 020, t=30):
#   SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_qd/spacetime_splash_v13" \
#   CANDIDATES="20 0" \
#   NAME_PREFIX=spacetime_splash_v13 \
#   bash scripts/campaigns/splash/hq_run.sh

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=shell
export GRTRESNA_MATTER_COUPLING=canonical
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-critical_collapse}"
export SPLASH_MODE="${SPLASH_MODE:-discovery}"
export GRTECLYN_EVOLVING_GEODESIC=0
export GRTECLYN_CENTRAL_TIMESERIES=1
export GRTECLYN_CENTRAL_BALL=1
export GRTECLYN_CENTRAL_RADIAL=1
export GRTECLYN_INCREMENTAL_SCORE=1
# HQ runs to stop_time; do not cut short on matter dispersion.
export GRTECLYN_SPLASH_EARLY_TERM="${GRTECLYN_SPLASH_EARLY_TERM:-0}"
export FRAMES_FIELDS="${FRAMES_FIELDS:-chi K rho_req local_speed weyl4}"
export PROJECTION_FIELDS="${PROJECTION_FIELDS:-chi rho_req}"
# Match search-stage AMR guard (0.02 can refine the full domain on splash ICs).
export REGRID_THRESHOLD="${REGRID_THRESHOLD:-0.1}"
export STOP_TIME="${STOP_TIME:-30}"

exec bash "${SCRIPT_DIR}/../hq/run_batch.sh"
