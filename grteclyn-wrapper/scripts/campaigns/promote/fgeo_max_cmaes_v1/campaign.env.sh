#!/usr/bin/env bash
# Campaign-local config for the f_geo_max CMA-ES → HQ / Richardson matrix
# (FMAX-*). Cloned from bicomplex_cmaes_v1 per GPU_RUN_PLAN.md §12.2 with the
# post-bicomplex conventions baked in — every value below that differs from
# the template is a convention, not a preference; see the plan section cited
# on each line. Sourced by freeze.sh / run.sh / launch.sh (not executed
# alone).
# shellcheck shell=bash

CAMPAIGN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROMOTE_ROOT="$(cd -- "${CAMPAIGN_DIR}/.." && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${PROMOTE_ROOT}/.." && pwd)"
SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
# shellcheck source=../../../lib/env.sh
source "${SCRIPTS_ROOT}/lib/env.sh"

export CAMPAIGN_NAME="fgeo_max_cmaes_v1"
export MANIFEST="${MANIFEST:-${CAMPAIGN_DIR}/manifest.json}"
export VALIDATION_LAUNCH_LOG_DIR="${VALIDATION_LAUNCH_LOG_DIR:-${GRTECLYN_ROOT}/research/neuralspacetime/validation/fgeo_max_cmaes/launches}"

export LIVE_RUN="${LIVE_RUN:-${GRTECLYN_ROOT}/runs/neuralspacetime/search/cma_es/qball_traj_fgeo_max_cmaes_v1}"
export FREEZE_ROOT="${FREEZE_ROOT:-${GRTECLYN_ROOT}/runs/neuralspacetime/hq/sources/qball_traj_fgeo_max_cmaes_v1}"

# Physics (must match the f_geo_max search, §12.3).
# Pump on for the ENTIRE run (§12.1): -1 writes no rl_pump_stop_time key and
# the evolution default is never-stop. The manifest carries no pump key
# either — the pump-convention validator refuses it if one sneaks back in.
export RL_PUMP_STOP_TIME="${RL_PUMP_STOP_TIME:--1}"
# A negative pump value erases the scorer's fallback emission floor — pin it
# (§12.1; the search used 4 via run.sh's default).
export GEODESIC_EMIT_MIN_TIME="${GEODESIC_EMIT_MIN_TIME:-4}"
export SCORE_PUMP_ENERGY_WEIGHT="${SCORE_PUMP_ENERGY_WEIGHT:-40}"
# f_geo_max has no exotic penalty (run_fgeo.sh pins 0; §12.3).
export SCORE_EXOTIC_PENALTY_WEIGHT="${SCORE_EXOTIC_PENALTY_WEIGHT:-0}"
export GRTRESNA_ALLOW_SIGN_MISMATCH="${GRTRESNA_ALLOW_SIGN_MISMATCH:-0}"
export GRTRESNA_MATTER_SECTOR="${GRTRESNA_MATTER_SECTOR:-boson_star}"
export GRTRESNA_MATTER_MODEL="${GRTRESNA_MATTER_MODEL:-grtresna_bicomplex_scalar}"
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-f_geo_max}"

# Solve: mpirun segfaults on this node — single-rank singleton only (§12.2;
# promote_common.sh would otherwise default to 8). Gate at the 5% the search
# accepted the champion under, not the promote default 10% (§12.2).
export GRTRESNA_RANKS="${GRTRESNA_RANKS:-1}"
export GRTRESNA_MAX_HAM_PCT="${GRTRESNA_MAX_HAM_PCT:-5}"
export GRTRESNA_MAX_MOM_PCT="${GRTRESNA_MAX_MOM_PCT:-5}"

# Mandatory paper-tier replay env (§12.1) — the promote framework exports
# neither, so they must live here or every cell ships a 65^3 metric stack
# and no free-fall certificate.
export GRTECLYN_METRIC_STACK_N_SPACE="${GRTECLYN_METRIC_STACK_N_SPACE:-257}"
export GRTECLYN_FREEFALL_OBSERVER_TIMING="${GRTECLYN_FREEFALL_OBSERVER_TIMING:-1}"

# Movie defaults (all matrix cells keep frames=0 via manifest).
export FRAMES_FIELDS="${FRAMES_FIELDS:-scalar_activity phi Pi phi_lump0 Pi_lump0 phi_lump1 Pi_lump1 phi_lump2 Pi_lump2 chi chi_minus_1 local_speed shift1 rho_req Weyl4_Re Weyl4_Im Weyl4_Mag}"
export PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity phi}"
export GRTECLYN_FRAMES_ZOOM="${GRTECLYN_FRAMES_ZOOM:-none}"
export GRTECLYN_FRAMES_STABLE_MOVIE="${GRTECLYN_FRAMES_STABLE_MOVIE:-1}"
unset GRTECLYN_FRAMES_AUTO_ZLIM 2>/dev/null || true
export GRTECLYN_FRAMES_DPI="${GRTECLYN_FRAMES_DPI:-110}"
export GRTECLYN_FRAMES_BUFF_CAP="${GRTECLYN_FRAMES_BUFF_CAP:-768}"
export GRTECLYN_FRAMES_SAMPLES_PER_CELL="${GRTECLYN_FRAMES_SAMPLES_PER_CELL:-1}"

export PROMOTE_LIB="${PROMOTE_ROOT}/lib"
export HQ_ENGINE="${CAMPAIGNS_ROOT}/hq"

campaign_resolve_source() {
  if [[ -n "${SOURCE_RUN:-}" ]]; then
    return 0
  fi
  if [[ -f "${FREEZE_ROOT}/CHAMPION.json" ]]; then
    export SOURCE_RUN="${FREEZE_ROOT}"
    if [[ -z "${SOURCE_EVAL_ID:-}" ]]; then
      SOURCE_EVAL_ID="$(python3 -c "import json; print(json.load(open('${FREEZE_ROOT}/CHAMPION.json'))['eval_id'])")"
      export SOURCE_EVAL_ID
    fi
  else
    export SOURCE_RUN="${LIVE_RUN}"
  fi
}
