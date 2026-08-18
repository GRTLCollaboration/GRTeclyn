#!/usr/bin/env bash
# Campaign-local config for bicomplex CMA-ES → HQ / Richardson.
# Sourced by freeze.sh / run.sh / launch.sh (not executed alone).
# shellcheck shell=bash

CAMPAIGN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROMOTE_ROOT="$(cd -- "${CAMPAIGN_DIR}/.." && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${PROMOTE_ROOT}/.." && pwd)"
SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
# shellcheck source=../../../lib/env.sh
source "${SCRIPTS_ROOT}/lib/env.sh"

export CAMPAIGN_NAME="bicomplex_cmaes_v1"
export MANIFEST="${MANIFEST:-${CAMPAIGN_DIR}/manifest.json}"
export VALIDATION_LAUNCH_LOG_DIR="${VALIDATION_LAUNCH_LOG_DIR:-${GRTECLYN_ROOT}/research/neuralspacetime/validation/bicomplex_cmaes/launches}"

export LIVE_RUN="${LIVE_RUN:-${GRTECLYN_ROOT}/runs/neuralspacetime/search/cma_es/qball_traj_bicomplex_cmaes_v1}"
export FREEZE_ROOT="${FREEZE_ROOT:-${GRTECLYN_ROOT}/runs/neuralspacetime/hq/sources/qball_traj_bicomplex_cmaes_v1}"

# Physics (must match QD/CMA search).
export RL_PUMP_STOP_TIME="${RL_PUMP_STOP_TIME:-4}"
export SCORE_PUMP_ENERGY_WEIGHT="${SCORE_PUMP_ENERGY_WEIGHT:-40}"
export SCORE_EXOTIC_PENALTY_WEIGHT="${SCORE_EXOTIC_PENALTY_WEIGHT:-0.2}"
export GRTRESNA_ALLOW_SIGN_MISMATCH="${GRTRESNA_ALLOW_SIGN_MISMATCH:-0}"
export GRTRESNA_MATTER_SECTOR="${GRTRESNA_MATTER_SECTOR:-boson_star}"
export GRTRESNA_MATTER_MODEL="${GRTRESNA_MATTER_MODEL:-grtresna_bicomplex_scalar}"
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-general_ftl}"

# RF movie defaults (RM/RC keep frames=0 via manifest).
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

bicomplex_resolve_source() {
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
