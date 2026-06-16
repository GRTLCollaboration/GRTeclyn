#!/usr/bin/env bash
# Shared HQ promotion defaults (stage 2) for campaigns/hq/run_batch.sh and replay_eval.py.
# Source after env.sh; expects GRTECLYN_ROOT and WRAPPER_ROOT.
set -euo pipefail

RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_promote}"
N_FULL="${N_FULL:-256}"
L_FULL="${L_FULL:-128}"
GRTRESNA_DOMAIN_L="${GRTRESNA_DOMAIN_L:-${L_FULL}}"
MAX_LEVEL="${MAX_LEVEL:-3}"
REGRID_THRESHOLD="${REGRID_THRESHOLD:-0.02}"
STOP_TIME="${STOP_TIME:-30}"
PLOT_INTERVAL="${PLOT_INTERVAL:-24}"
GRTRESNA_MAX_HAM_PCT="${GRTRESNA_MAX_HAM_PCT:-10}"
GRTRESNA_MAX_MOM_PCT="${GRTRESNA_MAX_MOM_PCT:-10}"
GRTRESNA_TIMEOUT="${GRTRESNA_TIMEOUT:-7200}"
GRTRESNA_ITERATIONS="${GRTRESNA_ITERATIONS:-30}"
GRTRESNA_RANKS="${GRTRESNA_RANKS:-8}"
# Keep >=3 plotfiles for evolved/geodesic FTL + effective EC scoring.
CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-3}"
FORCE="${FORCE:-0}"
INCLUDE_STOP_TIME_IN_NAME="${INCLUDE_STOP_TIME_IN_NAME:-1}"

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-1}"
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-lump_activity scalar_activity phi_lump_sum Pi_lump_sum chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_FRAMES_ZOOM="${FRAMES_ZOOM:-none}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-lump_activity scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${PROJECTION_METHOD:-mip}"
export GRTECLYN_EVOLVING_GEODESIC="${GRTECLYN_EVOLVING_GEODESIC:-1}"
export GRTECLYN_EVOLVING_GEODESIC_MODE="${GRTECLYN_EVOLVING_GEODESIC_MODE:-hq}"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
    PYTHON_BIN="${WRAPPER_ROOT}/.venv/bin/python"
  elif command -v uv >/dev/null 2>&1; then
    PYTHON_BIN="uv run --project ${WRAPPER_ROOT} python"
  else
    PYTHON_BIN="python"
  fi
fi

promote_build_name() {
  local eval_padded="$1"
  local base="l${L_FULL}n${N_FULL}"
  if [[ "${INCLUDE_STOP_TIME_IN_NAME}" == "1" ]]; then
    base="${base}t${STOP_TIME}"
  fi
  if [[ -n "${NAME_PREFIX:-}" ]]; then
    echo "${base}_${NAME_PREFIX}_qd_eval${eval_padded}"
  else
    echo "${base}_qd_eval${eval_padded}"
  fi
}
