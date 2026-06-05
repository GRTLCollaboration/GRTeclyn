#!/usr/bin/env bash
# High-resolution long validation of GRTresna shell QD winner eval_000057.
# Re-solves GRTresna at promoted gridinit resolution, then evolves on GPU with
# streaming frames (plotfiles deleted, keep-last-2).
#
# Usage:
#   ./run_tier2_grtresna_qd_eval057.sh [GPU] [NAME]
#
# Env overrides:
#   N_FULL=128 STOP_TIME=16 SOURCE_EVAL=.../eval_000057
#   L_FULL=256  (larger box; defaults to N_FULL when unset)
#   GRIDINIT_SOURCE=.../val16hq2_qd_eval057/initial_data.gridinit  (skip GRTresna)
set -euo pipefail

SEARCH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SEARCH_DIR}/../lib/env.sh"

GPU="${1:-0}"
NAME="${2:-val16hq_qd_eval057}"
N_FULL="${N_FULL:-128}"
L_FULL="${L_FULL:-${N_FULL}}"
GRTRESNA_DOMAIN_L="${GRTRESNA_DOMAIN_L:-${L_FULL}}"
MAX_LEVEL="${MAX_LEVEL:-3}"
REGRID_THRESHOLD="${REGRID_THRESHOLD:-0.02}"
STOP_TIME="${STOP_TIME:-16}"
PLOT_INTERVAL="${PLOT_INTERVAL:-48}"
SOURCE_EVAL="${SOURCE_EVAL:-${GRTECLYN_ROOT}/runs/grtresna_qd/qd_20260605T062448Z/eval_000057}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_promote}"
LOG_PATH="${RUNS_DIR}/${NAME}.log"

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-phi Pi scalar_activity chi chi_minus_1 local_speed shift1 rho_req}"
# Match slice-frame width to the evolved box (plot_consumer defaults to L_full).
export GRTECLYN_FRAMES_ZOOM="${FRAMES_ZOOM:-}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${GRTECLYN_PROJECTION_METHOD:-mip}"

mkdir -p "${RUNS_DIR}"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if command -v uv >/dev/null 2>&1; then
    PYTHON_BIN="uv run --project ${WRAPPER_ROOT} python"
  elif [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
    PYTHON_BIN="${WRAPPER_ROOT}/.venv/bin/python"
  else
    PYTHON_BIN="python"
  fi
fi

if [[ -z "${GRTECLYN_FRAMES_ZOOM}" ]]; then
  export GRTECLYN_FRAMES_ZOOM="${L_FULL}"
fi
GRIDINIT_ARGS=()
if [[ -n "${GRIDINIT_SOURCE:-}" ]]; then
  GRIDINIT_ARGS=(--gridinit "${GRIDINIT_SOURCE}")
  echo "[launch] ${NAME} on GPU ${GPU} (L=${L_FULL}, N=${N_FULL}, max_level=${MAX_LEVEL}, t=${STOP_TIME}, gpu-only, gridinit=${GRIDINIT_SOURCE})"
else
  echo "[launch] ${NAME} on GPU ${GPU} (L=${L_FULL}, N=${N_FULL}, max_level=${MAX_LEVEL}, t=${STOP_TIME}, frames_zoom=${GRTECLYN_FRAMES_ZOOM}, source=${SOURCE_EVAL})"
fi
echo "[log] ${LOG_PATH}"

# shellcheck disable=SC2086
${PYTHON_BIN} "${SEARCH_DIR}/replay_grtresna_eval.py" \
  "${SOURCE_EVAL}" \
  --name "${NAME}" \
  --runs-dir "${RUNS_DIR}" \
  --gpu "${GPU}" \
  --n-full "${N_FULL}" \
  --l-full "${L_FULL}" \
  --grtresna-domain-l "${GRTRESNA_DOMAIN_L}" \
  --max-level "${MAX_LEVEL}" \
  --regrid-threshold "${REGRID_THRESHOLD}" \
  --stop-time "${STOP_TIME}" \
  --plot-interval "${PLOT_INTERVAL}" \
  "${GRIDINIT_ARGS[@]}" \
  > "${LOG_PATH}" 2>&1

rc=$?
echo "[done] ${NAME} (exit ${rc})"
if [[ -d "${RUNS_DIR}/${NAME}" ]]; then
  du -sh "${RUNS_DIR}/${NAME}"
fi
exit "${rc}"
