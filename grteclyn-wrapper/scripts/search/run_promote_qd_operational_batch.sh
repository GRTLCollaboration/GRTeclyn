#!/usr/bin/env bash
# Parallel high-res promotion of operational-FTL QD elites.
# Fresh GRTresna solve per candidate at N>L (finer dx = L/N), one GPU each.
#
# Usage:
#   bash run_promote_qd_operational_batch.sh
#
# Env overrides:
#   QD_RUN=.../qd_20260605T155951Z
#   N_FULL=256 L_FULL=128 STOP_TIME=50 PLOT_INTERVAL=24
#   GRTRESNA_MAX_HAM_PCT=10 GRTRESNA_MAX_MOM_PCT=10
set -euo pipefail

SEARCH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SEARCH_DIR}/../lib/env.sh"

QD_RUN="${QD_RUN:-${GRTECLYN_ROOT}/runs/grtresna_qd/qd_20260605T155951Z}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_promote}"
N_FULL="${N_FULL:-256}"
L_FULL="${L_FULL:-128}"
GRTRESNA_DOMAIN_L="${GRTRESNA_DOMAIN_L:-${L_FULL}}"
MAX_LEVEL="${MAX_LEVEL:-3}"
REGRID_THRESHOLD="${REGRID_THRESHOLD:-0.02}"
STOP_TIME="${STOP_TIME:-50}"
PLOT_INTERVAL="${PLOT_INTERVAL:-24}"
GRTRESNA_MAX_HAM_PCT="${GRTRESNA_MAX_HAM_PCT:-10}"
GRTRESNA_MAX_MOM_PCT="${GRTRESNA_MAX_MOM_PCT:-10}"
GRTRESNA_TIMEOUT="${GRTRESNA_TIMEOUT:-7200}"

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-phi Pi scalar_activity chi chi_minus_1 local_speed shift1 rho_req}"
# Full-domain slice frames (no --frames-zoom crop). Set FRAMES_ZOOM=<code units> to override.
export GRTECLYN_FRAMES_ZOOM="${FRAMES_ZOOM:-none}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${PROJECTION_METHOD:-mip}"
CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-2}"

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

mkdir -p "${RUNS_DIR}"

# eval_id gpu_id (operational_ftl >= 0.03, search score order)
CANDIDATES=(
  "011 0"
  "121 1"
  "106 2"
  "077 3"
  "117 4"
  "094 5"
  "058 6"
  "016 7"
)

echo "== l${L_FULL}n${N_FULL} operational FTL promotion batch =="
echo "QD source     : ${QD_RUN}"
echo "Runs dir      : ${RUNS_DIR}"
echo "Domain        : L=${L_FULL} N=${N_FULL} max_level=${MAX_LEVEL} t=${STOP_TIME}"
echo "GRTresna gate : Ham<=${GRTRESNA_MAX_HAM_PCT}% Mom<=${GRTRESNA_MAX_MOM_PCT}%"
echo "Plotting      : consume+delete, frames_zoom=${GRTECLYN_FRAMES_ZOOM}, projections=${GRTECLYN_PROJECTION_FIELDS}"
echo

for entry in "${CANDIDATES[@]}"; do
  read -r EVAL_ID GPU_ID <<< "${entry}"
  NAME="l${L_FULL}n${N_FULL}_qd_eval${EVAL_ID}"
  SOURCE="${QD_RUN}/eval_$(printf '%06d' "$((10#${EVAL_ID}))")"
  LOG="${RUNS_DIR}/${NAME}.log"
  if [[ ! -d "${SOURCE}" ]]; then
    echo "[skip] missing ${SOURCE}"
    continue
  fi
  echo "[launch] ${NAME} GPU=${GPU_ID} source=eval_${EVAL_ID} log=${LOG}"
  # shellcheck disable=SC2086
  nohup ${PYTHON_BIN} "${SEARCH_DIR}/replay_grtresna_eval.py" \
    "${SOURCE}" \
    --name "${NAME}" \
    --runs-dir "${RUNS_DIR}" \
    --gpu "${GPU_ID}" \
    --n-full "${N_FULL}" \
    --l-full "${L_FULL}" \
    --grtresna-domain-l "${GRTRESNA_DOMAIN_L}" \
    --max-level "${MAX_LEVEL}" \
    --regrid-threshold "${REGRID_THRESHOLD}" \
    --stop-time "${STOP_TIME}" \
    --plot-interval "${PLOT_INTERVAL}" \
    --grtresna-ranks 8 \
    --grtresna-timeout "${GRTRESNA_TIMEOUT}" \
    --grtresna-max-ham-pct "${GRTRESNA_MAX_HAM_PCT}" \
    --grtresna-max-mom-pct "${GRTRESNA_MAX_MOM_PCT}" \
    --consumer-keep-last "${CONSUMER_KEEP_LAST}" \
    > "${LOG}" 2>&1 &
  echo "  pid=$!"
done

echo
echo "Launched ${#CANDIDATES[@]} promotions (background). Monitor:"
echo "  tail -f ${RUNS_DIR}/l${L_FULL}n${N_FULL}_qd_eval*.log"
