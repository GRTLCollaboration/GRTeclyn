#!/usr/bin/env bash
# Promote the top 8 center_fixed_search QD candidates with fresh HQ GRTresna solves.
#
# Defaults:
#   - source: qd_20260609T094634Z top-8 successful trajectory candidates
#   - grid:   L=128, N=256 (dx=0.5)
#   - time:   t=30
#   - output: runs/grtresna_promote/l128n256t30_center_fixed_qd_eval*/
#
# Frame extraction is on during the main GPU evolution. Processed plotfiles are
# deleted by the consumer (`--consumer-keep-last 0`) so only frames/small_data
# and compact metadata are retained.
set -euo pipefail

SEARCH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SEARCH_DIR}/../lib/env.sh"

QD_RUN="${QD_RUN:-${GRTECLYN_ROOT}/runs/grtresna_qd/center_fixed_search/qd_20260609T094634Z}"
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
CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-0}"
FORCE="${FORCE:-0}"

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-lump_activity scalar_activity phi Pi chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_FRAMES_ZOOM="${FRAMES_ZOOM:-none}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-lump_activity scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${PROJECTION_METHOD:-mip}"

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

# Top 8 successful GPU candidates from qd_20260609T094634Z trajectory score order.
CANDIDATES=(
  "324 0"
  "114 1"
  "146 2"
  "169 3"
  "358 4"
  "314 5"
  "136 6"
  "228 7"
)

echo "== center_fixed_search top-8 HQ promotion =="
echo "QD source     : ${QD_RUN}"
echo "Runs dir      : ${RUNS_DIR}"
echo "Domain        : L=${L_FULL} N=${N_FULL} dx=$(python3 - <<'PY' "${L_FULL}" "${N_FULL}"
import sys
print(float(sys.argv[1]) / int(sys.argv[2]))
PY
) max_level=${MAX_LEVEL} t=${STOP_TIME}"
echo "GRTresna gate : Ham<=${GRTRESNA_MAX_HAM_PCT}% Mom<=${GRTRESNA_MAX_MOM_PCT}% iterations=${GRTRESNA_ITERATIONS}"
echo "Plotting      : frames on the fly, delete plotfiles, keep_last=${CONSUMER_KEEP_LAST}"
echo "Frames fields : ${GRTECLYN_FRAMES_FIELDS}"
echo "Projections   : ${GRTECLYN_PROJECTION_FIELDS} axes=${GRTECLYN_PROJECTION_AXES}"
echo

launched=0
for entry in "${CANDIDATES[@]}"; do
  read -r EVAL_ID GPU_ID <<< "${entry}"
  eval_num="$((10#${EVAL_ID}))"
  eval_padded="$(printf "%06d" "${eval_num}")"
  name="l${L_FULL}n${N_FULL}t${STOP_TIME}_center_fixed_qd_eval${eval_padded}"
  source="${QD_RUN}/eval_${eval_padded}"
  out="${RUNS_DIR}/${name}"
  log="${RUNS_DIR}/${name}.log"

  if [[ ! -d "${source}" ]]; then
    echo "[skip] missing source ${source}"
    continue
  fi
  if [[ -e "${out}" && "${FORCE}" != "1" ]]; then
    echo "[skip] output exists ${out} (set FORCE=1 to replace manually after cleanup)"
    continue
  fi

  echo "[launch] ${name} GPU=${GPU_ID} source=eval_${eval_padded} log=${log}"
  # shellcheck disable=SC2086
  nohup ${PYTHON_BIN} "${SEARCH_DIR}/replay_grtresna_eval.py" \
    "${source}" \
    --name "${name}" \
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
    --grtresna-iterations "${GRTRESNA_ITERATIONS}" \
    --grtresna-timeout "${GRTRESNA_TIMEOUT}" \
    --grtresna-max-ham-pct "${GRTRESNA_MAX_HAM_PCT}" \
    --grtresna-max-mom-pct "${GRTRESNA_MAX_MOM_PCT}" \
    --consumer-keep-last "${CONSUMER_KEEP_LAST}" \
    > "${log}" 2>&1 &
  echo "  pid=$!"
  launched=$((launched + 1))
done

echo
echo "Launched ${launched} promotions (background). Monitor:"
echo "  tail -f ${RUNS_DIR}/l${L_FULL}n${N_FULL}t${STOP_TIME}_center_fixed_qd_eval*.log"
