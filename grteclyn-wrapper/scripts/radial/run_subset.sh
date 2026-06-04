#!/usr/bin/env bash
set -uo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"
SMOKE="${SCRIPT_DIR}/run_radialrecipe_gpu_smoke.sh"
STOP_TIME="${STOP_TIME:-2.0}"
N_FULL="${N_FULL:-64}"
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUNS_DIR="${RUNS_DIR:?set RUNS_DIR}"
CONSUME_PLOTFILES="${CONSUME_PLOTFILES:-1}"
PLOT_INTERVAL_ENV="${PLOT_INTERVAL_ENV:-}"
CUDA_START="${CUDA_START:-0}"
mkdir -p "${RUNS_DIR}"
# ENTRIES space-separated: label:kind(seed|nonsph):id
read -r -a ENTRY_ARR <<< "${ENTRIES}"
echo "== subset: stop_time=${STOP_TIME} consume=${CONSUME_PLOTFILES} stamp=${RUN_STAMP} dir=${RUNS_DIR} =="
PIDS=(); GPU="${CUDA_START}"
for entry in "${ENTRY_ARR[@]}"; do
  IFS=":" read -r label kind cid <<< "${entry}"
  log="${RUNS_DIR}/${label}.log"
  if [[ "${kind}" == "seed" ]]; then ID_ENV="SEED_NAME=${cid}"; else ID_ENV="NONSPHERICAL_ID=${cid}"; fi
  echo "GPU ${GPU}: ${label} (${kind}=${cid})"
  (
    env ${ID_ENV} BUILD=0 STOP_TIME="${STOP_TIME}" N_FULL="${N_FULL}" \
      CONSUME_PLOTFILES="${CONSUME_PLOTFILES}" CONSUMER_DELETE=1 CONSUMER_RADII="4 8" \
      ${PLOT_INTERVAL_ENV:+PLOT_INTERVAL=${PLOT_INTERVAL_ENV}} \
      FTL_L=8.0 ENABLE_FTL_SCORING=1 \
      RUN_STAMP="${RUN_STAMP}" RUNS_DIR="${RUNS_DIR}" GPU_NAME="${label}" \
      CUDA_VISIBLE_DEVICES_OVERRIDE="${GPU}" \
      bash "${SMOKE}" > "${log}" 2>&1
    echo "EXIT=$? ${label}" >> "${log}"
  ) &
  PIDS+=($!); GPU=$((GPU + 1))
done
for pid in "${PIDS[@]}"; do wait "${pid}"; done
echo "SUBSET_DONE ${RUNS_DIR}"
