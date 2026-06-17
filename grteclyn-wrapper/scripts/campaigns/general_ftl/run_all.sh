#!/usr/bin/env bash
# v20 general_ftl orchestrator — wormhole / ring / spin MAP-Elites (stage 0).
#
# Launches scripts/campaigns/qd/run.sh with matter-class pins from lib/general_ftl_pins.sh.
# All physics/grid/pipeline knobs come from lib/search_common.sh via qd/run.sh.
#
# Usage:
#   MODE=par bash scripts/campaigns/general_ftl/run_all.sh
#
# Single branch:
#   BRANCH=wormhole GPU_IDS="0 1 2 3" QD_TARGET_EVALS=80 bash scripts/campaigns/general_ftl/run_all.sh
#
# Single-GPU multi-slot test (2 concurrent evolutions on GPU 0):
#   BRANCH=wormhole PIPELINE_MONITOR=1 \
#     QD_TARGET_EVALS=8 GPU_IDS="0" GPU_SLOTS_PER_DEVICE=2 BATCH_SIZE=4 \
#     STOP_TIME=4.0 PLOT_INTERVAL=40 QD_ITERATIONS=4 SKIP_QD_PREFLIGHT_TESTS=1 \
#     bash scripts/campaigns/general_ftl/run_all.sh
#
# v21 production throughput (5 evols/GPU, 8 GPUs, pipelined):
#   BRANCH=wormhole PIPELINE_MONITOR=1 QD_NAME=general_ftl_wormhole_v21 \
#     QD_TARGET_EVALS=80 GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=5 BATCH_SIZE=40 \
#     bash scripts/campaigns/general_ftl/run_all.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGNS_LIB="$(cd -- "${SCRIPT_DIR}/../lib" && pwd)"
QD_RUN="${SCRIPT_DIR}/../qd/run.sh"
# shellcheck source=../lib/bootstrap.sh
source "${CAMPAIGNS_LIB}/bootstrap.sh"
_campaign_bootstrap "${SCRIPT_DIR}"
# shellcheck source=../lib/general_ftl_pins.sh
source "${CAMPAIGNS_LIB}/general_ftl_pins.sh"
# shellcheck source=../lib/pipeline_monitor.sh
source "${CAMPAIGNS_LIB}/pipeline_monitor.sh"

export OBJECTIVE_MODE=general_ftl
export DESCRIPTOR_MODE=ftl_lifetime
export GRTECLYN_GEO_DIRECTIONS="x y z"
export GRTRESNA_ANSATZ=shell
export STOP_TIME="${STOP_TIME:-16.0}"
export GRTECLYN_FRAMES=0
export GRTECLYN_EVOLVING_GEODESIC=1
export GRTECLYN_EVOLVING_GEODESIC_MODE="${GRTECLYN_EVOLVING_GEODESIC_MODE:-search}"

WORMHOLE_PINS="$(ftl_general_ftl_wormhole_pins)"
RING_PINS="$(ftl_general_ftl_ring_pins)"
SPIN_PINS="$(ftl_general_ftl_spin_pins)"

run_one() {  # name  pins  gpu_ids
  local name="$1" pins="$2" gpus="$3"
  local qd_name="${QD_NAME:-general_ftl_${name}}"
  local n_gpus batch_default
  n_gpus=$(wc -w <<< "${gpus}")
  batch_default=$((n_gpus * ${GPU_SLOTS_PER_DEVICE:-1}))

  if [[ "${PIPELINE_MONITOR:-0}" == "1" ]]; then
    ftl_pipeline_monitor_begin "${qd_name}" "${gpus}"
    RUNS_DIR="${GRTECLYN_ROOT}/runs/grtresna_qd" \
      QD_NAME="${qd_name}" \
      PIN_DIMS="${pins}" \
      GPU_IDS="${gpus}" \
      BATCH_SIZE="${BATCH_SIZE:-${batch_default}}" \
      QD_ITERATIONS="${QD_ITERATIONS:-30}" \
      bash "${QD_RUN}" 2>&1 | tee "${FTL_PIPELINE_LOG}"
    ftl_pipeline_monitor_end
  else
    RUNS_DIR="${GRTECLYN_ROOT}/runs/grtresna_qd" \
      QD_NAME="${qd_name}" \
      PIN_DIMS="${pins}" \
      GPU_IDS="${gpus}" \
      BATCH_SIZE="${BATCH_SIZE:-${batch_default}}" \
      QD_ITERATIONS="${QD_ITERATIONS:-30}" \
      bash "${QD_RUN}"
  fi
}

MODE="${MODE:-seq}"
BRANCH="${BRANCH:-all}"

if [[ "${MODE}" == "par" && "${BRANCH}" != "all" ]]; then
  echo "BRANCH=${BRANCH} with MODE=par is unsupported; use MODE=seq or BRANCH=all." >&2
  exit 2
fi

if [[ "${PIPELINE_MONITOR:-0}" == "1" ]]; then
  _campaign_resolve_python
fi

if [[ "${MODE}" == "par" ]]; then
  export CLUSTER_CPU_FRACTION="${CLUSTER_CPU_FRACTION:-0.10}"
  export PIPELINE_CPU_SHARE="${PIPELINE_CPU_SHARE:-0.333}"
  if [[ "${SKIP_QD_PREFLIGHT_TESTS:-0}" != "1" ]]; then
    _campaign_resolve_python
    ftl_search_common_preflight_tests
  fi
  export SKIP_QD_PREFLIGHT_TESTS=1
  run_one wormhole "${WORMHOLE_PINS}" "0 1 2" &
  run_one ring     "${RING_PINS}"     "3 4 5" &
  run_one spin     "${SPIN_PINS}"     "6 7"   &
  wait
elif [[ "${BRANCH}" == "wormhole" ]]; then
  export CLUSTER_CPU_FRACTION="${CLUSTER_CPU_FRACTION:-0.30}"
  export PIPELINE_CPU_SHARE="${PIPELINE_CPU_SHARE:-1.0}"
  run_one wormhole "${WORMHOLE_PINS}" "${GPU_IDS:-0 1 2 3 4 5 6 7}"
elif [[ "${BRANCH}" == "ring" ]]; then
  export CLUSTER_CPU_FRACTION="${CLUSTER_CPU_FRACTION:-0.30}"
  export PIPELINE_CPU_SHARE="${PIPELINE_CPU_SHARE:-1.0}"
  run_one ring "${RING_PINS}" "${GPU_IDS:-0 1 2 3 4 5 6 7}"
elif [[ "${BRANCH}" == "spin" ]]; then
  export CLUSTER_CPU_FRACTION="${CLUSTER_CPU_FRACTION:-0.30}"
  export PIPELINE_CPU_SHARE="${PIPELINE_CPU_SHARE:-1.0}"
  run_one spin "${SPIN_PINS}" "${GPU_IDS:-0 1 2 3 4 5 6 7}"
else
  export CLUSTER_CPU_FRACTION="${CLUSTER_CPU_FRACTION:-0.30}"
  export PIPELINE_CPU_SHARE="${PIPELINE_CPU_SHARE:-1.0}"
  run_one wormhole "${WORMHOLE_PINS}" "0 1 2 3 4 5 6 7"
  run_one ring     "${RING_PINS}"     "0 1 2 3 4 5 6 7"
  run_one spin     "${SPIN_PINS}"     "0 1 2 3 4 5 6 7"
fi
