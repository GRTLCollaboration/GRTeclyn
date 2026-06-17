#!/usr/bin/env bash
# GPU-only replay of each class champion at QD search resolution, reusing
# initial_data.gridinit and writing PNG frames via consume_plotfiles.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
# shellcheck source=../../lib/env.sh
source "${SCRIPTS_ROOT}/lib/env.sh"
# shellcheck source=../lib/promote_common.sh
source "${CAMPAIGNS_ROOT}/lib/promote_common.sh"

export OBJECTIVE_MODE=general_ftl
export GRTECLYN_EVOLVING_GEODESIC=1
export GRTECLYN_EVOLVING_GEODESIC_MODE=search
export GRTECLYN_GEO_DIRECTIONS="x y z"

# QD search grid (same as general_ftl MAP-Elites campaigns)
N_FULL=128
L_FULL=64
MAX_LEVEL=2
STOP_TIME=16.0
PLOT_INTERVAL=320
REGRID_THRESHOLD=0.01
CONSUMER_KEEP_LAST=3
FTL_L=8.0

GPU_IDS="${GPU_IDS:-0 1 2}"
# shellcheck disable=SC2206
GPUS=(${GPU_IDS})
REPLAY="${CAMPAIGNS_ROOT}/hq/replay_eval.py"
LOG_DIR="${GRTECLYN_ROOT}/runs/_logs"
mkdir -p "${LOG_DIR}"

declare -a JOBS=(
  "wormhole 31"
  "ring 43"
  "spin 49"
)

launched=0
for entry in "${JOBS[@]}"; do
  read -r class eval_id <<< "${entry}"
  gpu="${GPUS[$launched]:-${GPUS[0]}}"
  eval_padded="$(printf "%06d" "${eval_id}")"
  campaign="${GRTECLYN_ROOT}/runs/grtresna_qd/general_ftl_${class}"
  source_eval="${campaign}/eval_${eval_padded}"
  gridinit="${source_eval}/initial_data.gridinit"
  name="eval_${eval_padded}_frames"
  log="${LOG_DIR}/general_ftl_${class}_${name}.log"

  if [[ ! -f "${gridinit}" ]]; then
    echo "[skip] missing gridinit ${gridinit}" >&2
    continue
  fi
  if [[ -d "${campaign}/${name}/frames" && "${FORCE}" != "1" ]]; then
    echo "[skip] ${campaign}/${name} exists (FORCE=1 to replace)"
    continue
  fi

  echo "[launch] general_ftl_${class} eval ${eval_id} GPU=${gpu} -> ${campaign}/${name}"
  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    launched=$((launched + 1))
    continue
  fi

  # shellcheck disable=SC2086
  nohup ${PYTHON_BIN} "${REPLAY}" \
    "${source_eval}" \
    --name "${name}" \
    --runs-dir "${campaign}" \
    --gpu "${gpu}" \
    --gridinit "${gridinit}" \
    --n-full "${N_FULL}" \
    --l-full "${L_FULL}" \
    --max-level "${MAX_LEVEL}" \
    --regrid-threshold "${REGRID_THRESHOLD}" \
    --stop-time "${STOP_TIME}" \
    --plot-interval "${PLOT_INTERVAL}" \
    --ftl-L "${FTL_L}" \
    --consumer-keep-last "${CONSUMER_KEEP_LAST}" \
    --evolving-geodesic \
    --objective-mode general_ftl \
    > "${log}" 2>&1 &
  echo "  pid=$! log=${log}"
  launched=$((launched + 1))
done

echo "Launched ${launched} framed QD replays."
