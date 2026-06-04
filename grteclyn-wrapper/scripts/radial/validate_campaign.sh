#!/usr/bin/env bash
# One-off validation campaign: paper seeds + random nonspherical candidates,
# launched in parallel across GPUs 0-7 via the RadialRecipe smoke pipeline.
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"
SMOKE="${SCRIPT_DIR}/run_radialrecipe_gpu_smoke.sh"
STOP_TIME="${STOP_TIME:-16.0}"
N_FULL="${N_FULL:-64}"
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/validation_${RUN_STAMP}}"
mkdir -p "${RUNS_DIR}"

# label : kind(seed|nonsph) : id
ENTRIES=(
  "flat:seed:flat_minkowski"
  "ellis:seed:ellis_bronnikov"
  "alcubierre:seed:alcubierre_warp"
  "schwarzschild:seed:schwarzschild_puncture"
  "dipole:nonsph:dipole_lopsided_000"
  "quadrupole:nonsph:quadrupole_bubble_001"
  "mixed:nonsph:mixed_modes_002"
  "rand_angular:nonsph:random_angular_004"
)

echo "== Validation campaign =="
echo "stop_time=${STOP_TIME} N_full=${N_FULL} stamp=${RUN_STAMP}"
echo "runs_dir=${RUNS_DIR}"
echo

PIDS=()
GPU=0
for entry in "${ENTRIES[@]}"; do
  IFS=":" read -r label kind cid <<< "${entry}"
  log="${RUNS_DIR}/${label}.log"
  echo "GPU ${GPU}: ${label} (${kind}=${cid}) -> ${log}"
  if [[ "${kind}" == "seed" ]]; then
    ID_ENV="SEED_NAME=${cid}"
  else
    ID_ENV="NONSPHERICAL_ID=${cid}"
  fi
  (
    env ${ID_ENV} \
      BUILD=0 STOP_TIME="${STOP_TIME}" N_FULL="${N_FULL}" \
      CONSUME_PLOTFILES=1 CONSUMER_DELETE=1 CONSUMER_RADII="4 8" \
      FTL_L=8.0 ENABLE_FTL_SCORING=1 \
      RUN_STAMP="${RUN_STAMP}" RUNS_DIR="${RUNS_DIR}" \
      GPU_NAME="${label}" \
      CUDA_VISIBLE_DEVICES_OVERRIDE="${GPU}" \
      bash "${SMOKE}" > "${log}" 2>&1
    echo "EXIT=$? ${label}" >> "${log}"
  ) &
  PIDS+=($!)
  GPU=$((GPU + 1))
done

for pid in "${PIDS[@]}"; do wait "${pid}"; done
echo "ALL_DONE runs_dir=${RUNS_DIR}"
