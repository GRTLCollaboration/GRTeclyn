#!/usr/bin/env bash
# Freeze → RM (1 GPU). Optional Richardson companions RC (1 GPU) + RF (2 GPUs).
#
#   bash scripts/campaigns/promote/bicomplex_cmaes_v1/launch.sh
#   LAUNCH_LADDER=1 SKIP_FREEZE=1 bash scripts/campaigns/promote/bicomplex_cmaes_v1/launch.sh
set -euo pipefail

CAMPAIGN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=campaign.env.sh
source "${CAMPAIGN_DIR}/campaign.env.sh"

DRY_RUN="${DRY_RUN:-0}"
LAUNCH_LADDER="${LAUNCH_LADDER:-0}"
FORCE="${FORCE:-0}"
FOREGROUND="${FOREGROUND:-0}"
SKIP_FREEZE="${SKIP_FREEZE:-0}"
# shellcheck disable=SC2206
GPU_IDS=(${GPU_IDS:-0 1 2 3 4 5 6 7})

RUN="${CAMPAIGN_DIR}/run.sh"

if [[ "${SKIP_FREEZE}" != "1" ]]; then
  echo "== Freeze champion =="
  bash "${CAMPAIGN_DIR}/freeze.sh"
  echo
fi

bicomplex_resolve_source
export SOURCE_RUN SOURCE_EVAL_ID FORCE FOREGROUND DRY_RUN

rm_gpu="${GPU_IDS[0]:-0}"
echo "== BCMA-RM (HQ / Richardson mid, 1 GPU=${rm_gpu}) =="
GPU_ID="${rm_gpu}" bash "${RUN}" BCMA-RM

if [[ "${LAUNCH_LADDER}" != "1" ]]; then
  echo
  echo "After RM gate, Richardson companions:"
  echo "  LAUNCH_LADDER=1 SKIP_FREEZE=1 bash $0"
  echo "  GPU_ID=1 bash ${RUN} BCMA-RC    # 1 GPU"
  echo "  GPU_ID=2 bash ${RUN} BCMA-RF    # 3 GPUs → 2,3,4"
  exit 0
fi

rc_gpu="${GPU_IDS[1]:-1}"
rf_start="${GPU_IDS[2]:-2}"
echo
echo "== Richardson companions =="
echo "-- BCMA-RC GPU=${rc_gpu} (1 GPU)"
GPU_ID="${rc_gpu}" bash "${RUN}" BCMA-RC
echo "-- BCMA-RF GPU=${rf_start} (3 GPUs, auto-expanded by run_matrix.sh)"
GPU_ID="${rf_start}" bash "${RUN}" BCMA-RF

echo
echo "Series = RC + RM + RF under ${CAMPAIGN_NAME}"
echo "  ls -lt ${GRTECLYN_ROOT}/runs/neuralspacetime/hq/bcma_*_hq_eval* | head"
