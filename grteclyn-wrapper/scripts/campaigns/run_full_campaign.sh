#!/usr/bin/env bash
# End-to-end general_ftl campaign: QD (200 evals) -> CMA-ES (top-1 seed) -> HQ (top-1).
#
# All artifacts live under a single campaign root:
#   runs/campaigns/<CAMPAIGN_NAME>/
#     qd/           MAP-Elites (pruned to top-3 + FTL champions)
#     cmaes/        CMA-ES refinement (pruned to top-3 + FTL champions)
#     promote/      HQ replay of the CMA-ES champion
#     logs/         stage launch logs
#     campaign_state.json
#     campaign.json
#
# Usage:
#   cd grteclyn-wrapper
#   CAMPAIGN_NAME=general_ftl_wormhole_v22 \
#     bash scripts/campaigns/run_full_campaign.sh
#
# Partial / resume:
#   RESUME=1 STAGE=cmaes CAMPAIGN_NAME=... bash scripts/campaigns/run_full_campaign.sh
#
# Dry-run (no GPU jobs; validates wiring):
#   DRY_RUN=1 SKIP_QD_PREFLIGHT_TESTS=1 SKIP_CMA_PREFLIGHT_TESTS=1 \
#     CAMPAIGN_NAME=smoke_$(date +%s) bash scripts/campaigns/run_full_campaign.sh
set -euo pipefail

ORCH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/bootstrap.sh
source "${ORCH_DIR}/lib/bootstrap.sh"
_campaign_bootstrap "${ORCH_DIR}/qd"
CAMPAIGNS_ROOT="${ORCH_DIR}"
_campaign_resolve_python
# shellcheck source=lib/general_ftl_pins.sh
source "${CAMPAIGNS_ROOT}/lib/general_ftl_pins.sh"
# shellcheck source=lib/campaign_common.sh
source "${CAMPAIGNS_ROOT}/lib/campaign_common.sh"

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"

BRANCH="${BRANCH:-wormhole}"
STAGE="${STAGE:-all}"
RESUME="${RESUME:-0}"
DRY_RUN="${DRY_RUN:-0}"
export DRY_RUN

if [[ "${DRY_RUN}" == "1" ]]; then
  # Use DRY_RUN_* knobs so a pre-set QD_TARGET_EVALS=200 in the shell cannot leak in.
  QD_TARGET_EVALS="${DRY_RUN_QD_TARGET_EVALS:-1}"
  CMAES_TARGET_EVALS="${DRY_RUN_CMAES_TARGET_EVALS:-2}"
  CMAES_MAX_GENERATIONS="${DRY_RUN_CMAES_MAX_GENERATIONS:-1}"
  POPULATION="${DRY_RUN_POPULATION:-2}"
  if [[ "$(wc -w <<< "${GPU_IDS}")" -lt 2 ]]; then
    GPU_IDS="${GPU_IDS} 1"
  fi
  BATCH_SIZE="${BATCH_SIZE:-2}"
fi

QD_TARGET_EVALS="${QD_TARGET_EVALS:-200}"
QD_KEEP_TOP_EVAL_DIRS="${QD_KEEP_TOP_EVAL_DIRS:-3}"
QD_FTL_RETENTION="${QD_FTL_RETENTION:-1}"
CMAES_TARGET_EVALS="${CMAES_TARGET_EVALS:-150}"
CMAES_MAX_GENERATIONS="${CMAES_MAX_GENERATIONS:-50}"
CMAES_KEEP_TOP_EVAL_DIRS="${CMAES_KEEP_TOP_EVAL_DIRS:-3}"
CMAES_WARM_START_TOP_K="${CMAES_WARM_START_TOP_K:-1}"
CMAES_SIGMA0="${CMAES_SIGMA0:-0.05}"
CMAES_WARM_START_JITTER="${CMAES_WARM_START_JITTER:-0.05}"

GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
GPU_SLOTS_PER_DEVICE="${GPU_SLOTS_PER_DEVICE:-1}"
MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-3}"
BATCH_SIZE="${BATCH_SIZE:-$(($(wc -w <<< "${GPU_IDS}") * GPU_SLOTS_PER_DEVICE))}"

# search_common (via bootstrap) defaults to ftl_first; this pipeline uses general_ftl.
OBJECTIVE_MODE="${CAMPAIGN_OBJECTIVE_MODE:-general_ftl}"
export OBJECTIVE_MODE
export GRTECLYN_GEO_DIRECTIONS="${GRTECLYN_GEO_DIRECTIONS:-x y z}"

ftl_campaign_init_layout
ftl_campaign_write_manifest

PIN_DIMS=""
case "${BRANCH}" in
  wormhole) PIN_DIMS="$(ftl_general_ftl_wormhole_pins)" ;;
  ring) PIN_DIMS="$(ftl_general_ftl_ring_pins)" ;;
  spin) PIN_DIMS="$(ftl_general_ftl_spin_pins)" ;;
  *)
    echo "[campaign] unknown BRANCH=${BRANCH} (wormhole|ring|spin)" >&2
    exit 2
    ;;
esac

echo "== Full campaign: ${CAMPAIGN_NAME} =="
echo "Root          : ${CAMPAIGN_ROOT}"
echo "Branch        : ${BRANCH}"
echo "Stage         : ${STAGE} (RESUME=${RESUME})"
echo "QD target     : ${QD_TARGET_EVALS} evals (keep top ${QD_KEEP_TOP_EVAL_DIRS} + FTL champions)"
echo "CMA-ES target : ${CMAES_TARGET_EVALS} evals (warm-start top ${CMAES_WARM_START_TOP_K})"
echo "HQ            : top-1 CMA-ES eval"
echo "GPUs          : ${GPU_IDS}"
echo "Dry-run       : ${DRY_RUN}"
echo

# ---- Stage 0: MAP-Elites ------------------------------------------------------
if ftl_campaign_should_run_stage qd; then
  echo "== Stage 0: MAP-Elites (QD) =="
  RUNS_DIR="${CAMPAIGN_ROOT}" \
    QD_NAME=qd \
    QD_TARGET_EVALS="${QD_TARGET_EVALS}" \
    QD_KEEP_TOP_EVAL_DIRS="${QD_KEEP_TOP_EVAL_DIRS}" \
    QD_FTL_RETENTION="${QD_FTL_RETENTION}" \
    DESCRIPTOR_MODE=ftl_lifetime \
    GPU_IDS="${GPU_IDS}" \
    GPU_SLOTS_PER_DEVICE="${GPU_SLOTS_PER_DEVICE}" \
    MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA}" \
    BATCH_SIZE="${BATCH_SIZE}" \
    PIN_DIMS="${PIN_DIMS}" \
    DRY_RUN="${DRY_RUN}" \
    bash "${CAMPAIGNS_ROOT}/qd/run.sh" \
    2>&1 | tee "${CAMPAIGN_LOG_DIR}/qd.log"

  if [[ "${DRY_RUN}" != "1" ]]; then
    QD_TOP_EVAL="$(ftl_campaign_pick_top_eval "${QD_DIR}/trajectory.jsonl")"
    ftl_campaign_mark_stage qd "{\"top_eval\": ${QD_TOP_EVAL}, \"target_evals\": ${QD_TARGET_EVALS}}"
    echo "[campaign] QD complete. Top eval: ${QD_TOP_EVAL}"
  else
    ftl_campaign_mark_stage qd '{"dry_run": true}'
  fi
fi

# ---- Stage 1: CMA-ES ----------------------------------------------------------
if ftl_campaign_should_run_stage cmaes; then
  echo "== Stage 1: CMA-ES =="
  QD_TRAJECTORY="${QD_DIR}/trajectory.jsonl"
  if [[ "${DRY_RUN}" == "1" && ! -f "${QD_TRAJECTORY}" ]]; then
    mkdir -p "${QD_DIR}"
    echo '{"eval": 1, "status": "gpu_ok", "score": 100.0, "components": {}}' > "${QD_TRAJECTORY}"
  fi
  if [[ ! -f "${QD_TRAJECTORY}" ]]; then
    echo "[campaign] missing QD trajectory: ${QD_TRAJECTORY}" >&2
    exit 2
  fi

  RUNS_DIR="${CAMPAIGN_ROOT}" \
    RUN_NAME=cmaes \
    WARM_START_TRAJECTORY="${QD_TRAJECTORY}" \
    WARM_START_TOP_K="${CMAES_WARM_START_TOP_K}" \
    WARM_START_JITTER="${CMAES_WARM_START_JITTER}" \
    SIGMA0="${CMAES_SIGMA0}" \
    TARGET_EVALS="${CMAES_TARGET_EVALS}" \
    MAX_GENERATIONS="${CMAES_MAX_GENERATIONS}" \
    KEEP_TOP_EVAL_DIRS="${CMAES_KEEP_TOP_EVAL_DIRS}" \
    FTL_RETENTION="${QD_FTL_RETENTION}" \
    GPU_IDS="${GPU_IDS}" \
    GPU_SLOTS_PER_DEVICE="${GPU_SLOTS_PER_DEVICE}" \
    MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA}" \
    PIN_DIMS="${PIN_DIMS}" \
    DRY_RUN="${DRY_RUN}" \
    bash "${CAMPAIGNS_ROOT}/cmaes/run.sh" \
    2>&1 | tee "${CAMPAIGN_LOG_DIR}/cmaes.log"

  if [[ "${DRY_RUN}" != "1" ]]; then
    CMAES_TOP_EVAL="$(ftl_campaign_pick_top_eval "${CMAES_DIR}/trajectory.jsonl")"
    ftl_campaign_mark_stage cmaes "{\"top_eval\": ${CMAES_TOP_EVAL}, \"target_evals\": ${CMAES_TARGET_EVALS}}"
    echo "[campaign] CMA-ES complete. Top eval: ${CMAES_TOP_EVAL}"
  else
    ftl_campaign_mark_stage cmaes '{"dry_run": true}'
  fi
fi

# ---- Stage 2: HQ promotion ----------------------------------------------------
if ftl_campaign_should_run_stage hq; then
  echo "== Stage 2: HQ promotion =="
  CMAES_SOURCE="${CMAES_DIR}"
  CMAES_TRAJECTORY="${CMAES_DIR}/trajectory.jsonl"
  if [[ "${DRY_RUN}" == "1" && ! -f "${CMAES_TRAJECTORY}" ]]; then
    mkdir -p "${CMAES_DIR}"
    echo '{"eval": 1, "status": "gpu_ok", "score": 120.0, "components": {}}' > "${CMAES_TRAJECTORY}"
    mkdir -p "${CMAES_DIR}/eval_000001"
  fi
  if [[ ! -d "${CMAES_SOURCE}" || ! -f "${CMAES_TRAJECTORY}" ]]; then
    echo "[campaign] missing CMA-ES run: ${CMAES_SOURCE}" >&2
    exit 2
  fi

  HQ_GPU="$(ftl_campaign_first_gpu)"
  HQ_CANDIDATES="$(ftl_campaign_pick_top_pairs "${CMAES_TRAJECTORY}" "${HQ_GPU}")"
  HQ_EVAL_ID="${HQ_CANDIDATES%% *}"

  SOURCE_RUN="${CMAES_SOURCE}" \
    RUNS_DIR="${HQ_RUNS_DIR}" \
    NAME_PREFIX="${CAMPAIGN_NAME}" \
    OBJECTIVE_MODE="${OBJECTIVE_MODE}" \
    CANDIDATES="${HQ_CANDIDATES}" \
    FOREGROUND=1 \
    DRY_RUN="${DRY_RUN}" \
    bash "${CAMPAIGNS_ROOT}/hq/run_batch.sh" \
    2>&1 | tee "${CAMPAIGN_LOG_DIR}/hq.log"

  if [[ "${DRY_RUN}" != "1" ]]; then
    ftl_campaign_mark_stage hq "{\"eval_id\": ${HQ_EVAL_ID}, \"gpu\": ${HQ_GPU}}"
    echo "[campaign] HQ complete for eval ${HQ_EVAL_ID} on GPU ${HQ_GPU}"
  else
    ftl_campaign_mark_stage hq '{"dry_run": true}'
  fi
fi

ftl_campaign_write_manifest
echo
echo "[campaign] Done. Root: ${CAMPAIGN_ROOT}"
echo "  QD     : ${QD_DIR}"
echo "  CMA-ES : ${CMAES_DIR}"
echo "  HQ     : ${HQ_RUNS_DIR}"
echo "  State  : ${CAMPAIGN_STATE}"
