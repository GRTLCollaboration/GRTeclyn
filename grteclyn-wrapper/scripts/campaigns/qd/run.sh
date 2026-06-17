#!/usr/bin/env bash
# Stage 0 — MAP-Elites quality-diversity survey over the 23-D shell search space.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/bootstrap.sh
source "${SCRIPT_DIR}/../lib/bootstrap.sh"
_campaign_bootstrap "${SCRIPT_DIR}"
_campaign_resolve_python

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"

RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_qd}"
QD_ITERATIONS="${QD_ITERATIONS:-10}"
QD_TARGET_EVALS="${QD_TARGET_EVALS:-}"
QD_RESUME="${QD_RESUME:-0}"
QD_NAME="${QD_NAME:-}"
DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-ftl_lifetime}"
SEED_EVAL_DIRS="${SEED_EVAL_DIRS:-}"
QD_KEEP_TOP_EVAL_DIRS="${QD_KEEP_TOP_EVAL_DIRS:-10}"
QD_FTL_RETENTION="${QD_FTL_RETENTION:-1}"
BINS="${BINS:-8}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
SEED="${SEED:-7}"

mkdir -p "${RUNS_DIR}"

NAME_ARGS=()
[[ -n "${QD_NAME}" ]] && NAME_ARGS+=(--name "${QD_NAME}")
RESUME_ARGS=()
[[ "${QD_RESUME}" == "1" ]] && RESUME_ARGS+=(--resume)
TARGET_EVALS_ARGS=()
[[ -n "${QD_TARGET_EVALS}" ]] && TARGET_EVALS_ARGS+=(--target-evals "${QD_TARGET_EVALS}")
SEED_ARGS=()
if [[ -n "${SEED_EVAL_DIRS}" ]]; then
  # shellcheck disable=SC2206
  SEED_ARGS+=(--seed-eval-dirs ${SEED_EVAL_DIRS})
fi

FTL_RETENTION_ARGS=(--ftl-retention)
[[ "${QD_FTL_RETENTION}" == "0" ]] && FTL_RETENTION_ARGS=(--no-ftl-retention)

DRY_RUN_ARGS=()
[[ "${DRY_RUN:-0}" == "1" ]] && DRY_RUN_ARGS=(--dry-run)

SCORE_WEIGHT_ARGS=()
if [[ -n "${SCORE_WEIGHTS:-}" ]]; then
  # shellcheck disable=SC2206
  for pair in ${SCORE_WEIGHTS}; do
    SCORE_WEIGHT_ARGS+=(--score-weight "${pair}")
  done
fi

ftl_search_common_domain_args
ftl_search_common_grtresna_args
ftl_search_common_global_args
ftl_search_common_pipeline_args

if [[ "${SKIP_QD_PREFLIGHT_TESTS:-0}" != "1" ]]; then
  ftl_search_common_preflight_tests
fi

PIN_ARGS=()
if [[ -n "${PIN_DIMS:-}" ]]; then
  for kv in ${PIN_DIMS}; do PIN_ARGS+=(--pin-dimension "${kv}"); done
fi

# shellcheck disable=SC2086
exec ${PYTHON_BIN} -m grteclyn_wrapper \
  --runs-dir "${RUNS_DIR}" \
  "${NAME_ARGS[@]}" \
  "${FTL_GLOBAL_ARGS[@]}" \
  "${DRY_RUN_ARGS[@]}" \
  "${SCORE_WEIGHT_ARGS[@]}" \
  qd \
  --descriptor-mode "${DESCRIPTOR_MODE}" \
  --objective-mode "${OBJECTIVE_MODE}" \
  --iterations "${QD_ITERATIONS}" \
  --keep-top-eval-dirs "${QD_KEEP_TOP_EVAL_DIRS}" \
  "${FTL_RETENTION_ARGS[@]}" \
  "${TARGET_EVALS_ARGS[@]}" \
  "${RESUME_ARGS[@]}" \
  "${SEED_ARGS[@]}" \
  --batch-size "${BATCH_SIZE}" \
  --bins "${BINS}" \
  --seed "${SEED}" \
  --gpu-ids ${GPU_IDS} \
  "${FTL_PIPELINE_ARGS[@]}" \
  "${PIN_ARGS[@]}" \
  "${FTL_GRTRESNA_ARGS[@]}"
