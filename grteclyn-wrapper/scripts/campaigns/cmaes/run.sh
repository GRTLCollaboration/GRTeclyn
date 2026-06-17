#!/usr/bin/env bash
# Stage 1 — CMA-ES local refinement around QD (or prior CMA-ES) elites.
# Physics matches stage 0 via campaigns/lib/search_common.sh.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/bootstrap.sh
source "${SCRIPT_DIR}/../lib/bootstrap.sh"
_campaign_bootstrap "${SCRIPT_DIR}"
_campaign_resolve_python

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"

RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_cmaes}"
MAX_GENERATIONS="${MAX_GENERATIONS:-25}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
POPULATION="${POPULATION:-$(wc -w <<< "${GPU_IDS}")}"
SEED="${SEED:-7}"
SIGMA0="${SIGMA0:-0.08}"
RANDOM_INJECTION_FRACTION="${RANDOM_INJECTION_FRACTION:-0.1}"
EXOTIC_INJECTION_FRACTION="${EXOTIC_INJECTION_FRACTION:-0.1}"
WARM_START_TOP_K="${WARM_START_TOP_K:-8}"
WARM_START_JITTER="${WARM_START_JITTER:-0.05}"
KEEP_TOP_EVAL_DIRS="${KEEP_TOP_EVAL_DIRS:-10}"
FTL_RETENTION="${FTL_RETENTION:-1}"

CONSUME_ARGS=()
[[ "${NO_CONSUME:-0}" == "1" ]] && CONSUME_ARGS=(--no-consume-plotfiles)

ftl_search_common_domain_args
ftl_search_common_grtresna_args
ftl_search_common_global_args

PRE_ARGS=(--runs-dir "${RUNS_DIR}" "${FTL_GLOBAL_ARGS[@]}")
[[ -n "${RUN_NAME:-}" ]] && PRE_ARGS+=(--name "${RUN_NAME}")
[[ "${DRY_RUN:-0}" == "1" ]] && PRE_ARGS+=(--dry-run)

SCORE_WEIGHT_ARGS=()
if [[ -n "${SCORE_WEIGHTS:-}" ]]; then
  # shellcheck disable=SC2206
  for pair in ${SCORE_WEIGHTS}; do
    SCORE_WEIGHT_ARGS+=(--score-weight "${pair}")
  done
fi
PRE_ARGS+=("${SCORE_WEIGHT_ARGS[@]}")

mkdir -p "${RUNS_DIR}"

echo "== Stage 1: CMA-ES (QD-matched setup) =="
echo "Runs dir      : ${RUNS_DIR}"
echo "Evolution     : L=${GRTRESNA_EVOLUTION_L_FULL} N=${GRTRESNA_EVOLUTION_N_FULL} ml=${GRTRESNA_EVOLUTION_MAX_LEVEL} t=${STOP_TIME}"
echo "4D geodesic   : mode=${GRTECLYN_EVOLVING_GEODESIC_MODE} stack_n=${GRTECLYN_METRIC_STACK_N_SPACE}"
echo "Objective     : ${OBJECTIVE_MODE}"
echo "CMA-ES        : gen=${MAX_GENERATIONS} pop=${POPULATION} sigma0=${SIGMA0}"
echo "Warm start    : ${WARM_START_TRAJECTORY:-<none>}"
echo

WARM_START_ARGS=()
if [[ -n "${WARM_START_TRAJECTORY:-}" ]]; then
  IFS=',' read -r -a WARM_START_PATHS <<< "${WARM_START_TRAJECTORY}"
  for traj in "${WARM_START_PATHS[@]}"; do
    WARM_START_ARGS+=(--warm-start-trajectory "${traj}")
  done
fi

if [[ "${SKIP_CMA_PREFLIGHT_TESTS:-0}" != "1" ]]; then
  ftl_search_common_preflight_tests
fi

PIN_ARGS=()
if [[ -n "${PIN_DIMS:-}" ]]; then
  for kv in ${PIN_DIMS}; do PIN_ARGS+=(--pin-dimension "${kv}"); done
fi

# shellcheck disable=SC2086
exec ${PYTHON_BIN} -m grteclyn_wrapper \
  "${PRE_ARGS[@]}" \
  optimize \
  --objective-mode "${OBJECTIVE_MODE}" \
  --random-injection-fraction "${RANDOM_INJECTION_FRACTION}" \
  --exotic-injection-fraction "${EXOTIC_INJECTION_FRACTION}" \
  --warm-start-top-k "${WARM_START_TOP_K}" \
  --warm-start-jitter "${WARM_START_JITTER}" \
  --keep-top-eval-dirs "${KEEP_TOP_EVAL_DIRS}" \
  "$([[ "${FTL_RETENTION}" == "1" ]] && echo --ftl-retention || echo --no-ftl-retention)" \
  "${WARM_START_ARGS[@]}" \
  "${FTL_GRTRESNA_ARGS[@]}" \
  --max-generations "${MAX_GENERATIONS}" \
  --population-size "${POPULATION}" \
  --sigma0 "${SIGMA0}" \
  --seed "${SEED}" \
  --gpu-ids ${GPU_IDS} \
  "${PIN_ARGS[@]}" \
  "${CONSUME_ARGS[@]}"
