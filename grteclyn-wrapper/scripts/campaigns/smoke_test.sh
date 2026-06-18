#!/usr/bin/env bash
# Smoke-test QD, CMA-ES, and HQ campaign launchers (dry-run + preflight).
set -euo pipefail

ROOT="$(cd -- "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
STAMP="$(date +%Y%m%d_%H%M%S)"
SMOKE_RUNS="${GRTECLYN_ROOT:-$(cd .. && pwd)}/runs/_campaign_smoke_${STAMP}"
mkdir -p "${SMOKE_RUNS}"

pass=0
fail=0

ok() { echo "[PASS] $*"; pass=$((pass + 1)); }
bad() { echo "[FAIL] $*"; fail=$((fail + 1)); }

run_step() {
  local name="$1"
  shift
  echo
  echo "=== ${name} ==="
  if "$@"; then
    ok "${name}"
  else
    bad "${name}"
  fi
}

# ---- 1. Preflight pytest (shared QD + CMA-ES gate) ---------------------------
run_step "preflight pytest" bash -c '
  cd "'"${ROOT}"'"
  WRAPPER_ROOT="'"${ROOT}"'" PYTHON_BIN="uv run python" \
    bash -c "source scripts/campaigns/lib/search_common.sh && ftl_search_common_preflight_tests"
'

# ---- 2. QD dry-run -----------------------------------------------------------
QD_SMOKE="${SMOKE_RUNS}/qd"
run_step "QD dry-run" env \
  RUNS_DIR="${QD_SMOKE}" \
  QD_NAME="smoke_${STAMP}" \
  QD_ITERATIONS=1 \
  QD_TARGET_EVALS=1 \
  GPU_IDS="0" \
  BATCH_SIZE=1 \
  SKIP_QD_PREFLIGHT_TESTS=1 \
  DRY_RUN=1 \
  bash scripts/campaigns/qd/run.sh

# ---- 3. CMA-ES dry-run -------------------------------------------------------
CMA_SMOKE="${SMOKE_RUNS}/cmaes"
run_step "CMA-ES dry-run" env \
  RUNS_DIR="${CMA_SMOKE}" \
  RUN_NAME="smoke_${STAMP}" \
  MAX_GENERATIONS=1 \
  POPULATION=2 \
  GPU_IDS="0 1" \
  SKIP_CMA_PREFLIGHT_TESTS=1 \
  DRY_RUN=1 \
  bash scripts/campaigns/cmaes/run.sh

# ---- 4. CMA-ES warm-start wiring (parse only) --------------------------------
FTL4D_RUN="${GRTECLYN_ROOT:-$(cd .. && pwd)}/runs/grtresna_qd/ftl_4d/ftl_4d_v1"
if [[ -f "${FTL4D_RUN}/trajectory.jsonl" ]]; then
  run_step "CMA-ES warm-start from ftl_4d_v1 trajectory" bash -c '
    out=$(env \
      RUNS_DIR="'"${CMA_SMOKE}"'" \
      RUN_NAME="smoke_ws_'"${STAMP}"'" \
      WARM_START_TRAJECTORY="'"${FTL4D_RUN}/trajectory.jsonl"'" \
      WARM_START_TOP_K=2 \
      MAX_GENERATIONS=1 \
      POPULATION=2 \
      GPU_IDS="0 1" \
      SKIP_CMA_PREFLIGHT_TESTS=1 \
      DRY_RUN=1 \
      bash scripts/campaigns/cmaes/run.sh 2>&1)
    echo "${out}" | grep -q "loaded.*warm-start"
  '
else
  echo "[SKIP] CMA-ES warm-start (no ftl_4d_v1 trajectory)"
fi

# ---- 5. HQ replay_eval.py CLI ------------------------------------------------
run_step "HQ replay_eval.py --help" uv run python scripts/campaigns/hq/replay_eval.py --help

# ---- 6. HQ batch dry-run -----------------------------------------------------
SOURCE_RUN="${FTL4D_RUN}"
if [[ -d "${SOURCE_RUN}" ]]; then
  run_step "HQ batch dry-run (explicit candidates)" env \
    RUNS_DIR="${SMOKE_RUNS}/hq" \
    SOURCE_RUN="${SOURCE_RUN}" \
    CANDIDATES="90 0" \
    NAME_PREFIX=smoke \
    DRY_RUN=1 \
    bash scripts/campaigns/hq/run_batch.sh

  run_step "HQ batch TOP_K auto-pick" bash -c '
    out=$(env \
      RUNS_DIR="'"${SMOKE_RUNS}/hq"'" \
      SOURCE_RUN="'"${SOURCE_RUN}"'" \
      TOP_K=2 \
      MIN_FTL_GEO_EVOL=0.01 \
      NAME_PREFIX=smoke \
      DRY_RUN=1 \
      bash scripts/campaigns/hq/run_batch.sh 2>&1)
    echo "${out}" | grep -q "eval_"
  '
else
  bad "HQ batch (missing SOURCE_RUN ${SOURCE_RUN})"
fi

# ---- 7. Config parity QD vs CMA-ES -----------------------------------------
run_step "QD/CMA-ES share search_common.sh" bash -c '
  diff -q scripts/campaigns/lib/search_common.sh scripts/lib/ftl_search_common.sh >/dev/null 2>&1 || \
  grep -q campaigns/lib/search_common scripts/lib/ftl_search_common.sh
'

# ---- 8. Full campaign orchestrator dry-run (QD stage) ------------------------
run_step "full campaign orchestrator dry-run (QD)" env \
  CAMPAIGN_NAME="smoke_full_${STAMP}" \
  STAGE=qd \
  DRY_RUN=1 \
  QD_TARGET_EVALS=1 \
  SKIP_QD_PREFLIGHT_TESTS=1 \
  GPU_IDS="0" \
  BATCH_SIZE=1 \
  bash scripts/campaigns/run_full_campaign.sh

echo
echo "============================================"
echo "Campaign smoke: ${pass} passed, ${fail} failed"
echo "Artifacts (if any): ${SMOKE_RUNS}"
echo "============================================"
[[ "${fail}" -eq 0 ]]
