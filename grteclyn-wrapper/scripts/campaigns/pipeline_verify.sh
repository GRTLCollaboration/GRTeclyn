#!/usr/bin/env bash
# Short pipelined ring QD — delegates to general_ftl/run_all.sh (no duplicate launcher).
#
# Usage:
#   GPU_IDS="0 1 2 3" QD_TARGET_EVALS=20 bash scripts/campaigns/pipeline_verify.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export BRANCH=ring
export PIPELINE_MONITOR=1
export OBJECTIVE_MODE=general_ftl
export STOP_TIME="${STOP_TIME:-4.0}"
export PLOT_INTERVAL="${PLOT_INTERVAL:-40}"
export GPU_IDS="${GPU_IDS:-0 1 2 3}"
export BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
export QD_TARGET_EVALS="${QD_TARGET_EVALS:-20}"
export QD_ITERATIONS="${QD_ITERATIONS:-10}"
export USE_PIPELINE="${USE_PIPELINE:-1}"
export SKIP_QD_PREFLIGHT_TESTS="${SKIP_QD_PREFLIGHT_TESTS:-1}"

STAMP="$(date +%Y%m%d_%H%M%S)"
export QD_NAME="${QD_NAME:-pipeline_verify_ring_${STAMP}}"

exec bash "${SCRIPT_DIR}/general_ftl/run_all.sh"
