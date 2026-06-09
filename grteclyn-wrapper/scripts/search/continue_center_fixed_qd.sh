#!/usr/bin/env bash
# Resume the center_fixed_search MAP-Elites campaign to a larger eval budget.
#
# Default: wait for the current grteclyn_wrapper qd process (if WAIT_PID set),
# then resume qd_20260609T094634Z from archive + trajectory toward 400 evals.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

RUNS_DIR="${RUNS_DIR:-/home/jovyan/nachevsky/test/simulation/GRTeclyn/runs/grtresna_qd/center_fixed_search}"
QD_NAME="${QD_NAME:-qd_20260609T094634Z}"
QD_TARGET_EVALS="${QD_TARGET_EVALS:-400}"
QD_RESUME=1
BATCH_SIZE="${BATCH_SIZE:-8}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
RANKS="${RANKS:-8}"
SEED="${SEED:-7}"
BINS="${BINS:-8}"

# Match the original center_fixed_search campaign settings.
STOP_TIME="${STOP_TIME:-2.0}"
PLOT_INTERVAL="${PLOT_INTERVAL:-10}"
LUMPS="${LUMPS:-5}"
ITERATIONS="${ITERATIONS:-30}"
GRTRESNA_TIMEOUT="${GRTRESNA_TIMEOUT:-900}"
FRAMES_FIELDS="${FRAMES_FIELDS:-lump_activity scalar_activity phi Pi chi chi_minus_1 local_speed shift1 rho_req}"
PROJECTION_FIELDS="${PROJECTION_FIELDS:-lump_activity scalar_activity}"
PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
CONSUMER_RADII="${CONSUMER_RADII:-4 8}"

WAIT_PID="${WAIT_PID:-}"

if [[ -n "${WAIT_PID}" ]]; then
  echo "[continue] waiting for PID ${WAIT_PID} to finish..."
  while kill -0 "${WAIT_PID}" 2>/dev/null; do
    sleep 30
  done
  echo "[continue] PID ${WAIT_PID} finished."
fi

CAMPAIGN_DIR="${RUNS_DIR}/${QD_NAME}"
if [[ ! -f "${CAMPAIGN_DIR}/trajectory.jsonl" ]]; then
  echo "[continue] missing trajectory: ${CAMPAIGN_DIR}/trajectory.jsonl" >&2
  exit 1
fi

COMPLETED="$(python3 - <<'PY' "${CAMPAIGN_DIR}/trajectory.jsonl"
import json, sys
path = sys.argv[1]
last = 0
with open(path, encoding="utf-8") as fh:
    for line in fh:
        if line.strip():
            last = int(json.loads(line)["eval"])
print(last)
PY
)"
REMAINING=$((QD_TARGET_EVALS - COMPLETED))
if [[ "${REMAINING}" -le 0 ]]; then
  echo "[continue] already at ${COMPLETED} evals (target ${QD_TARGET_EVALS}); nothing to do."
  exit 0
fi

PLANNED_BATCHES=$(( (REMAINING + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "[continue] resuming ${QD_NAME}: ${COMPLETED} -> ${QD_TARGET_EVALS} evals (${PLANNED_BATCHES} batches of ${BATCH_SIZE})"

export RUNS_DIR QD_NAME QD_TARGET_EVALS QD_RESUME BATCH_SIZE GPU_IDS RANKS SEED BINS
export STOP_TIME PLOT_INTERVAL LUMPS ITERATIONS GRTRESNA_TIMEOUT
export FRAMES_FIELDS PROJECTION_FIELDS PROJECTION_AXES CONSUMER_RADII

exec bash "${SCRIPT_DIR}/run_grtresna_qd_search.sh"
