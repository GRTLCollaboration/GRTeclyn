#!/usr/bin/env bash
# Wait for boson_shell_ftl_rl_v1 to finish, then launch scalar_shell_ftl_rl_v1.
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
GRTECLYN_ROOT="$(cd -- "${ROOT}/.." && pwd)"
BOSON_RUN="${GRTECLYN_ROOT}/runs/grtresna_qd/boson_shell_ftl_rl_v1"
SCALAR_LOG="${GRTECLYN_ROOT}/runs/scalar_shell_ftl_rl_v1.launch.log"
WATCH_LOG="${GRTECLYN_ROOT}/runs/scalar_shell_ftl_rl_v1.watch.log"
TARGET="${QD_TARGET_EVALS:-200}"
POLL_SEC="${POLL_SEC:-120}"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" | tee -a "${WATCH_LOG}"; }

boson_running() {
  pgrep -f "grteclyn_wrapper.*--name boson_shell_ftl_rl_v1.* qd " >/dev/null 2>&1
}

boson_eval_counter() {
  python3 - <<'PY' "${BOSON_RUN}/metadata.json"
import json, sys
from pathlib import Path
p = Path(sys.argv[1])
if not p.exists():
    print(0)
    raise SystemExit
m = json.loads(p.read_text())
print(int(m.get("last_eval_counter") or 0))
PY
}

log "Watching boson_shell_ftl_rl_v1 (target_evals=${TARGET}, poll=${POLL_SEC}s)"

while true; do
  counter="$(boson_eval_counter)"
  if boson_running; then
    log "boson still running — last_eval_counter=${counter}/${TARGET}"
  else
    if (( counter >= TARGET )); then
      log "boson finished — last_eval_counter=${counter}"
      break
    fi
    log "boson process gone but counter=${counter}<${TARGET}; waiting for metadata catch-up"
  fi
  sleep "${POLL_SEC}"
done

if [[ -d "${GRTECLYN_ROOT}/runs/grtresna_qd/scalar_shell_ftl_rl_v1" ]]; then
  log "scalar_shell_ftl_rl_v1 dir already exists — not launching (remove dir to rerun)"
  exit 0
fi

log "Launching scalar_shell_ftl_rl_v1"
cd "${ROOT}"
uv sync >> "${WATCH_LOG}" 2>&1

QD_NAME=scalar_shell_ftl_rl_v1 \
QD_TARGET_EVALS="${TARGET}" \
QD_ITERATIONS=30 \
STOP_TIME=16.0 \
PLOT_INTERVAL=320 \
GRTECLYN_FRAMES=1 \
SKIP_QD_PREFLIGHT_TESTS="${SKIP_QD_PREFLIGHT_TESTS:-1}" \
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}" \
GPU_SLOTS_PER_DEVICE="${GPU_SLOTS_PER_DEVICE:-1}" \
MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-5}" \
BATCH_SIZE="${BATCH_SIZE:-8}" \
  nohup bash scripts/campaigns/general_ftl/scalar_shell_ftl_run.sh \
  >> "${SCALAR_LOG}" 2>&1 &

log "scalar launch pid=$! log=${SCALAR_LOG}"
