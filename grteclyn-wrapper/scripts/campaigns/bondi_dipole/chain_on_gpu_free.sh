#!/usr/bin/env bash
# Launch the next bondi cell as soon as the predecessor RELEASES ITS GPU.
#
#   bash chain_on_gpu_free.sh <wait-for-runs-dir> <new-tag> <launcher.sh> [VAR=VAL ...]
#
# Why this exists alongside chain_next.sh (2026-08-17): chain_next.sh waits for
# the whole predecessor launcher to exit, and the launcher outlives the
# evolution by the time its consumer needs to turn the leftover plotfiles into
# frames -- CPU-only work that ran ~1 h after the N=256 control reached t=60.
# During that hour the GPU sits idle.  This variant watches only the evolution
# binary (main3d), so the successor starts the moment the device is free while
# the predecessor's post-processing finishes in parallel on the CPUs.
#
# chain_next.sh is left untouched and stays the default: waiting for the whole
# launcher is the conservative choice for campaigns whose successor needs the
# predecessor's *processed* output.  Bondi cells are independent, so releasing
# on the GPU is safe here.
#
# Double-launch is impossible even if an older chain_next.sh watcher is still
# parked on the same successor: run_pair_selfgrav.sh exits without touching
# anything when its output directory already exists, and this script checks the
# same directory before launching.
# --- SUPERSEDED 2026-08-22 -------------------------------------------------
# This script belongs to a finished campaign and reads/writes runs/bondi/rerun/
# which no longer exists.  Running it now would either fail halfway or
# recreate a dead tree and mix superseded cells into the live campaign.
# It is also an auto-chaining launcher, and the campaign no longer works that
# way: cells are launched by hand, one at a time, because an orchestrator
# outlives the session that made it and will relaunch cells behind your back.
# What replaced it: runs/bondi/staging/ -- plan and launch commands in
# research/bondi_dipole/docs/GPU_RUN_PAPER.md.
# Kept for the reasoning in the header above.  Remove this block only if you
# deliberately want the old tree back.
printf '%s: superseded campaign; %s no longer exists.\n' "$(basename "$0")" "runs/bondi/rerun/" >&2
printf '  see research/bondi_dipole/docs/GPU_RUN_PAPER.md for the live campaign\n' >&2
exit 2
# ---------------------------------------------------------------------------

set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "${REPO_ROOT}"

WAIT_DIR="$1"; shift
TAG="$1"; shift
LAUNCHER="$1"; shift

CAMP="grteclyn-wrapper/scripts/campaigns/bondi_dipole"
OUT="${REPO_ROOT}/runs/bondi/rerun"

# Is the predecessor's evolution binary still holding the GPU?
#
# `pgrep -f "${WAIT_DIR}"` matches this script's own argv and every other
# watcher, so filter on the command being the evolution binary itself -- the
# consumer, launcher and replay driver deliberately do NOT count here.
gpu_held() {
  local pid cmd
  for pid in $(pgrep -f "${WAIT_DIR}" 2>/dev/null); do
    [ "${pid}" = "$$" ] && continue
    cmd="$(tr '\0' ' ' < "/proc/${pid}/cmdline" 2>/dev/null || true)"
    case "${cmd}" in
      *main3d*) return 0 ;;
      *) continue ;;
    esac
  done
  return 1
}

# The predecessor does a CPU-only constraint solve (up to ~70 min at N=256)
# BEFORE it ever touches the GPU.  Waiting only for main3d to be ABSENT would
# fire immediately during that solve -- which is exactly what happened on
# 2026-08-17, launching a whole successor ladder hours early.  So wait for the
# evolution to appear first, and only then for it to finish.  If the
# predecessor's launcher dies without ever reaching the GPU (failed solve),
# stop waiting rather than hanging forever.
launcher_alive() {
  local pid cmd
  for pid in $(pgrep -f "${WAIT_DIR}" 2>/dev/null); do
    [ "${pid}" = "$$" ] && continue
    cmd="$(tr '\0' ' ' < "/proc/${pid}/cmdline" 2>/dev/null || true)"
    case "${cmd}" in
      *chain_next.sh*|*chain_on_gpu_free.sh*) continue ;;
      *replay_eval.py*|*run_pair_selfgrav.sh*|*run_single_selfgrav.sh*) return 0 ;;
    esac
  done
  return 1
}

echo "[chain-gpu] waiting for '${WAIT_DIR}' to reach its GPU before starting ${TAG}"
while ! gpu_held; do
  if ! launcher_alive; then
    echo "[chain-gpu] '${WAIT_DIR}' exited without reaching the GPU -- not waiting further"
    break
  fi
  sleep 30
done

echo "[chain-gpu] waiting for '${WAIT_DIR}' to release its GPU"
while gpu_held; do
  sleep 30
done

if [[ -d "${OUT}/${TAG}" ]]; then
  echo "[chain-gpu] ${TAG} already exists -- another watcher got there first; nothing to do"
  exit 0
fi

echo "[chain-gpu] GPU released; launching ${TAG}"
setsid nohup /usr/bin/env "$@" BONDI_RUNS_DIR="${OUT}/${TAG}" \
  bash "${CAMP}/${LAUNCHER}" \
  > "${OUT}/${TAG}.launch.log" 2>&1 < /dev/null &
disown

sleep 10
echo "[chain-gpu] ${TAG} launched"
