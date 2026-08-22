#!/usr/bin/env bash
# Wait for a running bondi cell to exit, then launch the next one on its GPU.
#
#   bash chain_next.sh <wait-for-runs-dir> <new-tag> <launcher.sh> [VAR=VAL ...]
#
# Keeps the GPU pipeline full without hand-holding: the successor is launched
# only once the predecessor's launcher process is gone, so the two never share
# a device.  Note the explicit /usr/bin/env -- see wave1_rerun.sh for why a bare
# `env` silently launches nothing on this machine.
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

# Wait for the predecessor to exit.
#
# `pgrep -f "${WAIT_DIR}"` on its own is WRONG here: this script's own argv
# contains WAIT_DIR, so it matches itself and waits forever.  Filter out our own
# PID and any other chain_next.sh instance before deciding the predecessor is
# still alive.
# Only a real launcher counts as "still running".  Matching WAIT_DIR alone is
# not enough: any bystander whose command line happens to mention the run
# directory -- a `tail -F` over the launch logs, a du, an editor -- would pin the
# chain open forever.  That actually happened: a log monitor tailing
# runs/bondi/rerun/*.launch.log held every chain shut and left a GPU idle.
# So require the process to be the launcher script or its replay_eval child.
predecessor_alive() {
  local pid cmd
  for pid in $(pgrep -f "${WAIT_DIR}" 2>/dev/null); do
    [ "${pid}" = "$$" ] && continue
    cmd="$(tr '\0' ' ' < "/proc/${pid}/cmdline" 2>/dev/null || true)"
    case "${cmd}" in
      *chain_next.sh*) continue ;;
      *replay_eval.py*|*run_pair_selfgrav.sh*|*run_single_selfgrav.sh*) return 0 ;;
      *) continue ;;
    esac
  done
  return 1
}

echo "[chain] waiting for '${WAIT_DIR}' to finish before starting ${TAG}"
while predecessor_alive; do
  sleep 30
done
echo "[chain] predecessor gone; launching ${TAG}"

setsid nohup /usr/bin/env "$@" BONDI_RUNS_DIR="${OUT}/${TAG}" \
  bash "${CAMP}/${LAUNCHER}" \
  > "${OUT}/${TAG}.launch.log" 2>&1 < /dev/null &
disown

sleep 10
echo "[chain] ${TAG} launched"
