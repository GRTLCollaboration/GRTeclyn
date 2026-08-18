#!/usr/bin/env bash
# gpu_queue.sh — press-and-forget work-queue runner for a single node.
#
# One detached runner per pool keeps every slot busy: each worker claims the
# next pending job the moment its slot frees, so GPUs are re-fed continuously
# instead of being chained pairwise. Jobs may enqueue follow-up jobs (e.g. a
# CPU solve job that, on success, drops the matching GPU evolve job into a
# GPU pool) — the runner keeps polling until told to stop, so late-enqueued
# work is picked up.
#
#   bash gpu_queue.sh <queue-dir> <slot> [<slot>...]
#
#   <slot>       a CUDA device id (0..3) for a GPU pool, or any label
#                (s1, s2, ...) for a CPU pool. One worker per slot.
#   <queue-dir>  created on first run with this layout:
#                  pending/   *.job files (plain bash), dispatched in
#                             lexicographic order — prefix 010_, 020_ to order
#                  running/   claimed jobs (suffixed .slotN while running)
#                  done/ failed/  finished jobs, kept for the record
#                  logs/      one log per job + queue.log event stream
#                  STOP       sentinel: touch to stop dispatch (running jobs
#                             finish; runner exits when all slots idle)
#
# Each job runs with QUEUE_SLOT=<label> and QUEUE_GPU=<label> exported; GPU
# jobs should pass --gpu "${QUEUE_GPU}" / GPU_ID="${QUEUE_GPU}" themselves.
# Jobs must run in the FOREGROUND (no setsid/nohup inside a job) — the runner
# is what gets detached, once.
#
# Env knobs:
#   QUEUE_POLL_SEC        poll interval while idle/gated (default 30)
#   QUEUE_GPU_MEM_MAX_MB  if set and the slot is numeric, a job is dispatched
#                         only when nvidia-smi reports used memory on that GPU
#                         at or below this — lets the runner start while
#                         another campaign still owns the GPUs and take each
#                         device over the moment it clears
#   QUEUE_IDLE_EXIT=1     exit when pending and running are both empty
#                         (default: keep polling until STOP appears)
#
# A lock file guards against two runners on the same queue dir. Stale claims
# from a previously killed runner are re-queued at startup (safe: workers die
# with the runner, so a leftover running/ entry cannot still be executing).
#
# Failure policy: a job with nonzero exit moves to failed/ and the worker
# moves on. Nothing is retried automatically — inspect logs/<job>.log, fix,
# and `mv failed/X pending/X` to requeue.
set -euo pipefail

if (( $# < 2 )); then
  echo "usage: bash gpu_queue.sh <queue-dir> <slot> [<slot>...]" >&2
  exit 2
fi

QDIR="$(mkdir -p "$1" && cd -- "$1" && pwd)"
shift
SLOTS=("$@")
POLL="${QUEUE_POLL_SEC:-30}"

mkdir -p "${QDIR}/pending" "${QDIR}/running" "${QDIR}/done" "${QDIR}/failed" "${QDIR}/logs"

note() {
  printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >> "${QDIR}/queue.log"
}

# Single-runner lock (per queue dir).
exec 9> "${QDIR}/runner.lock"
if ! flock -n 9; then
  echo "another runner already owns ${QDIR} (runner.lock held)" >&2
  exit 1
fi
echo "$$" > "${QDIR}/runner.pid"

# Re-queue stale claims from a killed runner.
for stale in "${QDIR}/running"/*; do
  [[ -e "${stale}" ]] || continue
  base="$(basename "${stale}")"
  base="${base%.slot*}"
  mv -- "${stale}" "${QDIR}/pending/${base}"
  note "requeued stale claim ${base}"
done

gpu_busy() {
  # 0 (busy) when the slot is a CUDA id, a gate is configured, and the device
  # still holds more memory than the gate allows.
  local slot="$1" used
  [[ "${slot}" =~ ^[0-9]+$ ]] || return 1
  [[ -n "${QUEUE_GPU_MEM_MAX_MB:-}" ]] || return 1
  used="$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "${slot}" 2>/dev/null | tr -dc '0-9')"
  [[ -n "${used}" ]] || return 1
  (( used > QUEUE_GPU_MEM_MAX_MB ))
}

worker() {
  local slot="$1" job claimed
  while :; do
    [[ -e "${QDIR}/STOP" ]] && { note "slot ${slot}: STOP seen, exiting"; return 0; }
    if gpu_busy "${slot}"; then
      sleep "${POLL}"; continue
    fi
    job="$(ls -1 "${QDIR}/pending" 2>/dev/null | LC_ALL=C sort | head -1 || true)"
    if [[ -z "${job}" ]]; then
      if [[ "${QUEUE_IDLE_EXIT:-0}" == "1" ]] && ! ls -1A "${QDIR}/running" 2>/dev/null | grep -q .; then
        note "slot ${slot}: queue drained, exiting"
        return 0
      fi
      sleep "${POLL}"; continue
    fi
    claimed="${QDIR}/running/${job}.slot${slot}"
    mv -- "${QDIR}/pending/${job}" "${claimed}" 2>/dev/null || continue  # lost the claim race
    note "slot ${slot}: start ${job}"
    if QUEUE_SLOT="${slot}" QUEUE_GPU="${slot}" bash "${claimed}" \
        >> "${QDIR}/logs/${job}.log" 2>&1; then
      mv -- "${claimed}" "${QDIR}/done/${job}"
      note "slot ${slot}: done ${job}"
    else
      mv -- "${claimed}" "${QDIR}/failed/${job}"
      note "slot ${slot}: FAILED ${job} — see logs/${job}.log"
    fi
  done
}

note "runner start pid=$$ slots=${SLOTS[*]} gate=${QUEUE_GPU_MEM_MAX_MB:-none}"
pids=()
for slot in "${SLOTS[@]}"; do
  worker "${slot}" &
  pids+=("$!")
done
trap 'note "runner interrupted"; kill "${pids[@]}" 2>/dev/null || true' INT TERM
wait
note "runner exit"
