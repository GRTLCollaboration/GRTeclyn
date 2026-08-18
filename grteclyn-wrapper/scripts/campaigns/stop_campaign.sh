#!/usr/bin/env bash
# Universal detached-campaign stopper.  Works for ANY campaign in this repo:
# QD / MAP-Elites, CMA-ES, HQ replays, Bondi matrices, one-off ladders.
#
# Usage:
#   bash scripts/campaigns/stop_campaign.sh [--dry-run] <runs_dir | campaign_name> ...
#
#   --dry-run    print what would be killed, kill nothing
#   arguments    absolute runs dir, or a name resolved under
#                runs/neuralspacetime/search/map_elites/<name> then runs/<name>
#
# Why this exists (2026-08-05 post-mortem, bondi_dipole_v1): stopping a
# detached campaign by the pid captured at launch, or by pattern-killing its
# workers, does not work --
#   1. `$!` after `setsid nohup ... &` is the dead setsid parent, not the
#      session-leader launcher.
#   2. Killing workers (evolution binary, consumer) by path pattern looks like
#      a FINISHED STEP to the orchestrator, which then starts the next eval --
#      the campaign appears to refuse to die.
#   3. GRTresna solvers detach into their own session/pgid and survive
#      group-kills of the launcher.
#
# Working order, implemented here: freeze the queue (orchestrators + their
# shell ancestors FIRST), then sweep workers, then VERIFY and escalate.
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"

DRY=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY=1
  shift
fi
if [[ $# -eq 0 ]]; then
  echo "usage: stop_campaign.sh [--dry-run] <runs_dir | campaign_name> ..." >&2
  exit 2
fi

resolve_dir() {
  local arg="$1"
  if [[ -d "$arg" ]]; then (cd -- "$arg" && pwd); return 0; fi
  local base
  for base in "${REPO_ROOT}/runs/neuralspacetime/search/map_elites" "${REPO_ROOT}/runs/neuralspacetime/search/cma_es" "${REPO_ROOT}/runs/neuralspacetime/hq" "${REPO_ROOT}/runs"; do
    if [[ -d "${base}/${arg}" ]]; then echo "${base}/${arg}"; return 0; fi
  done
  return 1
}

# Never touch: this script, its callers, and Claude/interactive harness shells.
is_protected() {
  local pid="$1"
  [[ "$pid" == "$$" || "$pid" == "$PPID" ]] && return 0
  local cmd
  cmd="$(ps -o cmd= -p "$pid" 2>/dev/null)" || return 0
  case "$cmd" in
    *stop_campaign.sh*|*shell-snapshots*) return 0 ;;
  esac
  return 1
}

# Shell ancestors of a pid that are themselves campaign scripts (bash x.sh).
# Freezing these matters: a surviving launcher loop treats dead children as
# completed steps and starts the next one.
shell_ancestors() {
  local pid="$1" pp cmd out=""
  while :; do
    pp="$(ps -o ppid= -p "$pid" 2>/dev/null | tr -d ' ')" || break
    [[ -z "$pp" || "$pp" -le 1 ]] && break
    cmd="$(ps -o cmd= -p "$pp" 2>/dev/null)"
    case "$cmd" in
      *shell-snapshots*) break ;;
      *bash*.sh*) out+=" $pp"; pid="$pp" ;;
      *) break ;;
    esac
  done
  echo "$out"
}

live_filter() {
  local pid out=""
  for pid in "$@"; do
    [[ -n "$pid" ]] || continue
    kill -0 "$pid" 2>/dev/null || continue
    is_protected "$pid" && continue
    out+=" $pid"
  done
  echo "$out"
}

pgrep_safe() { pgrep -f "$1" 2>/dev/null | grep -vw "$$" || true; }

do_kill() {
  local sig="$1"; shift
  local pids
  pids="$(live_filter "$@" | tr ' ' '\n' | sort -un | tr '\n' ' ')"
  [[ -z "${pids// /}" ]] && return 0
  if [[ "$DRY" == 1 ]]; then
    local pid
    for pid in $pids; do
      echo "  [dry-run] would kill -$sig $(ps -o pid=,cmd= -p "$pid" 2>/dev/null | cut -c1-120)"
    done
  else
    kill "-$sig" $pids 2>/dev/null
  fi
}

overall_rc=0
for arg in "$@"; do
  dir="$(resolve_dir "$arg")" || { echo "[stop] no runs dir for '$arg'" >&2; overall_rc=1; continue; }
  name="$(basename "$dir")"
  echo "[stop] campaign: $name ($dir)"

  # ---- 1. Orchestrators: pid file, drivers by --name/--runs-dir, ancestors --
  orch=""
  if [[ -f "$dir/launcher.pid" ]]; then
    orch+=" $(cat "$dir/launcher.pid" 2>/dev/null)"
  fi
  orch+=" $(pgrep_safe "grteclyn_wrapper .*--name ${name}( |$)")"      # QD/CMA-ES driver
  orch+=" $(pgrep_safe "replay_eval.py .*--runs-dir ${dir}")"          # replay/matrix cells
  ancestors=""
  for pid in $(live_filter $orch); do
    ancestors+=" $(shell_ancestors "$pid")"
  done
  # One KILL batch: launcher loops + drivers together, so nothing survives to
  # treat a dead sibling as a finished step and advance the queue.
  do_kill KILL $ancestors $orch

  # ---- 2. Workers: anything carrying the runs dir or its scratch in argv ----
  workers=" $(pgrep_safe "${dir}/")"                    # evolution, GRTresna, consumer
  workers+=" $(pgrep_safe "grteclyn_scratch/${name}_")" # scratch-watching consumers
  do_kill TERM $workers

  # ---- 3. Verify; escalate ---------------------------------------------------
  if [[ "$DRY" == 1 ]]; then
    echo "  [dry-run] verification skipped"
    continue
  fi
  survivors=""
  for _pass in 1 2 3; do
    survivors="$(live_filter $(pgrep_safe "${dir}/|--name ${name}( |$)|grteclyn_scratch/${name}_"))"
    [[ -z "${survivors// /}" ]] && break
    kill -KILL $survivors 2>/dev/null
  done
  if [[ -n "${survivors// /}" ]]; then
    echo "[stop] WARNING: survivors for $name:"
    ps -o pid=,cmd= -p $survivors 2>/dev/null
    overall_rc=1
  else
    echo "[stop] $name fully stopped."
  fi
done
exit "$overall_rc"
