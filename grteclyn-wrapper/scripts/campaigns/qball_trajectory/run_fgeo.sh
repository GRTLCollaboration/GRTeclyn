#!/usr/bin/env bash
# Q-ball trajectory MAP-Elites — pure evolving-shortcut hunt (f_geo_max).
#
# Same physics, search space, probes and pipeline as run.sh; only the scoring
# differs:
#   * objective f_geo_max — the evolving-geodesic shortcut is the only
#     first-order reward (1% shortcut = 100 points).  Frozen-probe shortcut,
#     persistence and health are second-order shaping so the archive keeps a
#     gradient before the first live ray lands.
#   * NO exotic penalty — phantom lumps are free fuel.  f_geo_max has no
#     exotic term at all; the weight is zeroed here too so any
#     general_ftl-style rescore of this campaign's evals stays penalty-free.
# The graded horizon penalty and the pump-energy tax stay on so a collapsing
# or pump-inflated configuration cannot fake a shortcut.
#
# Usage (after the current campaign frees the GPUs):
#   cd grteclyn-wrapper
#   bash scripts/campaigns/qball_trajectory/run_fgeo.sh
#
# Overrides: same as run.sh (QD_NAME, QD_TARGET_EVALS, GPU_IDS, QD_RESUME=1, ...).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export QD_NAME="${QD_NAME:-qball_traj_fgeo_v1}"
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-f_geo_max}"
export SCORE_EXOTIC_PENALTY_WEIGHT="${SCORE_EXOTIC_PENALTY_WEIGHT:-0}"

# With the pump on for the whole run the geometry matures late (ray emission
# lands around t~10, vs t~6 in the pump-then-release era) and a flat-speed
# crossing takes 14.4 code units, so at stop_time=20 almost no evolving ray can
# arrive before the run ends (seed replays: frozen f_geo 0.27, evolving 0.0,
# n_reached 0).  26 gives every ray emitted up to t~11.6 room to land.
export STOP_TIME="${STOP_TIME:-26}"

# Warm-start from pump_v2's retained elites (mostly-phantom shortcut-makers,
# evolving f_geo up to 0.22) so the search starts in the productive corner of
# the 39-D space instead of rediscovering it from random throws.  Each seed is
# re-evaluated under the current physics and scored with f_geo_max.
PUMP_V2_DIR="$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)/runs/grtresna_qd/qball_traj_pump_v2"
if [[ -z "${SEED_EVAL_DIRS:-}" && -d "${PUMP_V2_DIR}" ]]; then
  SEED_EVAL_DIRS="$(ls -d "${PUMP_V2_DIR}"/eval_* 2>/dev/null | tr '\n' ' ')"
fi
export SEED_EVAL_DIRS

# This node: 4 GPUs, and mpirun segfaults cluster-wide — GRTresna must run as
# single-rank singleton solves (see runner.py).  Do not raise RANKS here.
export GPU_IDS="${GPU_IDS:-0 1 2 3}"
export RANKS="${RANKS:-1}"

# Stop handle for scripts/campaigns/stop_campaign.sh.  exec preserves the PID,
# so registering here records the pid run.sh will actually run under.
source "${SCRIPT_DIR}/../lib/launcher_common.sh"
campaign_register_launcher "$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)/runs/grtresna_qd/${QD_NAME}"

exec bash "${SCRIPT_DIR}/run.sh"
