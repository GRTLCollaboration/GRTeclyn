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

# This node: 4 GPUs, and mpirun segfaults cluster-wide — GRTresna must run as
# single-rank singleton solves (see runner.py).  Do not raise RANKS here.
export GPU_IDS="${GPU_IDS:-0 1 2 3}"
export RANKS="${RANKS:-1}"

exec bash "${SCRIPT_DIR}/run.sh"
