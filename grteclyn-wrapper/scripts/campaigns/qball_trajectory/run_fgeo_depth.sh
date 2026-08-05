#!/usr/bin/env bash
# Q-ball trajectory MAP-Elites — pure DEPTH hunt (f_geo_depth).
#
# Motivation (see runs/grtresna_qd/qball_traj_fgeo_v1/CAMPAIGN_RESULTS.md §3.2):
# under f_geo_max the depth signal saturates at 20% and is multiplied by
# matter retention, so the v1 campaign optimized retention, not depth — the
# deepest corridor ever seen (eval 370, f_geo=0.383) scored 4x below the
# champion.  This campaign flips that:
#   * objective f_geo_depth — the RAW evolving-geodesic depth is the only
#     first-order reward (1% path saving = 100 points), uncapped and NOT
#     multiplied by survival.  No stability / confinement / constraint-health
#     shaping at all; a maximally deep transient beats a shallow survivor.
#   * Honesty gates stay upstream: depth is zero unless the ray bundle is
#     complete, constraint drift stayed small and emission is post-pump-floor.
#     The pump-energy tax and the graded horizon penalty stay on.
#   * NO exotic penalty — phantom lumps are free fuel (weight zeroed for any
#     general_ftl-style rescore too, as in run_fgeo.sh).
#
# Physics window: v1 showed depth still RISING at both the last emission time
# sampled (t=12) and the final frozen-probe frame (t=26).  Give the implosion
# room: stop_time 26 -> 32 and emission sweep 7 -> 10 launches (t=0..18,
# capped at the last metric-stack slice; late emits only land if the shortcut
# is deep — exactly the gradient we want).
#
# Usage:
#   cd grteclyn-wrapper
#   bash scripts/campaigns/qball_trajectory/run_fgeo_depth.sh
#
# Overrides: same as run.sh (QD_NAME, QD_TARGET_EVALS, GPU_IDS, QD_RESUME=1, ...).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export QD_NAME="${QD_NAME:-qball_traj_fgeo_depth_v1}"
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-f_geo_depth}"
export SCORE_EXOTIC_PENALTY_WEIGHT="${SCORE_EXOTIC_PENALTY_WEIGHT:-0}"

export STOP_TIME="${STOP_TIME:-32}"
export GRTECLYN_GEO_MAX_EMISSIONS="${GRTECLYN_GEO_MAX_EMISSIONS:-10}"

# Warm-start from the fgeo_v1 elites: the depth record-holders (eval 370
# f_geo=0.383, eval 179 f_geo=0.333) and the retention champions are the
# productive corner of the 39-D space for a depth hunt.  Each seed is
# re-evaluated under the current physics (stop_time=32) and scored with
# f_geo_depth.
FGEO_V1_DIR="$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)/runs/grtresna_qd/qball_traj_fgeo_v1"
if [[ -z "${SEED_EVAL_DIRS:-}" && -d "${FGEO_V1_DIR}" ]]; then
  SEED_EVAL_DIRS="$(ls -d "${FGEO_V1_DIR}"/eval_* 2>/dev/null | tr '\n' ' ')"
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
