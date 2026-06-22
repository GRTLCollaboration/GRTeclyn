#!/usr/bin/env bash
# Gravitational splash MAP-Elites (canonical bosonic shell, critical_collapse).
# Superposes shell lump geometry into one U(1) complex scalar (BosonStarBH).
# Search space = scalar shell minus exotic dims + boson mass/lambda/omega.
#
# Example (smoke):
#   cd grteclyn-wrapper
#   QD_NAME=boson_shell_gpu_v1 QD_TARGET_EVALS=12 \
#     GPU_IDS="0 1 2 3" GPU_SLOTS_PER_DEVICE=1 \
#     bash scripts/campaigns/splash/run.sh
#
# Production (v22 pipeline: EvalPipeline + continuous GRTresna + pre-GPU learning):
#   QD_NAME=boson_shell_gpu_v7 QD_TARGET_EVALS=100 QD_ITERATIONS=30 \
#     GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=1 BATCH_SIZE=8 \
#     bash scripts/campaigns/splash/run.sh

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=shell
export GRTRESNA_MATTER_COUPLING=canonical
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-critical_collapse}"
export DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-wave_focusing}"
export SPLASH_MODE="${SPLASH_MODE:-discovery}"
export GRTRESNA_FULL_Z=1
# QD scoring uses central_timeseries / small_data — not slice PNGs.  Saves ~240
# PNGs/eval (11 z-fields + 6 projections × ~15 dumps).  Re-enable for viz:
#   GRTECLYN_FRAMES=1 FRAMES_FIELDS="chi K rho_req" bash scripts/campaigns/splash/run.sh
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-0}"
export GRTECLYN_FRAMES_AUTO_ZLIM="${GRTECLYN_FRAMES_AUTO_ZLIM:-1}"
export GRTECLYN_FRAMES_ZOOM="${GRTECLYN_FRAMES_ZOOM:-28}"
export FRAMES_FIELDS="${FRAMES_FIELDS:-chi K rho_req local_speed}"
export PROJECTION_FIELDS="${PROJECTION_FIELDS:-}"
# Plotfile dumps for central-ball / GW-proxy timeseries.  160 steps @ dt=0.01 →
# ~8 samples over a 16 s run (fewer with splash early-termination).
export PLOT_INTERVAL="${PLOT_INTERVAL:-160}"
# Stop once matter disperses after peak (typically t≈10–12 s); set 0 for full 16 s.
export GRTECLYN_SPLASH_EARLY_TERM="${GRTECLYN_SPLASH_EARLY_TERM:-1}"
# Belt-and-suspenders: search space pins sign; reject phantom coupling at launch.
# Static v13 default pins shell_static=1 (no velocity dims).  Moving test:
#   SPLASH_MOVING=1 PIN_DIMS="grtresna_scalar_sign=1 grtresna_shell_static=0" bash ...
if [[ "${SPLASH_MOVING:-0}" == "1" ]]; then
  export PIN_DIMS="${PIN_DIMS:-grtresna_scalar_sign=1 grtresna_shell_static=0}"
else
  export PIN_DIMS="${PIN_DIMS:-grtresna_scalar_sign=1 grtresna_shell_static=1}"
fi

# v22-style pipelined QD (MapElites.md): continuous GRTresna solves, pre-GPU learning auto with --grtresna
export USE_PIPELINE="${USE_PIPELINE:-1}"
export MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-5}"
# Relaxed from 2e-2 to 3e-2: recovers ~15% of postload-rejected configs whose
# Ham L2 sits in the 0.02-0.03 band (grid interpolation noise, not bad physics).
export POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-3e-2}"
export ITERATIONS="${ITERATIONS:-80}"

exec bash "${SCRIPT_DIR}/../qd/run.sh"
