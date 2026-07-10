#!/usr/bin/env bash
# HQ promotion for gw_beam v5 — replay a QD elite at full HQ resolution.
#
# Replays the best gw_beam_v5 eval at N=256, L=128, max_level=3 with:
#   - Multi-radius Psi4 extraction (12 16 20 24)
#   - Sponge zone (carried from source eval metadata)
#   - gw_beam v5 scoring (log(P)×beaming_gain×gate)
#   - No evolving geodesic (GW-only, no FTL)
#   - No live frames (metrics-only; post-run drain optional)
#   - Plotfile deletion after processing
#
# Usage:
#   # Promote the single best eval:
#   SOURCE_RUN=runs/grtresna_qd/gw_beam_v5 \
#   CANDIDATES="46 0" \
#   bash scripts/campaigns/gw_beam/promote_hq.sh
#
#   # Promote top 3 by score:
#   SOURCE_RUN=runs/grtresna_qd/gw_beam_v5 \
#   TOP_K=3 GPU_IDS="0 1 2" \
#   bash scripts/campaigns/gw_beam/promote_hq.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
HQ_RUN="${CAMPAIGNS_ROOT}/hq/run_batch.sh"

# --- GW beam HQ defaults ---
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-gw_beam}"
export NAME_PREFIX="${NAME_PREFIX:-gw_beam_v5_hq}"

# HQ resolution (same as FTL HQ promotions)
export N_FULL="${N_FULL:-256}"
export L_FULL="${L_FULL:-128}"
export MAX_LEVEL="${MAX_LEVEL:-3}"
export STOP_TIME="${STOP_TIME:-40}"
export PLOT_INTERVAL="${PLOT_INTERVAL:-24}"

# Multi-radius Psi4 extraction (wider for HQ box)
export CONSUMER_RADII="${CONSUMER_RADII:-12 16 20 24}"
export GRTECLYN_PSI4="${GRTECLYN_PSI4:-1}"
export GRTECLYN_PSI4_N_POINTS="${GRTECLYN_PSI4_N_POINTS:-128}"

# No evolving geodesic — GW-only promotion
export GRTECLYN_EVOLVING_GEODESIC="${GRTECLYN_EVOLVING_GEODESIC:-0}"
export GRTECLYN_GEO_DIRECTIONS=""

# Live frames on (render PNG movies during evolution; plotfiles kept for post-run drain)
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-1}"
export CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-3}"

# Sponge zone params (also carried from source eval metadata, but set as fallback)
export SPONGE_ENABLED="${SPONGE_ENABLED:-1}"
export SPONGE_INNER_RADIUS="${SPONGE_INNER_RADIUS:-56}"
export SPONGE_OUTER_RADIUS="${SPONGE_OUTER_RADIUS:-64}"
export SPONGE_STRENGTH="${SPONGE_STRENGTH:-4.0}"
export SPONGE_RAMP_POWER="${SPONGE_RAMP_POWER:-4}"

# GRTresna HQ solve settings
export GRTRESNA_TIMEOUT="${GRTRESNA_TIMEOUT:-7200}"
export GRTRESNA_ITERATIONS="${GRTRESNA_ITERATIONS:-30}"
export GRTRESNA_RANKS="${GRTRESNA_RANKS:-8}"
export GRTRESNA_MAX_HAM_PCT="${GRTRESNA_MAX_HAM_PCT:-10}"
export GRTRESNA_MAX_MOM_PCT="${GRTRESNA_MAX_MOM_PCT:-10}"

# Extra overrides: sponge + Weyl4 (passed to run_batch.sh as EXTRA_OVERRIDE)
export EXTRA_OVERRIDE="${EXTRA_OVERRIDE:-\
sponge_enabled=${SPONGE_ENABLED} \
sponge_inner_radius=${SPONGE_INNER_RADIUS} \
sponge_outer_radius=${SPONGE_OUTER_RADIUS} \
sponge_strength=${SPONGE_STRENGTH} \
sponge_ramp_power=${SPONGE_RAMP_POWER} \
amr.derive_plot_vars=Weyl4}"

echo "== GW beam v5 HQ promotion =="
echo "   Objective: ${OBJECTIVE_MODE}"
echo "   Domain: L=${L_FULL} N=${N_FULL} max_level=${MAX_LEVEL} t=${STOP_TIME}"
echo "   Psi4 radii: ${CONSUMER_RADII}  n_points=${GRTECLYN_PSI4_N_POINTS}"
echo "   Sponge: inner=${SPONGE_INNER_RADIUS} outer=${SPONGE_OUTER_RADIUS} strength=${SPONGE_STRENGTH}"
echo "   Geodesic: off (GW-only)"
echo "   Frames: on (PNG movies + post-run drain)"
echo

# Delegate to the shared HQ run_batch.sh
exec bash "${HQ_RUN}"
