#!/usr/bin/env bash
# Q-ball trajectory MAP-Elites — compact solitons on retrograde orbits.
#
# Uses the trajectory ansatz (5 independent per-lump orbits) with boson-star
# matter, compact Q-ball preset (m=1, λ=640, μ=85333, ω=0.8), ODE profile
# seeding, equilibrium amplitude cap, sub-luminal speed cap, and all-retrograde
# orbit enforcement.  Multi-ray emission sweep (Δt=2, 7 launches) maps the FTL
# channel lifetime via the evolving geodesic probe.
#
# Objective: general_ftl (gauge-invariant null-geodesic shortcut).
# Descriptor: ftl_lifetime (8×8 archive).
#
# Usage:
#   cd grteclyn-wrapper
#   QD_NAME=qball_traj_compact_v1 QD_TARGET_EVALS=200 \
#     GPU_IDS="0 1 2 3 4 5 6 7" \
#     bash scripts/campaigns/qball_trajectory/run.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# --- Matter model: boson star on trajectory orbits ---
export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=trajectory
export GRTRESNA_MATTER_COUPLING="${GRTRESNA_MATTER_COUPLING:-canonical}"
export GRTRESNA_FULL_Z=1

# --- Objective ---
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-general_ftl}"
export DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-ftl_lifetime}"

# --- Compact Q-ball physics (pinned, not searched) ---
# m=1, λ=640, μ=85333, ω=0.8 — thick-wall regime, localized soliton R~7.
# ODE profile seeding + equilibrium amplitude cap.
# All-retrograde orbits (removes 50% search space; HQ-validated).
export PIN_DIMS="${PIN_DIMS:-\
grtresna_scalar_mass=1.0 \
grtresna_scalar_lambda=640 \
grtresna_scalar_mu=85333 \
grtresna_bs_omega=0.8 \
grtresna_qball_ode_profile=1 \
grtresna_qball_equilibrium_amplitude=1 \
trajectory_retrograde_only=1}"

# --- Multi-ray emission sweep (7 launches, Δt=2 code units) ---
export GRTECLYN_GEO_EMIT_INTERVAL="${GRTECLYN_GEO_EMIT_INTERVAL:-2}"
export GRTECLYN_GEO_MAX_EMISSIONS="${GRTECLYN_GEO_MAX_EMISSIONS:-7}"

# --- Evolution ---
export STOP_TIME="${STOP_TIME:-16.0}"
export PLOT_INTERVAL="${PLOT_INTERVAL:-320}"

# --- Probe / frames ---
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-0}"
export GRTECLYN_EVOLVING_GEODESIC=1
export GRTECLYN_EVOLVING_GEODESIC_MODE="${GRTECLYN_EVOLVING_GEODESIC_MODE:-search}"
export GRTECLYN_GEO_DIRECTIONS="x y z"

# --- Pipeline ---
export USE_PIPELINE="${USE_PIPELINE:-1}"
export MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-5}"
export POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-3e-2}"

exec bash "${SCRIPT_DIR}/../qd/run.sh"
