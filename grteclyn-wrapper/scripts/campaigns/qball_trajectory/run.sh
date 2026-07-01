#!/usr/bin/env bash
# Q-ball trajectory MAP-Elites — compact solitons on retrograde orbits.
#
# Uses the trajectory ansatz (5 independent per-lump orbits) with boson-star
# matter, compact Q-ball preset (m=1, lambda=640, mu=85333, omega=0.8), ODE profile
# seeding, equilibrium amplitude cap, sub-luminal speed cap, and all-retrograde
# orbit enforcement.  Multi-ray emission sweep (dt=2, 7 launches) maps the FTL
# channel lifetime via the evolving geodesic probe.
#
# Objective: general_ftl (gauge-invariant null-geodesic shortcut).
# Descriptor: ftl_lifetime (8x8 archive).
#
# --- MAP-Elites search space (~40 continuous dimensions) ---
#
# Pinned (not searched):
#   Q-ball physics: m=1, lambda=640, mu=85333, omega=0.8 (thick-wall soliton).
#   ODE radial profile, equilibrium amplitude cap, retrograde-only orbits.
#
# Per lump (5 lumps, 7 dims each = 35 dims):
#   R0            — orbital radius [1.5, 8.0]
#   omega_rot     — angular velocity (negative = retrograde, speed-capped).
#                   Tangential speed v_t = R0 * |omega_rot| is capped at 0.3c
#                   jointly with v_rad so the total velocity stays sub-luminal.
#   phase0        — initial orbital phase [0, 2pi]
#   tilt_theta    — orbital plane polar tilt [0, pi]
#   tilt_phi      — orbital plane azimuthal tilt [0, 2pi]
#   v_rad         — radial drift speed for spiral orbits [-0.3, 0.3].
#                   v_rad > 0 spirals outward, v_rad < 0 spirals inward.
#   exotic        — continuous [0,1], rounded to binary 0 or 1 in config:
#                   0 = canonical Q-ball (+rho), 1 = phantom Q-ball (-rho).
#                   Each lump is purely one or the other, no mixture.
#   (well_depth is pinned to 0.15, which is clamped to the equilibrium
#    core amplitude 0.075 by QBallCouplings.cap_well_depth; this gives all
#    lumps the same, strongest on-attractor pump.)
#
# Global (shared, 4 dims):
#   A_breath      — radial breathing amplitude
#   omega_breath  — breathing frequency
#   z_amp         — vertical oscillation amplitude
#   omega_z       — vertical oscillation frequency
#   (well_width is not searched: the boson-star width is fixed by the
#    physics via bound_width(m, omega) = 1/sqrt(m^2 - omega^2) = 1.667.)
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

# Exotic matter is the fuel of this FTL engine, not a failure.  Reduce the
# scoring penalty so the optimizer explores the high-FTL exotic-rich region
# instead of hiding in long-lived but inert canonical configurations.
export SCORE_EXOTIC_PENALTY_WEIGHT="${SCORE_EXOTIC_PENALTY_WEIGHT:-0.2}"

# --- Compact Q-ball physics ---
# Search-space dimensions pinned (these exist in the trajectory-boson search space):
#   m=1, lambda=640, omega=0.8 — thick-wall regime, localized soliton R~7.
#   Per-lump well_depth pinned to 0.15; with equilibrium_amplitude=1 this is
#   clamped to the core_amplitude = sqrt(3*lambda/(4*mu)) = 0.075, so every lump
#   receives the same maximum on-attractor pump and the optimizer cannot weaken it.
export PIN_DIMS="${PIN_DIMS:-\
grtresna_scalar_mass=1.0 \
grtresna_scalar_lambda=640 \
grtresna_bs_omega=0.8 \
trajectory_lump0_well_depth=0.15 \
trajectory_lump1_well_depth=0.15 \
trajectory_lump2_well_depth=0.15 \
trajectory_lump3_well_depth=0.15 \
trajectory_lump4_well_depth=0.15 \
trajectory_well_width=1.667}"

# Non-search-space overrides passed as --set via EXTRA_SETS (base overrides).
# These keys are read by the config builder but are not optimizer dimensions.
export EXTRA_SETS="${EXTRA_SETS:-\
grtresna_scalar_mu=85333 \
grtresna_qball_ode_profile=1 \
grtresna_qball_equilibrium_amplitude=1 \
trajectory_retrograde_only=1}"

# --- Multi-ray emission sweep (7 launches, dt=2 code units) ---
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
