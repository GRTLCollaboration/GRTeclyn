#!/usr/bin/env bash
# GW beam v5 MAP-Elites — directional gravitational-wave laser campaign.
#
# v5 improvements over v4:
#   - Variable lump count: well_depth lower bound = 0 (search decides active lumps)
#   - Higher v_max = 0.5 (tighter, faster binaries possible)
#   - Sponge zone enabled: stop_time=40 without boundary reflections
#   - Multi-radius extraction (R=12,16,20) with 1/r wave-zone validity check
#   - Constrained scoring: log10(P) × beaming_gain × gate (no survival/stability bias)
#   - Junk-radiation truncation: samples before R_ext+margin excluded
#   - Descriptor y-axis: beaming_gain (replaces proxy beam_ratio)
#
# Physics choices:
#   - L=64, N=128 (dx=0.5) for QD cost.
#   - Sponge zone: inner_radius=24, outer_radius=32, strength=4.0 (quartic ramp).
#   - Multi-radius extraction R=12,16,20 for 1/r validation.
#   - stop_time=40: light-travel from orbit (R~5) to R=20 is ~15; leaves ~20 units
#     of clean post-junk observing window.  Sponge absorbs reflections.
#   - Canonical matter only (no exotic bullshit).
#
# Objective: gw_beam (v5: log10(P) × beaming_gain × gate).
# Descriptor: gw_beam (x=log power, y=beaming_gain).
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
QD_RUN="${SCRIPT_DIR}/../qd/run.sh"

# --- QD run identity / scale (override via env) ---
export QD_NAME="${QD_NAME:-gw_beam_v5}"
export QD_TARGET_EVALS="${QD_TARGET_EVALS:-200}"
export GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
export BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
export QD_ITERATIONS="${QD_ITERATIONS:-25}"

# --- Matter model: compact canonical Q-balls on trajectory orbits ---
export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=trajectory
export GRTRESNA_MATTER_COUPLING=canonical
export GRTRESNA_FULL_Z=1

# --- Objective: directional gravitational-wave emission (v5) ---
export OBJECTIVE_MODE=gw_beam
export DESCRIPTOR_MODE=gw_beam

# --- No exotic matter: every lump is canonical (+rho) ---
export SCORE_EXOTIC_PENALTY_WEIGHT=0.0

# --- v5: Variable lump count — well_depth can go to 0 (lumps switch off) ---
# Unpin well_depth from PIN_DIMS and lower bound in search space to 0.0.
# The search decides optimal number of active lumps.
export PIN_DIMS="${PIN_DIMS:-\
grtresna_scalar_mass=1.0 \
grtresna_scalar_lambda=640 \
grtresna_bs_omega=0.8 \
trajectory_lump0_exotic=0 \
trajectory_lump1_exotic=0 \
trajectory_lump2_exotic=0 \
trajectory_lump3_exotic=0 \
trajectory_lump4_exotic=0}"

# v5: Raise v_max to 0.5 (from default 0.3) — faster binaries.
export TRAJECTORY_V_MAX="${TRAJECTORY_V_MAX:-0.5}"

# --- Evolution: sponge zone enabled, longer stop_time ---
export STOP_TIME=40.0
# Multi-radius GW extraction for 1/r wave-zone validation.
export CONSUMER_RADII="${CONSUMER_RADII:-12 16 20}"
# Plot cadence: ~20 plotfiles over 40 units → interval ~960 coarse steps.
# At dx=0.5, dt_mult=0.02 → dt=0.01, so 960 steps ≈ 9.6 time units.
# Better: aim for ~15 Psi4 samples in the post-junk window (t>20 to t=40).
# That's 20 time units / 15 ≈ 1.3 units per sample → interval = 130 steps.
export PLOT_INTERVAL=640

# Sponge zone parameters (injected as EXTRA_SETS → params.txt overrides).
SPONGE_PARAMS="sponge_enabled=1 sponge_inner_radius=24 sponge_outer_radius=32 sponge_strength=4.0 sponge_ramp_power=4"

# --- GRTresna speed ---
export ITERATIONS=25
export GRTRESNA_MAX_LEVEL=1
export GRTRESNA_TIMEOUT=300
export GRTRESNA_NL_EXIT_TOLERANCE=2.0

# Non-search-space overrides.
export EXTRA_SETS="${EXTRA_SETS:-\
grtresna_scalar_mu=85333 \
grtresna_qball_ode_profile=1 \
grtresna_qball_equilibrium_amplitude=1 \
trajectory_retrograde_only=0 \
trajectory_r_min=0.1 \
trajectory_v_max=${TRAJECTORY_V_MAX} \
trajectory_well_depth_min=0.0 \
stop_time=${STOP_TIME} \
plot_interval=${PLOT_INTERVAL} \
amr.derive_plot_vars=Weyl4 \
${SPONGE_PARAMS}}"

# --- GW extraction ---
export GRTECLYN_PSI4=1
export GRTECLYN_PSI4_N_POINTS=64
export GRTECLYN_CONSUMER_DRAIN_MINIMAL=1
export CONSUMER_KEEP_LAST=0
export CONSUMER_JOBS=4
export GRTECLYN_EVOLVING_GEODESIC=0
export GRTECLYN_GEO_DIRECTIONS=""
export GRTECLYN_GEO_EMIT_INTERVAL=""
export GRTECLYN_GEO_MAX_EMISSIONS=""
export GRTECLYN_FRAMES=0
# Delete plotfiles after Psi4 extraction (frames disabled, so need this override).
export GRTECLYN_DELETE_WITHOUT_FRAMES=1

# --- Pipeline ---
export USE_PIPELINE="${USE_PIPELINE:-1}"
_grtresna_cap=6
if (( BATCH_SIZE < _grtresna_cap )); then _grtresna_cap="${BATCH_SIZE}"; fi
export MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-${_grtresna_cap}}"
export POSTLOAD_GATE=0
export POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-3e-2}"

echo "== GW beam v5 QD: ${QD_NAME} =="
echo "   GPUs: ${GPU_IDS} (batch=${BATCH_SIZE})  target_evals=${QD_TARGET_EVALS}"
echo "   Box: L=64 N=128  R_ext=${CONSUMER_RADII}  t_stop=${STOP_TIME}"
echo "   Sponge: inner=24 outer=32 strength=4.0"
echo "   v_max=${TRAJECTORY_V_MAX}  variable lumps (well_depth_min=0)"
echo "   Objective: ${OBJECTIVE_MODE} (v5: log(P)×beaming_gain×gate)"
echo "   Pipeline: max_grtresna=${MAX_CONCURRENT_GRTRESNA}"

if [[ "${PIPELINE_MONITOR:-1}" == "1" ]]; then
  CAMPAIGNS_LIB="${SCRIPT_DIR}/../lib"
  # shellcheck source=../lib/bootstrap.sh
  source "${CAMPAIGNS_LIB}/bootstrap.sh"
  _campaign_bootstrap "${SCRIPT_DIR}"
  # shellcheck source=../lib/pipeline_monitor.sh
  source "${CAMPAIGNS_LIB}/pipeline_monitor.sh"
  ftl_pipeline_monitor_begin "${QD_NAME}" "${GPU_IDS}"
  bash "${QD_RUN}" 2>&1 | tee "${FTL_PIPELINE_LOG}"
  ftl_pipeline_monitor_end
else
  exec bash "${QD_RUN}"
fi
