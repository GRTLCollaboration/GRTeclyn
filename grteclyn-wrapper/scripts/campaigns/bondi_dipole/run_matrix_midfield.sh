#!/usr/bin/env bash
# Bondi dipole runaway -- MID-FIELD matrix: the rung between the two ladders.
#
# Post-mortem of weakfield v2 (bondi_dipole_weakfield_v2, 2026-08-05): the
# corrected seed knobs (bs_phi_c=0.0375, bs_profile_width=2.0) did NOT cure
# dispersal.  The single_p calibration cell -- one lump, nothing to overlap
# with -- breathed in to rms 4.65 by t=8, then dispersed on v1's exact curve
# (rms 7.46 by t=20 vs v1's 8.7 by t=19).  So overlap and seed shape were
# never the root cause: at core amp 0.0375 the lump is simply too light for
# self-gravity to confine it once the pump is off.  (The strong-field ladder,
# core amp 0.075, had the opposite failure: collapse to min_chi=0.009 by
# t~25-30.)
#
# This rung sits at the geometric mean of the two failed amplitudes:
#   lambda = 1280, mu = 341333  (keeps lambda^2/mu = 4.8 -> omega_min = 0.316)
#   core amp sqrt(3*lambda/4mu) = 0.053   (strong 0.075 / weak 0.0375)
#   omega = 0.55 unchanged (45% binding), width 2.0 unchanged.
# single_p runs FIRST as the live calibration:
#   PASS: rms_radius_canon roughly flat to t=100 AND min_chi stays above ~0.3
#         -> the matrix proceeds to the pair cells on its own.
#   FAIL disperse: rms climbs like v1/v2 -> next rung up (lambda 960?).
#   FAIL collapse: min_chi dives toward 0 -> next rung down.
#   Stop with scripts/campaigns/stop_campaign.sh either way.
#
# Frames: plot every 40 steps (0.4 time units, ~250 frames per run) for
# smooth movies.  Physics, matrix cells, and diagnostics otherwise identical
# to the weak-field matrices -- see run_matrix_weakfield.sh for the Bondi
# (1957) background.
#
# Usage:
#   bash scripts/campaigns/bondi_dipole/run_matrix_midfield.sh
# Overrides: BONDI_GPU (default 1), BONDI_STOP_TIME (default 100),
#            BONDI_RUNS_DIR, BONDI_SEP (default 8).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi/dipole_midfield_v1}"
GPU="${BONDI_GPU:-1}"
STOP_TIME="${BONDI_STOP_TIME:-100}"
SEP="${BONDI_SEP:-8}"
R0="$(python3 -c "print(${SEP}/2)")"
PI="3.141592653589793"

mkdir -p "${RUNS_DIR}"

# Stop handle for scripts/campaigns/stop_campaign.sh.
source "${SCRIPT_DIR}/../lib/launcher_common.sh"
campaign_register_launcher "${RUNS_DIR}"

# Live per-sector centroid stream + radiation extraction.
export GRTECLYN_SECTOR_BARYCENTERS=1
export GRTECLYN_PSI4=1

# This node: mpirun segfaults cluster-wide -- GRTresna and the evolution both
# run single-rank.  Do not raise.
GRTRESNA_RANKS=1

# Everything shared by all five runs.  Lumps at rest: no orbital speed, no
# radial speed, no breathing, no z-motion, no initial-data boost, and the
# trajectory pump is off for the whole evolution (t >= 0).
common_overrides=(
  --extra-override grtresna_scalar_lambda=1280
  --extra-override grtresna_scalar_mu=341333
  --extra-override grtresna_bs_omega=0.55
  --extra-override grtresna_bs_phi_c=0.053
  --extra-override grtresna_bs_profile_width=2.0
  --extra-override trajectory_well_width=1.2
  --extra-override rl_pump_stop_time=0
  --extra-override grtresna_boost_lumps=0
  --extra-override trajectory_A_breath=0
  --extra-override trajectory_z_amp=0
  --extra-override trajectory_omega_z=0
  --extra-override trajectory_lump0_R0="${R0}"
  --extra-override trajectory_lump0_phase0=0
  --extra-override trajectory_lump0_tilt_theta=0
  --extra-override trajectory_lump0_tilt_phi=0
  --extra-override trajectory_lump0_v_rad=0
  --extra-override trajectory_lump0_omega_rot=0
  --extra-override trajectory_lump0_well_depth=0.15
  --extra-override trajectory_lump1_R0="${R0}"
  --extra-override trajectory_lump1_phase0="${PI}"
  --extra-override trajectory_lump1_tilt_theta=0
  --extra-override trajectory_lump1_tilt_phi=0
  --extra-override trajectory_lump1_v_rad=0
  --extra-override trajectory_lump1_omega_rot=0
  --extra-override trajectory_lump1_well_depth=0.15
)

# name | num_lumps | lump0_exotic | lump1_exotic (ignored for singles)
# single_p runs FIRST: it is the live calibration for this binding rung.
matrix=(
  "single_p 1 0 0"
  "pair_pm  2 0 1"
  "pair_pp  2 0 0"
  "pair_mm  2 1 1"
  "single_m 1 1 0"
)

for spec in "${matrix[@]}"; do
  read -r name num_lumps exotic0 exotic1 <<<"${spec}"
  out_name="bondi_mf_${name}"
  if [[ -d "${RUNS_DIR}/${out_name}" ]]; then
    echo "[bondi] ${out_name} already exists -- skipping"
    continue
  fi
  echo "[bondi] === ${out_name} (lumps=${num_lumps} s0=${exotic0} s1=${exotic1}) ==="
  PYTHONPATH="${WRAPPER_DIR}/src" "${WRAPPER_DIR}/.venv/bin/python" \
    "${WRAPPER_DIR}/scripts/campaigns/hq/replay_eval.py" \
    "${SOURCE_EVAL}" \
    --name "${out_name}" \
    --runs-dir "${RUNS_DIR}" \
    --gpu "${GPU}" \
    --n-full 128 --l-full 64 \
    --max-level 1 --regrid-threshold 0.02 \
    --stop-time "${STOP_TIME}" --plot-interval 40 \
    --ftl-L 8.0 \
    --grtresna-ranks "${GRTRESNA_RANKS}" \
    --grtresna-iterations 50 \
    --grtresna-max-level 3 \
    --grtresna-domain-l 128 \
    --grtresna-timeout 7200 \
    --consumer-radii 8 16 \
    --consumer-keep-last 2 \
    --objective-mode weighted \
    --extra-override trajectory_num_lumps="${num_lumps}" \
    --extra-override trajectory_lump0_exotic="${exotic0}" \
    --extra-override trajectory_lump1_exotic="${exotic1}" \
    "${common_overrides[@]}" \
    || echo "[bondi] WARNING: ${out_name} exited nonzero -- continuing with next cell"
done

echo "[bondi] matrix complete: ${RUNS_DIR}"
