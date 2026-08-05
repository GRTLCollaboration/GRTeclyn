#!/usr/bin/env bash
# Bondi dipole runaway -- weak-field matrix v2: CORRECTED SEEDS, DENSE FRAMES.
#
# Post-mortem of v1 (bondi_dipole_weakfield_v1, 2026-08-05): every lump
# dispersed (rms 5 -> 30 by t~50) because the GRTresna seed knobs were never
# overridden and fell back to defaults:
#   grtresna_bs_profile_width = 8.0  -- as wide as the whole pair separation,
#     so the two "lumps" overlapped from t=0;
#   grtresna_bs_phi_c = 0.08         -- 2.1x the Q-ball plateau amplitude
#     sqrt(3*lambda/4mu) = 0.0375 for lambda=2560, mu=1365333.
# (trajectory_well_width only shapes the PUMP wells, and the pump is off.)
#
# v2 seeds what the soliton actually wants:
#   grtresna_bs_phi_c         = 0.0375  (thin-wall plateau)
#   grtresna_bs_profile_width = 2.0     (natural bound-state width is 1.20;
#     2.0 keeps ~8 finest cells across the core and ~4.6x the charge of a
#     width-1.2 blob; sep/width = 4 so the pair starts cleanly separated)
# and runs single_p FIRST as a live calibration: if the lone lump still
# disperses, stop the matrix (scripts/campaigns/stop_campaign.sh) and pick a
# binding rung between the two ladders instead of burning the pair cells.
#
# Frames: plot every 40 steps (0.4 time units, ~250 frames per run) so the
# movies are smooth -- v1's 160-step cadence gave only 64 frames per run.
#
# Physics, matrix cells, and diagnostics otherwise identical to v1 -- see
# run_matrix_weakfield.sh header for the Bondi (1957) background.
#
# Usage:
#   bash scripts/campaigns/bondi_dipole/run_matrix_weakfield2.sh
# Overrides: BONDI_GPU (default 1), BONDI_STOP_TIME (default 100),
#            BONDI_RUNS_DIR, BONDI_SEP (default 8).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi_dipole_weakfield_v2}"
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
  --extra-override grtresna_scalar_lambda=2560
  --extra-override grtresna_scalar_mu=1365333
  --extra-override grtresna_bs_omega=0.55
  --extra-override grtresna_bs_phi_c=0.0375
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
# single_p runs FIRST: it is the live calibration for the corrected seed.
matrix=(
  "single_p 1 0 0"
  "pair_pm  2 0 1"
  "pair_pp  2 0 0"
  "pair_mm  2 1 1"
  "single_m 1 1 0"
)

for spec in "${matrix[@]}"; do
  read -r name num_lumps exotic0 exotic1 <<<"${spec}"
  out_name="bondi_wf2_${name}"
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
