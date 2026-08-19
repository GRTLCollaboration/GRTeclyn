#!/usr/bin/env bash
# Bondi dipole -- exact-amplitude calibration: ONE canonical lump, weak rung.
#
# Root cause of all four 2026-08-05 failures (bondi_dipole_debug.md §7): the
# painted lump was the exact Q-ball eigenstate at 95.45% of its amplitude --
# cap_well_depth clamps the amp to the thin-wall estimate sqrt(3*lambda/4mu),
# 4.55% under the ODE table's true phi_c, and both painters rescale the table
# by amp/phi_c.  grtresna_qball_exact_amplitude=1 paints at phi_c exactly.
#
# This is the §7.3 Step-1 gate: lone canonical lump at the WEAK rung (lowest
# compactness 0.14, so the remaining flat-space-profile approximation is
# smallest), stop_time 40.  Pass gate: rms_radius flat (+-10%) to t=40 AND
# min_chi plateau above 0.3.  Only then burn pair cells.
#
# Launch-time verification (t=0 row of small_data/sector_barycenters.dat):
#   total_canon ~= 36.30  -> exact amplitude painted (fix active)
#   total_canon ~= 34.65  -> old clamped seed (fix did NOT reach the solver)
# Ground truth even earlier: grtresna/params.txt lump0_amp must read
# 0.03928682593468194 (the ODE phi_c), not 0.03750000457763756 (the clamp).
#
# The old grtresna_bs_phi_c / grtresna_bs_profile_width overrides are dropped:
# they are dead knobs on the lump path (written to params.txt, read by nobody).
#
# Usage:
#   bash scripts/campaigns/bondi_dipole/run_single_exact.sh
# Overrides: BONDI_GPU (default 1), BONDI_STOP_TIME (default 40),
#            BONDI_RUNS_DIR, BONDI_SEP (default 8).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi_dipole_exact_v1}"
GPU="${BONDI_GPU:-1}"
STOP_TIME="${BONDI_STOP_TIME:-40}"
SEP="${BONDI_SEP:-8}"
R0="$(python3 -c "print(${SEP}/2)")"

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

out_name="bondi_ex_single_p"
if [[ -d "${RUNS_DIR}/${out_name}" ]]; then
  echo "[bondi] ${out_name} already exists -- delete it or set BONDI_RUNS_DIR"
  exit 1
fi

echo "[bondi] === ${out_name} (weak rung, exact eigenstate amplitude) ==="
PYTHONPATH="${WRAPPER_DIR}/src" "${WRAPPER_DIR}/.venv/bin/python" \
  "${WRAPPER_DIR}/scripts/campaigns/hq/replay_eval.py" \
  "${SOURCE_EVAL}" \
  --name "${out_name}" \
  --runs-dir "${RUNS_DIR}" \
  --gpu "${GPU}" \
  --n-full 128 --l-full 64 \
  --max-level 1 --regrid-threshold 0.02 \
  --stop-time "${STOP_TIME}" --plot-interval 80 \
  --ftl-L 8.0 \
  --grtresna-ranks "${GRTRESNA_RANKS}" \
  --grtresna-iterations 50 \
  --grtresna-max-level 3 \
  --grtresna-domain-l 128 \
  --grtresna-timeout 7200 \
  --consumer-radii 8 16 \
  --consumer-keep-last 2 \
  --objective-mode weighted \
  --extra-override grtresna_scalar_lambda=2560 \
  --extra-override grtresna_scalar_mu=1365333 \
  --extra-override grtresna_bs_omega=0.55 \
  --extra-override grtresna_qball_exact_amplitude=1 \
  --extra-override trajectory_well_width=1.2 \
  --extra-override rl_pump_stop_time=0 \
  --extra-override grtresna_boost_lumps=0 \
  --extra-override trajectory_A_breath=0 \
  --extra-override trajectory_z_amp=0 \
  --extra-override trajectory_omega_z=0 \
  --extra-override trajectory_num_lumps=1 \
  --extra-override trajectory_lump0_R0="${R0}" \
  --extra-override trajectory_lump0_phase0=0 \
  --extra-override trajectory_lump0_tilt_theta=0 \
  --extra-override trajectory_lump0_tilt_phi=0 \
  --extra-override trajectory_lump0_v_rad=0 \
  --extra-override trajectory_lump0_omega_rot=0 \
  --extra-override trajectory_lump0_well_depth=0.15 \
  --extra-override trajectory_lump0_exotic=0

echo "[bondi] calibration cell complete: ${RUNS_DIR}/${out_name}"
