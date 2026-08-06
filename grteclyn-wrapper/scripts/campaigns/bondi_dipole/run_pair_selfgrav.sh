#!/usr/bin/env bash
# Bondi dipole -- gravity-dressed PAIR cell, ultraweak rung.
#
# Two self-gravitating omega=0.55 sextic stars at rest, separation SEP along x
# (lump0 at +R0, lump1 at -R0).  Sector per lump via BONDI_S0/BONDI_S1
# (0=canonical, 1=phantom):
#   S0=0 S1=1 -> pair_pm  (the runaway cell: canonical at +x, phantom at -x;
#                          negative mass chases positive -> both drift +x)
#   S0=0 S1=0 -> pair_pp  (control: mutual attraction, no net drift)
#   S0=1 S1=1 -> pair_mm  (control: mutual repulsion, no net drift)
#
# Star data (bondi_dipole_debug.md §8): canonical phi_c=0.019695 ADM=+0.0640,
# phantom phi_c=0.019589 ADM=-0.0770.  |M-| > |M+| by ~20% -- expect imperfect
# momentum cancellation in pair_pm.
#
# Launch-time verification (t=0 row of small_data/sector_barycenters.dat):
#   canonical lump present: total_canon ~= 15.92, rms_canon ~= 5.05
#   phantom  lump present: total_phantom ~= 20.99, rms_phantom ~= 5.43
# Runaway signal: canon x (col 3) and phantom x (col 8) both drifting the
# same direction, canonical side leading.
#
# Default stop 60, not 100: the massive-field radiation bath cannot exit the
# massless-wave boundaries, and chi drifted 0.99 -> 0.86 by t=40 in the
# singles -- beyond t~60 the bath starts to poison the drift measurement.
#
# Usage:
#   BONDI_S0=0 BONDI_S1=1 BONDI_GPU=3 bash scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
# Overrides: BONDI_GPU (default 3), BONDI_STOP_TIME (default 60),
#            BONDI_RUNS_DIR, BONDI_SEP (default 8), BONDI_S0/BONDI_S1.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
GPU="${BONDI_GPU:-3}"
STOP_TIME="${BONDI_STOP_TIME:-60}"
SEP="${BONDI_SEP:-8}"
S0="${BONDI_S0:-0}"
S1="${BONDI_S1:-1}"
R0="$(python3 -c "print(${SEP}/2)")"
PI="3.141592653589793"

sector_tag() { if [[ "$1" == "1" ]]; then echo m; else echo p; fi; }
suffix="$(sector_tag "${S0}")$(sector_tag "${S1}")"
out_name="bondi_sg_pair_${suffix}"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi_dipole_selfgrav_${suffix}}"

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

if [[ -d "${RUNS_DIR}/${out_name}" ]]; then
  echo "[bondi] ${out_name} already exists -- delete it or set BONDI_RUNS_DIR"
  exit 1
fi

echo "[bondi] === ${out_name} (ultraweak rung, gravity-dressed pair s0=${S0} s1=${S1}) ==="
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
  --grtresna-nl-exit-tolerance 0.1 \
  --grtresna-nl-stall-tolerance 0.002 \
  --grtresna-max-level 3 \
  --grtresna-domain-l 128 \
  --grtresna-timeout 7200 \
  --consumer-radii 8 16 \
  --consumer-keep-last 2 \
  --objective-mode weighted \
  --extra-override grtresna_scalar_lambda=10240 \
  --extra-override grtresna_scalar_mu=21845333 \
  --extra-override grtresna_bs_omega=0.55 \
  --extra-override grtresna_bs_selfgrav=1 \
  --extra-override trajectory_well_width=1.2 \
  --extra-override rl_pump_stop_time=0 \
  --extra-override grtresna_boost_lumps=0 \
  --extra-override trajectory_A_breath=0 \
  --extra-override trajectory_z_amp=0 \
  --extra-override trajectory_omega_z=0 \
  --extra-override trajectory_num_lumps=2 \
  --extra-override trajectory_lump0_R0="${R0}" \
  --extra-override trajectory_lump0_phase0=0 \
  --extra-override trajectory_lump0_tilt_theta=0 \
  --extra-override trajectory_lump0_tilt_phi=0 \
  --extra-override trajectory_lump0_v_rad=0 \
  --extra-override trajectory_lump0_omega_rot=0 \
  --extra-override trajectory_lump0_well_depth=0.15 \
  --extra-override trajectory_lump0_exotic="${S0}" \
  --extra-override trajectory_lump1_R0="${R0}" \
  --extra-override trajectory_lump1_phase0="${PI}" \
  --extra-override trajectory_lump1_tilt_theta=0 \
  --extra-override trajectory_lump1_tilt_phi=0 \
  --extra-override trajectory_lump1_v_rad=0 \
  --extra-override trajectory_lump1_omega_rot=0 \
  --extra-override trajectory_lump1_well_depth=0.15 \
  --extra-override trajectory_lump1_exotic="${S1}"

echo "[bondi] pair cell complete: ${RUNS_DIR}/${out_name}"
