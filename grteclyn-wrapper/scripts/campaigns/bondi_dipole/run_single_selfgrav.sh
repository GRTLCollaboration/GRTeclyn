#!/usr/bin/env bash
# Bondi dipole -- gravity-dressed calibration: ONE canonical lump, ultraweak rung.
#
# Why this rung and this seed (bondi_dipole_debug.md §7, 2026-08-06 follow-up):
# the exact-amplitude FLAT Q-ball at the weak rung (lam=2560) still blew off
# (rms 5.22 -> 6.6 by t=20, gate broken) because at that rung's compactness
# (2E/R ~ 0.14) the gravitationally dressed heavy-branch star DOES NOT EXIST --
# the dressed family maps to feather-light stars (ADM ~ 0.05) or ultra-compact
# unstable ones (ADM ~ 2.3), nothing near the flat ball's E=0.28.  The lump had
# no equilibrium to settle into.
#
# Two rungs up (lam=10240, mu=21845333; lam^2/mu = 4.8 so omega_min = 0.316
# unchanged), the dressed omega=0.55 star exists: phi_c = 0.019695, ADM =
# 0.0640, alpha(0) = 0.977.  This run seeds THAT star via the fixed-frequency
# self-grav ODE solve (grtresna_bs_selfgrav=1; phi_c is a solved output, the
# 3-column table carries the lapse, and the painters' amp/phi_c rescale is
# exactly 1).  Expect breathing ~1%, not the +-10-18% of the flat seeds.
#
# Launch-time verification (t=0 row of small_data/sector_barycenters.dat):
#   total_canon ~= 15.92, rms ~= 5.05  -> dressed star painted
# Ground truth earlier: grtresna/params.txt lump0_amp ~= 0.019695 and
# bs_omega ~= 0.55; qball_profile.dat header says "self-gravitating".
#
# Pass gate unchanged: rms flat (+-10%) to t=40 AND min_chi plateau above 0.3.
#
# Usage:
#   bash scripts/campaigns/bondi_dipole/run_single_selfgrav.sh
# Overrides: BONDI_GPU (default 1), BONDI_STOP_TIME (default 40),
#            BONDI_RUNS_DIR, BONDI_SEP (default 8), BONDI_EXOTIC=1 (lone
#            phantom star single_m -- reads the PHANTOM columns downstream).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi_dipole_selfgrav_v1}"
GPU="${BONDI_GPU:-1}"
STOP_TIME="${BONDI_STOP_TIME:-40}"
SEP="${BONDI_SEP:-8}"
EXOTIC="${BONDI_EXOTIC:-0}"
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

out_name="bondi_sg_single_p"
if [[ "${EXOTIC}" == "1" ]]; then out_name="bondi_sg_single_m"; fi
if [[ -d "${RUNS_DIR}/${out_name}" ]]; then
  echo "[bondi] ${out_name} already exists -- delete it or set BONDI_RUNS_DIR"
  exit 1
fi

echo "[bondi] === ${out_name} (ultraweak rung, gravity-dressed star) ==="
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
  --extra-override trajectory_num_lumps=1 \
  --extra-override trajectory_lump0_R0="${R0}" \
  --extra-override trajectory_lump0_phase0=0 \
  --extra-override trajectory_lump0_tilt_theta=0 \
  --extra-override trajectory_lump0_tilt_phi=0 \
  --extra-override trajectory_lump0_v_rad=0 \
  --extra-override trajectory_lump0_omega_rot=0 \
  --extra-override trajectory_lump0_well_depth=0.15 \
  --extra-override trajectory_lump0_exotic="${EXOTIC}"

echo "[bondi] dressed-star calibration complete: ${RUNS_DIR}/${out_name}"
