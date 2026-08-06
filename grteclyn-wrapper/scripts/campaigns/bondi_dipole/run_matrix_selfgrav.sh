#!/usr/bin/env bash
# Bondi dipole runaway -- GRAVITY-DRESSED matrix, ultraweak rung, dense frames.
#
# The five-cell control matrix (bondi_dipole_debug.md §1) on true equilibria:
# each lump is a self-gravitating sextic star solved at fixed frequency
# omega=0.55 (bondi_dipole_debug.md §8).  Canonical star: phi_c=0.019695,
# ADM=+0.0640, alpha(0)=0.977.  Phantom star (repulsive self-gravity, its own
# table): phi_c=0.019589, ADM=-0.0770, alpha(0)=1.025.  Note |M-| > |M+| by
# ~20% -- the dressed sectors are NOT mirror images; expect imperfect
# momentum cancellation in the (+,-) drift.
#
# Rung: lam=10240, mu=21845333 (lam^2/mu = 4.8, omega_min=0.316) -- the first
# rung where the dressed omega=0.55 star exists at all (§8.4).  Two tables are
# emitted per mixed cell (qball_profile.dat + qball_profile_exotic.dat) with
# per-lump profile_path params; both sectors share omega, so the painters'
# global-omega momentum stays exact.
#
# Cell order: single_p FIRST as live calibration (it passed the stop-40 gate
# as bondi_sg_single_p before this matrix existed), then the pair cells, then
# single_m -- the never-yet-run lone phantom.
#
# Launch-time verification (t=0 row of sector_barycenters.dat):
#   single_p / pair cells: total_canon ~= 15.92, rms_canon ~= 5.05
#   cells with a phantom lump: total_phantom ~= 15.8, rms_phantom ~= 5.1
#   (phantom slightly puffed; both far from the flat-seed 34.65/36.30 numbers)
#
# Usage:
#   bash scripts/campaigns/bondi_dipole/run_matrix_selfgrav.sh
# Overrides: BONDI_GPU (default 1), BONDI_STOP_TIME (default 100),
#            BONDI_RUNS_DIR, BONDI_SEP (default 8).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi_dipole_selfgrav_matrix_v1}"
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

# Lumps at rest: no orbital speed, no radial speed, no breathing, no z-motion,
# no initial-data boost; trajectory pump off for the whole evolution (t >= 0).
common_overrides=(
  --extra-override grtresna_scalar_lambda=10240
  --extra-override grtresna_scalar_mu=21845333
  --extra-override grtresna_bs_omega=0.55
  --extra-override grtresna_bs_selfgrav=1
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
matrix=(
  "single_p 1 0 0"
  "pair_pm  2 0 1"
  "pair_pp  2 0 0"
  "pair_mm  2 1 1"
  "single_m 1 1 0"
)

for spec in "${matrix[@]}"; do
  read -r name num_lumps exotic0 exotic1 <<<"${spec}"
  out_name="bondi_sgm_${name}"
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
  --grtresna-nl-exit-tolerance 0.1 \
  --grtresna-nl-stall-tolerance 0.002 \
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
