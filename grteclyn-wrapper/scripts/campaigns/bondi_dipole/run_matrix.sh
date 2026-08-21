#!/usr/bin/env bash
# Bondi dipole runaway in full NR -- the 2x2 control matrix (N1, Appendix B).
#
# Bondi (1957): a positive/negative-mass pair self-accelerates -- the phantom
# lump falls toward the canonical lump's well while the canonical lump rolls
# off the phantom's hill, so both accelerate the SAME way, phantom chasing
# canonical, with P_ADM ~ 0.  Never evolved in full 3+1 NR with dynamical,
# constraint-solved matter.  The bicomplex model realizes the required sign
# structure exactly: both sectors obey the same Klein-Gordon equation
# (positive inertial/passive mass); the sign flip sits only in the Einstein
# source (negative active mass).
#
# The matrix (each cell a falsifiable prediction):
#   pair_pm  (+,-) : runs away along +x (phantom at -3 chases canonical at +3)
#   pair_pp  (+,+) : attracts -- merges or orbits, no net drift
#   pair_mm  (-,-) : mutually repels (each digs a hill the other rolls off)
#   single_p / single_m : drift nowhere; calibrate per-sector dispersal rates
#
# Setup per Appendix B: two lumps AT REST on the x axis at separation d=6
# (R0=3, phase0=0 vs pi), evolution pump DISABLED from t=0
# (rl_pump_stop_time=0 -- well_depth stays nonzero only as the GRTresna seed
# amplitude), no breathing / z-motion / boosts.  Pure self-gravity.
#
# Diagnostic: per-sector barycentres, streamed live to sector_barycenters.dat
# (GRTECLYN_SECTOR_BARYCENTERS=1).  The aggregate barycentre cancels for the
# mixed pair, and plotfiles are purged after consumption -- this cannot be
# recovered post hoc.  Psi4 stays on: what an accelerating dipole radiates is
# one of the open questions.
#
# Seed scaffolding: the git-tracked v1 champion eval (stable, never pruned);
# every physics knob that matters is overridden below.  Sequential on ONE GPU
# so the depth campaign keeps three slots clean.
#
# Usage:
#   bash scripts/campaigns/bondi_dipole/run_matrix.sh
# Overrides: BONDI_GPU (default 0), BONDI_STOP_TIME (default 60),
#            BONDI_RUNS_DIR, BONDI_SEP (default 6).
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi/dipole_v1}"
GPU="${BONDI_GPU:-0}"
STOP_TIME="${BONDI_STOP_TIME:-60}"
SEP="${BONDI_SEP:-6}"
R0="$(python3 -c "print(${SEP}/2)")"
PI="3.141592653589793"

mkdir -p "${RUNS_DIR}"

# Stop handle for scripts/campaigns/stop_campaign.sh ($! at launch is the dead
# setsid parent -- only the launcher itself knows its true pid).
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
  "pair_pm  2 0 1"
  "pair_pp  2 0 0"
  "pair_mm  2 1 1"
  "single_p 1 0 0"
  "single_m 1 1 0"
)

for spec in "${matrix[@]}"; do
  read -r name num_lumps exotic0 exotic1 <<<"${spec}"
  out_name="bondi_${name}"
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
    --max-level 1 --regrid-threshold 0.1 \
    --stop-time "${STOP_TIME}" --plot-interval 160 \
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
