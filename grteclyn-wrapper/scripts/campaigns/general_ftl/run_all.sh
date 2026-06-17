#!/usr/bin/env bash
# Orchestrator: launch the general-FTL discovery campaigns (wormhole / ring /
# spinning) over 8 GPUs.  MODE=seq (default) runs each on all 8 GPUs in turn;
# MODE=par splits the 8 GPUs across the three classes.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
QD_RUN="${SCRIPT_DIR}/../qd/run.sh"

export OBJECTIVE_MODE=general_ftl
export DESCRIPTOR_MODE=ftl_lifetime
export GRTECLYN_GEO_DIRECTIONS="x y z"
export GRTRESNA_ANSATZ=shell
export STOP_TIME="${STOP_TIME:-16.0}"

# --- config matrix: NAME | PIN_DIMS ---------------------------------------
# Shared static-matter lock: shell_static=1 already zeros every current, so the
# toroidal/poloidal/radial/omega dims are inert and pinned out to avoid wasting
# optimizer variance on knobs with no physical effect.
STATIC_INERT="grtresna_shell_toroidal_velocity=0 grtresna_shell_poloidal_velocity=0 grtresna_shell_radial_velocity=0 grtresna_shell_omega=0"
# Wormhole: two mouths on the axis, axis aligned to +x (axis_theta=pi/2,
# axis_phi=0 -> a=(1,0,0)); static matter.
WORMHOLE_PINS="grtresna_matter_layout=2 grtresna_shell_axis_theta=1.5708 grtresna_shell_axis_phi=0 grtresna_shell_static=1 ${STATIC_INERT}"
# Toroidal waveguide: ring lies in the plane orthogonal to the (searched)
# polar axis; static matter.  The x/y/z probe scan covers any orientation.
RING_PINS="grtresna_matter_layout=3 grtresna_shell_static=1 ${STATIC_INERT}"
# Spinning frame-drag: NOT static — translation velocities + shift pinned to 0,
# spin (shell_omega) left free.  shell_static is pinned to 0 so the optimizer
# cannot flip the static toggle and silently zero omega.
SPIN_PINS="grtresna_matter_layout=0 grtresna_shell_static=0 grtresna_shell_toroidal_velocity=0 grtresna_shell_poloidal_velocity=0 grtresna_shell_radial_velocity=0 grtresna_shift_seed=0"

run_one () {  # name  pins  gpu_ids
  local name="$1" pins="$2" gpus="$3"
  RUNS_DIR="${GRTECLYN_ROOT:-$PWD}/runs/general_ftl_${name}" \
  QD_NAME="general_ftl_${name}" \
  PIN_DIMS="${pins}" \
  GPU_IDS="${gpus}" \
  QD_ITERATIONS="${QD_ITERATIONS:-30}" \
  bash "${QD_RUN}"
}

MODE="${MODE:-seq}"
if [[ "${MODE}" == "par" ]]; then
  # Run the preflight pytest gate once here, then skip it inside each concurrent
  # campaign so the three launches don't redundantly re-run the same suite.
  if [[ "${SKIP_QD_PREFLIGHT_TESTS:-0}" != "1" ]]; then
    source "${SCRIPT_DIR}/../lib/bootstrap.sh"
    _campaign_bootstrap "${SCRIPT_DIR}"
    _campaign_resolve_python
    ftl_search_common_preflight_tests
  fi
  export SKIP_QD_PREFLIGHT_TESTS=1
  run_one wormhole "${WORMHOLE_PINS}" "0 1 2" &
  run_one ring     "${RING_PINS}"     "3 4 5" &
  run_one spin     "${SPIN_PINS}"     "6 7"   &
  wait
else
  run_one wormhole "${WORMHOLE_PINS}" "0 1 2 3 4 5 6 7"
  run_one ring     "${RING_PINS}"     "0 1 2 3 4 5 6 7"
  run_one spin     "${SPIN_PINS}"     "0 1 2 3 4 5 6 7"
fi
