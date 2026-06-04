#!/usr/bin/env bash
# GRTresna-in-the-loop CMA-ES search.
#
# Each evaluation:
#   CMA-ES proposes grtresna_lump{k}_* matter params
#     -> GRTresna 3D elliptic solve (constraint-satisfying, momentum-carrying ID)
#     -> initial_data.gridinit  -> GRTeclyn loads it -> short GPU evolution
#     -> plotfiles streamed into radial profiles + frames (raw plotfiles deleted)
#     -> score  -> back to CMA-ES
#
# The plotfile consumer is ON by default (frames are produced and the heavy raw
# plotfiles are deleted as they are consumed), so disk stays bounded across a
# long conveyor search. Pass NO_CONSUME=1 to keep raw plotfiles instead.
#
# Each evaluation runs GRTresna under MPI (default RANKS=8). mpirun is resolved
# from the conda/micromamba env that built the .MPI. executable
# (grtresna/solver.py::_resolve_mpirun). Override with GRTRESNA_MPIRUN /
# GRTRESNA_ENV / GRTRESNA_ENV_NAME, or set RANKS=1 for the serial binary only.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/env.sh"

# Resolve a Python that has the wrapper's deps (cma, h5py). env.sh cd's to
# GRTECLYN_ROOT, so a bare `uv run` would miss the grteclyn-wrapper project;
# prefer the wrapper's own venv, then uv with an explicit --project.
PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
    PYTHON_BIN="${WRAPPER_ROOT}/.venv/bin/python"
  elif command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
    PYTHON_BIN="uv run --project ${WRAPPER_ROOT} python"
  else
    PYTHON_BIN="python"
  fi
fi

# GRTresna source tree (holds Examples/ScalarFieldBH/Main_*.ex).
export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"

# ---- Search configuration (override via environment) -----------------------
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_search}"
LUMPS="${LUMPS:-5}"                       # scalar lumps in the matter basis (11 dims each)
GRTRESNA_ANSATZ="${GRTRESNA_ANSATZ:-ring}" # free=11*LUMPS dims, ring=14 template dims
RANKS="${RANKS:-8}"                       # MPI ranks per GRTresna solve (>1 => .MPI. executable)
ITERATIONS="${ITERATIONS:-30}"            # max non-linear iterations per solve
if [[ -z "${GRTRESNA_FULL_Z+x}" ]]; then
  if [[ "${GRTRESNA_ANSATZ}" == "shell" ]]; then
    GRTRESNA_FULL_Z=1
  else
    GRTRESNA_FULL_Z=0
  fi
fi
GRTRESNA_MAX_LEVEL="${GRTRESNA_MAX_LEVEL:-3}"
GRTRESNA_REFINE_THRESHOLD="${GRTRESNA_REFINE_THRESHOLD:-0.5}"
GRTRESNA_REGRID_RADIUS="${GRTRESNA_REGRID_RADIUS:-0}"
GRTRESNA_JACOBIAN_CAP="${GRTRESNA_JACOBIAN_CAP:-25.0}"
GRTRESNA_EVOLUTION_L_FULL="${GRTRESNA_EVOLUTION_L_FULL:-64.0}"
GRTRESNA_EVOLUTION_N_FULL="${GRTRESNA_EVOLUTION_N_FULL:-64}"
GRTRESNA_DOMAIN_L="${GRTRESNA_DOMAIN_L:-128.0}"
GRTRESNA_DOMAIN_NX="${GRTRESNA_DOMAIN_NX:-64}"
GRTRESNA_DOMAIN_NY="${GRTRESNA_DOMAIN_NY:-64}"
GRTRESNA_DOMAIN_NZ="${GRTRESNA_DOMAIN_NZ:-}"
GRTRESNA_GRIDINIT_NX="${GRTRESNA_GRIDINIT_NX:-64}"
GRTRESNA_GRIDINIT_NY="${GRTRESNA_GRIDINIT_NY:-64}"
GRTRESNA_GRIDINIT_NZ="${GRTRESNA_GRIDINIT_NZ:-64}"
MAX_GENERATIONS="${MAX_GENERATIONS:-50}"
GPU_IDS="${GPU_IDS:-0 1 2 3}"
POPULATION="${POPULATION:-$(wc -w <<< "${GPU_IDS}")}"
SEED="${SEED:-7}"
SIGMA0="${SIGMA0:-0.3}"
OBJECTIVE_MODE="${OBJECTIVE_MODE:-ftl_first}"
RANDOM_INJECTION_FRACTION="${RANDOM_INJECTION_FRACTION:-0.25}"
EXOTIC_INJECTION_FRACTION="${EXOTIC_INJECTION_FRACTION:-0.25}"
WARM_START_TOP_K="${WARM_START_TOP_K:-8}"
WARM_START_JITTER="${WARM_START_JITTER:-0.08}"

# Pre-GPU solved-geometry FTL gate policy.  These are intentionally centralized
# here so exploratory/strict runs do not require code edits.
SOLVED_FTL_F_OP_FLOOR="${SOLVED_FTL_F_OP_FLOOR:-1.0e-4}"
SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR="${SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR:-0.99}"
SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR="${SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR:-1.01}"
SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR="${SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR:-0.02}"
SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED="${SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED:-8.0}"
SOLVED_FTL_MAX_PHYSICAL_F_OP="${SOLVED_FTL_MAX_PHYSICAL_F_OP:-0.85}"
SOLVED_FTL_REJECTION_SPEED_TARGET="${SOLVED_FTL_REJECTION_SPEED_TARGET:-1.01}"

# Evolution / scoring knobs.
CONSUMER_RADII="${CONSUMER_RADII:-4 8}"
FTL_L="${FTL_L:-8.0}"
ENABLE_FTL_SCORING="${ENABLE_FTL_SCORING:-1}"

# Slice-frame movie fields. chi/lapse/K are the trace/gauge quantities and stay
# near-trivial (and look identical across candidates) for weak momentum-carrying
# scalar matter; the cloud + momentum show up in phi/Pi (scalar field & its
# conjugate momentum), scalar_activity, shift1 (frame dragging), rho_req
# (energy density) and local_speed (FTL precursor map).
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-phi Pi scalar_activity chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${PROJECTION_METHOD:-mip}"

CONSUME_ARGS=()
if [[ "${NO_CONSUME:-0}" == "1" ]]; then
  CONSUME_ARGS=(--no-consume-plotfiles)
fi

# Global (pre-subcommand) args. NOTE on ordering: --consumer-radii is greedy
# (nargs="+"), so it must not be the last global arg or it would swallow the
# "optimize" subcommand token. We always terminate the global section with the
# single-valued --ftl-L so the subcommand boundary is unambiguous.
PRE_ARGS=(--runs-dir "${RUNS_DIR}" --example RadialRecipe --consumer-radii ${CONSUMER_RADII})
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  PRE_ARGS+=(--dry-run)
fi
# FTL-focused reweighting. Two deliberate changes make the FTL objective
# actually reachable for this matter family:
#   * ftl_precursor (5.0): a continuous shaping gradient that ramps up as the
#     light cones tilt superluminal locally -- BEFORE operational_ftl's hard
#     end-to-end-channel gate fires -- so CMA-ES has a slope to climb out of the
#     flat-space basin instead of staring at a flat 0.
#   * exotic_penalty (2.0, down from the 8.0 default): we now ALLOW exotic
#     matter (per-lump phantom lumps source genuine negative energy) and accept
#     the resulting NEC violation as the price of a persistent FTL channel, so
#     the penalty must no longer dominate the budget and veto every exotic
#     geometry. operational_ftl (9.0) still rules: exotic that buys no FTL loses.
if [[ "${ENABLE_FTL_SCORING}" == "1" ]]; then
  PRE_ARGS+=(--score-weight ftl_shortcut=5.0 \
             --score-weight expansion_asymmetry=2.0 \
             --score-weight nonflat_geometry=1.0 \
             --score-weight comoving_stability=2.5 \
             --score-weight ftl_precursor="${FTL_PRECURSOR_WEIGHT:-5.0}" \
             --score-weight exotic_penalty="${EXOTIC_PENALTY_WEIGHT:-2.0}")
fi
PRE_ARGS+=(--ftl-L "${FTL_L:-8.0}")

mkdir -p "${RUNS_DIR}"

echo "== GRTresna-in-the-loop search =="
echo "GRTeclyn root : ${GRTECLYN_ROOT}"
echo "GRTresna root : ${GRTRESNA_ROOT}"
echo "Runs dir      : ${RUNS_DIR}"
if [[ "${GRTRESNA_ANSATZ}" == "ring" ]]; then
  echo "Lumps/ansatz  : ${LUMPS} generated from ring template (=> 14 search dims, planar)"
elif [[ "${GRTRESNA_ANSATZ}" == "shell" ]]; then
  echo "Lumps/ansatz  : ${LUMPS} generated from full-sphere shell template (=> 16 search dims, 3D discovery)"
else
  echo "Lumps/ansatz  : ${LUMPS} free lumps (=> $((LUMPS * 11)) search dims)"
fi
echo "Solve         : RANKS=${RANKS}  ITERATIONS=${ITERATIONS}  max_level=${GRTRESNA_MAX_LEVEL}"
echo "Domain        : $([[ "${GRTRESNA_FULL_Z}" == "1" ]] && echo "full-z (no reflective z=0 plane)" || echo "half-z reflective")"
echo "Domain sizes  : evolution L=${GRTRESNA_EVOLUTION_L_FULL} N=${GRTRESNA_EVOLUTION_N_FULL}; GRTresna L=${GRTRESNA_DOMAIN_L} N=${GRTRESNA_DOMAIN_NX},${GRTRESNA_DOMAIN_NY},${GRTRESNA_DOMAIN_NZ:-auto}; gridinit=${GRTRESNA_GRIDINIT_NX},${GRTRESNA_GRIDINIT_NY},${GRTRESNA_GRIDINIT_NZ}"
echo "AMR           : refine_threshold=${GRTRESNA_REFINE_THRESHOLD} regrid_radius=${GRTRESNA_REGRID_RADIUS}"
echo "CMA-ES        : generations=${MAX_GENERATIONS} population=${POPULATION} sigma0=${SIGMA0} seed=${SEED}"
echo "GPUs          : ${GPU_IDS}"
echo "Objective     : ${OBJECTIVE_MODE}  random_injection=${RANDOM_INJECTION_FRACTION} exotic_injection=${EXOTIC_INJECTION_FRACTION}"
echo "Warm start    : ${WARM_START_TRAJECTORY:-<none>}"
echo "Consumer      : $([[ "${NO_CONSUME:-0}" == "1" ]] && echo DISABLED || echo "frames+delete radii=${CONSUMER_RADII}")"
echo "Frame fields  : ${GRTECLYN_FRAMES_FIELDS}"
echo "Projections   : ${GRTECLYN_PROJECTION_FIELDS} axes=${GRTECLYN_PROJECTION_AXES} method=${GRTECLYN_PROJECTION_METHOD}"
echo "FTL scoring   : ENABLE_FTL_SCORING=${ENABLE_FTL_SCORING}  FTL_L=${FTL_L}"
echo "Solved gate   : F_op>${SOLVED_FTL_F_OP_FLOOR} near_max_c>=${SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR}(<1) super_max_c>=${SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR} frac>=${SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR}"
echo "Degenerate    : max_c>${SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED} or F_op>${SOLVED_FTL_MAX_PHYSICAL_F_OP}"
echo

WARM_START_ARGS=()
if [[ -n "${WARM_START_TRAJECTORY:-}" ]]; then
  IFS=',' read -r -a WARM_START_PATHS <<< "${WARM_START_TRAJECTORY}"
  for traj in "${WARM_START_PATHS[@]}"; do
    WARM_START_ARGS+=(--warm-start-trajectory "${traj}")
  done
fi

DOMAIN_ARGS=(
  --grtresna-evolution-l-full "${GRTRESNA_EVOLUTION_L_FULL}"
  --grtresna-evolution-n-full "${GRTRESNA_EVOLUTION_N_FULL}"
  --grtresna-domain-l "${GRTRESNA_DOMAIN_L}"
  --grtresna-domain-nx "${GRTRESNA_DOMAIN_NX}"
  --grtresna-domain-ny "${GRTRESNA_DOMAIN_NY}"
  --grtresna-gridinit-nx "${GRTRESNA_GRIDINIT_NX}"
  --grtresna-gridinit-ny "${GRTRESNA_GRIDINIT_NY}"
  --grtresna-gridinit-nz "${GRTRESNA_GRIDINIT_NZ}"
)
if [[ -n "${GRTRESNA_DOMAIN_NZ}" ]]; then
  DOMAIN_ARGS+=(--grtresna-domain-nz "${GRTRESNA_DOMAIN_NZ}")
fi

# shellcheck disable=SC2086
exec ${PYTHON_BIN} -m grteclyn_wrapper \
  "${PRE_ARGS[@]}" \
  optimize \
  --objective-mode "${OBJECTIVE_MODE}" \
  --random-injection-fraction "${RANDOM_INJECTION_FRACTION}" \
  --exotic-injection-fraction "${EXOTIC_INJECTION_FRACTION}" \
  --warm-start-top-k "${WARM_START_TOP_K}" \
  --warm-start-jitter "${WARM_START_JITTER}" \
  "${WARM_START_ARGS[@]}" \
  --grtresna \
  --grtresna-ansatz "${GRTRESNA_ANSATZ}" \
  --grtresna-lumps "${LUMPS}" \
  "$([[ "${GRTRESNA_FULL_Z}" == "1" ]] && echo --grtresna-full-z || echo --no-grtresna-full-z)" \
  "${DOMAIN_ARGS[@]}" \
  --grtresna-ranks "${RANKS}" \
  --grtresna-iterations "${ITERATIONS}" \
  --grtresna-max-level "${GRTRESNA_MAX_LEVEL}" \
  --grtresna-refine-threshold "${GRTRESNA_REFINE_THRESHOLD}" \
  --grtresna-regrid-radius "${GRTRESNA_REGRID_RADIUS}" \
  --grtresna-jacobian-cap "${GRTRESNA_JACOBIAN_CAP}" \
  --solved-ftl-f-op-floor "${SOLVED_FTL_F_OP_FLOOR}" \
  --solved-ftl-near-luminal-speed-floor "${SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR}" \
  --solved-ftl-superluminal-speed-floor "${SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR}" \
  --solved-ftl-superluminal-fraction-floor "${SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR}" \
  --solved-ftl-max-physical-coord-speed "${SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED}" \
  --solved-ftl-max-physical-f-op "${SOLVED_FTL_MAX_PHYSICAL_F_OP}" \
  --solved-ftl-rejection-speed-target "${SOLVED_FTL_REJECTION_SPEED_TARGET}" \
  --max-generations "${MAX_GENERATIONS}" \
  --population-size "${POPULATION}" \
  --sigma0 "${SIGMA0}" \
  --seed "${SEED}" \
  --gpu-ids ${GPU_IDS} \
  "${CONSUME_ARGS[@]}"
