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
# mpirun for the GRTresna solve is auto-resolved from the conda/micromamba env
# that built GRTresna (see grtresna/solver.py::_resolve_mpirun); override with
# GRTRESNA_MPIRUN / GRTRESNA_ENV / GRTRESNA_ENV_NAME if needed.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

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
LUMPS="${LUMPS:-3}"                       # scalar lumps in the matter basis (10 dims each)
RANKS="${RANKS:-8}"                       # MPI ranks per GRTresna solve
ITERATIONS="${ITERATIONS:-30}"            # max non-linear iterations per solve
MAX_GENERATIONS="${MAX_GENERATIONS:-10}"
POPULATION="${POPULATION:-8}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
SEED="${SEED:-7}"
SIGMA0="${SIGMA0:-0.3}"

# Evolution / scoring knobs.
CONSUMER_RADII="${CONSUMER_RADII:-4 8}"
FTL_L="${FTL_L:-8.0}"
ENABLE_FTL_SCORING="${ENABLE_FTL_SCORING:-1}"

# Slice-frame movie fields. chi/lapse/K are the trace/gauge quantities and stay
# near-trivial (and look identical across candidates) for weak momentum-carrying
# scalar matter; the cloud + momentum show up in phi/Pi (scalar field & its
# conjugate momentum), shift1 (frame dragging) and rho_req (energy density).
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-phi Pi chi shift1 rho_req}"

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
if [[ "${ENABLE_FTL_SCORING}" == "1" ]]; then
  PRE_ARGS+=(--score-weight ftl_shortcut=5.0 \
             --score-weight expansion_asymmetry=2.0 \
             --score-weight nonflat_geometry=1.0 \
             --score-weight comoving_stability=2.5)
fi
PRE_ARGS+=(--ftl-L "${FTL_L:-8.0}")

mkdir -p "${RUNS_DIR}"

echo "== GRTresna-in-the-loop search =="
echo "GRTeclyn root : ${GRTECLYN_ROOT}"
echo "GRTresna root : ${GRTRESNA_ROOT}"
echo "Runs dir      : ${RUNS_DIR}"
echo "Lumps         : ${LUMPS}  (=> $((LUMPS * 10)) search dims)"
echo "Solve         : RANKS=${RANKS}  ITERATIONS=${ITERATIONS}"
echo "CMA-ES        : generations=${MAX_GENERATIONS} population=${POPULATION} sigma0=${SIGMA0} seed=${SEED}"
echo "GPUs          : ${GPU_IDS}"
echo "Consumer      : $([[ "${NO_CONSUME:-0}" == "1" ]] && echo DISABLED || echo "frames+delete radii=${CONSUMER_RADII}")"
echo "Frame fields  : ${GRTECLYN_FRAMES_FIELDS}"
echo "FTL scoring   : ENABLE_FTL_SCORING=${ENABLE_FTL_SCORING}  FTL_L=${FTL_L}"
echo

# shellcheck disable=SC2086
exec ${PYTHON_BIN} -m grteclyn_wrapper \
  "${PRE_ARGS[@]}" \
  optimize \
  --grtresna \
  --grtresna-lumps "${LUMPS}" \
  --grtresna-ranks "${RANKS}" \
  --grtresna-iterations "${ITERATIONS}" \
  --max-generations "${MAX_GENERATIONS}" \
  --population-size "${POPULATION}" \
  --sigma0 "${SIGMA0}" \
  --seed "${SEED}" \
  --gpu-ids ${GPU_IDS} \
  "${CONSUME_ARGS[@]}"
