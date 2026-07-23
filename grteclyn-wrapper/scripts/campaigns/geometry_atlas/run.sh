#!/usr/bin/env bash
# Pure-geometry MAP-Elites atlas — Stage-1 inverse-design scout.
#
# Searches smooth stationary asymptotically flat 4-metrics (no matter, no
# GRTresna, no GRTeclyn evolution).  Scores frozen null-geodesic f_geo and
# stationary free-fall f_ff, archiving shortcut strength vs exotic-energy cost.
#
# Usage (CPU smoke):
#   bash scripts/campaigns/geometry_atlas/run.sh
#
# Probe calibration (Alcubierre + optional cand.146):
#   MODE=calibrate bash scripts/campaigns/geometry_atlas/run.sh
#
# Focused CMA-ES refine (recommended before a wide MAP-Elites hunt):
#   MODE=cmaes TARGET_EVALS=40 N=32 NO_FF=1 bash scripts/campaigns/geometry_atlas/run.sh
#
# Larger MAP-Elites campaign:
#   MODE=map_elites TARGET_EVALS=200 N=32 BINS=8 bash scripts/campaigns/geometry_atlas/run.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/bootstrap.sh
source "${SCRIPT_DIR}/../lib/bootstrap.sh"
_campaign_bootstrap "${SCRIPT_DIR}"
_campaign_resolve_python

RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/geometry_atlas}"
MODE="${MODE:-map_elites}"
NAME="${NAME:-geometry_atlas_${MODE}_$(date -u +%Y%m%dT%H%M%SZ)}"
TARGET_EVALS="${TARGET_EVALS:-16}"
BINS="${BINS:-8}"
SEED="${SEED:-7}"
BATCH_SIZE="${BATCH_SIZE:-2}"
N="${N:-24}"
L="${L:-64.0}"
N_CENTERS="${N_CENTERS:-7}"
SUPPORT_RADIUS="${SUPPORT_RADIUS:-12.0}"
N_RAYS="${N_RAYS:-3}"
RESUME="${RESUME:-0}"
NO_FF="${NO_FF:-0}"
NO_ALC="${NO_ALC:-0}"
FULLBOX_PROBE="${FULLBOX_PROBE:-0}"
ALC_ONLY="${ALC_ONLY:-0}"
ALC_VELOCITY_MAX="${ALC_VELOCITY_MAX:-2.0}"
CMA_OBJECTIVE="${CMA_OBJECTIVE:-f_geo}"
CMA_SIGMA0="${CMA_SIGMA0:-0.25}"
CMA_POPSIZE="${CMA_POPSIZE:-}"
SEED_GENOME="${SEED_GENOME:-}"

mkdir -p "${RUNS_DIR}"

EXTRA_ARGS=()
[[ "${RESUME}" == "1" ]] && EXTRA_ARGS+=(--resume)
[[ "${NO_FF}" == "1" ]] && EXTRA_ARGS+=(--no-ff)
[[ "${NO_ALC}" == "1" ]] && EXTRA_ARGS+=(--no-alcubierre)
[[ "${FULLBOX_PROBE}" == "1" ]] && EXTRA_ARGS+=(--fullbox-probe)
[[ "${ALC_ONLY}" == "1" ]] && EXTRA_ARGS+=(--alc-only)
EXTRA_ARGS+=(--alc-velocity-max "${ALC_VELOCITY_MAX}")
EXTRA_ARGS+=(--cma-objective "${CMA_OBJECTIVE}")
EXTRA_ARGS+=(--cma-sigma0 "${CMA_SIGMA0}")
[[ -n "${CMA_POPSIZE}" ]] && EXTRA_ARGS+=(--cma-popsize "${CMA_POPSIZE}")
[[ -n "${SEED_GENOME}" ]] && EXTRA_ARGS+=(--seed-genome "${SEED_GENOME}")

# shellcheck disable=SC2086
exec ${PYTHON_BIN} -m grteclyn_wrapper \
  --runs-dir "${RUNS_DIR}" \
  --name "${NAME}" \
  geometry_atlas \
  --mode "${MODE}" \
  --target-evals "${TARGET_EVALS}" \
  --bins "${BINS}" \
  --seed "${SEED}" \
  --batch-size "${BATCH_SIZE}" \
  --n "${N}" \
  --L "${L}" \
  --n-centers "${N_CENTERS}" \
  --support-radius "${SUPPORT_RADIUS}" \
  --n-rays "${N_RAYS}" \
  "${EXTRA_ARGS[@]}"
