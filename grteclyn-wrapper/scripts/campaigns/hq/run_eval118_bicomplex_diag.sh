#!/usr/bin/env bash
# Phase-5 diagnostic: re-solve+replay eval-118 genome under bicomplex.
#
# Two variants at search resolution (L=64, N=128, t=16, max_level=1):
#   1) pump-on   — search-tier sanity (pump stays on for full stop_time)
#   2) pump-off  — transient igniter (RL_PUMP_STOP_TIME); f_geo from post-stop rays
#
# Uses the MPI+CUDA RadialRecipe binary (main3d.gnu.MPI.CUDA.ex) via
# resolve_executable preference; single-rank unless EVOLUTION_MPI_RANKS>1.
#
# Usage:
#   bash scripts/campaigns/hq/run_eval118_bicomplex_diag.sh            # both
#   bash scripts/campaigns/hq/run_eval118_bicomplex_diag.sh pump-on
#   bash scripts/campaigns/hq/run_eval118_bicomplex_diag.sh pump-off
#   DRY_RUN=1 bash scripts/campaigns/hq/run_eval118_bicomplex_diag.sh
set -euo pipefail

# NOTE: env.sh overwrites SCRIPT_DIR; keep launcher path under HQ_DIR.
HQ_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${HQ_DIR}/.." && pwd)"
SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
# shellcheck source=../../lib/env.sh
source "${SCRIPTS_ROOT}/lib/env.sh"

VARIANT="${1:-both}"
SOURCE_EVAL="${SOURCE_EVAL:-${GRTECLYN_ROOT}/runs/grtresna_qd/qball_traj_spiral_v2/eval_000118}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_promote}"
GPU_ID="${GPU_ID:-0}"
EVOLUTION_MPI_RANKS="${EVOLUTION_MPI_RANKS:-1}"
N_FULL="${N_FULL:-128}"
L_FULL="${L_FULL:-64}"
STOP_TIME="${STOP_TIME:-16.0}"
MAX_LEVEL="${MAX_LEVEL:-1}"
REGRID_THRESHOLD="${REGRID_THRESHOLD:-0.1}"
PLOT_INTERVAL="${PLOT_INTERVAL:-320}"
RL_PUMP_STOP_TIME="${RL_PUMP_STOP_TIME:-4.0}"
DRY_RUN="${DRY_RUN:-0}"

if [[ ! -d "${SOURCE_EVAL}" ]]; then
  echo "[diag] missing source eval: ${SOURCE_EVAL}" >&2
  exit 2
fi

REPLAY="${HQ_DIR}/replay_eval.py"
COMMON_ARGS=(
  "${SOURCE_EVAL}"
  --runs-dir "${RUNS_DIR}"
  --gpu "${GPU_ID}"
  --evolution-mpi-ranks "${EVOLUTION_MPI_RANKS}"
  --n-full "${N_FULL}"
  --l-full "${L_FULL}"
  --grtresna-domain-l "${L_FULL}"
  --stop-time "${STOP_TIME}"
  --max-level "${MAX_LEVEL}"
  --regrid-threshold "${REGRID_THRESHOLD}"
  --plot-interval "${PLOT_INTERVAL}"
  --objective-mode general_ftl
  --evolving-geodesic
  --qball-ode-profile
  --qball-equilibrium-amplitude
  --extra-override grtresna_matter_model=grtresna_bicomplex_scalar
  --extra-override grtresna_matter_coupling=canonical
)

export GRTRESNA_ALLOW_SIGN_MISMATCH=0
export OBJECTIVE_MODE=general_ftl
export GRTECLYN_EVOLVING_GEODESIC=1
export GRTECLYN_GEO_DIRECTIONS="${GRTECLYN_GEO_DIRECTIONS:-x y z}"
export GRTECLYN_GEO_EMIT_INTERVAL="${GRTECLYN_GEO_EMIT_INTERVAL:-2}"
export GRTECLYN_GEO_MAX_EMISSIONS="${GRTECLYN_GEO_MAX_EMISSIONS:-7}"
export SCORE_PUMP_ENERGY_WEIGHT="${SCORE_PUMP_ENERGY_WEIGHT:-40}"
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-0}"
unset GRTECLYN_KEEP_PLOTFILES 2>/dev/null || true

run_one() {
  local name="$1"
  shift
  local log="${RUNS_DIR}/${name}.log"
  mkdir -p "${RUNS_DIR}"
  echo "[diag] === ${name} (mpi_ranks=${EVOLUTION_MPI_RANKS} gpu=${GPU_ID}) ===" | tee -a "${log}"
  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "[diag] dry-run: uv run python ${REPLAY} ${COMMON_ARGS[*]} --name ${name} $*" | tee -a "${log}"
    return 0
  fi
  uv run python "${REPLAY}" "${COMMON_ARGS[@]}" --name "${name}" "$@" \
    2>&1 | tee -a "${log}"
}

case "${VARIANT}" in
  pump-on|both)
    unset RL_PUMP_STOP_TIME || true
    export GRTECLYN_EVOLVING_GEODESIC_MODE=search
    run_one "e118_bicomplex_pumpon_L${L_FULL%.*}_N${N_FULL}_t${STOP_TIME%.*}"
    ;;&
  pump-off|both)
    export RL_PUMP_STOP_TIME
    export GRTECLYN_EVOLVING_GEODESIC_MODE=hq
    run_one \
      "e118_bicomplex_pumpoff_tstop${RL_PUMP_STOP_TIME%.*}_L${L_FULL%.*}_N${N_FULL}_t${STOP_TIME%.*}" \
      --extra-override "rl_pump_stop_time=${RL_PUMP_STOP_TIME}"
    ;;&
  pump-on|pump-off|both)
    ;;
  *)
    echo "usage: $0 [pump-on|pump-off|both]" >&2
    exit 2
    ;;
esac

echo "[diag] done (variant=${VARIANT}). Inspect metadata.json f_geo / pump_energy in ${RUNS_DIR}/e118_bicomplex_*"
