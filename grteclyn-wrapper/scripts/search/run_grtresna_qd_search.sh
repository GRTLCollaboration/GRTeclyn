#!/usr/bin/env bash
# GRTresna-in-the-loop MAP-Elites search.
#
# This is the quality-diversity counterpart to run_grtresna_search.sh: it keeps
# a diverse archive of shell candidates instead of collapsing immediately onto
# the highest scalar CMA-ES score.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/env.sh"

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

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"

RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_qd}"
LUMPS="${LUMPS:-5}"
SHELL_PROFILE="${SHELL_PROFILE:-compact}"
RANKS="${RANKS:-8}"
# Raised from 30: oscillating near-misses (Ham~7%/Mom~15%) need more iterations
# to settle below the 5% gate instead of being cut off early.
ITERATIONS="${ITERATIONS:-50}"
# GRTresna early-exit: stop once Ham/Mom are good (%%) or improvement stalls.
# Stall tolerance tightened from 0.02 so near-converging solves are not
# abandoned prematurely on a small residual plateau.
GRTRESNA_NL_EXIT_TOLERANCE="${GRTRESNA_NL_EXIT_TOLERANCE:-1.0}"
GRTRESNA_NL_STALL_TOLERANCE="${GRTRESNA_NL_STALL_TOLERANCE:-0.005}"
# Parallel Chombo→gridinit conversion (0 = auto, min(32, cpu_count)).
GRTRESNA_GRIDINIT_WORKERS="${GRTRESNA_GRIDINIT_WORKERS:-0}"
GRTRESNA_MAX_LEVEL="${GRTRESNA_MAX_LEVEL:-3}"
GRTRESNA_REFINE_THRESHOLD="${GRTRESNA_REFINE_THRESHOLD:-0.5}"
GRTRESNA_REGRID_RADIUS="${GRTRESNA_REGRID_RADIUS:-0}"
GRTRESNA_JACOBIAN_CAP="${GRTRESNA_JACOBIAN_CAP:-25.0}"
GRTRESNA_TIMEOUT="${GRTRESNA_TIMEOUT:-900}"
# Evolution base resolution: dx = L_full / N_full = 64 / 128 = 0.5.  Raised from
# N_full=64 (dx=1.0): the HQ-promotion campaign repeatedly showed QD "winners"
# whose superluminal coordinate channel is an *under-resolved* artifact of the
# dx=1.0 base grid -- it evaporates the moment the same data is evolved at dx=0.5
# (eval 231 f_geo 5.30%->0, eval 489 3.57%->0).  Evolving the QD loop itself at
# dx=0.5 stops the search from banking those artifacts, at ~16x compute/eval
# (8x cells + 2x timesteps from CFL).  Domain L=64 is kept so the box boundary
# stays at r=32, far from the r=8 FTL probe (shrinking L would risk boundary
# contamination).  max_level stays 2 (finest dx 0.125): the control experiments
# showed ml=3 vs ml=2 does NOT change the shortcut (eval 231: f_geo 4.02% vs
# 4.01% at t=16) -- the real variable behind the QD->HQ collapse is *time*
# (diffusion by t~30), which the time-averaged scoring now captures directly.
# ml=2 is therefore the right call: same physics, far cheaper compute, and
# lighter plotfiles -> much less NFS I/O for the in-flight FTL probe.
GRTRESNA_EVOLUTION_L_FULL="${GRTRESNA_EVOLUTION_L_FULL:-64.0}"
GRTRESNA_EVOLUTION_N_FULL="${GRTRESNA_EVOLUTION_N_FULL:-128}"
GRTRESNA_EVOLUTION_MAX_LEVEL="${GRTRESNA_EVOLUTION_MAX_LEVEL:-2}"
GRTRESNA_DOMAIN_L="${GRTRESNA_DOMAIN_L:-128.0}"
GRTRESNA_DOMAIN_NX="${GRTRESNA_DOMAIN_NX:-64}"
GRTRESNA_DOMAIN_NY="${GRTRESNA_DOMAIN_NY:-64}"
GRTRESNA_DOMAIN_NZ="${GRTRESNA_DOMAIN_NZ:-}"
# Gridinit (initial-data) grid bumped 64->128 to match the dx=0.5 evolution base
# grid; otherwise the solved (phi, Pi) would be loaded at dx=1.0 and interpolated
# up, throwing away the fidelity the finer evolution is meant to provide.
GRTRESNA_GRIDINIT_NX="${GRTRESNA_GRIDINIT_NX:-128}"
GRTRESNA_GRIDINIT_NY="${GRTRESNA_GRIDINIT_NY:-128}"
GRTRESNA_GRIDINIT_NZ="${GRTRESNA_GRIDINIT_NZ:-128}"
QD_ITERATIONS="${QD_ITERATIONS:-10}"
QD_TARGET_EVALS="${QD_TARGET_EVALS:-}"
QD_RESUME="${QD_RESUME:-0}"
QD_NAME="${QD_NAME:-}"
# MAP-Elites behaviour grid: ftl_lifetime (peak gauge-invariant strength x
# FTL-lifetime fraction -- the time-resolved grid that separates a transient,
# Alcubierre-like shortcut from a sustained one, now that the per-frame FTL
# stream is scored), channel (path x mechanism, needs shift>0), speed_horizon
# (cone-tilt x horizon-free), speed_super (cone-tilt x superluminal_fraction),
# legacy.
DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-ftl_lifetime}"
# Optional: warm-start the initial population from prior eval dirs (survivors).
SEED_EVAL_DIRS="${SEED_EVAL_DIRS:-}"
# Keep disk bounded: retain top-N by total score plus one champion dir per FTL
# peak metric (f_geo, speed, lifetime, ...); see ftl_retention.jsonl.
QD_KEEP_TOP_EVAL_DIRS="${QD_KEEP_TOP_EVAL_DIRS:-10}"
QD_FTL_RETENTION="${QD_FTL_RETENTION:-1}"
BINS="${BINS:-8}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
SEED="${SEED:-7}"
# Long enough for structural dissipation AND late-time instability to show: at
# t=8 a survivor has plateaued vs a dissipator, but the HQ campaign revealed
# QD-resolution "survivors" that still go unstable later (areal-radius drift,
# chi drop) -- the QD window was too short to catch them, so they were promoted
# only to fail at HQ.  Extending to t=16 lets the persistence/stability metrics
# observe the second half of the evolution and reject those late-unstable
# candidates *in the QD loop* instead of wasting an HQ promotion on them.
STOP_TIME="${STOP_TIME:-16.0}"
# Plotfile cadence in coarse steps.  dt = dt_multiplier * dx = 0.02 * 0.5 = 0.01
# at the dx=0.5 base grid, so STOP_TIME=16 is now ~1600 coarse steps (was 800 at
# dx=1.0).  PLOT_INTERVAL bumped 160->320 to keep ~6 plotfiles (t=0,3.2,..,16) ->
# ~6 frames/field, the same post-processing budget as before (late-time AMR
# refinement makes each yt slice frame expensive and they render in a serial
# burst at the batch boundary with the GPUs idle, so frame *count*, not
# stop_time, sets that gap).  Still >= the 3 plotfiles the evolved/geodesic FTL
# probes need.
PLOT_INTERVAL="${PLOT_INTERVAL:-320}"

SOLVED_FTL_F_OP_FLOOR="${SOLVED_FTL_F_OP_FLOOR:-1.0e-4}"
SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR="${SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR:-0.95}"
SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR="${SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR:-1.01}"
SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR="${SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR:-0.02}"
SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED="${SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED:-8.0}"
SOLVED_FTL_MAX_PHYSICAL_F_OP="${SOLVED_FTL_MAX_PHYSICAL_F_OP:-0.85}"
SOLVED_FTL_REJECTION_SPEED_TARGET="${SOLVED_FTL_REJECTION_SPEED_TARGET:-1.01}"
GRTRESNA_MAX_HAM_PCT="${GRTRESNA_MAX_HAM_PCT:-5.0}"
GRTRESNA_MAX_MOM_PCT="${GRTRESNA_MAX_MOM_PCT:-5.0}"

POSTLOAD_GATE="${POSTLOAD_GATE:-1}"
POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-1e-2}"
POSTLOAD_MAX_MOM_L2="${POSTLOAD_MAX_MOM_L2:-1e-2}"
export POSTLOAD_GATE

# Keep the GRTresna Chombo HDF5 + workdir per eval (conversion validation/debug).
GRTRESNA_KEEP_SOURCE="${GRTRESNA_KEEP_SOURCE:-0}"
KEEP_SOURCE_ARGS=()
if [[ "${GRTRESNA_KEEP_SOURCE}" == "1" ]]; then
  KEEP_SOURCE_ARGS+=(--grtresna-keep-source)
fi

CONSUMER_RADII="${CONSUMER_RADII:-4 8}"
# Retain the last few plotfiles so evolved/geodesic FTL and effective energy
# conditions can be scored (>=3 required); the rest are still consumed+deleted.
CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-3}"
FTL_L="${FTL_L:-8.0}"

# Frame/movie rendering OFF in the QD loop.  At max_level=3 each plotfile needs
# ~12 yt renders (9 slices + 3 MIP projections integrating the full refined 3D
# volume) at 10-70s each; with 8 evals consuming concurrently the render backlog
# cannot keep up during evolution and the per-eval consumers stall (observed:
# only the first finished eval ever produced frames).  QD scoring needs only the
# FTL time-series + scalar metrics (~6.5s/plotfile), so we skip frames here and
# render movies only for promoted candidates (run_promote_qd_batch.sh / manual
# make_movies.sh).  Set GRTECLYN_FRAMES=1 to re-enable for debugging.
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-0}"
# 4D evolving null-geodesic trace in the QD loop (fast search profile).
export GRTECLYN_EVOLVING_GEODESIC="${GRTECLYN_EVOLVING_GEODESIC:-1}"
export GRTECLYN_EVOLVING_GEODESIC_MODE="${GRTECLYN_EVOLVING_GEODESIC_MODE:-search}"
export GRTECLYN_METRIC_STACK_N_SPACE="${GRTECLYN_METRIC_STACK_N_SPACE:-33}"
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-lump_activity scalar_activity phi_lump_sum Pi_lump_sum chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${GRTECLYN_PROJECTION_METHOD:-mip}"

mkdir -p "${RUNS_DIR}"

NAME_ARGS=()
if [[ -n "${QD_NAME}" ]]; then
  NAME_ARGS+=(--name "${QD_NAME}")
fi
RESUME_ARGS=()
if [[ "${QD_RESUME}" == "1" ]]; then
  RESUME_ARGS+=(--resume)
fi
TARGET_EVALS_ARGS=()
if [[ -n "${QD_TARGET_EVALS}" ]]; then
  TARGET_EVALS_ARGS+=(--target-evals "${QD_TARGET_EVALS}")
fi
SEED_ARGS=()
if [[ -n "${SEED_EVAL_DIRS}" ]]; then
  # shellcheck disable=SC2206
  SEED_ARGS+=(--seed-eval-dirs ${SEED_EVAL_DIRS})
fi

FTL_RETENTION_ARGS=(--ftl-retention)
if [[ "${QD_FTL_RETENTION}" == "0" ]]; then
  FTL_RETENTION_ARGS=(--no-ftl-retention)
fi

DRY_RUN_ARGS=()
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  DRY_RUN_ARGS=(--dry-run)
fi

OBJECTIVE_MODE="${OBJECTIVE_MODE:-ftl_first}"
SCORE_WEIGHT_ARGS=()
if [[ -n "${SCORE_WEIGHTS:-}" ]]; then
  # shellcheck disable=SC2206
  for pair in ${SCORE_WEIGHTS}; do
    SCORE_WEIGHT_ARGS+=(--score-weight "${pair}")
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

# Fail fast -- do not burn cluster GPUs on untested v13 matter/geometry code.
if [[ "${SKIP_QD_PREFLIGHT_TESTS:-0}" != "1" ]]; then
  echo "Running QD preflight pytest gate..."
  ${PYTHON_BIN} -m pytest \
    "${WRAPPER_ROOT}/tests/grtresna/test_scalar_lambda_potential.py" \
    "${WRAPPER_ROOT}/tests/grtresna/test_grtresna_shell_ansatz.py" \
    "${WRAPPER_ROOT}/tests/grtresna/test_matter_geometry_consistency.py" \
    "${WRAPPER_ROOT}/tests/search/test_ftl_retention.py" \
    "${WRAPPER_ROOT}/tests/search/test_descriptors_4d.py" \
    "${WRAPPER_ROOT}/tests/search/test_qd_4d_smoke.py" \
    "${WRAPPER_ROOT}/tests/metrics/ftl/test_ftl_peak_metrics.py" \
    "${WRAPPER_ROOT}/tests/metrics/ftl/test_ftl_peak_metrics_4d.py" \
    "${WRAPPER_ROOT}/tests/metrics/ftl/test_evolving_geodesic_search_mode.py" \
    "${WRAPPER_ROOT}/tests/metrics/score/test_ftl_4d_gate.py" \
    -q --tb=short
fi

# shellcheck disable=SC2086
exec ${PYTHON_BIN} -m grteclyn_wrapper \
  --runs-dir "${RUNS_DIR}" \
  "${NAME_ARGS[@]}" \
  --example RadialRecipe \
  --set stop_time="${STOP_TIME}" \
  --set plot_interval="${PLOT_INTERVAL}" \
  --set max_level="${GRTRESNA_EVOLUTION_MAX_LEVEL}" \
  "${DRY_RUN_ARGS[@]}" \
  --consume-plotfiles \
  --consumer-delete \
  --consumer-keep-last "${CONSUMER_KEEP_LAST}" \
  --consumer-radii ${CONSUMER_RADII} \
  --ftl-L "${FTL_L}" \
  "${SCORE_WEIGHT_ARGS[@]}" \
  qd \
  --descriptor-mode "${DESCRIPTOR_MODE}" \
  --objective-mode "${OBJECTIVE_MODE}" \
  --iterations "${QD_ITERATIONS}" \
  --keep-top-eval-dirs "${QD_KEEP_TOP_EVAL_DIRS}" \
  "${FTL_RETENTION_ARGS[@]}" \
  "${TARGET_EVALS_ARGS[@]}" \
  "${RESUME_ARGS[@]}" \
  "${SEED_ARGS[@]}" \
  --batch-size "${BATCH_SIZE}" \
  --bins "${BINS}" \
  --seed "${SEED}" \
  --gpu-ids ${GPU_IDS} \
  --grtresna \
  --grtresna-ansatz shell \
  --grtresna-shell-profile "${SHELL_PROFILE}" \
  --grtresna-lumps "${LUMPS}" \
  --grtresna-full-z \
  "${DOMAIN_ARGS[@]}" \
  --grtresna-ranks "${RANKS}" \
  --grtresna-iterations "${ITERATIONS}" \
  --grtresna-nl-exit-tolerance "${GRTRESNA_NL_EXIT_TOLERANCE}" \
  --grtresna-nl-stall-tolerance "${GRTRESNA_NL_STALL_TOLERANCE}" \
  --grtresna-gridinit-workers "${GRTRESNA_GRIDINIT_WORKERS}" \
  --grtresna-timeout "${GRTRESNA_TIMEOUT}" \
  --grtresna-max-level "${GRTRESNA_MAX_LEVEL}" \
  --grtresna-refine-threshold "${GRTRESNA_REFINE_THRESHOLD}" \
  --grtresna-regrid-radius "${GRTRESNA_REGRID_RADIUS}" \
  --grtresna-jacobian-cap "${GRTRESNA_JACOBIAN_CAP}" \
  --grtresna-max-ham-pct "${GRTRESNA_MAX_HAM_PCT}" \
  --grtresna-max-mom-pct "${GRTRESNA_MAX_MOM_PCT}" \
  --solved-ftl-f-op-floor "${SOLVED_FTL_F_OP_FLOOR}" \
  --solved-ftl-near-luminal-speed-floor "${SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR}" \
  --solved-ftl-superluminal-speed-floor "${SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR}" \
  --solved-ftl-superluminal-fraction-floor "${SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR}" \
  --solved-ftl-max-physical-coord-speed "${SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED}" \
  --solved-ftl-max-physical-f-op "${SOLVED_FTL_MAX_PHYSICAL_F_OP}" \
  --solved-ftl-rejection-speed-target "${SOLVED_FTL_REJECTION_SPEED_TARGET}" \
  --grtresna-postload-gate \
  --postload-max-ham-l2 "${POSTLOAD_MAX_HAM_L2}" \
  --postload-max-mom-l2 "${POSTLOAD_MAX_MOM_L2}" \
  "${KEEP_SOURCE_ARGS[@]}"
