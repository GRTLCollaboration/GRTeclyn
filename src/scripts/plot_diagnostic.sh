#!/usr/bin/env bash
set -euo pipefail

# Plot all standard diagnostics for a run.
# All output goes to a single folder: src/visualisation/plots/
# The folder is wiped on each invocation so results are always fresh.
#
# Plots produced:
#   constraints_plot.*              — constraint norms
#   collapse_diagnostics_plot.*   — collapse diagnostics (+ areal radius + K-decay lifetime)
#   psi4_analysis.*               — combined 3x2 panel: waveforms, PSD, propagation speed,
#                                   strain PSD, and LIGO sensitivity overlay
#   evolution_K_z_panel.*         — K_z frame strip (color, frames 2 1002 3001 4002 4050)
#   evolution_embedding_4panels.* — embedding frame strip (--keep-title, frames 2 1002 3001 4002)
#
# Notes:
# - The ESD (energy spectral density) panel can be frequency-cut to hide the flat high-f tail.
#   Default: cut at f=20 (code units M^-1). Override with env var:
#     ESD_FMAX=30 ./src/scripts/plot_diagnostic.sh ...
#   Disable passing the cutoff entirely with:
#     ESD_FMAX=none ./src/scripts/plot_diagnostic.sh ...
# - The strain vs Advanced LIGO panel is plotted in the paper-style ASD form by default
#   (sqrt(S_h) and sqrt(S_n), units 1/sqrt(Hz), band 10–5000 Hz).
#   By default we use a "detectable" scaling (M=1000 Msun, D=0.002 Mpc ~ 2 kpc)
#   so the example signal appears in-band. Control scaling via:
#     MASS_MSUN=30 DISTANCE_MPC=10 ./src/scripts/plot_diagnostic.sh ...
#   Optional:
#     LIGO_QUANTITY=hchar  # use characteristic strain instead of ASD
# - Evolution panels need PNG frames under src/visualisation/visualize/K_z/frames and
#   src/visualisation/visualize/embedding/frames (e.g. from make_movies). Missing frame
#   indices cause that step to be skipped with a warning.
#
# Usage:
#   ./src/scripts/plot_diagnostic.sh [RUN_DIR] [RADIUS ...]
#
# Examples:
#   ./src/scripts/plot_diagnostic.sh
#   ./src/scripts/plot_diagnostic.sh "/home/jovyan/nachevsky/test/simulation/data_2gpu"
#   ./src/scripts/plot_diagnostic.sh "/home/jovyan/nachevsky/test/simulation/data" 10 14
#
# If no radii are given, plot_extracted_psi4.py will plot all radii found in
# psi4_mode_l2m0.dat.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIS_DIR="$(cd "${SCRIPT_DIR}/../visualisation" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SIM_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"
PLOTS_DIR="${VIS_DIR}/plots"

ESD_FMAX_DEFAULT="20"
ESD_FMAX="${ESD_FMAX:-$ESD_FMAX_DEFAULT}"

MASS_MSUN_DEFAULT="1000"
DISTANCE_MPC_DEFAULT="0.002"
LIGO_QUANTITY_DEFAULT="asd"

MASS_MSUN="${MASS_MSUN:-$MASS_MSUN_DEFAULT}"
DISTANCE_MPC="${DISTANCE_MPC:-$DISTANCE_MPC_DEFAULT}"
LIGO_QUANTITY="${LIGO_QUANTITY:-$LIGO_QUANTITY_DEFAULT}"

candidate_run_mtime() {
  local run_dir="$1"
  local constraint_file="${run_dir}/data/constraint_norms.dat"
  local collapse_file="${run_dir}/data/collapse_diagnostics.dat"
  local psi4_file="${run_dir}/small_data/psi4_mode_l2m0.dat"

  if [[ ! -f "${constraint_file}" || ! -f "${collapse_file}" || ! -f "${psi4_file}" ]]; then
    return 1
  fi

  local t1 t2 t3 newest
  t1=$(stat -c %Y "${constraint_file}")
  t2=$(stat -c %Y "${collapse_file}")
  t3=$(stat -c %Y "${psi4_file}")
  newest="${t1}"
  if (( t2 > newest )); then newest="${t2}"; fi
  if (( t3 > newest )); then newest="${t3}"; fi
  printf '%s\n' "${newest}"
}

choose_default_run_dir() {
  local candidates=("${SIM_ROOT}/data_2gpu" "${SIM_ROOT}/data_supported" "${SIM_ROOT}/data")
  local best_dir=""
  local best_mtime=-1
  local dir mtime

  for dir in "${candidates[@]}"; do
    if mtime=$(candidate_run_mtime "${dir}"); then
      if (( mtime > best_mtime )); then
        best_mtime="${mtime}"
        best_dir="${dir}"
      fi
    fi
  done

  if [[ -n "${best_dir}" ]]; then
    printf '%s\n' "${best_dir}"
    return 0
  fi

  if [[ -d "${SIM_ROOT}/data_2gpu" ]]; then
    printf '%s\n' "${SIM_ROOT}/data_2gpu"
  else
    printf '%s\n' "${SIM_ROOT}/data"
  fi
}

DEFAULT_RUN_DIR="$(choose_default_run_dir)"

RUN_DIR="${1:-$DEFAULT_RUN_DIR}"
if [[ $# -gt 0 ]]; then
  shift
fi

CONSTRAINT_FILE="${RUN_DIR}/data/constraint_norms.dat"
COLLAPSE_FILE="${RUN_DIR}/data/collapse_diagnostics.dat"
PSI4_FILE="${RUN_DIR}/small_data/psi4_mode_l2m0.dat"

if [[ ! -f "${CONSTRAINT_FILE}" ]]; then
  echo "Missing constraint norms file: ${CONSTRAINT_FILE}" >&2
  exit 1
fi
if [[ ! -f "${COLLAPSE_FILE}" ]]; then
  echo "Missing collapse diagnostics file: ${COLLAPSE_FILE}" >&2
  exit 1
fi
if [[ ! -f "${PSI4_FILE}" ]]; then
  echo "Missing Psi4 extracted file: ${PSI4_FILE}" >&2
  exit 1
fi

echo "Using run directory: ${RUN_DIR}"
echo "All plots -> ${PLOTS_DIR}"

rm -rf "${PLOTS_DIR}"
mkdir -p "${PLOTS_DIR}"

RADII_ARGS=()
if [[ $# -gt 0 ]]; then
  RADII_ARGS=(--radii "$@")
fi

ESD_ARGS=()
if [[ -n "${ESD_FMAX}" && "${ESD_FMAX}" != "none" ]]; then
  ESD_ARGS=(--esd-fmax "${ESD_FMAX}")
fi

LIGO_ARGS=(--ligo-quantity "${LIGO_QUANTITY}")

echo "[1/5] Plotting constraint norms..."
python3 -m src.visualisation.constraines \
  "${CONSTRAINT_FILE}" \
  -o "${PLOTS_DIR}/constraints_plot.eps"

echo "[2/5] Plotting collapse diagnostics (+ areal radius + K-decay lifetime)..."
python3 "${VIS_DIR}/diagnostic/diagnostic.py" \
  "${COLLAPSE_FILE}" \
  --data "${RUN_DIR}" \
  --out "${PLOTS_DIR}"

echo "[3/5] Plotting combined Psi4 analysis (waveforms + PSD + propagation + strain + LIGO)..."

# Define multiple combinations of Mass (M_sun) and Distance (Mpc) as "MASS:DISTANCE"
# You can add or modify the values in this array to plot different configurations!
CONFIGS=(
  "${MASS_MSUN}:${DISTANCE_MPC}"   # Keeps the default/env-variable configuration
  "30:10"                          # 30 M_sun at 10 Mpc
  "1000:0.002"                     # 1000 M_sun at 0.002 Mpc
  "1000:1"                          # 100 M_sun at 1 Mpc
)

# Remove duplicates if any
UNIQUE_CONFIGS=($(printf "%s\n" "${CONFIGS[@]}" | sort -u))

for CONFIG in "${UNIQUE_CONFIGS[@]}"; do
  M_VAL="${CONFIG%%:*}"
  D_VAL="${CONFIG##*:}"
  
  # Format output filename so plots don't overwrite each other
  OUT_NAME="psi4_analysis_M${M_VAL}_D${D_VAL}.eps"
  echo "  -> Generating ${OUT_NAME} for Mass=${M_VAL} M_sun, Distance=${D_VAL} Mpc"
  
  python3 -m src.visualisation.process_wave.plot_extracted_psi4 \
    "${PSI4_FILE}" \
    "${RADII_ARGS[@]}" \
    "${ESD_ARGS[@]}" \
    "${LIGO_ARGS[@]}" \
    --out "${PLOTS_DIR}" \
    --name "${OUT_NAME}" \
    --combined \
    --strain --mass-msun "${M_VAL}" --distance-mpc "${D_VAL}"
done

echo "[4/5] Plotting K_z evolution panel (color)..."
python3 "${VIS_DIR}/make_evolution_panel/make_evolution_panel.py" \
  --frame_dir "${VIS_DIR}/visualize/K_z/frames" \
  --out evolution_K_z_panel \
  --frames 20 1002 2001 3002 4500 \
  || echo "Warning: K_z evolution panel failed (missing frames under visualize/K_z/frames?)." >&2

echo "[5/5] Plotting embedding evolution panel..."
python3 "${VIS_DIR}/make_evolution_panel/make_evolution_panel.py" \
  --frame_dir "${VIS_DIR}/visualize/embedding/frames" \
  --mode embedding \
  --keep-title \
  --out evolution_embedding_4panels \
  --frames 20 1002 2001 3002 4500 \
  || echo "Warning: embedding evolution panel failed (missing frames under visualize/embedding/frames?)." >&2

echo ""
echo "All plots saved to: ${PLOTS_DIR}"
ls -1 "${PLOTS_DIR}"