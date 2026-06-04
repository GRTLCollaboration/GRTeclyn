#!/usr/bin/env bash
set -euo pipefail
# Writes to grteclyn-wrapper/src/grteclyn_wrapper/visualisation/plots/ (directory cleared each run).
# Usage: ./grteclyn-wrapper/scripts/plot/plot_diagnostic.sh [RUN_DIR] [RADIUS ...]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"
VIS_DIR="${WRAPPER_ROOT}/src/grteclyn_wrapper/visualisation"
SIM_ROOT="$(cd "${GRTECLYN_ROOT}/.." && pwd)"

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

PLOTS_DIR="${RUN_DIR}/plots"
FRAMES_DIR="${RUN_DIR}/frames"

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

PYTHON=(python3)
if command -v uv >/dev/null 2>&1 && [[ -f "${WRAPPER_ROOT}/pyproject.toml" ]]; then
  PYTHON=(uv run --directory "${WRAPPER_ROOT}" python)
fi

echo "[1/5] Plotting constraint norms..."
"${PYTHON[@]}" -m grteclyn_wrapper.visualisation.constraines \
  "${CONSTRAINT_FILE}" \
  -o "${PLOTS_DIR}/constraints_plot.eps"

echo "[2/5] Plotting collapse diagnostics (+ areal radius + K-decay lifetime)..."
"${PYTHON[@]}" "${VIS_DIR}/diagnostic/diagnostic.py" \
  "${COLLAPSE_FILE}" \
  --data "${RUN_DIR}" \
  --out "${PLOTS_DIR}"

echo "[3/5] Plotting combined Psi4 analysis (waveforms + PSD + propagation + strain + LIGO)..."

CONFIGS=(
  "${MASS_MSUN}:${DISTANCE_MPC}"
  "30:10"
  "1000:0.002"
  "1000:1"
)

UNIQUE_CONFIGS=($(printf "%s\n" "${CONFIGS[@]}" | sort -u))

for CONFIG in "${UNIQUE_CONFIGS[@]}"; do
  M_VAL="${CONFIG%%:*}"
  D_VAL="${CONFIG##*:}"

  OUT_NAME="psi4_analysis_M${M_VAL}_D${D_VAL}.eps"
  echo "  -> Generating ${OUT_NAME} for Mass=${M_VAL} M_sun, Distance=${D_VAL} Mpc"

  "${PYTHON[@]}" -m grteclyn_wrapper.visualisation.process_wave.plot_extracted_psi4 \
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
"${PYTHON[@]}" "${VIS_DIR}/make_evolution_panel/make_evolution_panel.py" \
  --frame_dir "${FRAMES_DIR}/K_z/frames" \
  --out evolution_K_z_panel \
  --frames 20 1002 2001 3002 4500 \
  || echo "Warning: K_z evolution panel failed (missing frames under ${FRAMES_DIR}/K_z/frames?)." >&2

echo "[5/5] Plotting embedding evolution panel..."
"${PYTHON[@]}" "${VIS_DIR}/make_evolution_panel/make_evolution_panel.py" \
  --frame_dir "${FRAMES_DIR}/embedding/frames" \
  --mode embedding \
  --keep-title \
  --out evolution_embedding_4panels \
  --frames 20 1002 2001 3002 4500 \
  || echo "Warning: embedding evolution panel failed (missing frames under ${FRAMES_DIR}/embedding/frames?)." >&2

echo ""
echo "All plots saved to: ${PLOTS_DIR}"
ls -1 "${PLOTS_DIR}"
