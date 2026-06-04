#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=env.sh
source "${SCRIPT_DIR}/env.sh"
WRAPPER_ROOT="${SCRIPT_DIR}/.."
SIM_ROOT="$(cd "${GRTECLYN_ROOT}/.." && pwd)"
VIS_ROOT="${WRAPPER_ROOT}/src/grteclyn_wrapper/visualisation"

DEFAULT_DATA_DIR=""
for candidate in "${SIM_ROOT}/data_2gpu" "${SIM_ROOT}/data_supported" "${SIM_ROOT}/data"; do
  if [[ -d "${candidate}" ]]; then
    DEFAULT_DATA_DIR="${candidate}"
    break
  fi
done
if [[ -z "${DEFAULT_DATA_DIR}" ]]; then
  DEFAULT_DATA_DIR="${SIM_ROOT}/data_2gpu"
fi

REMOVE_STALE=true
DATA_DIR=""
EXTRA_ARGS=""
JOBS=64

while [[ $# -gt 0 ]]; do
  case $1 in
    --not-remove)
      REMOVE_STALE=false
      EXTRA_ARGS="${EXTRA_ARGS} --keep-existing-frames"
      shift
      ;;
    -j|--jobs)
      if [[ -n "${2:-}" ]] && [[ "$2" != -* ]]; then
        JOBS="$2"
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    *)
      if [[ -z "$DATA_DIR" ]]; then
        DATA_DIR="$1"
      fi
      shift
      ;;
  esac
done

if [[ -z "$DATA_DIR" ]]; then
  DATA_DIR="${DEFAULT_DATA_DIR}"
fi
RUN_DATA_DIR="${DATA_DIR}/data"
SMALL_DATA_DIR="${DATA_DIR}/small_data"
FRAMES_DIR="${DATA_DIR}/frames"

if [[ ! -d "${DATA_DIR}" ]]; then
  echo "Run directory does not exist: ${DATA_DIR}" >&2
  exit 1
fi

if [ "$REMOVE_STALE" = true ]; then
  echo "=========================================="
  echo "Removing stale plotfiles from: ${DATA_DIR}"
  shopt -s nullglob
  STALE_PLOTS=(
    "${DATA_DIR}"/WormholePlt*
    "${DATA_DIR}"/SupportedWormholePlt*
    "${DATA_DIR}"/RotatingWormholePlt*
    "${DATA_DIR}"/plt*
  )
  if (( ${#STALE_PLOTS[@]} > 0 )); then
    rm -rf "${STALE_PLOTS[@]}"
    echo "Removed ${#STALE_PLOTS[@]} stale plotfile directories"
  else
    echo "No stale plotfiles found"
  fi
  shopt -u nullglob

  if [[ -d "${SMALL_DATA_DIR}" ]]; then
    rm -f "${SMALL_DATA_DIR}/psi4_mode_l2m0.dat" \
          "${SMALL_DATA_DIR}/areal_radius.dat" \
          "${SMALL_DATA_DIR}/consume_state.json"
    echo "Reset extracted small-data state in: ${SMALL_DATA_DIR}"
  fi

  if [[ -d "${RUN_DATA_DIR}" ]]; then
    rm -f "${RUN_DATA_DIR}/constraint_norms.dat" \
          "${RUN_DATA_DIR}/collapse_diagnostics.dat"
    echo "Reset run diagnostics in: ${RUN_DATA_DIR}"
  fi

  rm -f "${VIS_ROOT}/constraines/constraints_plot.png" \
        "${VIS_ROOT}/constraines/constraints_plot.eps" \
        "${VIS_ROOT}/constraines/constraints_plot.pdf" \
        "${VIS_ROOT}/diagnostic/collapse_diagnostics_plot.png" \
        "${VIS_ROOT}/diagnostic/collapse_diagnostics_plot.eps" \
        "${VIS_ROOT}/diagnostic/collapse_diagnostics_plot.pdf" \
        "${VIS_ROOT}/process_wave/psi4_extracted_R8_R12_R16.png" \
        "${VIS_ROOT}/process_wave/psi4_extracted_R8_R12_R16.eps" \
        "${VIS_ROOT}/process_wave/psi4_extracted_R8_R12_R16.pdf" \
        "${VIS_ROOT}/process_wave/psi4_extracted_simulation.png" \
        "${VIS_ROOT}/process_wave/psi4_extracted_simulation.eps" \
        "${VIS_ROOT}/process_wave/psi4_extracted_simulation.pdf" \
        "${VIS_ROOT}/process_wave/psi4_strain_analysis.png" \
        "${VIS_ROOT}/process_wave/psi4_strain_analysis.eps" \
        "${VIS_ROOT}/process_wave/psi4_strain_analysis.pdf" \
        "${VIS_ROOT}/process_wave/psi4_propagation_speed.png" \
        "${VIS_ROOT}/process_wave/psi4_propagation_speed.eps" \
        "${VIS_ROOT}/process_wave/psi4_propagation_speed.pdf"
  if [[ -d "${FRAMES_DIR}" ]]; then
    rm -rf "${FRAMES_DIR}"
    echo "Reset frames in: ${FRAMES_DIR}"
  fi
  echo "Removed shared generated plot images"
else
  echo "=========================================="
  echo "Keeping existing plotfiles and data in: ${DATA_DIR}"
fi

mkdir -p "${FRAMES_DIR}"
echo "=========================================="
echo "Watching plotfiles in: ${DATA_DIR}"
echo "Writing small-data to: ${SMALL_DATA_DIR}"
echo "Writing frames to:     ${FRAMES_DIR}"
echo "=========================================="

PYTHON=(python)
if command -v uv >/dev/null 2>&1 && [[ -f "${WRAPPER_ROOT}/pyproject.toml" ]]; then
  PYTHON=(uv run --directory "${WRAPPER_ROOT}" python)
fi

"${PYTHON[@]}" -m grteclyn_wrapper.visualisation.process_wave.consume_plotfiles \
  --data "${DATA_DIR}" \
  --out "${DATA_DIR}/small_data" \
  --radii 12 16 20 24 \
  --n-points 64 \
  --areal-radius \
  --embedding --embedding-rmax 5.0 \
  --frames-fields K Weyl4_Re \
  --frames-axis z \
  --frames-corner \
  --frames-out "${FRAMES_DIR}" \
  --watch --delete --keep-last 2 \
  --verbose -j "${JOBS}" ${EXTRA_ARGS}
