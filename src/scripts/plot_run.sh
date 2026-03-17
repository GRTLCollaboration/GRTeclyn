#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SIM_ROOT="$(cd "${REPO_DIR}/.." && pwd)"

DEFAULT_DATA_DIR="${SIM_ROOT}/data_2gpu"
if [[ ! -d "${DEFAULT_DATA_DIR}" && -d "${SIM_ROOT}/data" ]]; then
  DEFAULT_DATA_DIR="${SIM_ROOT}/data"
fi

DATA_DIR="${1:-${DEFAULT_DATA_DIR}}"
RUN_DATA_DIR="${DATA_DIR}/data"
SMALL_DATA_DIR="${DATA_DIR}/small_data"
VIS_ROOT="${REPO_DIR}/src/visualisation"

if [[ ! -d "${DATA_DIR}" ]]; then
  echo "Run directory does not exist: ${DATA_DIR}" >&2
  exit 1
fi

cd "${REPO_DIR}"

echo "=========================================="
echo "Removing stale plotfiles from: ${DATA_DIR}"
shopt -s nullglob
STALE_PLOTS=(
  "${DATA_DIR}"/WormholePlt*
  "${DATA_DIR}"/SupportedWormholePlt*
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
      "${VIS_ROOT}/process_wave/psi4_extracted_simulation.pdf"
echo "Removed shared generated plot images"

echo "=========================================="
echo "Watching plotfiles in: ${DATA_DIR}"
echo "Writing small-data to: ${DATA_DIR}/small_data"
echo "=========================================="

python -m src.visualisation.process_wave.consume_plotfiles \
  --data "${DATA_DIR}" \
  --out "${DATA_DIR}/small_data" \
  --radii 8 12 16 \
  --n-points 32 \
  --frames-fields chi K Weyl4_Re Weyl4_Mag \
  --frames-axis z \
  --frames-corner \
  --frames-out "${REPO_DIR}/src/visualisation/visualize" \
  --watch --delete --keep-last 2 \
  --verbose
