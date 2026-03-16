#!/usr/bin/env bash
set -euo pipefail

# Plot all standard diagnostics for a run:
# - constraint norms
# - collapse diagnostics
# - extracted Psi4 waveform / PSD in retarded time
# - extracted Psi4 waveform / PSD in simulation time
#
# Usage:
#   ./src/scripts/plot_diagnostic.sh [RUN_DIR] [RADIUS ...]
#
# Examples:
#   ./src/scripts/plot_diagnostic.sh
#   ./src/scripts/plot_diagnostic.sh "/home/jovyan/nachevsky/test/simulation/data_supported"
#   ./src/scripts/plot_diagnostic.sh "/home/jovyan/nachevsky/test/simulation/data_supported" 10 14
#
# If no radii are given, plot_extracted_psi4.py will plot all radii found in
# psi4_mode_l2m0.dat.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIS_DIR="$(cd "${SCRIPT_DIR}/../visualisation" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_RUN_DIR="$(cd "${REPO_ROOT}/.." && pwd)/data_2gpu"

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

RADII_ARGS=()
if [[ $# -gt 0 ]]; then
  RADII_ARGS=(--radii "$@")
fi

echo "[1/3] Plotting constraint norms..."
python -m src.visualisation.constraines "${CONSTRAINT_FILE}"

echo "[2/3] Plotting collapse diagnostics..."
python -m src.visualisation.diagnostic.diagnostic \
  "${COLLAPSE_FILE}" \
  --out "${VIS_DIR}/diagnostic"

echo "[3/4] Plotting extracted Psi4 in retarded time..."
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  "${PSI4_FILE}" \
  "${RADII_ARGS[@]}" \
  --time-axis retarded \
  --out "${VIS_DIR}/process_wave" \
  --plot-psd

echo "[4/4] Plotting extracted Psi4 in simulation time..."
python -m src.visualisation.process_wave.plot_extracted_psi4 \
  "${PSI4_FILE}" \
  "${RADII_ARGS[@]}" \
  --time-axis simulation \
  --out "${VIS_DIR}/process_wave" \
  --name "psi4_extracted_simulation.png" \
  --plot-psd

echo "Done."