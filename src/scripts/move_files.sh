#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

BASE_VISUALISATION_DIR="${REPO_ROOT}/src/visualisation"

RUN_TYPE="${1:-SupportedWormholeCollapse}"

if [ "$RUN_TYPE" == "SupportedWormholeCollapse" ]; then
    DATA_DIR="$(cd "${REPO_ROOT}/.." && pwd)/data_supported"
    PARAMS_FILE="${REPO_ROOT}/Examples/SupportedWormholeCollapse/params_2gpu.txt"
    VIS_FRAMES_DIR="${BASE_VISUALISATION_DIR}/visualize"
else
    DATA_DIR="$(cd "${REPO_ROOT}/.." && pwd)/data_2gpu"
    PARAMS_FILE="${REPO_ROOT}/Examples/WormholeCollapse/params_2gpu.txt"
    VIS_FRAMES_DIR="${BASE_VISUALISATION_DIR}/visualize"
fi

if [ ! -f "$PARAMS_FILE" ]; then
    echo "Error: Parameters file not found at $PARAMS_FILE"
    exit 1
fi

RADIUS=$(awk -F'=' '/^[ \t]*wormhole_throat_radius[ \t]*=/ {print $2}' "$PARAMS_FILE" | tr -d ' ' | tr -d '\r')
A0=$(awk -F'=' '/^[ \t]*wormhole_k_monopole_amplitude[ \t]*=/ || /^[ \t]*wormhole_phi_monopole_amplitude[ \t]*=/ {print $2}' "$PARAMS_FILE" | tr -d ' ' | tr -d '\r')
A2=$(awk -F'=' '/^[ \t]*wormhole_k_quadrupole_amplitude[ \t]*=/ || /^[ \t]*wormhole_phi_perturbation_amplitude[ \t]*=/ {print $2}' "$PARAMS_FILE" | tr -d ' ' | tr -d '\r')
SIGMA=$(awk -F'=' '/^[ \t]*wormhole_k_width[ \t]*=/ || /^[ \t]*wormhole_phi_perturbation_width[ \t]*=/ {print $2}' "$PARAMS_FILE" | tr -d ' ' | tr -d '\r')

RADIUS=${RADIUS:-unknown}
A0=${A0:-unknown}
A2=${A2:-unknown}
SIGMA=${SIGMA:-unknown}

FOLDER_NAME="Run_R${RADIUS}_A0${A0}_A2${A2}_sigma${SIGMA}"
TARGET_DIR="${REPO_ROOT}/src/SimResults/$FOLDER_NAME"

echo "=========================================="
echo "Extracted Parameters:"
echo "Radius (R) = $RADIUS"
echo "A0 (Monopole) = $A0"
echo "A2 (Quadrupole) = $A2"
echo "Sigma (Width) = $SIGMA"
echo ""
echo "Creating target directory: $TARGET_DIR"
echo "=========================================="
mkdir -p "$TARGET_DIR"

echo "Copying frames (chi, K, Weyl4_Mag, Weyl4_Re, embedding)..."
cp -ur "$VIS_FRAMES_DIR/chi"* "$TARGET_DIR/" 2>/dev/null || true
cp -ur "$VIS_FRAMES_DIR/K"* "$TARGET_DIR/" 2>/dev/null || true
cp -ur "$VIS_FRAMES_DIR/Weyl4_Mag"* "$TARGET_DIR/" 2>/dev/null || true
cp -ur "$VIS_FRAMES_DIR/Weyl4_Re"* "$TARGET_DIR/" 2>/dev/null || true
cp -ur "$VIS_FRAMES_DIR/embedding"* "$TARGET_DIR/" 2>/dev/null || true

PLOTS_DIR="${BASE_VISUALISATION_DIR}/plots"
echo "Copying plot files from ${PLOTS_DIR} ..."
if [ -d "$PLOTS_DIR" ]; then
  cp -u "$PLOTS_DIR"/* "$TARGET_DIR/" 2>/dev/null || true
else
  echo "  Warning: plots directory not found at $PLOTS_DIR — run plot_diagnostic.sh first"
fi

echo "Copying data files from $DATA_DIR..."
cp -u "$DATA_DIR/data/collapse_diagnostics.dat" "$TARGET_DIR/" 2>/dev/null || true
cp -u "$DATA_DIR/data/constraint_norms.dat" "$TARGET_DIR/" 2>/dev/null || true
cp -u "$DATA_DIR/small_data/consume_state.json" "$TARGET_DIR/" 2>/dev/null || true
cp -u "$DATA_DIR/small_data/psi4_mode_l2m0.dat" "$TARGET_DIR/" 2>/dev/null || true
cp -u "$DATA_DIR/small_data/areal_radius.dat" "$TARGET_DIR/" 2>/dev/null || true

echo "Copying parameter file..."
cp -u "$PARAMS_FILE" "$TARGET_DIR/" 2>/dev/null || true

echo "=========================================="
echo "All files successfully copied to: $TARGET_DIR"
echo "=========================================="
