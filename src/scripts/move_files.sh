#!/bin/bash

# This script is used to move the simulation files to the new folder so the simulation results are 
# saved for future validation and verification.

# Base directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

BASE_VISUALISATION_DIR="${REPO_ROOT}/src/visualisation"

# Read run_type from argument or default to WormholeCollapse
RUN_TYPE="${1:-WormholeCollapse}"

# Automatically select the correct input directories based on the run type
if [ "$RUN_TYPE" == "SupportedWormholeCollapse" ]; then
    DATA_DIR="$(cd "${REPO_ROOT}/.." && pwd)/data_supported"
    PARAMS_FILE="${REPO_ROOT}/Examples/SupportedWormholeCollapse/params_2gpu.txt"
    VIS_FRAMES_DIR="${BASE_VISUALISATION_DIR}/visualize_supported"
else
    DATA_DIR="$(cd "${REPO_ROOT}/.." && pwd)/data_2gpu"
    PARAMS_FILE="${REPO_ROOT}/Examples/WormholeCollapse/params_2gpu.txt"
    VIS_FRAMES_DIR="${BASE_VISUALISATION_DIR}/visualize"
fi

if [ ! -f "$PARAMS_FILE" ]; then
    echo "Error: Parameters file not found at $PARAMS_FILE"
    exit 1
fi

# Extract parameters using awk
RADIUS=$(awk -F'=' '/^[ \t]*wormhole_throat_radius[ \t]*=/ {print $2}' "$PARAMS_FILE" | tr -d ' ' | tr -d '\r')
A0=$(awk -F'=' '/^[ \t]*wormhole_k_monopole_amplitude[ \t]*=/ {print $2}' "$PARAMS_FILE" | tr -d ' ' | tr -d '\r')
A2=$(awk -F'=' '/^[ \t]*wormhole_k_quadrupole_amplitude[ \t]*=/ {print $2}' "$PARAMS_FILE" | tr -d ' ' | tr -d '\r')
SIGMA=$(awk -F'=' '/^[ \t]*wormhole_k_width[ \t]*=/ {print $2}' "$PARAMS_FILE" | tr -d ' ' | tr -d '\r')

# Fallback values if missing
RADIUS=${RADIUS:-unknown}
A0=${A0:-unknown}
A2=${A2:-unknown}
SIGMA=${SIGMA:-unknown}

# Construct folder name automatically
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

echo "Copying frames (chi, K, Weyl4_Mag, Weyl4_Re)..."
cp -ur "$VIS_FRAMES_DIR/chi"* "$TARGET_DIR/" 2>/dev/null || true
cp -ur "$VIS_FRAMES_DIR/K"* "$TARGET_DIR/" 2>/dev/null || true
cp -ur "$VIS_FRAMES_DIR/Weyl4_Mag"* "$TARGET_DIR/" 2>/dev/null || true
cp -ur "$VIS_FRAMES_DIR/Weyl4_Re"* "$TARGET_DIR/" 2>/dev/null || true

echo "Copying plot files (PDF, EPS, PNG)..."
cp -u "$BASE_VISUALISATION_DIR/constraines/constraints_plot".* "$TARGET_DIR/" 2>/dev/null || true
cp -u "$BASE_VISUALISATION_DIR/diagnostic/collapse_diagnostics_plot".* "$TARGET_DIR/" 2>/dev/null || true
cp -u "$BASE_VISUALISATION_DIR/process_wave/psi4_extracted_R"*.png "$TARGET_DIR/" 2>/dev/null || true
cp -u "$BASE_VISUALISATION_DIR/process_wave/psi4_extracted_R"*.pdf "$TARGET_DIR/" 2>/dev/null || true
cp -u "$BASE_VISUALISATION_DIR/process_wave/psi4_extracted_R"*.eps "$TARGET_DIR/" 2>/dev/null || true
cp -u "$BASE_VISUALISATION_DIR/process_wave/psi4_extracted_simulation".* "$TARGET_DIR/" 2>/dev/null || true

echo "Copying data files from $DATA_DIR..."
cp -u "$DATA_DIR/data/collapse_diagnostics.dat" "$TARGET_DIR/" 2>/dev/null || true
cp -u "$DATA_DIR/data/constraint_norms.dat" "$TARGET_DIR/" 2>/dev/null || true
cp -u "$DATA_DIR/small_data/consume_state.json" "$TARGET_DIR/" 2>/dev/null || true
cp -u "$DATA_DIR/small_data/psi4_mode_l2m0.dat" "$TARGET_DIR/" 2>/dev/null || true

echo "Copying parameter file..."
cp -u "$PARAMS_FILE" "$TARGET_DIR/" 2>/dev/null || true

echo "=========================================="
echo "All files successfully copied to: $TARGET_DIR"
echo "=========================================="
