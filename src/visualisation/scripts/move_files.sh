#!/bin/bash

# This script is used to move the simulation files to the new folder so the simulation results are 
# saved for future validation and verification.

# Check if a folder name argument is provided
if [ -z "$1" ]; then
    echo "Error: Please provide a folder name for the simulation results."
    echo "Usage: $0 <folder_name>"
    exit 1
fi

FOLDER_NAME="$1"

# Base directories
BASE_VISUALISATION_DIR="/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation"
TARGET_DIR="/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/SimResults/$FOLDER_NAME"
DATA_DIR="/home/jovyan/nachevsky/test/simulation/data_2gpu"
PARAMS_FILE="/home/jovyan/nachevsky/test/simulation/GRTeclyn/Examples/WormholeCollapse/params_2gpu.txt"

echo "=========================================="
echo "Creating target directory: $TARGET_DIR"
echo "=========================================="
mkdir -p "$TARGET_DIR"

echo "Copying frames (chi_z, K_z, Weyl4_Mag_z, Weyl4_Re_z)..."
# Using cp -ur to copy directories and update only if source is newer
cp -ur "$BASE_VISUALISATION_DIR/visualize/chi_z" "$TARGET_DIR/"
cp -ur "$BASE_VISUALISATION_DIR/visualize/K_z" "$TARGET_DIR/"
cp -ur "$BASE_VISUALISATION_DIR/visualize/Weyl4_Mag_z" "$TARGET_DIR/"
cp -ur "$BASE_VISUALISATION_DIR/visualize/Weyl4_Re_z" "$TARGET_DIR/"

echo "Copying plot files..."
# Using cp -u to copy files and update only if source is newer
cp -u "$BASE_VISUALISATION_DIR/constraines/constraints_plot.png" "$TARGET_DIR/"
cp -u "$BASE_VISUALISATION_DIR/diagnostic/collapse_diagnostics_plot.png" "$TARGET_DIR/"
cp -u "$BASE_VISUALISATION_DIR/process_wave/psi4_extracted_R8_R12_R16.png" "$TARGET_DIR/"
cp -u "$BASE_VISUALISATION_DIR/process_wave/psi4_extracted_simulation.png" "$TARGET_DIR/"

echo "Copying data files from $DATA_DIR..."
cp -u "$DATA_DIR/data/collapse_diagnostics.dat" "$TARGET_DIR/"
cp -u "$DATA_DIR/data/constraint_norms.dat" "$TARGET_DIR/"
cp -u "$DATA_DIR/small_data/consume_state.json" "$TARGET_DIR/"
cp -u "$DATA_DIR/small_data/psi4_mode_l2m0.dat" "$TARGET_DIR/"

echo "Copying parameter file..."
cp -u "$PARAMS_FILE" "$TARGET_DIR/"

echo "=========================================="
echo "All files successfully copied to: $TARGET_DIR"
echo "=========================================="
