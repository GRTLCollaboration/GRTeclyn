#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

DATA_DIR="/home/jovyan/nachevsky/test/simulation/data_2gpu"
REPO_DIR="/home/jovyan/nachevsky/test/simulation/GRTeclyn"

echo "=========================================="
echo "Cleaning up previous run data..."
echo "=========================================="
if [ -d "$DATA_DIR" ]; then
    # Delete all files and subdirectories inside the data directory
    rm -rf "${DATA_DIR:?}/"*
    echo "Deleted contents of $DATA_DIR"
else
    echo "Data directory $DATA_DIR does not exist. Skipping cleanup."
fi

# Change to the root of the repository so python -m src... works properly
cd "$REPO_DIR"

echo "=========================================="
echo "Starting plot processing..."
echo "=========================================="
python -m src.visualisation.process_wave.consume_plotfiles \
  --data "$DATA_DIR" \
  --out "$DATA_DIR/small_data" \
  --radii 8 12 16 \
  --n-points 32 \
  --frames-fields chi K Weyl4_Re Weyl4_Mag \
  --frames-axis z \
  --frames-corner \
  --frames-out "/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/visualize" \
  --watch --delete --keep-last 2 \
  --verbose
