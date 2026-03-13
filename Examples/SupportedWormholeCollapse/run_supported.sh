#!/bin/bash
set -euo pipefail

PARAMS_FILE="${1:-params_2gpu.txt}"
if [[ $# -gt 0 ]]; then
    shift
fi

if [[ -x ./main3d.gnu.MPI.CUDA.ex ]]; then
    EXE=./main3d.gnu.MPI.CUDA.ex
elif [[ -x ./main3d.gnu.MPI.ex ]]; then
    EXE=./main3d.gnu.MPI.ex
else
    echo "No supported executable found. Build the example first." >&2
    exit 1
fi

mpiexec -n 1 "$EXE" "$PARAMS_FILE" "$@" > run.out

if [[ -f /home/jovyan/nachevsky/test/simulation/data_supported/data/collapse_diagnostics.dat ]]; then
    python3 ../../src/visualisation/diagnostic/diagnostic.py \
        --data /home/jovyan/nachevsky/test/simulation/data_supported \
        --out . \
        --name supported_collapse_plot.png
    echo "Done! Plot saved to supported_collapse_plot.png"
else
    echo "Run completed, but collapse diagnostics were not written."
fi
