#!/usr/bin/env bash
# Run the GRTresna BosonStarBH solver on the rotating-wormhole winding test params.
set -euo pipefail
G=/home/jovyan/nachevsky/test/simulation/GRTresna
GRTRESNA_ENV="${GRTRESNA_ENV:-$HOME/.mlspace/envs/grtresna}"
export PATH="${GRTRESNA_ENV}/bin:${PATH}"
export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
EXE="$G/Examples/BosonStarBH/Main_BosonStarBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex"
cd "$G/Examples/BosonStarBH"
rm -rf Outputs_rotwh Ham_and_Mom_errors.txt
exec mpirun --oversubscribe -np "${1:-2}" "$EXE" params_rotating_wormhole_test.txt
