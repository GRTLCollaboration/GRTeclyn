#!/usr/bin/env bash
# Phase 2, middle rung of the resolution ladder.  Same physics as the archived
# runaway_pair_d10_L64_N128_lev0; the ONLY thing that changes is the cell size,
# and the solve is retuned with it: solve cells = 192*(128/64) = 384 keeps the
# spacing match (rule 1), and the exit tolerance scales as dx^4 (rule 2) --
# 0.002 * (128/192)^4 = 0.000395 -- so the ladder measures convergence, not the
# tolerance floor.  Frames off: figure cells are the N=256 rungs.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=1 BONDI_STOP_TIME=200 BONDI_NFULL=192 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=120 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.000395 BONDI_NL_STALL_TOL=0.0000079 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=384 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/runaway_pair_d10_L64_N192_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
