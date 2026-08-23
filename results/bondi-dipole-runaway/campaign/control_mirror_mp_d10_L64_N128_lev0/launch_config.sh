#!/usr/bin/env bash
# Phase 1 of the paper campaign -- the sign matrix.  sectors swapped: the runaway must reverse, equal and opposite to the archived d=10 cell.
#
# Byte-identical to runs/bondi/runaway_paper/runaway_pair_d10_L64_N128_lev0
# except for the two sector flags and their frequencies, so the comparison
# isolates the sign of the mass and nothing else.  Frames are off: this cell's
# product is numbers.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=2 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=1 BONDI_S1=0 BONDI_S0_OMEGA=0.7603 BONDI_S1_OMEGA=0.75 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/control_mirror_mp_d10_L64_N128_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
