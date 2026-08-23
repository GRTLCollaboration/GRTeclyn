#!/usr/bin/env bash
# Phase 1 of the paper campaign -- the sign matrix.  two canonical stars: the gap must SHRINK (mutual attraction) while the pair barycentre stays at the lone-star noise floor.
#
# Byte-identical to runs/bondi/runaway_paper/runaway_pair_d10_L64_N128_lev0
# except for the two sector flags and their frequencies, so the comparison
# isolates the sign of the mass and nothing else.  Frames are off: this cell's
# product is numbers.
#
# THE ONE EXTRA KNOB, and why only this cell needs it.  GRTresna switches
# maximal slicing on BY ITSELF whenever any lump carries negative energy,
# because the CTTK ansatz K = sign*sqrt(24 pi G rho) is imaginary for rho < 0.
# Every phantom-bearing cell therefore starts from K = 0.  This cell is the
# only one in the matrix with no phantom star, so without asking it would be
# born on an already-collapsing slice and could not be read against the mixed
# pair at all -- which is the entire job of a null control.  Forcing it also
# picks up the rest of the matched path (psi relaxation 0.6, psi floor 0.1,
# arithmetic coefficient averaging), so the two cells differ only in the sign
# of the mass.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=0 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=0 BONDI_S0_OMEGA=0.75 \
BONDI_GRTRESNA_MAXIMAL_SLICING=1 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/control_pair_pp_d10_L64_N128_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
