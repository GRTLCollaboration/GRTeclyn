#!/usr/bin/env bash
# Phase 2, PP null at the middle rung.  The barycentre residual of the
# two-canonical pair must SHRINK as the grid refines -- that is what proves the
# mixed pair's motion is not grid noise.  Carries the maximal-slicing flag
# because it is all-canonical (README rule 9: nothing else forces the flat
# K = 0 start every phantom-bearing cell gets automatically).
# Same-sign pair: the dx^4 tolerance ladder does NOT apply here.  The solve pins
# the conformal factor to exactly flat at its outer boundary, which is only true
# for a zero-net-mass configuration; a same-sign pair carries 2M, so the residual
# floors near 5.4e-04 % regardless of cell size.  The gate is therefore held flat
# at the N=128 value (0.002) across every rung -- see GPU_RUN_PAPER.md, Phase 2.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=3 BONDI_STOP_TIME=200 BONDI_NFULL=192 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=120 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=0 BONDI_S0_OMEGA=0.75 \
BONDI_GRTRESNA_MAXIMAL_SLICING=1 \
BONDI_GRTRESNA_N=384 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/control_pair_pp_d10_L64_N192_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
