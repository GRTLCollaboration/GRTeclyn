#!/usr/bin/env bash
# Wide-pair companion to the runaway campaign: identical to the mass-matched
# cell pm_eqm in every respect EXCEPT the separation, 12 instead of 10.
#
# WHY: at sep=10 the two lumps overlap heavily (r99 = 8.99 each), and the pair
# is not accelerating in lockstep -- the phantom runs ~3x the canonical, so the
# gap closes 9.999 -> 9.928 by t=31.  Ideal Bondi holds the separation fixed.
# If that closing is a finite-size overlap effect it must weaken at sep=12; if
# it is unchanged, the cause is something else.  This cell is the discriminator.
#
# Newtonian prediction at this separation: a = G M / d^2 = 0.014350/144
#   = 9.965e-05, so midpoint displacement = 1.993 at t=200 (vs 2.870 at sep=10).
set -euo pipefail
REPO="$GRTECLYN_ROOT"
cd "${REPO}"
BONDI_GPU=3 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=12 \
BONDI_GRTRESNA_RANKS=8 BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/runaway/pm_eqm_sep12" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
