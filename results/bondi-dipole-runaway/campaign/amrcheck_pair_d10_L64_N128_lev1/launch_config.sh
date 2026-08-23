#!/usr/bin/env bash
# Referee-proofing only: the headline d=10 cell with mesh refinement switched
# ON, against the whole campaign's max_level = 0 rule.
#
# The campaign runs uniform grids because the convergence ladder must mean
# exactly one thing -- the cell size -- and a tagging criterion re-evaluated per
# rung refines a different region on each.  A referee will still ask whether the
# result depends on that choice.  This cell answers it once, cheaply.
#
# PREDICTION: identical to lev0, because the tagger never fires.  Refinement
# triggers at |chi - 1| = 0.02 and the corrected runs peak near 0.005, so this
# cell should evolve on level 0 for its entire life and pay AMR bookkeeping for
# nothing.  If it DOES tag, that is the interesting outcome and the reason to
# run it -- it would mean the peak curvature is four times what the uniform runs
# report.
#
# GATE: drift and acceleration match the lev0 cell to the same tolerance the
# N=192 and N=256 rungs match each other (0.4 % / 0.9 %), and the AMR
# constraint columns show level 1 was never populated.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=2 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=1 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/amrcheck_pair_d10_L64_N128_lev1" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
