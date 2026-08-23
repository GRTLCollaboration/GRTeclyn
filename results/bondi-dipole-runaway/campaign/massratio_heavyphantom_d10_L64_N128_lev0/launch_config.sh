#!/usr/bin/env bash
# Phase 3b, cell 1 of 3 -- the reversed mass ordering.  The one worth running.
#
# Every uneven cell so far has the gap CLOSING, so a sceptic can argue the
# closing is the artefact and the drift is a consequence of it.  This cell
# reverses which star is heavier, and predicts the opposite: the gap OPENS.
# No artefact story predicts a sign flip that tracks the mass ordering.
#
# The ordering is reversed by LIGHTENING THE CANONICAL, not by fattening the
# phantom.  Reaching |M-| > M+ from the phantom side would need |M-| ~ 0.0187,
# which is off the end of the scanned branch; lightening the canonical gets the
# same ordering inside frequencies the stability survey already covers.
#
# Frequencies picked off both branches (CPU scan, 2026-08-23):
#   canonical omega 0.75   -> M = +0.014350   (the matched cell)
#   canonical omega 0.81   -> M = +0.010721   <-- this cell
#   phantom   omega 0.7603 -> M = -0.014350   <-- this cell, unchanged
# so |M-| / M+ = 1.338: the phantom is 34% heavier than its partner.
#
# GATE: the SIGN of d(separation)/dt, not its size.  The gap must open.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=0 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.81 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/massratio_heavyphantom_d10_L64_N128_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
