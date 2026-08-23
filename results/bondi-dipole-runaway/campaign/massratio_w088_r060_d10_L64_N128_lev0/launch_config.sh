#!/usr/bin/env bash
# Phase 3b, cell 2 of 3 -- third point on the mass-scaling law.
#
# Phase 3 leaves the scaling claim resting on two points, ratio 1.000 and
# 0.7995, and two points fit any monotone law through them.  This cell puts a
# third point at ratio 0.60 -- twice the lever arm of phase 3, which clears the
# ~3% systematics in the separation-corrected fit.
#
# Frequencies picked off the phantom branch (CPU scan, 2026-08-23):
#   phantom omega 0.7603 -> M = -0.014350   (matched; ratio 1.000)
#   phantom omega 0.8040 -> M = -0.011472   (phase 3;  ratio 0.7995)
#   phantom omega 0.8800 -> M = -0.008573   <-- this cell, ratio 0.5974
# The canonical partner is unchanged at omega 0.75, M = +0.014350.
#
# GATE: the canonical star's pull ratio from the separation-corrected fit lands
# within a few % of 0.5974, and the three ratios together are LINEAR in the
# mass ratio rather than merely monotone.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=1 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.88 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/massratio_w088_r060_d10_L64_N128_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
