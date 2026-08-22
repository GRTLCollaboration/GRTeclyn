#!/usr/bin/env bash
# CORRECTED INITIAL DATA -- aligned solve grid.  Launched 2026-08-21.
#
# WHAT IS BEING FIXED.  The solve ran on its own grid (box 128, 128 cells ->
# spacing 1.0) and was refined three times to spacing 0.125, then copied onto
# the evolution grid (spacing 0.5).  That copy is piecewise-constant and, where
# the solve cell is finer than the evolution cell, the LAST source cell to
# touch a target cell wins -- always the one above it.  The metric therefore
# arrives displaced downward by a fraction of a cell, identically on all three
# axes.  The matter does not move with it: it is repainted analytically on the
# target grid, exact to machine precision.  Every star was born sitting above
# the centre of its own gravitational well -- canonical falls back down it,
# phantom is pushed further up it, which is the diagonal drift and 90% of the
# gap closing.  Measured offset: -0.10 at N=128, -0.055 at N=192; it does not
# converge away.
#
# THE FIX.  Solve cells = 128 * (128/64) = 256 with no solve refinement, so the
# solve spacing equals the evolution spacing (0.5) and the transfer is a
# straight copy with zero bias.  Ranks raised 8 -> 32 because the solve grid is
# 8x larger; timeout raised to match.  EVERYTHING ELSE IS BYTE-IDENTICAL to the
# cell it is paired against, so the comparison isolates the initial data.
#
# PASS CRITERION.  Read initial_data.gridinit as soon as the solve lands: the
# centroid of chi minus the centroid of the matter, in y, must be ~0 (it was
# -0.10).  Then the drift in sector_dynamics.dat must stay flat.
set -euo pipefail
REPO="$GRTECLYN_ROOT"
cd "${REPO}"
# PAIRED AGAINST: pm_eqm_sep16 (same pair, spacing 16).  That cell was the
# worst hit by the artefact -- 30% of its measured gap closing was the bug --
# so it is the most valuable second point for the separation scan.  With d=10
# (fix_pm_eqm) already running, these two give the widest lever arm on the
# inverse-square law: the predicted common acceleration falls from 1.435e-4 at
# d=10 to 5.6e-5 here, a factor 2.56 that is far outside any plausible error.
BONDI_GPU=3 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=16 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/runaway/fix_pm_eqm_sep16" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
