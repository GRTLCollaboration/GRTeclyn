#!/usr/bin/env bash
# Separation-scan cell: mass-matched, d=16 -- the widest rung.
# Identical to the runaway baseline except the separation.
#
# WHY: extends the d = 8/10/12 scan to a point where the two lumps barely
# overlap (r99 = 8.99 each, so the tails only just touch at d=16).  Any
# finite-size or tail-mediated cause of the gap closing must be far weaker
# here; a separation-independent cause will look identical.
#
# Newtonian prediction: a = G M / d^2 = 0.014350/256 = 5.605e-05,
# midpoint displacement 1.121 at t=200.
set -euo pipefail
REPO="$GRTECLYN_ROOT"
cd "${REPO}"
BONDI_GPU=0 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=16 \
BONDI_GRTRESNA_RANKS=8 BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/runaway/pm_eqm_sep16" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
