#!/usr/bin/env bash
# Long-baseline rerun of the headline matched pair: identical to pm_eqm in
# every respect EXCEPT the stop time, 400 instead of 200.
#
# WHY: pm_eqm tracked the Newtonian prediction to 4.2% over t=200 with the
# midpoint displacement reaching 0.73.  Doubling the baseline quadruples the
# predicted displacement to 11.48, so the runaway signal grows far faster than
# the tracker's ~0.05 noise floor and than the birth-drift contamination
# (which is linear in t).  It is the cheapest way to sharpen the measurement
# without touching resolution.
#
# Newtonian prediction: a = G M / d^2 = 0.014350/100 = 1.435e-04,
# midpoint displacement 0.5*a*t^2 = 11.48 at t=400.
#
# CAVEAT: the gap closes at a fixed 3.72e-03 per unit t (a birth defect, not
# gravity), so by t=400 the separation is ~8.5 rather than 10.  Expect the
# measured/predicted ratio to run further above 1 than the 1.042 seen at
# t=200 -- that excess is the closing, not a failure of the runaway.
set -euo pipefail
REPO="$GRTECLYN_ROOT"
cd "${REPO}"
BONDI_GPU=1 BONDI_STOP_TIME=400 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_GRTRESNA_RANKS=8 BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/runaway/pm_eqm_t400" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
