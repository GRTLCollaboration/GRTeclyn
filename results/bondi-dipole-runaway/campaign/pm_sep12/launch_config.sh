#!/usr/bin/env bash
# Separation-scan cell: unmatched (phantom heavier), d=12
# Identical to the runaway campaign's baseline in every respect EXCEPT
# separation (12) and the lump1 frequency.
#
# WHY: the matched pair's gap closes by ~0.17 over t=57 at d=10, which ideal
# Bondi forbids.  Scanning d = 8/10/12 in both the matched and mass-mismatched
# flavours separates a finite-size overlap cause (must fall off fast with d)
# from a separation-independent one.
#
# Newtonian prediction here: a = G M / d^2 = 9.965e-05
set -euo pipefail
REPO="$GRTECLYN_ROOT"
cd "${REPO}"
BONDI_GPU=0 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=12 \
BONDI_GRTRESNA_RANKS=8 BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75  \
BONDI_RUNS_DIR="${REPO}/runs/bondi/runaway/pm_sep12" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
