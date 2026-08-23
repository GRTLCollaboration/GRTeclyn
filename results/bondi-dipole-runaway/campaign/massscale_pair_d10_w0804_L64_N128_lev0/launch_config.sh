#!/usr/bin/env bash
# Phase 3 of the paper campaign -- gravity scales with the source.
#
# The phantom partner is made 20% lighter and NOTHING else changes.  Each star
# is pulled by the OTHER star's mass, so the canonical star's acceleration must
# fall by exactly the mass ratio while the phantom's own acceleration must not
# move at all.  Total momentum is now expected to be non-zero, at a rate that
# is itself a prediction.
#
# The retuned frequency was picked off the phantom branch (CPU, no GPU):
#   omega 0.7603 -> M = -0.014295   (matched, the archived d=10 cell)
#   omega 0.8040 -> M = -0.011472   (79.95% of it)  <-- this cell
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=3 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.8040 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/massscale_pair_d10_w0804_L64_N128_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
