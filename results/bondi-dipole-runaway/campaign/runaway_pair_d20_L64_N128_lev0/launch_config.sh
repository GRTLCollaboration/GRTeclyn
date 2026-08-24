#!/usr/bin/env bash
# Fifth point of the separation series: d = 20.
#
# WHY.  The force-law excess a*d^2 - GM tightens from +4.4% at d=8 to ~1% at
# d = 12-16 and then sits there.  One wider point decides whether that ~1%
# floor is physics (residual envelope overlap, which must keep shrinking) or a
# resolution/systematic floor (which must stay put).  Star centres land at
# x = 22 and 42 -- 14 (>3 rms radii) clear of the sponge's inner edge, tighter
# than d=16's clearance but still comfortable; expected drift ~ +0.7 over
# t = 200, 400x the single-star noise floor.
#
# Identical to the packed d=8..16 series cells otherwise (same gate, same solve
# grid, same evolution grid), so it drops straight onto the same fit.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=0 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=20 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/runaway_pair_d20_L64_N128_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
