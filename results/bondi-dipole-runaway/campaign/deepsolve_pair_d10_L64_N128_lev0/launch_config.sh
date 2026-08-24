#!/usr/bin/env bash
# The headline cell with ONE knob turned: the solve gate 10x tighter.
#
# WHY.  The resolution ladder is non-monotone (N128 low, N192 slightly above
# N256) and the paper's explanation is that two error sources fight: the t=0
# initial-data noise RISES with evolution resolution (1/dx^2) while the
# evolution's own error FALLS.  If that is right, then giving the N128 cell
# cleaner initial data -- same evolution grid, same everything, solve gate
# 0.002 -> 0.0002 -- must move its drift TOWARD the fine-grid value
# (+2.8815 -> toward +3.00).  If the drift does not move, the explanation is
# wrong and the ladder's non-monotonicity needs another story.
# Either outcome is a paper sentence; that is the point of the cell.
#
# Identical to the archived runaway_pair_d10_L64_N128_lev0 in every other
# respect (mixed pair has zero net mass, so the same-sign solve floor does NOT
# apply -- the solver converges geometrically, ~2 extra NL iterations expected:
# the original exited at 8.6e-04 % after 12).  Stall tolerance tightened by the
# same factor so the stalled door cannot fire before the new gate is reached.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=2 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.0002 BONDI_NL_STALL_TOL=0.000004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/deepsolve_pair_d10_L64_N128_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
