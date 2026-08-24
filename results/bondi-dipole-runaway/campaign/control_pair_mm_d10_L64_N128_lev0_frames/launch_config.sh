#!/usr/bin/env bash
# The MM null at N=128, re-run for ONE reason: frames.
#
# WHY THIS CELL EXISTS.  For a same-sign pair the sector splitter assigns matter
# by field sign, so BOTH stars land in the same sector and the tracker reports a
# single core at the pair's midpoint (x = 32.0000) with coord_sep = nan.  Every
# same-sign number in the campaign is therefore the combined centroid: it shows
# the pair as a whole does not move, and says nothing about what the two stars
# do relative to each other.
#
# The PP cell at N=256 was launched with frames, so its two wells could be
# tracked directly in the images -- they close from 8.75 to zero by t ~ 35, a
# real gravitational infall.  The MM cells were launched with frames OFF, so the
# mirror measurement does not exist: whether two negative masses push APART, as
# they should, is currently a prediction in this campaign and not a result.
#
# This run answers that.  Same cell as control_pair_mm_d10_L64_N128_lev0 in
# every respect except GRTECLYN_FRAMES=1 and the slice cache, written to a new
# directory so the completed cell and its packed extract are untouched.
#
# The same-sign boundary limit still applies: the box floods after t ~ 32, so
# the two-well tracking is only trustworthy BEFORE that.  t=200 is kept anyway
# so this cell is a like-for-like twin of the original.
#
# Same-sign pair, so the dx^4 tolerance ladder does not apply and the solve gate
# is held flat at 0.002 -- see GPU_RUN_PAPER.md, Phase 2.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=1 GRTECLYN_FRAMES_CACHE_SLICES=1 \
BONDI_GPU=0 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=1 BONDI_S1=1 BONDI_S0_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/control_pair_mm_d10_L64_N128_lev0_frames" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
