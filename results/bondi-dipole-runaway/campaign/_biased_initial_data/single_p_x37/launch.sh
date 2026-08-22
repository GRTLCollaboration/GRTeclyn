#!/usr/bin/env bash
# Lone-star drift control: lone CANONICAL star at x = 37, the canonical's own position in pm_eqm
#
# WHY THIS EXISTS.  The coordinate gap between the two lumps closes by 7-9% over
# t = 200, and the closing is NOT a fixed birth drift as first recorded: the
# rate starts near zero, peaks around t = 50-150, then turns over, earliest at
# the smallest separation.  Two explanations survive -- the pair really is
# oscillating about a constant mean separation (which is what Bondi needs), or
# each star is being pulled somewhere by the construction, independently of its
# companion.  A star ALONE in the box separates them: with no companion there
# is nothing for a pair oscillation to be about, so any drift here is
# construction, not physics.
#
# DESIGN.  Byte-identical to the pm_eqm cell in grid, rung, potential, sponge,
# solve tolerance and cadence.  The only differences are one lump instead of
# two, and the placement below.  The lone canonical and lone phantom sit at the
# SAME off-centre position, so comparing them isolates the sector from the
# position: a pull toward the box centre moves both the same way, while the
# sign-follows-energy-density pattern seen in the transverse drift moves them
# oppositely.  The centred canonical cell is the baseline -- a star at x = 32
# has nowhere to be pulled toward and must sit still.
#
# MAXIMAL SLICING IS FORCED ON.  GRTresna picks maximal slicing automatically
# for exotic matter only, so a purely canonical cell silently gets K = 1e-01
# and is born on an already-collapsing slice, while every phantom-bearing cell
# starts at K = 0.  Without this flag the lone canonical cells would not be
# comparable with pm_eqm or with the lone phantom.
set -euo pipefail
REPO="$GRTECLYN_ROOT"
cd "${REPO}"
BONDI_GPU=0 BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_GRTRESNA_RANKS=8 BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_GRTRESNA_MAXIMAL_SLICING=1 \
BONDI_EXOTIC=0 BONDI_OMEGA=0.75 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/runaway/single_p_x37" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_single_selfgrav.sh"
