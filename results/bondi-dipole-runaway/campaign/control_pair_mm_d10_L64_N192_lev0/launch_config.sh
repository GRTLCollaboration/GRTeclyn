#!/usr/bin/env bash
# Phase 5 garnish -- MM null at the N=192 rung, completing the null-ladder
# figure alongside the PP rungs.
#
# Same-sign pair, so the dx^4 tolerance ladder does NOT apply: the solve pins
# the conformal factor to exactly flat at its outer boundary, which is only
# true at zero net ADM mass.  A same-sign pair carries 2M, so the residual
# floors near 5.4e-04 % regardless of cell size and the gate is held flat at
# the N=128 value (0.002) on every rung.
#
# KNOWN LIMIT, measured 2026-08-23 on the PP rungs and expected here.  The same
# wrong boundary condition emits a wave at t = 0 that reaches the box centre one
# light-crossing time later -- t ~ 32 for a box of half-width 32 -- and grows the
# scalar content of the box about sevenfold before draining out.  The stars
# survive it (peak field moves ~5 %, constraints stay bounded near 2e-05) but the
# diagnostics do not.  The growth factor was 7.1/7.0/6.9 at N=128/192/256 on the
# PP side: flat in resolution, which is what identifies it as a boundary effect
# rather than a grid effect.  THIS CELL IS QUOTABLE ONLY FOR t < 32.
# It is run to complete the figure, not to extend the null in time.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=2 BONDI_STOP_TIME=200 BONDI_NFULL=192 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=120 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=1 BONDI_S1=1 BONDI_S0_OMEGA=0.7603 \
BONDI_GRTRESNA_N=384 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/control_pair_mm_d10_L64_N192_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
