#!/usr/bin/env bash
# Phase 2, finest rung.  Solve cells = 256*(128/64) = 512; tolerance scales as
# dx^4: 0.002 * (128/256)^4 = 0.000125.  This is a figure cell, so frames stay
# ON and every slice is cached so the movies can be re-rendered on one fixed
# colour scale afterwards (rule 6).  The 512^3 solve is the long pole (~4 h on
# 32 ranks), which is exactly why this launches while the card is still busy:
# the evolution only asks for the GPU when the solve is done.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=1 GRTECLYN_FRAMES_CACHE_SLICES=1 \
BONDI_GPU=0 BONDI_STOP_TIME=200 BONDI_NFULL=256 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=160 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.000125 BONDI_NL_STALL_TOL=0.0000025 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=512 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=43200 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/runaway_pair_d10_L64_N256_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
