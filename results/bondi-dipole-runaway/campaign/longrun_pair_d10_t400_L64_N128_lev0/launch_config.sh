#!/usr/bin/env bash
# Phase 5, the sustained-acceleration money plot: the archived headline cell,
# byte-identical, run to t=400 instead of 200.  Exists for the picture -- the
# velocity still growing linearly long past the fit window -- so frames stay ON
# and every slice is cached for the fixed-colour-scale re-render (rule 6).
# Quantitative fits still come from t <= 200 only.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=1 GRTECLYN_FRAMES_CACHE_SLICES=1 \
BONDI_GPU=3 BONDI_STOP_TIME=400 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/longrun_pair_d10_t400_L64_N128_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
