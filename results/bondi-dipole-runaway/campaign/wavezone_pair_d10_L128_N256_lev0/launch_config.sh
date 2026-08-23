#!/usr/bin/env bash
# Phase 4, the wave zone.  Doubled box at the SAME cell size (dx = 0.5, so the
# baseline tolerance already matches -- dx^4 scaling changes nothing here); the
# solve box doubles with it (256, solved at 512 = 256*(256/128), rule 1 kept).
# The sponge is pushed out to 48->60 explicitly: its radii are absolute, and at
# the defaults (24->32) the dissipation band would run through three of the
# four extraction shells.  In-code Weyl extraction ON at 16 24 32 40 -- shells
# outside the stars (matter ends ~11.5) and inside the sponge.
set -euo pipefail
REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "${REPO}"
GRTECLYN_FRAMES=0 \
BONDI_GPU=2 BONDI_STOP_TIME=200 BONDI_NFULL=256 BONDI_LFULL=128 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_SPONGE_INNER=48 BONDI_SPONGE_OUTER=60 \
BONDI_RADII="16 24 32 40" BONDI_EXTRACTION_RADII="16 24 32 40" BONDI_PSI4_HIGHER_L=1 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_DOMAIN_L=256 BONDI_GRTRESNA_N=512 BONDI_GRTRESNA_MAXLEVEL=0 BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=43200 \
BONDI_RUNS_DIR="${REPO}/runs/bondi/staging/wavezone_pair_d10_L128_N256_lev0" \
  bash "${REPO}/grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh"
