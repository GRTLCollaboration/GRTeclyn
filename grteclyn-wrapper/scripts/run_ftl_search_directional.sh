#!/bin/bash
# Closed-loop CMA-ES FTL search -- DIRECTIONAL stage
# (campaign: ftl_search_directional_01).
#
# Builds on the non-spherical stage but removes the last constraint blocking a
# genuine fore/aft asymmetry:
#
#   * FULL z-domain, no reflection.  The earlier runs used a z=0 REFLECTIVE_BC,
#     which forces the geometry to be EVEN about z=0 -- so a true ell=1 dipole
#     (odd in z) could never survive.  Here lo_boundary z is switched to
#     Sommerfeld (1 1 1) and the center moved to the box middle (z = L/2), so
#     the full z-axis is simulated and dipole (directional) channels are
#     allowed.  This is SAFE: alpha/beta are gauge and the 1D constrained
#     (phi-from-chi) solve is spherically radial regardless of the symmetry.
#
#   * 21-D angular search.  --nonspherical now also searches the radial CENTER
#     (rc) and WIDTH (rw) of every Legendre mode on the lapse and shift, not
#     just the amplitudes -- the optimizer decides where along the radius each
#     directional lobe sits and how sharp it is.
#
# --no-phantom keeps the NORMAL-matter (rho>=0) constrained solve; the graded
# exotic penalty rewards FTL with the least exotic matter.  --consumer-keep-last
# retains final plotfiles so evolved F_op^ev (persistence) guides the search.
#
# NOTE: the full z-domain roughly doubles cost vs the reflective runs.
set -u
cd "$(dirname "$0")/.." || exit 1

NAME="${1:-ftl_search_directional_01}"
MAXGEN="${2:-20}"
SEED="${3:-23}"
POPSIZE="${4:-16}"   # 21-D space wants a larger population than 8

uv run python -m grteclyn_wrapper \
  --example RadialRecipe \
  --consume-plotfiles --consumer-radii 4.0 8.0 --consumer-keep-last 5 \
  --no-phantom \
  --ftl-L 8.0 \
  --set "lo_boundary=1 1 1" \
  --set "center=32.0 32.0 32.0" \
  --name "$NAME" \
  --runs-dir ../runs \
  optimize \
    --nonspherical \
    --max-generations "$MAXGEN" \
    --population-size "$POPSIZE" \
    --gpu-ids 0 1 2 3 4 5 6 7 \
    --seed "$SEED"
