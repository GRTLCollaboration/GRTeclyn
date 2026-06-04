#!/bin/bash
# Closed-loop CMA-ES FTL metric search -- NON-SPHERICAL stage
# (campaign: ftl_search_nonspherical_01).
#
# The spherically-symmetric radial recipe can only produce breathing shells and
# radial channels.  --nonspherical opens up DIRECTIONAL geometries by adding
# axisymmetric Legendre angular modes (ell=1 dipole, ell=2 quadrupole) on the
# LAPSE and SHIFT only.  Those fields are pure GAUGE -- they do not enter the
# t=0 Hamiltonian/momentum constraints -- so the angular freedom:
#   * reshapes/tilts the local light cones (the Alcubierre mechanism), giving a
#     genuinely directional operational-FTL channel the search can sculpt, and
#   * adds ZERO exotic matter and does NOT break the 1D constrained (phi-from-
#     chi) solve.
# (Angular modes on chi/K are intentionally excluded: they would violate the
#  spherical 1D constraint solve and surface as forced exotic content.)
#
# --no-phantom keeps the constrained solve on NORMAL matter (rho>=0): the goal
# is FTL with the LEAST exotic matter.  The graded exotic penalty now keeps a
# usable gradient through the int_neg ~ 0.5-2 band, so the search minimises
# exotic content among the geometries that actually sustain a channel.
#
# --consumer-keep-last 5 retains the final plotfiles so the EVOLVED operational
# FTL (F_op^ev) and effective energy conditions guide the search (persistence),
# not just the t=0 reconstructed slice.
#
# Output: runs/<name>/{trajectory.jsonl, metadata.json, result.json, eval_*/}
set -u
cd "$(dirname "$0")/../.." || exit 1

NAME="${1:-ftl_search_nonspherical_01}"
MAXGEN="${2:-25}"
SEED="${3:-11}"

uv run python -m grteclyn_wrapper \
  --example RadialRecipe \
  --consume-plotfiles --consumer-radii 4.0 8.0 --consumer-keep-last 5 \
  --no-phantom \
  --ftl-L 8.0 \
  --name "$NAME" \
  --runs-dir ../runs \
  optimize \
    --nonspherical \
    --max-generations "$MAXGEN" \
    --gpu-ids 0 1 2 3 4 5 6 7 \
    --seed "$SEED"
