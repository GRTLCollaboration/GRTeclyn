#!/bin/bash
# Closed-loop CMA-ES FTL metric search (campaign: ftl_search_cmaes_01).
#
# Runs the custom CMA-ES driver over the 9-D RadialRecipe coefficient space
# (chi/alpha/beta basis), 8-wide across all H100s.  --no-phantom forces a
# NORMAL-matter (rho>=0) constrained solve: the goal is FTL *without* exotic
# matter, so geometries that demand negative energy fail the Hamiltonian
# constraint (pruned by preflight) or are penalized by the graded exotic
# penalty.  Scored by the gated objective S (non-triviality gate + log-amplified
# operational FTL - graded exotic-matter penalty).
#
# --consumer-keep-last 5 retains the final plotfiles per episode so the
# evolved-frame operational FTL (F_op^ev) and effective energy conditions
# (T^eff = G/8pi) are computed at run time and actually guide the search,
# rather than only the t=0 reconstructed slice.
#
# Output: runs/<name>/{trajectory.jsonl, metadata.json, result.json, eval_*/}
set -u
cd "$(dirname "$0")/../.." || exit 1

NAME="${1:-ftl_search_cmaes_01}"
MAXGEN="${2:-25}"
SEED="${3:-7}"

uv run python -m grteclyn_wrapper \
  --example RadialRecipe \
  --consume-plotfiles --consumer-radii 4.0 8.0 --consumer-keep-last 5 \
  --no-phantom \
  --ftl-L 8.0 \
  --name "$NAME" \
  --runs-dir ../runs \
  optimize \
    --max-generations "$MAXGEN" \
    --gpu-ids 0 1 2 3 4 5 6 7 \
    --seed "$SEED"
