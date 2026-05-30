#!/bin/bash
# Tier-2 extended validation of the non-spherical FTL winner
# (ftl_search_nonspherical_01 / eval_000188, score 18.30).
#
# Higher resolution (N=96, 2 AMR levels), longer evolution (t=16), dense
# plotting, and ALL plotfiles retained (no --consume-plotfiles) so we can:
#   * rebuild an x-z movie of the directional channel, and
#   * re-measure the evolved operational FTL (F_op^ev) and effective energy
#     conditions on the long-time slice to confirm the channel persists.
#
# Geometry == the exact 13-parameter winner (radial chi/alpha/beta basis +
# axisymmetric Legendre angular modes on lapse and shift).  --constrained with
# NO --phantom reproduces the search's NORMAL-matter (rho>=0) constrained solve.
# z=0 stays a reflection plane (same symmetry as the search) so this validates
# the discovered candidate rather than changing the experiment.
set -u
cd /home/jovyan/nachevsky/test/simulation/GRTeclyn/grteclyn-wrapper || exit 1

GPU="${1:-0}"
NAME="${2:-val16_nonsph_eval188}"

# --- winner geometry (eval_000188) ---
WINNER="--set recipe_chi_coeff_0=0.08289991897087551 \
 --set recipe_chi_coeff_1=0.0009085209432527861 \
 --set recipe_chi_coeff_2=-0.009425178652544067 \
 --set recipe_chi_coeff_3=0.010783822534279716 \
 --set recipe_basis_width=2.6910380503377973 \
 --set recipe_alpha_coeff_0=0.14341573215332026 \
 --set recipe_alpha_coeff_1=0.1686591284641204 \
 --set recipe_beta_coeff_0=-0.4124718353373008 \
 --set recipe_beta_coeff_1=-0.4381748383223989 \
 --set recipe_lapse_mode_amp_0=0.1828751882998318 \
 --set recipe_lapse_mode_amp_1=0.09114723946808133 \
 --set recipe_beta_mode_amp_0=-0.11870389880218904 \
 --set recipe_beta_mode_amp_1=-0.23861505704348432"

# --- angular-mode activation (fixed config used during the search) ---
ANGULAR="--set recipe_num_lapse_angular_modes=2 \
 --set recipe_lapse_mode_ell_0=1 --set recipe_lapse_mode_rc_0=2.5 --set recipe_lapse_mode_rw_0=2.0 \
 --set recipe_lapse_mode_ell_1=2 --set recipe_lapse_mode_rc_1=2.5 --set recipe_lapse_mode_rw_1=2.0 \
 --set recipe_num_beta_angular_modes=2 \
 --set recipe_beta_mode_ell_0=1 --set recipe_beta_mode_rc_0=2.5 --set recipe_beta_mode_rw_0=2.0 \
 --set recipe_beta_mode_ell_1=2 --set recipe_beta_mode_rc_1=2.5 --set recipe_beta_mode_rw_1=2.0"

# --- higher-resolution / longer-time / dense-plot validation grid ---
GRID="--set N_full=96 --set max_level=2 --set stop_time=16 \
 --set plot_interval=60 --set checkpoint_interval=-1 \
 --set dt_multiplier=0.02 --set regrid_threshold=0.01"

echo "[launch] $NAME on GPU $GPU (N=96, t=16, all plotfiles retained)"
nohup uv run python -m grteclyn_wrapper \
  --example RadialRecipe --constrained --ftl-L 8.0 --runs-dir ../runs \
  $WINNER $ANGULAR $GRID \
  --cuda-devices "$GPU" --name "$NAME" \
  reproduce > "../runs/${NAME}.log" 2>&1 &
echo "[pid] $!"
