#!/bin/bash
# High-quality validation of the non-spherical FTL winner (eval_000188), using
# the SAME streaming approach as the search runs: frames are extracted ON THE
# FLY by the plot consumer and the heavy plotfiles are deleted immediately
# (keep-last-2), so disk never fills up with multi-GB plotfiles.
#
# Quality (cf. Examples/SupportedWormholeCollapse/params_2gpu.txt):
#   * High base resolution N=192 (dx ~ 0.33) -- this is where the quality comes
#     from; the channel is well resolved even at the base level.
#   * Modest refinement max_level=2 with regrid_threshold=0.02.  (We deliberately
#     do NOT use max_level=5 here: the on-the-fly frames are yt SlicePlots, and
#     deep refinement makes them block-littered.  High base N gives the quality
#     without that artifact.)
#   * max_box_size=32 / min_box_size=16, 4th-order derivatives, nan_check.
#   * t=16, ~50 frames.
#
# Output kept: runs/<name>/{frames, small_data, score.json, params.txt, run.log}.
# Plotfiles are streamed through the consumer and removed (keep-last-2).
set -u
cd /home/jovyan/nachevsky/test/simulation/GRTeclyn/grteclyn-wrapper || exit 1

GPU="${1:-0}"
NAME="${2:-val16hq_nonsph_eval188}"
RUNS=/home/jovyan/nachevsky/test/simulation/GRTeclyn/runs

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

ANGULAR="--set recipe_num_lapse_angular_modes=2 \
 --set recipe_lapse_mode_ell_0=1 --set recipe_lapse_mode_rc_0=2.5 --set recipe_lapse_mode_rw_0=2.0 \
 --set recipe_lapse_mode_ell_1=2 --set recipe_lapse_mode_rc_1=2.5 --set recipe_lapse_mode_rw_1=2.0 \
 --set recipe_num_beta_angular_modes=2 \
 --set recipe_beta_mode_ell_0=1 --set recipe_beta_mode_rc_0=2.5 --set recipe_beta_mode_rw_0=2.0 \
 --set recipe_beta_mode_ell_1=2 --set recipe_beta_mode_rc_1=2.5 --set recipe_beta_mode_rw_1=2.0"

# High base resolution; modest refinement so streamed SlicePlot frames stay clean.
GRID="--set N_full=192 --set max_level=2 --set regrid_threshold=0.02 \
 --set max_box_size=32 --set min_box_size=16 \
 --set stop_time=16 --set plot_interval=48 \
 --set checkpoint_interval=-1 --set dt_multiplier=0.02 \
 --set max_spatial_derivative_order=4 --set nan_check=1"

echo "[launch] HQ streaming validation $NAME on GPU $GPU (N=192, frames on the fly, plotfiles deleted keep-last-2)"
uv run python -m grteclyn_wrapper \
  --example RadialRecipe --constrained --ftl-L 8.0 --runs-dir "$RUNS" \
  --consume-plotfiles --consumer-delete --consumer-keep-last 2 --consumer-radii 4.0 8.0 \
  $WINNER $ANGULAR $GRID \
  --cuda-devices "$GPU" --name "$NAME" \
  reproduce > "$RUNS/${NAME}.log" 2>&1
echo "[done] $NAME (exit $?)"
du -sh "$RUNS/$NAME"
