#!/bin/bash
# Tier-2 extended validation of the top warp candidates from ftl_search_cmaes_01.
# N=96 (higher precision), t=16 (longer), plotfiles retained for evolved FTL + effective EC.
set -u
cd /home/jovyan/nachevsky/test/simulation/GRTeclyn/grteclyn-wrapper

COMMON="--example RadialRecipe --constrained --phantom --runs-dir ../runs --ftl-L 8.0 \
 --set recipe_exotic_matter=1 --set N_full=96 --set max_level=2 --set stop_time=16 \
 --set plot_interval=300 --set checkpoint_interval=-1 --set dt_multiplier=0.02 --set regrid_threshold=0.01"

declare -A C
C[eval_000196]="--set recipe_chi_coeff_0=-0.3796813719 --set recipe_chi_coeff_1=0.2622340700 --set recipe_chi_coeff_2=0.1133712391 --set recipe_chi_coeff_3=0.0517640228 --set recipe_basis_width=2.1408702450 --set recipe_alpha_coeff_0=-0.1076117157 --set recipe_alpha_coeff_1=0.1485128015 --set recipe_beta_coeff_0=-0.4995935625 --set recipe_beta_coeff_1=-0.4948732316"
C[eval_000166]="--set recipe_chi_coeff_0=-0.3851429853 --set recipe_chi_coeff_1=0.2615323037 --set recipe_chi_coeff_2=0.1343336176 --set recipe_chi_coeff_3=0.0417995589 --set recipe_basis_width=2.1327158550 --set recipe_alpha_coeff_0=-0.0157981685 --set recipe_alpha_coeff_1=0.1800308706 --set recipe_beta_coeff_0=-0.4710164498 --set recipe_beta_coeff_1=-0.4993035157"
C[eval_000192]="--set recipe_chi_coeff_0=-0.3877858065 --set recipe_chi_coeff_1=0.2612222534 --set recipe_chi_coeff_2=0.1294223179 --set recipe_chi_coeff_3=0.0714402973 --set recipe_basis_width=2.1462688371 --set recipe_alpha_coeff_0=-0.0700023241 --set recipe_alpha_coeff_1=0.1413657927 --set recipe_beta_coeff_0=-0.4767922030 --set recipe_beta_coeff_1=-0.4829977174"
C[eval_000146]="--set recipe_chi_coeff_0=-0.3193747641 --set recipe_chi_coeff_1=0.2157708792 --set recipe_chi_coeff_2=0.1599065626 --set recipe_chi_coeff_3=0.0318797977 --set recipe_basis_width=1.9399682590 --set recipe_alpha_coeff_0=-0.0634993691 --set recipe_alpha_coeff_1=0.1644813739 --set recipe_beta_coeff_0=-0.3791197810 --set recipe_beta_coeff_1=-0.4999565643"

gpu=0
for ep in eval_000196 eval_000166 eval_000192 eval_000146; do
  name="val16_${ep}"
  echo "[launch] $name on GPU $gpu"
  nohup uv run python -m grteclyn_wrapper \
     $COMMON ${C[$ep]} --cuda-devices "$gpu" --name "$name" \
     reproduce > "../runs/${name}.log" 2>&1 &
  gpu=$((gpu+1))
done
wait
echo "[done] all validation runs finished"
