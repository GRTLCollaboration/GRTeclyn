#!/usr/bin/env bash
# Resolution + long-time promotion ladder for smoke survivors.
# Skips REJECTED_LONG_RUN candidates (e.g. random_000).
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
if command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
  PYTHON_BIN="uv run python"
fi

STOP_TIME="${STOP_TIME:-5.0}"
RESOLUTIONS="${RESOLUTIONS:-64 96 128}"
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/radialrecipe_gpu_promote}"
BUILD="${BUILD:-0}"
CUDA_DEVICE_START="${CUDA_DEVICE_START:-0}"
FTL_L="${FTL_L:-8.0}"
ENABLE_FTL_SCORING="${ENABLE_FTL_SCORING:-1}"

# Comma-separated override, else default promoted survivors.
if [[ -n "${PROMOTE_CANDIDATES:-}" ]]; then
  IFS=',' read -r -a CANDIDATES <<< "${PROMOTE_CANDIDATES}"
else
  CANDIDATES=(ellis_bronnikov bubble_wall_016 alcubierre_warp)
fi

REJECTED=(random_000)

echo "== RadialRecipe GPU promotion ladder =="
echo "stop_time=${STOP_TIME}  N_full in: ${RESOLUTIONS}"
echo "candidates: ${CANDIDATES[*]}"
echo "runs_dir: ${RUNS_DIR}"
echo "ftl: ENABLE_FTL_SCORING=${ENABLE_FTL_SCORING}  FTL_L=${FTL_L}"
echo

gpu="${CUDA_DEVICE_START}"
for candidate in "${CANDIDATES[@]}"; do
  for rejected in "${REJECTED[@]}"; do
    if [[ "${candidate}" == "${rejected}" ]]; then
      echo "SKIP ${candidate}: in REJECTED_LONG_RUN (${rejected})"
      continue 2
    fi
  done

  for n_full in ${RESOLUTIONS}; do
    echo "--- ${candidate} N_full=${n_full} t=${STOP_TIME} (GPU ${gpu}) ---"
    if [[ "${candidate}" == "ellis_bronnikov" || "${candidate}" == "flat_minkowski" || "${candidate}" == "alcubierre_warp" || "${candidate}" == "schwarzschild_puncture" ]]; then
      SEED_NAME="${candidate}" CANDIDATE_ID="" NONSPHERICAL_ID="" \
        STOP_TIME="${STOP_TIME}" N_FULL="${n_full}" BUILD="${BUILD}" \
        FTL_L="${FTL_L}" ENABLE_FTL_SCORING="${ENABLE_FTL_SCORING}" \
        RUN_STAMP="${RUN_STAMP}" RUNS_DIR="${RUNS_DIR}" \
        CUDA_VISIBLE_DEVICES_OVERRIDE="${gpu}" \
        bash "${SCRIPT_DIR}/run_radialrecipe_gpu_smoke.sh"
    else
      SEED_NAME="" CANDIDATE_ID="${candidate}" NONSPHERICAL_ID="" \
        STOP_TIME="${STOP_TIME}" N_FULL="${n_full}" BUILD="${BUILD}" \
        FTL_L="${FTL_L}" ENABLE_FTL_SCORING="${ENABLE_FTL_SCORING}" \
        RUN_STAMP="${RUN_STAMP}" RUNS_DIR="${RUNS_DIR}" \
        CUDA_VISIBLE_DEVICES_OVERRIDE="${gpu}" \
        bash "${SCRIPT_DIR}/run_radialrecipe_gpu_smoke.sh"
    fi
    gpu=$((gpu + 1))
    BUILD=0
  done
done

echo
echo "== Promotion summary (upgraded FTL scoring) =="
${PYTHON_BIN} - <<PY
import json
from pathlib import Path
from grteclyn_wrapper.metrics import read_episode_metrics
from grteclyn_wrapper.score import score_episode

runs = Path(${RUNS_DIR@Q})
stamp = ${RUN_STAMP@Q}
ftl_L = float(${FTL_L@Q})
rows = []
for ep in sorted(runs.glob(f"*_gpu_t*_{stamp}")):
    try:
        m = read_episode_metrics(ep, ftl_L=ftl_L)
        s = score_episode(m, target_stop_time=float(${STOP_TIME@Q}))
    except Exception as exc:
        rows.append({"episode": ep.name, "error": str(exc)})
        continue
    rows.append({
        "episode": ep.name,
        "S": round(s.total, 2),
        "F_log": round(s.components.get("ftl_shortcut", 0), 3),
        "F_asym": round(s.components.get("expansion_asymmetry", 0), 3),
        "s_comoving": round(s.components.get("comoving_stability", 0), 3),
        "s_euler": round(s.components.get("stability", 0), 3),
        "max_H_L2": m.constraints.max_hamiltonian_l2 if m.constraints else None,
        "min_lapse": m.collapse.min_lapse if m.collapse else None,
        "min_chi": m.collapse.min_chi if m.collapse else None,
        "stationary_fallback": bool(m.comoving.stationary) if m.comoving else False,
    })
rows.sort(key=lambda r: r.get("S", -1), reverse=True)
print(json.dumps(rows, indent=2))
print()
print(f"{'episode':<52} {'S':>7} {'F_log':>6} {'s_com':>6} {'s_eul':>6} {'stat':>5}")
print("-" * 82)
for r in rows:
    if "error" in r:
        print(f"{r['episode']:<52} ERROR")
        continue
    stat = "Y" if r.get("stationary_fallback") else "N"
    print(
        f"{r['episode']:<52} {r['S']:7.2f} {r['F_log']:6.3f} "
        f"{r['s_comoving']:6.3f} {r['s_euler']:6.3f} {stat:>5}"
    )
PY
