#!/usr/bin/env bash
# GPU smoke batch for accepted non-spherical RadialRecipe candidates.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"
SMOKE="${SCRIPT_DIR}/run_radialrecipe_gpu_smoke.sh"

RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/radialrecipe_nonspherical}"
STOP_TIME="${STOP_TIME:-2.0}"
N_FULL="${N_FULL:-64}"
BUILD="${BUILD:-0}"
CUDA_START="${CUDA_START:-0}"
CONSUME_PLOTFILES="${CONSUME_PLOTFILES:-1}"
CONSUMER_DELETE="${CONSUMER_DELETE:-1}"
CONSUMER_RADII="${CONSUMER_RADII:-4 8}"
FTL_L="${FTL_L:-8.0}"
ENABLE_FTL_SCORING="${ENABLE_FTL_SCORING:-1}"

# Accepted from nonspherical_ray_validation.csv (skip rejected_negative_chi).
CANDIDATES=(
  dipole_lopsided_000
  quadrupole_bubble_001
  mixed_modes_002
  random_angular_004
  random_angular_005
  random_angular_006
  random_angular_007
)

# Optional: include bad control for comparison (preflight may still fail on GPU).
if [[ "${INCLUDE_BAD_CONTROL:-0}" == "1" ]]; then
  CANDIDATES+=(strong_quadrupole_bad_003)
fi

echo "== Non-spherical GPU batch =="
echo "candidates: ${#CANDIDATES[@]}"
echo "runs_dir:   ${RUNS_DIR}"
echo "stop_time:  ${STOP_TIME}  N_full: ${N_FULL}"
echo "consumer:   CONSUME=${CONSUME_PLOTFILES} DELETE=${CONSUMER_DELETE} RADII=${CONSUMER_RADII}"
echo "ftl:        ENABLE_FTL_SCORING=${ENABLE_FTL_SCORING}  FTL_L=${FTL_L}"
echo

mkdir -p "${RUNS_DIR}"
PIDS=()
GPU="${CUDA_START}"

for cid in "${CANDIDATES[@]}"; do
  log="${RUNS_DIR}/${cid}_${RUN_STAMP}.log"
  echo "Launch ${cid} on GPU ${GPU} -> ${log}"
  (
    NONSPHERICAL_ID="${cid}" \
    SEED_NAME="" CANDIDATE_ID="" \
    BUILD=0 STOP_TIME="${STOP_TIME}" N_FULL="${N_FULL}" \
    CONSUME_PLOTFILES="${CONSUME_PLOTFILES}" \
    CONSUMER_DELETE="${CONSUMER_DELETE}" \
    CONSUMER_RADII="${CONSUMER_RADII}" \
    FTL_L="${FTL_L}" \
    ENABLE_FTL_SCORING="${ENABLE_FTL_SCORING}" \
    RUN_STAMP="${RUN_STAMP}" RUNS_DIR="${RUNS_DIR}" \
    CUDA_VISIBLE_DEVICES_OVERRIDE="${GPU}" \
    bash "${SMOKE}" > "${log}" 2>&1
    echo "DONE ${cid} (GPU ${GPU})" >> "${log}"
  ) &
  PIDS+=($!)
  GPU=$((GPU + 1))
done

fail=0
for pid in "${PIDS[@]}"; do
  wait "${pid}" || fail=1
done

echo
echo "== Batch summary (with FTL scores) =="
export PYTHONPATH="${SCRIPT_DIR}/../src:${PYTHONPATH:-}"
python3 - <<PY
import json
from pathlib import Path
from grteclyn_wrapper.metrics import read_episode_metrics
from grteclyn_wrapper.score import score_episode

runs = Path("${RUNS_DIR}")
stamp = "${RUN_STAMP}"
ftl_L = float("${FTL_L}")
rows = []
for ep in sorted(runs.glob(f"*_gpu_t*_{stamp}")):
    try:
        score_path = ep / "score.json"
        if score_path.exists():
            payload = json.loads(score_path.read_text())
            total = payload["score"]["total"]
            comps = payload["score"]["components"]
            ftl = payload.get("metrics", {}).get("ftl")
        else:
            m = read_episode_metrics(ep, ftl_L=ftl_L)
            s = score_episode(m, target_stop_time=2.0)
            total = s.total
            comps = s.components
            ftl = None if m.ftl is None else {
                "f_null": m.ftl.f_null,
                "f_portal": m.ftl.f_portal,
                "f_throat": m.ftl.f_throat,
                "s_nonflat": m.ftl.s_nonflat,
                "f_shortcut": m.ftl.f_shortcut,
            }
        rows.append({
            "episode": ep.name,
            "total_score": round(total, 3),
            "ftl_shortcut": round(comps.get("ftl_shortcut", 0.0), 4),
            "nonflat_geometry": round(comps.get("nonflat_geometry", 0.0), 4),
            "stability": round(comps.get("stability", 0.0), 4),
            "survival": round(comps.get("survival", 0.0), 4),
            "ftl": ftl,
        })
    except Exception as exc:
        rows.append({"episode": ep.name, "error": str(exc)})

rows.sort(key=lambda r: r.get("total_score", -1e9), reverse=True)
print(json.dumps(rows, indent=2))
print()
print(f"{'episode':<55} {'score':>7} {'F_FTL':>7} {'s_nf':>6} {'stab':>6}")
print("-" * 85)
for r in rows:
    if "error" in r:
        print(f"{r['episode']:<55} ERROR: {r['error']}")
        continue
    ftl = r.get("ftl") or {}
    print(
        f"{r['episode']:<55} {r['total_score']:7.2f} "
        f"{r['ftl_shortcut']:7.4f} {r['nonflat_geometry']:6.3f} {r['stability']:6.3f}"
    )
PY

exit "${fail}"
