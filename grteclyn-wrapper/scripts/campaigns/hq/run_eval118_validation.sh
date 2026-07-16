#!/usr/bin/env bash
# Eval 118 validation matrix launcher (Phases 2–5).
#
# Usage:
#   bash scripts/campaigns/hq/run_eval118_validation.sh E118-RM
#   DRY_RUN=1 bash scripts/campaigns/hq/run_eval118_validation.sh E118-RM
#   FORCE=1 FOREGROUND=1 bash scripts/campaigns/hq/run_eval118_validation.sh E118-RM
#
# Lists available run IDs when called with --list.
set -euo pipefail

# NOTE: env.sh overwrites SCRIPT_DIR; keep launcher path under HQ_DIR.
HQ_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${HQ_DIR}/.." && pwd)"
SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
# shellcheck source=../../lib/env.sh
source "${SCRIPTS_ROOT}/lib/env.sh"
# shellcheck source=../lib/promote_common.sh
source "${CAMPAIGNS_ROOT}/lib/promote_common.sh"

MANIFEST="${MANIFEST:-${HQ_DIR}/eval118_manifest.json}"
RUN_ID="${1:-}"
DRY_RUN="${DRY_RUN:-0}"
FORCE="${FORCE:-0}"
FOREGROUND="${FOREGROUND:-0}"
GPU_ID="${GPU_ID:-0}"
EVOLUTION_MPI_RANKS="${EVOLUTION_MPI_RANKS:-1}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "[eval118] missing manifest: ${MANIFEST}" >&2
  exit 2
fi

if [[ "${RUN_ID}" == "--list" || -z "${RUN_ID}" ]]; then
  python3 - <<'PY' "${MANIFEST}"
import json, sys
m = json.load(open(sys.argv[1], encoding="utf-8"))
print("Eval 118 validation runs:")
for r in m["runs"]:
    deps = ",".join(r.get("depends_on") or []) or "-"
    print(f"  {r['id']:10s}  phase={r['phase']}  N={r['n_full']} L={r['l_full']} t={r.get('stop_time', m['defaults']['stop_time'])}  deps={deps}")
PY
  exit 0
fi

# Resolve run config from manifest into shell variables.
eval "$(python3 - <<'PY' "${MANIFEST}" "${RUN_ID}"
import json, shlex, sys
m = json.load(open(sys.argv[1], encoding="utf-8"))
rid = sys.argv[2]
runs = {r["id"]: r for r in m["runs"]}
if rid not in runs:
    raise SystemExit(f"unknown run id {rid!r}; known={sorted(runs)}")
r = runs[rid]
d = m["defaults"]
source_run = m["source_run"]
eval_id = int(m["source_eval"])
n = int(r["n_full"])
l = float(r["l_full"])
stop = float(r.get("stop_time", d["stop_time"]))
plot = int(r.get("plot_interval", d["plot_interval"]))
radii = r.get("consumer_radii", d["consumer_radii"])
prefix = r["name_prefix"]
print(f"SOURCE_RUN={shlex.quote(source_run)}")
print(f"EVAL_ID={eval_id}")
print(f"N_FULL={n}")
print(f"L_FULL={l}")
print(f"GRTRESNA_DOMAIN_L={l}")
print(f"STOP_TIME={stop}")
print(f"PLOT_INTERVAL={plot}")
print(f"NAME_PREFIX={shlex.quote(prefix)}")
print(f"CONSUMER_RADII={shlex.quote(' '.join(str(x) for x in radii))}")
print(f"RUN_ROLE={shlex.quote(r.get('role',''))}")
print(f"PHASE={int(r['phase'])}")
print(f"OBJECTIVE_MODE={shlex.quote(m.get('objective_mode','general_ftl'))}")
print(f"GEO_EMIT_INTERVAL={float(d['geo_emit_interval'])}")
print(f"GEO_MAX_EMISSIONS={int(d['geo_max_emissions'])}")
print(f"GEO_DIRECTIONS={shlex.quote(d['geo_directions'])}")
print(f"PSI4_N_POINTS={int(d['grteclyn_psi4_n_points'])}")
print(f"FRAMES={int(d['grteclyn_frames'])}")
PY
)"

export SOURCE_RUN="${GRTECLYN_ROOT}/${SOURCE_RUN}"
export N_FULL L_FULL GRTRESNA_DOMAIN_L STOP_TIME PLOT_INTERVAL NAME_PREFIX
export OBJECTIVE_MODE
export MAX_LEVEL="${MAX_LEVEL:-3}"
export REGRID_THRESHOLD="${REGRID_THRESHOLD:-0.02}"
export GRTECLYN_FRAMES="${FRAMES}"
# Metrics-only validation: never keep multi-GB plotfiles for a later frame drain.
unset GRTECLYN_KEEP_PLOTFILES 2>/dev/null || true
export GRTECLYN_PSI4=1
export GRTECLYN_PSI4_N_POINTS="${PSI4_N_POINTS}"
export CONSUMER_RADII
export CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-3}"
export GRTECLYN_GEO_DIRECTIONS="${GEO_DIRECTIONS}"
export GRTECLYN_GEO_EMIT_INTERVAL="${GEO_EMIT_INTERVAL}"
export GRTECLYN_GEO_MAX_EMISSIONS="${GEO_MAX_EMISSIONS}"
export GRTECLYN_EVOLVING_GEODESIC_MODE=hq
export GRTECLYN_EVOLVING_GEODESIC=1
export FORCE FOREGROUND
export CANDIDATES="${EVAL_ID} ${GPU_ID}"

# Refuse overwrite unless FORCE=1.
OUT_DIR="${RUNS_DIR}/${NAME_PREFIX}_hq_eval$(printf '%06d' "${EVAL_ID}")"

if [[ -e "${OUT_DIR}" && "${FORCE}" != "1" ]]; then
  echo "[eval118] output exists: ${OUT_DIR}" >&2
  echo "  set FORCE=1 to replace, or pick another run id" >&2
  exit 3
fi

echo "== Eval 118 validation launch =="
echo "  run_id     : ${RUN_ID} (phase ${PHASE}, ${RUN_ROLE})"
echo "  source     : ${SOURCE_RUN}/eval_$(printf '%06d' "${EVAL_ID}")"
echo "  domain     : L=${L_FULL} N=${N_FULL} ml=${MAX_LEVEL} t=${STOP_TIME}"
echo "  h          : $(python3 -c "print(float('${L_FULL}')/float('${N_FULL}'))")"
echo "  geodesic   : emit_interval=${GEO_EMIT_INTERVAL} max_emissions=${GEO_MAX_EMISSIONS} dirs=${GEO_DIRECTIONS}"
echo "  psi4       : on radii=${CONSUMER_RADII} n_points=${PSI4_N_POINTS}"
echo "  frames     : ${GRTECLYN_FRAMES}"
echo "  output     : ${OUT_DIR}"
echo "  gpu        : ${GPU_ID}"
echo "  mpi ranks  : ${EVOLUTION_MPI_RANKS}"
echo "  dry_run    : ${DRY_RUN}"
echo

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[dry-run] would call:"
  echo "  SOURCE_RUN=${SOURCE_RUN} CANDIDATES=\"${CANDIDATES}\" NAME_PREFIX=${NAME_PREFIX} \\"
  echo "  N_FULL=${N_FULL} L_FULL=${L_FULL} GRTRESNA_DOMAIN_L=${GRTRESNA_DOMAIN_L} \\"
  echo "  STOP_TIME=${STOP_TIME} PLOT_INTERVAL=${PLOT_INTERVAL} MAX_LEVEL=${MAX_LEVEL} \\"
  echo "  GRTECLYN_FRAMES=${GRTECLYN_FRAMES} GRTECLYN_PSI4=1 CONSUMER_RADII=\"${CONSUMER_RADII}\" \\"
  echo "  GRTECLYN_GEO_EMIT_INTERVAL=${GEO_EMIT_INTERVAL} GRTECLYN_GEO_MAX_EMISSIONS=${GEO_MAX_EMISSIONS} \\"
  echo "  FOREGROUND=${FOREGROUND} FORCE=${FORCE} \\"
  echo "  EVOLUTION_MPI_RANKS=${EVOLUTION_MPI_RANKS} \\"
  echo "  bash ${HQ_DIR}/run_batch.sh"
  exit 0
fi

# Persist a launch record next to the eventual episode (or log dir).
LAUNCH_LOG_DIR="${GRTECLYN_ROOT}/research/neuralspacetime/validation/eval118/launches"
mkdir -p "${LAUNCH_LOG_DIR}"
LAUNCH_RECORD="${LAUNCH_LOG_DIR}/${RUN_ID}_$(date -u +%Y%m%dT%H%M%SZ).json"
python3 - <<PY
import json, os, subprocess, pathlib
rec = {
  "run_id": "${RUN_ID}",
  "phase": int("${PHASE}"),
  "role": "${RUN_ROLE}",
  "source_run": "${SOURCE_RUN}",
  "eval_id": int("${EVAL_ID}"),
  "n_full": int("${N_FULL}"),
  "l_full": float("${L_FULL}"),
  "h": float("${L_FULL}") / float("${N_FULL}"),
  "stop_time": float("${STOP_TIME}"),
  "plot_interval": int("${PLOT_INTERVAL}"),
  "max_level": int("${MAX_LEVEL}"),
  "regrid_threshold": float("${REGRID_THRESHOLD}"),
  "consumer_radii": "${CONSUMER_RADII}".split(),
  "geo_emit_interval": float("${GEO_EMIT_INTERVAL}"),
  "geo_max_emissions": int("${GEO_MAX_EMISSIONS}"),
  "gpu_id": "${GPU_ID}",
  "evolution_mpi_ranks": int("${EVOLUTION_MPI_RANKS}"),
  "out_dir": "${OUT_DIR}",
  "git_commit": subprocess.check_output(["git", "-C", "${GRTECLYN_ROOT}", "rev-parse", "HEAD"], text=True).strip(),
  "env": {
    "GRTECLYN_PSI4": os.environ.get("GRTECLYN_PSI4"),
    "GRTECLYN_GEO_EMIT_INTERVAL": os.environ.get("GRTECLYN_GEO_EMIT_INTERVAL"),
    "GRTECLYN_GEO_MAX_EMISSIONS": os.environ.get("GRTECLYN_GEO_MAX_EMISSIONS"),
    "GRTECLYN_FRAMES": os.environ.get("GRTECLYN_FRAMES"),
  },
}
pathlib.Path("${LAUNCH_RECORD}").write_text(json.dumps(rec, indent=2) + "\n")
print("Wrote ${LAUNCH_RECORD}")
PY

export CANDIDATES NAME_PREFIX SOURCE_RUN N_FULL L_FULL GRTRESNA_DOMAIN_L
export STOP_TIME PLOT_INTERVAL MAX_LEVEL REGRID_THRESHOLD OBJECTIVE_MODE
export CONSUMER_RADII GRTECLYN_FRAMES GRTECLYN_PSI4 GRTECLYN_PSI4_N_POINTS
export GRTECLYN_GEO_DIRECTIONS GRTECLYN_GEO_EMIT_INTERVAL GRTECLYN_GEO_MAX_EMISSIONS
export GRTECLYN_EVOLVING_GEODESIC_MODE GRTECLYN_EVOLVING_GEODESIC
export FORCE FOREGROUND
export EVOLUTION_MPI_RANKS

bash "${HQ_DIR}/run_batch.sh"
