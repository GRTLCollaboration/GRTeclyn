#!/usr/bin/env bash
# Generic manifest-driven HQ / Richardson matrix launcher.
#
#   MANIFEST=.../manifest.json bash scripts/campaigns/promote/lib/run_matrix.sh BCMA-RM
#   bash scripts/campaigns/promote/lib/run_matrix.sh --list
#
# Calls campaigns/hq/run_batch.sh (shared replay engine).
set -euo pipefail

PROMOTE_LIB="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${PROMOTE_LIB}/../.." && pwd)"
SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
HQ_DIR="${CAMPAIGNS_ROOT}/hq"
# shellcheck source=../../../lib/env.sh
source "${SCRIPTS_ROOT}/lib/env.sh"
# shellcheck source=../../lib/promote_common.sh
source "${CAMPAIGNS_ROOT}/lib/promote_common.sh"

MANIFEST="${MANIFEST:-}"
RUN_ID="${1:-}"
DRY_RUN="${DRY_RUN:-0}"
FORCE="${FORCE:-0}"
FOREGROUND="${FOREGROUND:-0}"
GPU_ID="${GPU_ID:-0}"
EVOLUTION_MPI_RANKS="${EVOLUTION_MPI_RANKS:-}"

if [[ -z "${MANIFEST}" || ! -f "${MANIFEST}" ]]; then
  echo "[matrix] set MANIFEST to a campaign manifest.json" >&2
  exit 2
fi

# Pump convention (GPU_RUN_PLAN.md §12.1): the pump runs for the entire
# simulation. Refuse any manifest/env that stops it mid-run unless the
# manifest declares itself a deliberate pump-off control
# ("pump_off_control": true), and require the emission floor on pump-on
# manifests. Runs on every --list and launch.
python3 "${PROMOTE_LIB}/validate_pump_convention.py" "${MANIFEST}"

if [[ "${RUN_ID}" == "--list" || -z "${RUN_ID}" ]]; then
  python3 - <<'PY' "${MANIFEST}"
import json, sys
m = json.load(open(sys.argv[1], encoding="utf-8"))
print(f"{m.get('campaign', 'matrix')} validation runs:")
print(f"  source_run : {m.get('source_run')}")
print(f"  source_eval: {m.get('source_eval')}")
for r in m["runs"]:
    deps = ",".join(r.get("depends_on") or []) or "-"
    ranks = int(r.get("evolution_mpi_ranks", m.get("defaults", {}).get("evolution_mpi_ranks", 1)))
    frames = int(r.get("grteclyn_frames", m["defaults"].get("grteclyn_frames", 0)))
    print(
        f"  {r['id']:10s}  phase={r['phase']}  N={r['n_full']} L={r['l_full']} "
        f"t={r.get('stop_time', m['defaults']['stop_time'])}  "
        f"gpus={ranks}  frames={frames}  deps={deps}"
    )
PY
  exit 0
fi

eval "$(
  python3 - <<'PY' "${MANIFEST}" "${RUN_ID}" "${GRTECLYN_ROOT}" "${SOURCE_EVAL_ID:-}" "${SOURCE_RUN:-}"
import json, shlex, sys
from pathlib import Path

m = json.load(open(sys.argv[1], encoding="utf-8"))
rid = sys.argv[2]
root = Path(sys.argv[3])
env_eval = sys.argv[4].strip()
env_source = sys.argv[5].strip()

runs = {r["id"]: r for r in m["runs"]}
if rid not in runs:
    raise SystemExit(f"unknown run id {rid!r}; known={sorted(runs)}")
r = runs[rid]
d = m["defaults"]

source_run = env_source or m["source_run"]
src_path = Path(source_run)
if not src_path.is_absolute():
    src_path = root / source_run

sev = m.get("source_eval", "auto")
if env_eval:
    eval_id = int(env_eval)
elif sev in (None, "auto", ""):
    traj = src_path / "trajectory.jsonl"
    if not traj.is_file():
        raise SystemExit(f"auto source_eval needs trajectory: {traj}")
    best = None
    with traj.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("status") != "gpu_ok":
                continue
            score = float(rec.get("score") or -1e300)
            if best is None or score > best[0]:
                best = (score, int(rec["eval"]))
    if best is None:
        raise SystemExit(f"no gpu_ok rows in {traj}")
    eval_id = best[1]
else:
    eval_id = int(sev)

n = int(r["n_full"])
l = float(r["l_full"])
stop = float(r.get("stop_time", d["stop_time"]))
plot = int(r.get("plot_interval", d["plot_interval"]))
radii = r.get("consumer_radii", d["consumer_radii"])
prefix = r["name_prefix"]
extras = m.get("extra_overrides") or []
frames = int(r.get("grteclyn_frames", d.get("grteclyn_frames", 0)))
ranks = int(r.get("evolution_mpi_ranks", d.get("evolution_mpi_ranks", 1)))
print(f"SOURCE_RUN={shlex.quote(str(src_path))}")
print(f"EVAL_ID={eval_id}")
print(f"N_FULL={n}")
print(f"L_FULL={l}")
print(f"GRTRESNA_DOMAIN_L={l}")
print(f"STOP_TIME={stop}")
print(f"PLOT_INTERVAL={plot}")
print(f"NAME_PREFIX={shlex.quote(prefix)}")
print(f"CONSUMER_RADII={shlex.quote(' '.join(str(x) for x in radii))}")
print(f"RUN_ROLE={shlex.quote(r.get('role', ''))}")
print(f"PHASE={int(r['phase'])}")
print(f"OBJECTIVE_MODE={shlex.quote(m.get('objective_mode', 'general_ftl'))}")
print(f"GEO_EMIT_INTERVAL={float(d['geo_emit_interval'])}")
print(f"GEO_MAX_EMISSIONS={int(d['geo_max_emissions'])}")
print(f"GEO_DIRECTIONS={shlex.quote(d['geo_directions'])}")
print(f"PSI4_N_POINTS={int(d['grteclyn_psi4_n_points'])}")
print(f"FRAMES={frames}")
print(f"MANIFEST_MPI_RANKS={ranks}")
print(f"CAMPAIGN_SLUG={shlex.quote(str(m.get('campaign', 'matrix')))}")
print(f"MANIFEST_EXTRA_OVERRIDE={shlex.quote(' '.join(str(x) for x in extras))}")
PY
)"

export SOURCE_RUN N_FULL L_FULL GRTRESNA_DOMAIN_L STOP_TIME PLOT_INTERVAL NAME_PREFIX
export OBJECTIVE_MODE
export MAX_LEVEL="${MAX_LEVEL:-3}"
export REGRID_THRESHOLD="${REGRID_THRESHOLD:-0.02}"
export GRTECLYN_FRAMES="${FRAMES}"
if [[ "${FRAMES}" == "0" ]]; then
  unset GRTECLYN_KEEP_PLOTFILES 2>/dev/null || true
else
  export CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-4}"
fi
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

if [[ -z "${EVOLUTION_MPI_RANKS}" ]]; then
  EVOLUTION_MPI_RANKS="${MANIFEST_MPI_RANKS:-1}"
fi
export EVOLUTION_MPI_RANKS
if (( EVOLUTION_MPI_RANKS > 1 )) && [[ "${GPU_ID}" != *,* ]]; then
  _start="${GPU_ID}"
  _ids=()
  for (( _i = 0; _i < EVOLUTION_MPI_RANKS; _i++ )); do
    _ids+=("$((_start + _i))")
  done
  GPU_ID="$(IFS=,; echo "${_ids[*]}")"
fi
export GPU_ID
export CANDIDATES="${EVAL_ID} ${GPU_ID}"

if [[ -n "${MANIFEST_EXTRA_OVERRIDE:-}" ]]; then
  EXTRA_OVERRIDE="${MANIFEST_EXTRA_OVERRIDE} ${EXTRA_OVERRIDE:-}"
fi
export EXTRA_OVERRIDE="${EXTRA_OVERRIDE:-}"

SOURCE_EVAL_DIR="${SOURCE_RUN}/eval_$(printf '%06d' "${EVAL_ID}")"
if [[ ! -d "${SOURCE_EVAL_DIR}" ]]; then
  echo "[matrix] missing source eval dir: ${SOURCE_EVAL_DIR}" >&2
  echo "  Freeze the champion first: ${FREEZE_HINT:-<campaign>/freeze.sh}" >&2
  exit 4
fi
if [[ ! -f "${SOURCE_EVAL_DIR}/metadata.json" ]]; then
  echo "[matrix] source eval missing metadata.json: ${SOURCE_EVAL_DIR}" >&2
  exit 4
fi

OUT_DIR="${RUNS_DIR}/${NAME_PREFIX}_hq_eval$(printf '%06d' "${EVAL_ID}")"
if [[ -e "${OUT_DIR}" && "${FORCE}" != "1" ]]; then
  echo "[matrix] output exists: ${OUT_DIR}" >&2
  echo "  set FORCE=1 to replace, or pick another run id" >&2
  exit 3
fi

echo "== Matrix validation launch =="
echo "  campaign   : ${CAMPAIGN_SLUG}"
echo "  run_id     : ${RUN_ID} (phase ${PHASE}, ${RUN_ROLE})"
echo "  source     : ${SOURCE_EVAL_DIR}"
echo "  domain     : L=${L_FULL} N=${N_FULL} ml=${MAX_LEVEL} t=${STOP_TIME}"
echo "  h          : $(python3 -c "print(float('${L_FULL}')/float('${N_FULL}'))")"
echo "  frames     : ${GRTECLYN_FRAMES}"
echo "  gpu        : ${GPU_ID}  (mpi ranks=${EVOLUTION_MPI_RANKS})"
echo "  output     : ${OUT_DIR}"
echo "  dry_run    : ${DRY_RUN}"
echo

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[dry-run] would call ${HQ_DIR}/run_batch.sh"
  echo "  CANDIDATES=\"${CANDIDATES}\" EXTRA_OVERRIDE=\"${EXTRA_OVERRIDE}\" N=${N_FULL} L=${L_FULL}"
  exit 0
fi

LAUNCH_LOG_DIR="${VALIDATION_LAUNCH_LOG_DIR:-${GRTECLYN_ROOT}/research/neuralspacetime/validation/${CAMPAIGN_SLUG}/launches}"
mkdir -p "${LAUNCH_LOG_DIR}"
LAUNCH_RECORD="${LAUNCH_LOG_DIR}/${RUN_ID}_$(date -u +%Y%m%dT%H%M%SZ).json"
python3 - <<PY
import json, os, subprocess, pathlib

_ROOT = "${GRTECLYN_ROOT}"

def _rel(path):
    """Repo-relative so launch records never carry host/user/home literals."""
    if not path:
        return path
    try:
        return os.path.relpath(path, _ROOT)
    except ValueError:
        return path

rec = {
  "run_id": "${RUN_ID}",
  "campaign": "${CAMPAIGN_SLUG}",
  "phase": int("${PHASE}"),
  "role": "${RUN_ROLE}",
  "manifest": _rel("${MANIFEST}"),
  "source_run": _rel("${SOURCE_RUN}"),
  "eval_id": int("${EVAL_ID}"),
  "n_full": int("${N_FULL}"),
  "l_full": float("${L_FULL}"),
  "h": float("${L_FULL}") / float("${N_FULL}"),
  "stop_time": float("${STOP_TIME}"),
  "plot_interval": int("${PLOT_INTERVAL}"),
  "max_level": int("${MAX_LEVEL}"),
  "gpu_id": "${GPU_ID}",
  "evolution_mpi_ranks": int("${EVOLUTION_MPI_RANKS}"),
  "out_dir": _rel("${OUT_DIR}"),
  "git_commit": subprocess.check_output(["git", "-C", "${GRTECLYN_ROOT}", "rev-parse", "HEAD"], text=True).strip(),
}
pathlib.Path("${LAUNCH_RECORD}").write_text(json.dumps(rec, indent=2) + "\n")
print("Wrote ${LAUNCH_RECORD}")
PY

export CANDIDATES NAME_PREFIX SOURCE_RUN N_FULL L_FULL GRTRESNA_DOMAIN_L
export STOP_TIME PLOT_INTERVAL MAX_LEVEL REGRID_THRESHOLD OBJECTIVE_MODE
export CONSUMER_RADII GRTECLYN_FRAMES GRTECLYN_PSI4 GRTECLYN_PSI4_N_POINTS
export GRTECLYN_GEO_DIRECTIONS GRTECLYN_GEO_EMIT_INTERVAL GRTECLYN_GEO_MAX_EMISSIONS
export GRTECLYN_EVOLVING_GEODESIC_MODE GRTECLYN_EVOLVING_GEODESIC
export FORCE FOREGROUND EVOLUTION_MPI_RANKS EXTRA_OVERRIDE

bash "${HQ_DIR}/run_batch.sh"
