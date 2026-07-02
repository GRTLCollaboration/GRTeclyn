#!/usr/bin/env bash
# Matter-geometry consistency smokes: canonical, exotic, mixed overlapping shell.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/env.sh"

GRTRESNA_ENV="${GRTRESNA_ENV:-/home/jovyan/.mlspace/envs/grtresna}"
CHOMBO_HOME="${CHOMBO_HOME:-${GRTECLYN_ROOT}/../Chombo/lib}"
GRTRESNA_ROOT="${GRTRESNA_ROOT:-${GRTECLYN_ROOT}/../GRTresna}"
SCALAR_BH="${GRTRESNA_ROOT}/Examples/ScalarFieldBH"
EXE="${SCALAR_BH}/Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex"
SMOKE_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/matter_geometry_smokes}/smoke_$(date +%Y%m%dT%H%M%SZ)"
RANKS="${RANKS:-8}"
CUDA_DEVICE="${CUDA_DEVICE:-0}"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
    PYTHON_BIN="${WRAPPER_ROOT}/.venv/bin/python"
  elif command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
    PYTHON_BIN="uv run --project ${WRAPPER_ROOT} python"
  else
    PYTHON_BIN="python"
  fi
fi

export PATH="${GRTRESNA_ENV}/bin:${PATH}"
export CONDA_PREFIX="${GRTRESNA_ENV}"
export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
export GRTRESNA_ROOT

mkdir -p "${SMOKE_DIR}"

if [[ ! -x "${EXE}" ]]; then
  echo "Building GRTresna MPI binary..." >&2
  make -C "${SCALAR_BH}" all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
fi

if [[ ! -x "${GRTECLYN_ROOT}/Examples/RadialRecipe/main3d.gnu.CUDA.ex" ]]; then
  echo "Building RadialRecipe GPU binary..." >&2
  BUILD=1 bash "${WRAPPER_ROOT}/scripts/radial/run_radialrecipe_gpu_smoke.sh" DRY_RUN=1
fi

echo "=== GRTresna AMR solver smokes ==="
for case in canonical exotic mixed_exotic; do
  echo "--- ${case} ---"
  (
    cd "${SCALAR_BH}"
    mpirun --oversubscribe -np "${RANKS}" "./$(basename "${EXE}")" "params_${case}_amr_test.txt"
  )
done

run_postload_smoke() {
  local label="$1"
  local exotic="$2"
  local out="${SMOKE_DIR}/${label}"
  mkdir -p "${out}"

  # shellcheck disable=SC2086
  ${PYTHON_BIN} - <<PY
from pathlib import Path
from grteclyn_wrapper.grtresna.solver import (
    GRTresnaConfig,
    apply_exotic_safe_solver,
    config_has_exotic_lump,
    solve,
)
from grteclyn_wrapper.grtresna.matter.wiring import evolution_overrides_from_config
from grteclyn_wrapper.projection.postload_gate import PostLoadGateConfig, run_postload_gate

out = Path("${out}")
cfg = GRTresnaConfig(
    mpi_ranks=${RANKS},
    max_NL_iterations=30,
    N=(32, 32, 16),
    L=64.0,
    max_level=2,
    phi_0=0.0,
    dphi=0.0,
    pi_0=0.0,
    dpi=0.0,
    scalar_mass=0.1,
    lumps=[
        {
            "amp": 0.15,
            "width": 4.0,
            "center": (0.0, 0.0, 0.0),
            "velocity": (0.12, 0.0, 0.0),
            "omega": 0.0,
            "mode": 0,
            "exotic": ${exotic},
        }
    ],
    gridinit_nx=32,
    gridinit_ny=32,
    gridinit_nz=32,
    target_center=(32.0, 32.0, 0.0),
    cleanup=False,
)
if config_has_exotic_lump(cfg):
    apply_exotic_safe_solver(cfg)
gridinit = solve(cfg, work_dir=out / "grtresna", gridinit_path=out / "initial_data.gridinit")
overrides = {
    "recipe_initial_data_file": str(gridinit),
    **evolution_overrides_from_config(cfg),
    "calculate_constraint_norms": 1,
    "stop_time": 0.01,
}
gate = run_postload_gate(
    gridinit,
    out_dir=out / "postload_gate",
    config=PostLoadGateConfig(max_hamiltonian_l2=1.0e-1, max_momentum_l2=1.0e-1),
    cuda_devices="${CUDA_DEVICE}",
    overrides=overrides,
)
print("${label}", "passed=", gate.passed, "ham=", gate.max_hamiltonian_l2, "mom=", gate.max_momentum_l2)
if not gate.passed:
    raise SystemExit(1)
PY
}

echo "=== Bridge + post-load smokes (independent scalar model) ==="
run_postload_smoke canonical 0
run_postload_smoke exotic 1

# Mixed overlapping shell: two lumps with opposite signs.
# shellcheck disable=SC2086
${PYTHON_BIN} - <<PY
from pathlib import Path
from grteclyn_wrapper.grtresna.solver import (
    GRTresnaConfig,
    apply_exotic_safe_solver,
    config_has_exotic_lump,
    solve,
)
from grteclyn_wrapper.grtresna.matter.wiring import evolution_overrides_from_config
from grteclyn_wrapper.projection.postload_gate import PostLoadGateConfig, run_postload_gate

out = Path("${SMOKE_DIR}/mixed_shell")
cfg = GRTresnaConfig(
    mpi_ranks=${RANKS},
    max_NL_iterations=30,
    N=(32, 32, 16),
    L=64.0,
    max_level=2,
    phi_0=0.0,
    dphi=0.0,
    pi_0=0.0,
    dpi=0.0,
    scalar_mass=0.1,
    lumps=[
        {
            "amp": 0.14,
            "width": 4.0,
            "center": (-2.0, 0.0, 0.0),
            "velocity": (0.10, 0.0, 0.0),
            "omega": 0.0,
            "mode": 0,
            "exotic": 0,
        },
        {
            "amp": 0.12,
            "width": 4.0,
            "center": (2.0, 0.0, 0.0),
            "velocity": (0.08, 0.0, 0.0),
            "omega": 0.0,
            "mode": 0,
            "exotic": 1,
        },
    ],
    gridinit_nx=32,
    gridinit_ny=32,
    gridinit_nz=32,
    target_center=(32.0, 32.0, 0.0),
    cleanup=False,
)
if config_has_exotic_lump(cfg):
    apply_exotic_safe_solver(cfg)
gridinit = solve(cfg, work_dir=out / "grtresna", gridinit_path=out / "initial_data.gridinit")
overrides = {
    "recipe_initial_data_file": str(gridinit),
    **evolution_overrides_from_config(cfg),
    "calculate_constraint_norms": 1,
    "stop_time": 0.01,
}
gate = run_postload_gate(
    gridinit,
    out_dir=out / "postload_gate",
    config=PostLoadGateConfig(max_hamiltonian_l2=1.0e-1, max_momentum_l2=1.0e-1),
    cuda_devices="${CUDA_DEVICE}",
    overrides=overrides,
)
print("mixed_shell", "passed=", gate.passed, "ham=", gate.max_hamiltonian_l2, "mom=", gate.max_momentum_l2)
if not gate.passed:
    raise SystemExit(1)
PY

echo "All matter-geometry smokes passed. Artifacts: ${SMOKE_DIR}"
