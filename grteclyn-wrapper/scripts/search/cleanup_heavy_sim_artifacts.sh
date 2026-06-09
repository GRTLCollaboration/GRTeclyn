#!/usr/bin/env bash
# Remove leftover AMReX plotfile/checkpoint trees after streaming extraction.
#
# Targets per eval_* folder:
#   - RadialRecipePlt* / RadialRecipeChk* at episode root (consumer keep-last)
#   - RadialRecipePlt* / RadialRecipeChk* under postload_gate/postload_gate/
#
# Keeps: frames/, small_data/, data/*.dat, metadata, score, gridinit, gate JSON.
set -euo pipefail

ROOT="${1:-}"
if [[ -z "${ROOT}" ]]; then
  echo "usage: $0 <qd_campaign_dir>" >&2
  exit 1
fi
if [[ ! -d "${ROOT}" ]]; then
  echo "not a directory: ${ROOT}" >&2
  exit 1
fi

python3 - "${ROOT}" <<'PY'
import shutil
import sys
from pathlib import Path

root = Path(sys.argv[1])
prefixes = ("RadialRecipePlt", "RadialRecipeChk")
removed_dirs = 0
freed = 0

def dir_bytes(path: Path) -> int:
    total = 0
    for f in path.rglob("*"):
        if f.is_file():
            try:
                total += f.stat().st_size
            except OSError:
                pass
    return total

def purge(parent: Path, label: str) -> None:
    global removed_dirs, freed
    if not parent.is_dir():
        return
    for child in parent.iterdir():
        if not child.is_dir():
            continue
        if not any(child.name.startswith(p) for p in prefixes):
            continue
        nbytes = dir_bytes(child)
        shutil.rmtree(child)
        removed_dirs += 1
        freed += nbytes
        print(f"removed {label}/{child.name} ({nbytes / 1e6:.1f} MB)")

for ev in sorted(root.glob("eval_*")):
    purge(ev, ev.name)
    purge(ev / "postload_gate" / "postload_gate", f"{ev.name}/postload_gate")

print(f"done: {removed_dirs} dirs, freed {freed / 1e9:.2f} GB")
PY
