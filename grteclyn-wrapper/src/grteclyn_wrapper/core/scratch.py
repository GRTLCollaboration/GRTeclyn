"""Node-local scratch for the plotfile stream.

A plotfile is written once, read once, and deleted.  At 256^3 ml=3 that is
~3.2 GB making the round trip every few minutes, per run.  Sending it over NFS
is what capped this node at two concurrent runs: consumers sat in
uninterruptible I/O wait, backlogs grew, and plotfiles accumulated faster than
they were extracted.  Writing them to the node's own NVMe instead drops shared
filesystem traffic to the ~1.4 MB of numbers each plotfile distils down to, and
the concurrency ceiling becomes the GPU count.

The promotion queues (``scripts/campaigns/rl/pump_ladder_queue.py``,
``pump_confine_queue.py``) each did this by hand, rewriting ``amr.plot_file``
in their own params.  This module is the same policy applied at the one place
every campaign passes through, so QD and CMA-ES inherit it without a launcher
change.

Only the transients move.  ``output_path`` -- ``data/*.dat``, ``small_data/``,
``run.log``, ``params.txt`` -- stays on the shared filesystem, because that is
the part anyone reads after the run.

Set ``GRTECLYN_SCRATCH=0`` to keep plotfiles beside the episode (the old
behaviour), or ``GRTECLYN_SCRATCH=/some/path`` to move the scratch root.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from pathlib import Path

__all__ = [
    "DEFAULT_SCRATCH_ROOT",
    "cache_env",
    "free_gb",
    "plotfile_dir",
    "purge_orphan_scratch",
    "purge_plotfile_scratch",
    "scratch_enabled",
    "scratch_root",
]

DEFAULT_SCRATCH_ROOT = Path("/tmp/grteclyn_scratch")

#: Below this much free space the scratch root is still used, but loudly.  A
#: run needs ``(keep_last + jobs) * plotfile_size``; six concurrent HQ runs is
#: ~70 GB.  Falling back to NFS instead would re-introduce the stall this
#: module exists to avoid, so a tight disk is reported, not worked around.
MIN_FREE_GB = 50.0

_DISABLED = {"0", "off", "no", "none", "false"}

#: Anything AMReX writes into the scratch dir: plotfiles and checkpoints.
_TRANSIENT_RE = re.compile(r"(Plt|Chk)\d+$|^(plt|chk)\d+$")
_PLOTFILE_RE = re.compile(r"Plt\d+$|^plt\d+$")

#: Written into each scratch directory so the mapping can be read backwards.
#: The directory name carries a digest, which does not invert -- without this
#: there is no way to tell an abandoned scratch from a live one.
_OWNER_FILE = ".episode"

_root_cache: dict[str, Path | None] = {}
_warned: set[str] = set()


def _warn_once(key: str, message: str) -> None:
    if key in _warned:
        return
    _warned.add(key)
    print(f"[scratch] {message}", flush=True)


def free_gb(path: Path) -> float:
    """Free space in GB on the filesystem holding ``path`` (0.0 if unknown)."""
    try:
        return shutil.disk_usage(path).free / 1e9
    except OSError:
        return 0.0


def _probe_writable(root: Path) -> bool:
    """Can we actually create files here?  A read-only or full mount cannot."""
    try:
        root.mkdir(parents=True, exist_ok=True)
        handle, name = tempfile.mkstemp(prefix=".probe-", dir=root)
    except OSError:
        return False
    os.close(handle)
    os.unlink(name)
    return True


def scratch_root() -> Path | None:
    """Where transients go, or ``None`` to leave them beside the episode.

    ``None`` means either the operator disabled scratch, or the configured root
    is unusable -- in which case we say so once and fall back rather than
    killing a campaign over a directory.
    """
    raw = os.environ.get("GRTECLYN_SCRATCH", "").strip()
    if raw in _root_cache:
        return _root_cache[raw]

    if raw.lower() in _DISABLED:
        _root_cache[raw] = None
        return None

    root = Path(raw).expanduser() if raw else DEFAULT_SCRATCH_ROOT
    if not _probe_writable(root):
        _warn_once(
            f"unusable:{root}",
            f"{root} is not writable; plotfiles will stay in the episode "
            f"directory. If that is on NFS, expect the consumer to stall. "
            f"Set GRTECLYN_SCRATCH to a node-local path.",
        )
        _root_cache[raw] = None
        return None

    available = free_gb(root)
    if available < MIN_FREE_GB:
        _warn_once(
            f"tight:{root}",
            f"only {available:.0f} GB free on {root}; a run needs roughly "
            f"(keep_last + jobs) x plotfile_size, ~9 GB for one HQ run",
        )
    _root_cache[raw] = root
    return root


def scratch_enabled() -> bool:
    """True when transients are being routed away from the episode directory."""
    return scratch_root() is not None


def _scratch_name(episode_dir: Path) -> str:
    """A directory name that is readable in ``ls`` and still unique.

    Two campaigns both have an ``eval_000007``, so the campaign name alone does
    not separate them and the full path is too long to be a directory name.
    The digest carries the uniqueness; the prefix is there so a human looking
    at the scratch root can tell whose it is.
    """
    resolved = episode_dir.expanduser().resolve()
    digest = hashlib.sha1(str(resolved).encode("utf-8")).hexdigest()[:8]
    return f"{resolved.parent.name}_{resolved.name}_{digest}"


def plotfile_dir(episode_dir: Path | str, *, create: bool = False) -> Path:
    """Where this episode's plotfiles and checkpoints live.

    Falls back to the episode directory itself whenever scratch is off or
    unusable, so every caller can use this unconditionally.
    """
    episode_dir = Path(episode_dir)
    root = scratch_root()
    if root is None:
        return episode_dir.expanduser().resolve()
    target = root / _scratch_name(episode_dir)
    if create:
        try:
            target.mkdir(parents=True, exist_ok=True)
            owner = target / _OWNER_FILE
            if not owner.exists():
                owner.write_text(
                    f"{episode_dir.expanduser().resolve()}\n", encoding="utf-8"
                )
        except OSError as exc:
            _warn_once(
                f"mkdir:{target}",
                f"cannot create {target} ({exc}); falling back to the episode "
                f"directory",
            )
            return episode_dir.expanduser().resolve()
    return target


#: Library caches that otherwise land in ``$HOME`` -- shared, and on this
#: cluster partly not even writable.  Pinning them into scratch keeps the whole
#: pipeline's write set to two places: node-local scratch and the run directory.
_CACHE_VARS = {
    "XDG_CACHE_HOME": "",
    "UV_CACHE_DIR": "uv",
    "MPLCONFIGDIR": "mpl",
    "TMPDIR": "tmp",
    "PYTHONPYCACHEPREFIX": "pyc",
}


def cache_env() -> dict[str, str]:
    """Cache/temp variables to add to a child process, if not already set.

    An operator who has pinned one of these deliberately keeps their value --
    we only fill in what is unset, and only when there is a scratch root to
    put it in.
    """
    root = scratch_root()
    if root is None:
        return {}
    base = root / "_cache"
    env: dict[str, str] = {}
    for name, sub in _CACHE_VARS.items():
        if os.environ.get(name):
            continue
        target = base / sub if sub else base
        try:
            target.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        env[name] = str(target)
    return env


def _load_ledger(episode_dir: Path) -> dict[str, bool]:
    """What ``consume_plotfiles`` has recorded as successfully extracted."""
    path = Path(episode_dir) / "small_data" / "consume_state.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(payload, dict):
        return {}
    return {str(key): bool(value) for key, value in payload.items()}


def purge_plotfile_scratch(
    episode_dir: Path | str, *, force: bool = False
) -> tuple[list[str], str | None]:
    """Reclaim an episode's scratch, keeping anything not yet extracted.

    Deleting from scratch is final -- on NFS a stalled consumer could be re-run
    against the plotfiles that survived, and here it cannot.  So the only
    plotfiles removed are the ones ``consume_state.json`` records as extracted.
    One that failed extraction is left behind and reported, rather than
    quietly costing us the sample it carried.

    ``force`` (or ``GRTECLYN_SCRATCH_FORCE_PURGE=1``) removes everything --
    used when the whole eval is being discarded anyway, so there is no sample
    left to protect.

    Returns the paths removed and, if anything was kept, why.
    """
    episode_dir = Path(episode_dir)
    target = plotfile_dir(episode_dir)
    if target == episode_dir.expanduser().resolve() or not target.is_dir():
        return [], None

    if force or os.environ.get("GRTECLYN_SCRATCH_FORCE_PURGE", "").strip() == "1":
        shutil.rmtree(target, ignore_errors=True)
        return [str(target)], None

    ledger = _load_ledger(episode_dir)
    removed: list[str] = []
    kept: list[str] = []
    for child in sorted(target.iterdir()):
        if not child.is_dir() or not _TRANSIENT_RE.search(child.name):
            continue
        # A checkpoint never enters the ledger (the consumer does not read
        # them), and search runs write none; it is a transient like the rest.
        if _PLOTFILE_RE.search(child.name) and not ledger.get(child.name, False):
            kept.append(child.name)
            continue
        shutil.rmtree(child, ignore_errors=True)
        removed.append(str(child))

    if kept:
        return removed, (
            f"kept {len(kept)} unextracted plotfile(s) in {target}: "
            f"{', '.join(kept[:3])}{'...' if len(kept) > 3 else ''}"
        )

    shutil.rmtree(target, ignore_errors=True)
    return [*removed, str(target)], None


def purge_orphan_scratch() -> list[str]:
    """Drop scratch directories whose episode is gone.

    A campaign killed mid-flight leaves plotfiles behind that nothing will come
    back for -- and unlike the shared filesystem, the local disk is small
    enough that the next campaign notices.  A directory is only removed once
    the episode it names no longer exists, so a live run is never touched, and
    a directory without an owner marker (``_cache``, anything hand-made) is
    left completely alone.
    """
    root = scratch_root()
    if root is None:
        return []
    removed: list[str] = []
    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        owner = child / _OWNER_FILE
        try:
            episode = Path(owner.read_text(encoding="utf-8").strip())
        except OSError:
            continue
        if episode.exists():
            continue
        shutil.rmtree(child, ignore_errors=True)
        removed.append(str(child))
    if removed:
        print(
            f"[scratch] reclaimed {len(removed)} orphaned scratch dir(s) in "
            f"{root}; {free_gb(root):.0f} GB free",
            flush=True,
        )
    return removed
