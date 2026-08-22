#!/usr/bin/env python3
"""List, and optionally kill, the MPI ranks a campaign has left running.

WHY THIS EXISTS (measured 2026-08-22).  ``kill -TERM`` on a launcher's process
group does NOT stop the work.  After three elliptic solves were "stopped" that
way, 133 GRTresna ranks were still in state ``R`` burning 4800% CPU, and the
four GPU evolutions sharing the node stayed starved at a fifth of their proper
speed until someone noticed.  The tell-tale is ``.nfs*`` files left behind in
the run directory: those are handles held by processes that are very much
alive.

AND THE USUAL TOOLS DO NOT WORK HERE.  ``/proc/stat`` and ``/proc/loadavg`` are
broken on this node ("Transport endpoint is not connected"), so ``ps``, ``top``,
``uptime``, ``pgrep``, ``pkill``, ``killall`` and ``free`` all fail or silently
return nothing, and ``nvidia-smi --query-compute-apps`` reports no PIDs.
Counting ``pout.N`` files is not a substitute -- they persist after the writer
dies, so a finished solve looks identical to a running one.

The per-PID files under ``/proc`` are fine, so this walks those directly.  That
is also why killing is done with ``os.kill`` per PID rather than by process
group: the group leader is usually the thing that already exited.

    # what is running, with real per-process CPU
    grteclyn-wrapper/scripts/ops/sweep_ranks.py

    # stop every elliptic solve, leave the GPU evolutions untouched
    grteclyn-wrapper/scripts/ops/sweep_ranks.py --kill solves

    # stop everything belonging to one cell
    grteclyn-wrapper/scripts/ops/sweep_ranks.py --kill all --match <cell-name>

Killing always re-checks afterwards and exits non-zero if anything survived, so
"it is stopped" is a verified statement rather than a hope.
"""

from __future__ import annotations

import argparse
import os
import signal
import sys
import time

# What each kind of process looks like in /proc/<pid>/cmdline.
SOLVE_MARKERS = ("GRTresna", "Main_BosonStarBH", "Main_ScalarField")
EVOLVE_MARKERS = ("RadialRecipe", "main3d.gnu")
CONSUMER_MARKERS = ("plot_consumer", "consume_plotfiles")


def _cmdline(pid: str) -> str | None:
    try:
        with open(f"/proc/{pid}/cmdline", "rb") as fh:
            return fh.read().replace(b"\0", b" ").decode(errors="replace").strip()
    except OSError:
        return None


def _cpu_ticks(pid: str) -> int | None:
    """utime + stime from /proc/<pid>/stat, skipping the parenthesised comm."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            fields = fh.read().rsplit(")", 1)[1].split()
        return int(fields[11]) + int(fields[12])
    except (OSError, IndexError, ValueError):
        return None


def _kind(cmd: str) -> str | None:
    if any(m in cmd for m in SOLVE_MARKERS):
        return "solve"
    if any(m in cmd for m in EVOLVE_MARKERS):
        return "evolve"
    if any(m in cmd for m in CONSUMER_MARKERS):
        return "consumer"
    return None


def scan(match: str | None) -> dict[str, tuple[str, str]]:
    """pid -> (kind, cmdline) for everything this campaign owns."""
    found = {}
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        cmd = _cmdline(pid)
        if not cmd:
            continue
        kind = _kind(cmd)
        if kind is None:
            continue
        if match and match not in cmd:
            continue
        found[pid] = (kind, cmd)
    return found


def measure(pids, window: float) -> dict[str, float]:
    """Per-process CPU percent over a short window -- the only load metric here."""
    hz = os.sysconf("SC_CLK_TCK")
    before = {p: _cpu_ticks(p) for p in pids}
    time.sleep(window)
    out = {}
    for p in pids:
        a, b = before.get(p), _cpu_ticks(p)
        if a is not None and b is not None:
            out[p] = (b - a) / hz / window * 100.0
    return out


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--kill",
        choices=("solves", "evolutions", "consumers", "all"),
        default=None,
        help="what to SIGKILL. Default is to list only and change nothing.",
    )
    ap.add_argument(
        "--match",
        default=None,
        help="only touch processes whose command line contains this (a cell name)",
    )
    ap.add_argument("--window", type=float, default=5.0, help="CPU sample seconds")
    args = ap.parse_args(argv)

    found = scan(args.match)
    if not found:
        print("nothing running" + (f" matching {args.match!r}" if args.match else ""))
        return 0

    cpu = measure(list(found), args.window) if args.kill is None else {}
    by_kind: dict[str, list[str]] = {}
    for pid, (kind, _) in found.items():
        by_kind.setdefault(kind, []).append(pid)

    print(f"{'kind':<10}{'count':>6}{'CPU%':>10}")
    for kind, pids in sorted(by_kind.items()):
        total = sum(cpu.get(p, 0.0) for p in pids)
        print(f"{kind:<10}{len(pids):>6}{total:>10.0f}")
    if cpu:
        cores = os.cpu_count() or 0
        print(f"\ntotal {sum(cpu.values()):.0f}% of {cores * 100}% ({cores} cores)")

    if args.kill is None:
        return 0

    wanted = (
        {"solve", "evolve", "consumer"}
        if args.kill == "all"
        else {"solves": "solve", "evolutions": "evolve", "consumers": "consumer"}[
            args.kill
        ]
    )
    if isinstance(wanted, str):
        wanted = {wanted}
    victims = [p for p, (kind, _) in found.items() if kind in wanted]
    kept = len(found) - len(victims)
    print(f"\nSIGKILL {len(victims)} process(es); leaving {kept} alone")
    for pid in victims:
        try:
            os.kill(int(pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass

    time.sleep(3)
    still = [p for p, (kind, _) in scan(args.match).items() if kind in wanted]
    if still:
        print(f"FAIL  {len(still)} survived -- rerun, they may be in uninterruptible I/O")
        return 1
    print("PASS  none of the targeted processes survive")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
