#!/usr/bin/env python3
"""Rewrite machine-specific absolute paths in a packed artefact file.

Same contract as the scrubber inside
``research/neuralspacetime/pack_publishable_results.sh``: identity tokens are
derived at RUNTIME from the environment (never hard-coded host/user literals),
so the packed results carry no machine identity into git.

Usage: scrub_paths.py FILE [FILE ...]
Env knobs (optional): ROOT, SIM_ROOT, GRTRESNA_ENV, OPENMPI_ROOT, CHOMBO_HOME.
"""

from __future__ import annotations

import os
import pathlib
import re
import socket
import sys

# Path segments that are product/campaign names and must survive scrubbing.
KEEP = {
    "grteclyn", "grtresna", "chombo", "grteclyn-wrapper", "amrex", "openmpi",
    "radialrecipe", "bondi", "neuralspacetime",
}
GENERIC = {
    "", "home", "usr", "tmp", "var", "opt", "etc", "lib", "bin", "local",
    "path", "to", "users", "user", "envs", "env", "conda", "mamba",
    "test", "simulation", "src", "lib64", "share", "include", "runs",
}
BINARY_SUFFIXES = {".png", ".jpg", ".jpeg", ".pdf", ".npz", ".bin", ".ex", ".so", ".mp4"}


def _identity_tokens() -> set[str]:
    tokens: set[str] = set()

    def add(value: str) -> None:
        value = (value or "").strip().strip("/.")
        if len(value) < 3 or value.lower() in GENERIC or value.lower() in KEEP:
            return
        if re.split(r"[-_]", value, maxsplit=1)[0].lower() in KEEP:
            return  # versioned product dir, e.g. openmpi-5.0.8
        tokens.add(value)

    add(os.environ.get("USER", ""))
    add(os.environ.get("LOGNAME", ""))
    add(pathlib.Path(os.environ.get("HOME", "")).name)
    host = socket.gethostname().split(".")[0]
    add(socket.gethostname())
    add(host)
    for part in re.split(r"[-_.]", host):
        if len(part) >= 4:
            add(part)
    for key in ("ROOT", "SIM_ROOT", "GRTRESNA_ENV", "OPENMPI_ROOT", "CHOMBO_HOME", "HOME"):
        for part in pathlib.Path(os.environ.get(key) or "").parts:
            add(part)
    return tokens


def scrub(path: pathlib.Path) -> None:
    if path.suffix.lower() in BINARY_SUFFIXES:
        return
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return

    root = os.path.abspath(os.environ.get("ROOT", ""))
    sim = os.path.abspath(os.environ.get("SIM_ROOT", ""))

    prefixes = []
    for key in ("ROOT", "SIM_ROOT", "GRTRESNA_ENV", "OPENMPI_ROOT", "CHOMBO_HOME", "HOME"):
        value = os.environ.get(key) or ""
        if value and os.path.isabs(value):
            prefixes.append(os.path.abspath(value))
    for prefix in sorted(set(prefixes), key=len, reverse=True):
        if prefix == root:
            text = text.replace(prefix + "/", "").replace(prefix, "$GRTECLYN_ROOT")
        elif prefix == sim:
            text = text.replace(prefix + "/", "").replace(prefix, "$SIM_ROOT")
        else:
            text = text.replace(prefix + "/", "$SITE/").replace(prefix, "$SITE")

    # Any remaining absolute home-style path (pattern built at runtime).
    home = "/" + "home" + "/"
    text = re.sub(rf"(?i){home}[^/\s\"']+/[^\s\"']*", r"$HOME/<redacted>", text)
    text = re.sub(rf"(?i){home}[^/\s\"']+", r"$HOME/<redacted>", text)

    # Word-boundary replace so campaign names (bondi_sg_pair_pm) stay intact.
    for token in sorted(_identity_tokens(), key=len, reverse=True):
        text = re.sub(
            rf"(?<![A-Za-z0-9_]){re.escape(token)}(?![A-Za-z0-9_])",
            "<redacted>",
            text,
            flags=re.IGNORECASE,
        )

    path.write_text(text, encoding="utf-8")


def main(argv: list[str]) -> int:
    for arg in argv:
        scrub(pathlib.Path(arg))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
