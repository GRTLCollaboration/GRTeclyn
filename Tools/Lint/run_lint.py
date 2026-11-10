#!/usr/bin/env python3
"""Run the clang-tidy and clang-format checks locally, replicating the lint
GitHub action (.github/workflows/lint.yml) that runs in CI.

Usage:
    Tools/Lint/run_lint.py                          # lint files changed vs HEAD
    Tools/Lint/run_lint.py --base origin/main       # lint files changed on this branch
    Tools/Lint/run_lint.py --all-files             # lint the whole repo
    Tools/Lint/run_lint.py Source/CCZ4/foo.hpp     # lint specific files

Requires the Python packages listed in Tools/Lint/requirements.txt,
installed into a virtual environment. The venv is located via $VIRTUAL_ENV
if one is activated, else Tools/Lint/.venv next to the script.
The compile_commands.json database is generated the same way as in CI and
cached at the repo root; pass --skip-db to reuse it. AMReX is expected next
to the repo (../amrex) unless AMREX_HOME is set.
"""

import argparse
import fnmatch
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from importlib import metadata
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# The directories whose build commands are captured into the compilation
# database, in the same order as in lint.yml
LINT_DIRS = [
    "Examples/BinaryBH",
    "Examples/KleinGordon",
    "Examples/ScalarField",
    "Tests",
]

def python_deps():
    """The required Python packages, read from requirements.txt so the list
    lives in one place. Comments and version specifiers are stripped."""
    deps = []
    for line in (REPO_ROOT / "Tools" / "Lint" / "requirements.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            deps.append(re.split(r"[<>=!~]", line, maxsplit=1)[0].strip())
    return deps

COMPILATION_DATABASE = REPO_ROOT / "compile_commands.json"

DEFAULT_BUILD_ARGS = "COMP=llvm USE_MPI=TRUE"

# The bin directory of the lint virtual environment, resolved per call:
# an activated venv ($VIRTUAL_ENV) wins, else the default Tools/Lint/.venv
# next to the script. Its bin directory is
# prepended to PATH so that versioned tools like clang-tidy-23 installed
# there are found even when the venv is not activated.
DEFAULT_VENV_BIN = REPO_ROOT / "Tools" / "Lint" / ".venv" / "bin"


def lint_venv_bin():
    """The bin directory of the lint virtual environment: $VIRTUAL_ENV if a
    venv is activated, else the default Tools/Lint/.venv next to the script."""
    if os.environ.get("VIRTUAL_ENV"):
        return Path(os.environ["VIRTUAL_ENV"]) / "bin"
    return DEFAULT_VENV_BIN


def tool_path(name):
    """Look up an executable, preferring the lint venv over the system PATH."""
    venv_bin = lint_venv_bin()
    if venv_bin.is_dir():
        candidate = venv_bin / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return shutil.which(name)


def parse_llvm_version():
    """Read LLVM_VERSION from the lint workflow so the local run always uses
    the same clang tools version as CI."""
    lint_yml = REPO_ROOT / ".github" / "workflows" / "lint.yml"
    match = re.search(r"LLVM_VERSION:\s*'(\d+)'", lint_yml.read_text())
    if not match:
        sys.exit(
            "error: could not parse LLVM_VERSION from .github/workflows/lint.yml"
        )
    return match.group(1)


def get_amrex_home():
    home = Path(os.environ.get("AMREX_HOME", REPO_ROOT.parent / "amrex"))
    if not (home / "Src").is_dir():
        sys.exit(
            f"error: AMReX not found at {home}. Clone AMReX next to GRTeclyn "
            "or set the AMREX_HOME environment variable."
        )
    return home


def installed_version(package):
    """Version of an installed Python package in the lint venv (or the
    current environment), via importlib.metadata."""
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def latest_version(package):
    """Latest version of a package on PyPI. CI installs these unpinned, so it
    always uses the latest release."""
    url = f"https://pypi.org/pypi/{package}/json"
    with urllib.request.urlopen(url, timeout=10) as response:
        return json.loads(response.read())["info"]["version"]


def check_dependencies(skip_check):
    deps = python_deps()
    missing = [pkg for pkg in deps if tool_path(pkg) is None]
    if missing:
        sys.exit(
            "error: missing Python dependencies: "
            + ", ".join(missing)
            + "\nInstall them into a venv next to the script:\n"
            "  python3 -m venv Tools/Lint/.venv\n"
            "  source Tools/Lint/.venv/bin/activate\n"
            "  pip install -r Tools/Lint/requirements.txt"
        )
    if skip_check:
        return
    for package in deps:
        installed = installed_version(package)
        try:
            latest = latest_version(package)
        except Exception:
            print(
                f"note: could not reach PyPI to check the latest version of "
                f"{package}; skipping the version check"
            )
            continue
        if installed and installed != latest:
            print(
                f"warning: {package} {installed} is installed but {latest} is "
                "the latest on PyPI. CI installs the latest, so its results "
                "may differ from yours."
            )
            print(
                "  upgrade with: Tools/Lint/.venv/bin/pip install --upgrade "
                "-r Tools/Lint/requirements.txt"
            )


def venv_env():
    """A copy of the environment with the lint venv's bin directory prepended
    to PATH, so the venv's compiledb/compdb/cpp-linter/clang-tools are used
    even when the venv is not activated."""
    env = dict(os.environ)
    venv_bin = lint_venv_bin()
    if venv_bin.is_dir():
        env["PATH"] = f"{venv_bin}{os.pathsep}{env.get('PATH', '')}"
    return env


def generate_compilation_database(build_args, amrex_home):
    """Generate compile_commands.json exactly as the lint workflow does:
    dry-run make with compiledb in each example/test directory, restore the
    compiler executable compiledb drops, then add header entries with compdb."""
    env = dict(venv_env(), AMREX_HOME=str(amrex_home))
    args = build_args.split()
    for i, rel_dir in enumerate(LINT_DIRS):
        directory = REPO_ROOT / rel_dir
        if i == 0:
            subprocess.run(
                ["make", *args, "AMReX_Config.H"], cwd=directory, env=env, check=True
            )
        print(f"capturing compile commands from {rel_dir} ...")
        subprocess.run(
            ["compiledb", "-o", str(COMPILATION_DATABASE), "-n", "make", *args],
            cwd=directory,
            env=env,
            check=True,
        )

    # compiledb drops the compiler executable from the captured commands, so
    # clang-tidy would consume the first flag (e.g. -std=c++20) as the compiler
    # name. Restore a compiler so the flags are not lost (as in lint.yml).
    database = json.loads(COMPILATION_DATABASE.read_text())
    for entry in database:
        arguments = entry.get("arguments")
        if arguments and arguments[0].startswith("-"):
            arguments.insert(0, "clang++")
    COMPILATION_DATABASE.write_text(json.dumps(database, indent=1) + "\n")

    # compdb adds entries for headers to the database
    with_headers = REPO_ROOT / "compile_commands_with_headers.json"
    with open(with_headers, "w") as handle:
        subprocess.run(
            ["compdb", "-p", str(REPO_ROOT), "list"],
            cwd=REPO_ROOT,
            env=env,
            stdout=handle,
            check=True,
        )
    os.replace(with_headers, COMPILATION_DATABASE)
    print(f"compilation database written to {COMPILATION_DATABASE}")


def build_ignore_tidy():
    """Build the ignore-tidy pattern list from .lint-ignore using the same
    filtering as the lint workflow: drop comments, blank lines and the
    lint-only: prefix, then join with '|'."""
    patterns = []
    for line in (REPO_ROOT / ".lint-ignore").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        patterns.append(line.removeprefix("lint-only:").strip())
    return "|".join(patterns)


def resolve_clang_version_arg(llvm_version, ignore_version_check):
    """The value to pass to cpp-linter's --version: the version number, provided
    the versioned executables clang-tidy-<version>/clang-format-<version> are
    found (on PATH, or in the lint venv's bin directory, where the clang-tools
    pip package installs them via 'clang-tools install clang-tidy clang-format
    --version <version>'). Exits with an error if they are missing, unless
    --ignore-llvm-version was given."""
    venv_bin = lint_venv_bin()
    if all(
        tool_path(f"{tool}-{llvm_version}")
        for tool in ("clang-tidy", "clang-format")
    ):
        return llvm_version
    if not ignore_version_check:
        sys.exit(
            f"error: clang-tidy-{llvm_version} and/or clang-format-{llvm_version} "
            f"not found (searched PATH and {venv_bin}).\n"
            f"CI lints with LLVM {llvm_version} (from .github/workflows/lint.yml), "
            "so a different local version may report different results.\n"
            "Install the matching tools, e.g.:\n"
            "  python3 -m venv Tools/Lint/.venv\n"
            "  Tools/Lint/.venv/bin/pip install -r Tools/Lint/requirements.txt\n"
            f"  source Tools/Lint/.venv/bin/activate\n"
            f"  clang-tools install clang-format clang-tidy --version {llvm_version} --directory $VIRTUAL_ENV/bin\n"
            "Re-run with --ignore-llvm-version to skip this check."
        )
    print(
        f"warning: clang-tidy-{llvm_version}/clang-format-{llvm_version} not "
        "found; cpp-linter will fall back to whatever clang tools it finds, "
        "which may report different results than CI (--ignore-llvm-version "
        "was given)."
    )
    if venv_bin.is_dir() and all(
        (venv_bin / tool).is_file() for tool in ("clang-tidy", "clang-format")
    ):
        return str(venv_bin)
    return llvm_version


def run_linter(version_arg, ignore_tidy, files_changed, base, jobs, files):
    command = [
        "cpp-linter",
        "--verbosity=info",
        "--style=file",  # clang-format style from .clang-format
        "--tidy-checks=",  # use the checks in the .clang-tidy file
        "--repo-root",
        str(REPO_ROOT),
        "--ignore=External",
        "--ignore-tidy",
        ignore_tidy,
        "--version",
        version_arg,  # clang tools version, resolved to match CI
        "--database",
        str(REPO_ROOT),  # compilation database directory
        "--lines-changed-only=false",
        "--no-lgtm=true",  # no LGTM comments outside GitHub
        "--thread-comments=false",  # no PR thread comments
        "--step-summary=false",
        "--jobs",
        str(jobs),
    ]
    if files_changed:
        command += ["--files-changed-only=true", "--diff-base", base]
    else:
        command += ["--files-changed-only=false"]
    if files:
        # cpp-linter's CLI declares 'files' as a positional with nargs='*'
        # alongside a 'version' subcommand, and its argparse setup rejects
        # trailing positional tokens, so explicit files cannot be passed to
        # cpp-linter directly. Run clang-tidy and clang-format on them
        # directly with the same configuration instead.
        return lint_specific_files(version_arg, ignore_tidy, files)
    return _run_cpp_linter(command)


def _run_cpp_linter(command):
    """Run cpp-linter and return 1 if it reported any failed checks (like the
    'Fail job if checks fail' step in lint.yml), 0 otherwise."""
    print("running:", " ".join(command))
    env = venv_env()
    result = subprocess.run(
        command, cwd=REPO_ROOT, env=env, capture_output=True, text=True
    )
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    match = re.search(r"INFO:CPP Linter:(\d+) checks-failed", result.stdout + result.stderr)
    if match and int(match.group(1)) > 0:
        return 1
    return result.returncode


def lint_specific_files(version_arg, ignore_tidy, files):
    """Lint specific files by invoking clang-tidy and clang-format directly
    with the same configuration cpp-linter uses (see run_linter). version_arg
    is the resolved --version value: a major version number or a bin
    directory."""
    ignore_patterns = ignore_tidy.split("|") if ignore_tidy else []

    def is_ignored(path):
        return any(fnmatch.fnmatch(path, pat) for pat in ignore_patterns)

    if Path(version_arg).is_dir():
        tidy = tool_path("clang-tidy")
        format = tool_path("clang-format")
    else:
        tidy = tool_path(f"clang-tidy-{version_arg}") or tool_path("clang-tidy")
        format = (
            tool_path(f"clang-format-{version_arg}")
            or tool_path("clang-format")
        )
    returncode = 0
    for path in files:
        if is_ignored(path):
            print(f"skipping {path} (in .lint-ignore)")
            continue
        print(f"linting {path}")
        result = subprocess.run(
            # --use-color keeps the ANSI colours even though the output is
            # piped through the filter below
            [tidy, f"-p={REPO_ROOT}", "--config-file=.clang-tidy", "--use-color", path],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        # clang-tidy prints its diagnostics on stdout and reports warnings
        # from all included headers too (see HeaderFilterRegex in
        # .clang-tidy); cpp-linter filters these down to the file being
        # checked, so do the same here
        output = result.stdout + result.stderr
        found = False
        # ANSI colour escapes may appear anywhere in the line when
        # --use-color is on; strip them for matching, keep them for printing
        ansi = re.compile(r"\x1b\[[0-9;]*m")
        diagnostic = re.compile(r"^(?P<file>.*?):(?P<line>\d+):(?P<col>\d+): warning:")
        lines = output.splitlines()
        for i, line in enumerate(lines):
            match = diagnostic.match(ansi.sub("", line))
            if match:
                file_part = match.group("file")
                # clang-tidy prints a bare filename for warnings in the main
                # file when invoked with a relative path, and absolute paths
                # for headers
                if file_part in (path, str(REPO_ROOT / path), Path(path).name) or file_part.endswith(f"/{path}"):
                    # print the warning and its context (source excerpt and
                    # notes) up to the next diagnostic line
                    print(line)
                    found = True
                    # the warning's context (source excerpt and notes, which
                    # may reference other files) runs until the next
                    # warning/error or a summary line
                    for context_line in lines[i + 1:]:
                        stripped = ansi.sub("", context_line)
                        if re.match(r"^.*:\d+:\d+: (warning|error):", stripped):
                            break
                        if re.match(r"^(Suppressed|.* warnings? generated|.* warnings? and|Use -|Error while|Found compiler)", stripped):
                            break
                        print(context_line)
        if not found and result.returncode != 0:
            print(output)
        rc = 1 if found else result.returncode
        returncode = max(returncode, rc)
        rc = subprocess.run(
            [format, "--style=file", "--dry-run", "--Werror", path],
            cwd=REPO_ROOT,
        ).returncode
        returncode = max(returncode, rc)
    return returncode


def main():
    parser = argparse.ArgumentParser(
        description="Run the clang-tidy and clang-format checks locally, "
        "replicating .github/workflows/lint.yml."
    )
    parser.add_argument(
        "--base",
        default="HEAD",
        metavar="BASE",
        help="git ref to diff against when linting only changed files "
        "(default: HEAD; use e.g. origin/main to lint everything on your "
        "branch)",
    )
    parser.add_argument(
        "--all-files",
        action="store_true",
        help="lint the whole repository instead of only changed files. This "
        "matches what CI does when files-changed-only is false; note that a "
        "clean changed-files run does not guarantee a clean all-files run.",
    )
    parser.add_argument(
        "--skip-db",
        action="store_true",
        help="reuse the existing compile_commands.json instead of regenerating it",
    )
    parser.add_argument(
        "--build-args",
        default=os.environ.get("BUILD_ARGS", DEFAULT_BUILD_ARGS),
        help="make arguments used to generate the compilation database "
        f"(default: '{DEFAULT_BUILD_ARGS}'; env var BUILD_ARGS is also read)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(4, os.cpu_count() or 1),
        help="number of parallel clang-tidy jobs (default: min(4, CPU count))",
    )
    parser.add_argument(
        "--skip-dep-check",
        action="store_true",
        help="skip the PyPI check that the installed compiledb/compdb/cpp-linter "
        "versions match the latest (useful offline)",
    )
    parser.add_argument(
        "--ignore-llvm-version",
        action="store_true",
        help="skip the check that clang-tidy/clang-format match the LLVM version "
        "used in CI (warns instead of failing)",
    )
    parser.add_argument(
        "files",
        nargs="*",
        help="specific files to lint, relative to the repo root, e.g. "
        "'run_lint.py --skip-db Source/CCZ4/CCZ4.hpp'",
    )
    args = parser.parse_args()

    check_dependencies(args.skip_dep_check)
    if not COMPILATION_DATABASE.exists() or not args.skip_db:
        generate_compilation_database(args.build_args, get_amrex_home())
    ignore_tidy = build_ignore_tidy()
    llvm_version = parse_llvm_version()
    version_arg = resolve_clang_version_arg(llvm_version, args.ignore_llvm_version)
    sys.exit(
        run_linter(
            version_arg,
            ignore_tidy,
            not args.all_files,
            args.base,
            args.jobs,
            args.files,
        )
    )


if __name__ == "__main__":
    main()
