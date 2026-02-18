#!/usr/bin/env python3
"""
Plot constraint norms from constraint_norms.dat (SmallDataIO format).

Expected data format (whitespace separated):
  t  L2_Ham  L2_Mom

Usage:
  python plot_constraints.py [path/to/constraint_norms.dat]
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Use non-interactive backend by default
matplotlib.use("Agg")

# Set style to match scientific publication standards (similar to hamiltonian/__main__.py)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
    "mathtext.fontset": "stix",
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.linewidth": 1.2,
    "grid.alpha": 0.5,
    "legend.fontsize": 12,
    "figure.figsize": (10, 8),
})


def parse_constraint_file(filepath: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parses the constraint_norms.dat file.
    Returns arrays: time, L2_Ham, L2_Mom.
    Handles headers and restart lines.
    """
    times = []
    l2_ham = []
    l2_mom = []

    with open(filepath, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            # Header check: if line contains non-numeric characters, skip
            try:
                # Try parsing the first token as float
                float(parts[0])
            except ValueError:
                # Likely a header line like "L2_Ham L2_Mom"
                continue

            if len(parts) < 3:
                print(f"Warning: Line {line_num} has fewer than 3 columns. Skipping.")
                continue

            try:
                t = float(parts[0])
                ham = float(parts[1])
                mom = float(parts[2])

                times.append(t)
                l2_ham.append(ham)
                l2_mom.append(mom)
            except ValueError as e:
                print(f"Warning: Could not parse line {line_num}: {e}")
                continue

    if not times:
        raise ValueError(f"No valid data found in {filepath}")

    # Convert to numpy arrays
    times = np.array(times)
    l2_ham = np.array(l2_ham)
    l2_mom = np.array(l2_mom)

    # Sort by time (in case of restarts that appended out of order, though typically SmallDataIO appends correctly)
    sort_idx = np.argsort(times)

    # Remove duplicates (if any) - simply keep unique times if exact duplicates exist
    # For simplicity, just sorting here. A more robust approach might filter exact duplicates.

    return times[sort_idx], l2_ham[sort_idx], l2_mom[sort_idx]


def plot_constraints(times, l2_ham, l2_mom, output_path):
    """
    Plots the constraint norms.
    """
    fig, axes = plt.subplots(2, 1, sharex=True)

    # Hamiltonian Plot
    ax1 = axes[0]
    ax1.semilogy(times, l2_ham, label=r'$\|\mathcal{H}\|_{L^2}$', color='blue', linewidth=1.5)
    ax1.set_ylabel(r'$\|\mathcal{H}\|_{L^2}$')
    ax1.set_title(r'Hamiltonian Constraint: $\mathcal{H}$')
    ax1.grid(True, which="both", ls="--", alpha=0.6)
    ax1.legend(loc='upper right', frameon=True, framealpha=0.9)
    ax1.tick_params(axis='both', which='major', direction='in', top=True, right=True)

    # Momentum Plot
    ax2 = axes[1]
    ax2.semilogy(times, l2_mom, label=r'$\|\mathcal{M}\|_{L^2}$', color='red', linewidth=1.5)
    ax2.set_xlabel(r'$t$ $[M]$')
    ax2.set_ylabel(r'$\|\mathcal{M}\|_{L^2}$')
    ax2.set_title(r'Momentum Constraint: $\mathcal{M}^i$')
    ax2.grid(True, which="both", ls="--", alpha=0.6)
    ax2.legend(loc='upper right', frameon=True, framealpha=0.9)
    ax2.tick_params(axis='both', which='major', direction='in', top=True, right=True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot constraint norms from simulation data.")

    # Default path setup
    # Default path:
    # We assume we are in /home/jovyan/nachevsky/test/simulation/GRTeclyn
    # and data is in ../data/data/constraint_norms.dat
    script_dir = Path(__file__).resolve().parent

    # script_dir is src/visualisation/constraines
    # .parent -> src/visualisation
    # .parent.parent -> src
    # .parent.parent.parent -> GRTeclyn (project root)
    project_root = script_dir.parent.parent.parent

    # Path relative to project root (since user likely runs from project root)
    default_data_path = project_root.parent / "data" / "data" / "constraint_norms.dat"

    # Fallback to current directory if default doesn't exist
    if not default_data_path.exists():
        default_data_path = Path("constraint_norms.dat")

    parser = argparse.ArgumentParser(description="Plot constraint norms from simulation data.")
    parser.add_argument("input_file", nargs="?", default=str(default_data_path),
                        help=f"Path to constraint_norms.dat file (default: {default_data_path})")
    parser.add_argument("-o", "--output", default="constraints_plot.png",
                        help="Output filename (default: constraints_plot.png)")

    args = parser.parse_args()

    input_path = Path(args.input_file)
    output_path = Path(args.output)

    # If output path is just a filename, put it in the script directory
    if output_path.name == str(output_path):
        output_path = script_dir / output_path

    if not input_path.exists():
        print(f"Error: Input file not found at {input_path}")
        print("Please provide the path to constraint_norms.dat")
        sys.exit(1)

    print(f"Reading data from: {input_path}")

    try:
        times, l2_ham, l2_mom = parse_constraint_file(input_path)

        print(f"Data range: t = [{times.min():.4f}, {times.max():.4f}]")
        print(f"Hamiltonian Norm: min={l2_ham.min():.4e}, max={l2_ham.max():.4e}")
        print(f"Momentum Norm:    min={l2_mom.min():.4e}, max={l2_mom.max():.4e}")

        plot_constraints(times, l2_ham, l2_mom, output_path)

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
