#!/usr/bin/env python3
"""
Plot constraint convergence using constraint_norms.dat from different runs.

Expected data format (whitespace separated):
  t  L2_Ham  L2_Mom

Usage:
  python convergence.py [options]
"""

import argparse
import sys
from pathlib import Path
from typing import Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from grteclyn_wrapper.core.config import SIM_RESULTS_DIR

matplotlib.use("Agg")

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
    "figure.figsize": (10, 6),
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
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            try:
                float(parts[0])
            except ValueError:
                continue

            if len(parts) < 3:
                continue

            try:
                t = float(parts[0])
                ham = float(parts[1])
                mom = float(parts[2])

                times.append(t)
                l2_ham.append(ham)
                l2_mom.append(mom)
            except ValueError:
                continue

    if not times:
        raise ValueError(f"No valid data found in {filepath}")

    times = np.array(times)
    l2_ham = np.array(l2_ham)
    l2_mom = np.array(l2_mom)

    sort_idx = np.argsort(times)

    return times[sort_idx], l2_ham[sort_idx], l2_mom[sort_idx]


def main():
    parser = argparse.ArgumentParser(description="Plot constraint convergence from simulation data.")

    default_med_path = str(
        SIM_RESULTS_DIR / "Run_R0.5_A00.0_A20.0_sigma0.5_unperturbed/constraint_norms.dat"
    )
    default_high_path = str(
        SIM_RESULTS_DIR / "Run_R0.5_A00.0_A20.02_sigma0.5_perturbed/constraint_norms.dat"
    )

    parser.add_argument("--med", default=default_med_path,
                        help=f"Path to Medium resolution constraint_norms.dat (default: {default_med_path})")
    parser.add_argument("--high", default=default_high_path,
                        help=f"Path to High resolution constraint_norms.dat (default: {default_high_path})")
    parser.add_argument("-o", "--output", default="convergence_plot",
                        help="Output filename base (default: convergence_plot)")
    parser.add_argument("-t", "--time", type=float, default=None,
                        help="Maximum time to plot (e.g., 2 or 2.0)")

    args = parser.parse_args()

    med_path = Path(args.med)
    high_path = Path(args.high)
    output_path = Path(args.output)
    max_time = args.time

    script_dir = Path(__file__).resolve().parent

    if output_path.name == str(output_path):
        output_path = script_dir / output_path

    if not med_path.exists():
        print(f"Error: Medium resolution file not found at {med_path}")
        sys.exit(1)

    if not high_path.exists():
        print(f"Error: High resolution file not found at {high_path}")
        sys.exit(1)

    print(f"Reading Medium resolution data from: {med_path}")
    print(f"Reading High resolution data from: {high_path}")

    try:
        t_med, ham_med, mom_med = parse_constraint_file(med_path)
        t_high, ham_high, mom_high = parse_constraint_file(high_path)

        if max_time is not None:
            mask_med = t_med <= max_time
            t_med = t_med[mask_med]
            ham_med = ham_med[mask_med]
            mom_med = mom_med[mask_med]

            mask_high = t_high <= max_time
            t_high = t_high[mask_high]
            ham_high = ham_high[mask_high]
            mom_high = mom_high[mask_high]
            
            print(f"Filtered data to t <= {max_time}")

        ham_high_scaled_4th = ham_high * 3.16
        mom_high_scaled_4th = mom_high * 3.16
        
        ham_high_scaled_2nd = ham_high * 1.77
        mom_high_scaled_2nd = mom_high * 1.77

        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 8))

        ax1 = axes[0]
        ax1.semilogy(t_med, ham_med, label=r'Medium ($N=192$)', color='black', linestyle='-', linewidth=1.5)
        ax1.semilogy(t_high, ham_high_scaled_4th, label=r'High ($N=256$) $\times 3.16$ (4th order)', color='black', linestyle='--', linewidth=1.5)
        ax1.semilogy(t_high, ham_high_scaled_2nd, label=r'High ($N=256$) $\times 1.77$ (2nd order)', color='black', linestyle='-.', linewidth=1.5)
        ax1.semilogy(t_high, ham_high, label=r'High ($N=256$) (unscaled)', color='black', linestyle=':', linewidth=1.5)
        ax1.set_ylabel(r'$\|\mathcal{H}\|_{L^2}$')
        ax1.set_title(r'Hamiltonian Constraint Convergence: $\mathcal{H}$')
        ax1.grid(True, which="both", ls="--", alpha=0.6)
        ax1.legend(loc='best', frameon=True, framealpha=0.9)
        ax1.tick_params(axis='both', which='major', direction='in', top=True, right=True)

        ax2 = axes[1]
        ax2.semilogy(t_med, mom_med, label=r'Medium ($N=192$)', color='black', linestyle='-', linewidth=1.5)
        ax2.semilogy(t_high, mom_high_scaled_4th, label=r'High ($N=256$) $\times 3.16$ (4th order)', color='black', linestyle='--', linewidth=1.5)
        ax2.semilogy(t_high, mom_high_scaled_2nd, label=r'High ($N=256$) $\times 1.77$ (2nd order)', color='black', linestyle='-.', linewidth=1.5)
        ax2.semilogy(t_high, mom_high, label=r'High ($N=256$) (unscaled)', color='black', linestyle=':', linewidth=1.5)
        ax2.set_xlabel(r'$t$ $[M]$')
        ax2.set_ylabel(r'$\|\mathcal{M}\|_{L^2}$')
        ax2.set_title(r'Momentum Constraint Convergence: $\mathcal{M}^i$')
        ax2.grid(True, which="both", ls="--", alpha=0.6)
        ax2.legend(loc='best', frameon=True, framealpha=0.9)
        ax2.tick_params(axis='both', which='major', direction='in', top=True, right=True)

        plt.tight_layout()
        
        out_png = output_path.with_suffix(".png")
        out_eps = output_path.with_suffix(".eps")
        out_pdf = output_path.with_suffix(".pdf")
        
        plt.savefig(out_png, dpi=600, bbox_inches='tight')
        plt.savefig(out_eps, dpi=600, bbox_inches='tight')
        plt.savefig(out_pdf, dpi=600, bbox_inches='tight')
        
        print(f"Plot saved to: {out_png}, .eps, and .pdf")
        plt.close(fig)

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
