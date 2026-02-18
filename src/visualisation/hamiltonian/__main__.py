import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import argparse

def main():
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _project_root = os.path.dirname(os.path.dirname(os.path.dirname(_script_dir)))
    _default_data = os.path.abspath(os.path.join(_project_root, "..", "data"))

    parser = argparse.ArgumentParser(description="Plot Hamiltonian and Momentum Constraint Violations")
    parser.add_argument("--data", default=_default_data, help="Directory containing WormholePlt* / plt*")
    parser.add_argument("--out", default=_script_dir, help="Output directory for constraints_plot.png")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.data, "WormholePlt*")))
    if not files:
        print("No plotfiles found!")
        return

    times = []
    ham_l2 = []
    mom_l2 = []

    print(f"Found {len(files)} files. Computing constraint norms...")

    for f in files:
        try:
            ds = yt.load(f)
            t = float(ds.current_time)

            # Load all data (this might be memory intensive for huge grids,
            # consider iterating over grids if needed)
            # ad = ds.all_data()

            # Check if fields exist
            if ('boxlib', 'Ham') not in ds.field_list:
                print(f"Skipping {f}: 'Ham' field not found")
                continue

            # Use chunking to avoid loading all data at once
            integral_ham2 = 0.0
            total_vol = 0.0
            integral_mom2 = 0.0

            has_mom = all(('boxlib', f'Mom{i}') in ds.field_list for i in [1,2,3])

            # Iterate over grids to save memory
            for grid in ds.index.grids:
                vol = grid['boxlib', 'cell_volume']
                ham = grid['boxlib', 'Ham']

                integral_ham2 += np.sum(ham**2 * vol)
                total_vol += np.sum(vol)

                if has_mom:
                    mom1 = grid['boxlib', 'Mom1']
                    mom2 = grid['boxlib', 'Mom2']
                    mom3 = grid['boxlib', 'Mom3']
                    mom_sq = mom1**2 + mom2**2 + mom3**2
                    integral_mom2 += np.sum(mom_sq * vol)

            l2_ham = np.sqrt(integral_ham2 / total_vol)
            l2_mom = np.sqrt(integral_mom2 / total_vol) if has_mom else np.nan

            times.append(t)
            ham_l2.append(l2_ham)
            mom_l2.append(l2_mom)

            print(f"Processed t={t:.2f}: Ham_L2={l2_ham:.2e}")

        except Exception as e:
            print(f"Skipping {f}: {e}")

    # Plotting
    if not times:
        print("No data to plot.")
        return

    # Convert to numpy arrays for sorting
    times = np.array(times)
    ham_l2 = np.array(ham_l2)
    mom_l2 = np.array(mom_l2)

    # Sort by time
    idx = np.argsort(times)
    times = times[idx]
    ham_l2 = ham_l2[idx]
    mom_l2 = mom_l2[idx]

    # LaTeX scientific paper style
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "stix",
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.linewidth": 1.0,
        "grid.alpha": 0.5,
    })
    nrows = 2 if not np.all(np.isnan(mom_l2)) else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(10, 8), sharex=(nrows == 2))
    if nrows == 1:
        axes = [axes]

    # Hamiltonian Plot
    ax1 = axes[0]
    ax1.semilogy(times, ham_l2, label=r'$\|\mathcal{H}\|_{L^2}$', color='blue', linewidth=1.5)
    ax1.set_ylabel(r'$\|\mathcal{H}\|_{L^2}$')
    ax1.set_title(r'Hamiltonian constraint: $\mathcal{H}$')
    ax1.grid(True, which="both", ls="--", alpha=0.6)
    ax1.legend(loc='upper right')
    ax1.tick_params(axis='both', which='major', direction='in')
    if nrows == 1:
        ax1.set_xlabel(r'$t$')

    # Momentum Plot
    if nrows == 2:
        ax2 = axes[1]
        ax2.semilogy(times, mom_l2, label=r'$\|\mathcal{M}\|_{L^2}$', color='red', linewidth=1.5)
        ax2.set_xlabel(r'$t$')
        ax2.set_ylabel(r'$\|\mathcal{M}\|_{L^2}$')
        ax2.set_title(r'Momentum constraint: $\mathcal{M}^i$')
        ax2.grid(True, which="both", ls="--", alpha=0.6)
        ax2.legend(loc='upper right')
        ax2.tick_params(axis='both', which='major', direction='in')

    plt.tight_layout()
    out_file = os.path.join(args.out, "constraints_plot.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {out_file}")

if __name__ == "__main__":
    main()
