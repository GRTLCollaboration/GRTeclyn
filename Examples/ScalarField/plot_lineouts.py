#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Edit these plotting options as needed.
data_path = Path("data")
dt = 0.25
output_file = "lineouts.png"

xlim = (0, 50)
phi_ylim = None  # For example: (-0.1, 0.1)
rho_ylim = None  # For example: (0, 0.01)

phi_files = sorted(data_path.glob("phi_profile_*.dat"))
rho_files = sorted(data_path.glob("rho_profile_*.dat"))

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
colours = plt.cm.viridis(np.linspace(0, 1, len(phi_files)))

for phi_file, rho_file, colour in zip(phi_files, rho_files, colours):
    phi_data = np.loadtxt(phi_file)
    rho_data = np.loadtxt(rho_file)

    # Distance along the diagonal from the first extraction point.
    radius = np.linalg.norm(phi_data[:, :3] - phi_data[0, :3], axis=1)
    step = int(phi_file.stem.split("_")[-1])

    axes[0].plot(radius, phi_data[:, 3], color=colour, label=f"t = {step * dt:g}")
    axes[1].plot(radius, rho_data[:, 3], color=colour)

axes[0].set(xlabel=r"$r$", ylabel=r"$\phi$", xlim=xlim, ylim=phi_ylim)
axes[1].set(xlabel=r"$r$", ylabel=r"$\rho$", xlim=xlim, ylim=rho_ylim)

for axis in axes:
    axis.grid(alpha=0.25)

#axes[0].legend()
fig.subplots_adjust(
    left=0.1, right=0.97, bottom=0.15, top=0.95, wspace=0.3)

fig.savefig(output_file, dpi=200)
plt.show()

