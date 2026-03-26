import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def wormhole_embedding(r, b0, z_shift=0, scale=1.0, throat_pinch=1.0):
    """
    Parametric embedding of a wormhole
    r: radial coordinate
    b0: throat radius
    throat_pinch: factor to simulate throat collapse/expansion (1.0 = open, <1.0 = pinch, >1.0 = expansion)
    """
    # Modify the embedding to show pinch-off or expansion
    if throat_pinch < 1.0:
        # As throat_pinch -> 0, the throat gets narrower and the funnels stretch
        r_safe = np.maximum(r, b0 * throat_pinch + 1e-5)
        z = b0 * np.arccosh(r_safe / (b0 * throat_pinch)) / scale
    else:
        # For expansion, we effectively increase the throat radius to b0 * throat_pinch
        B = b0 * throat_pinch
        r_safe = np.maximum(r, B + 1e-5)
        # To make the asymptotic regions move closer (flatten), we scale down z by throat_pinch
        z = B * np.arccosh(r_safe / B) / scale / throat_pinch
        
    return z + z_shift

def plot_wormhole_stages(output_dir):
    """
    Creates a panel plot showing the stages of wormhole collapse
    """
    # Use standard text but configure fonts nicely, no LaTeX dependency
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 18,        # Increased from 14
        "axes.titlesize": 22,   # Increased from 16
    })
    
    # We will create two rows: one for collapse, one for expansion
    fig = plt.figure(figsize=(20, 12))  # Increased height for two rows
    
    # Parameters
    b0 = 1.0
    r_max = 5.0
    # Reduce theta resolution significantly
    theta = np.linspace(0, 2*np.pi, 16)
    # Reduce radial resolution significantly
    num_r = 8
    
    # Use pure white for the surface, with only black lines
    collapse_stages = [
        {"title": "(a) Initial Traversable Wormhole", "pinch": 1.0, "with_waves": False},
        {"title": "(b) Dynamical Constriction", "pinch": 0.4, "with_waves": False},
        {"title": "(c) Pinch-off & Black Holes", "pinch": 0.05, "with_waves": True} 
    ]
    
    expansion_stages = [
        {"title": "(d) Initial Traversable Wormhole", "pinch": 1.0, "with_waves": False},
        {"title": "(e) Dynamical Expansion", "pinch": 1.5, "with_waves": False},
        {"title": "(f) Expanded Wormhole", "pinch": 2.5, "with_waves": False} 
    ]
    
    all_stages = [collapse_stages, expansion_stages]
    
    for row_idx, row_stages in enumerate(all_stages):
        
        for col_idx, stage in enumerate(row_stages):
            # 2 rows, 3 columns
            ax = fig.add_subplot(2, 3, row_idx * 3 + col_idx + 1, projection='3d')
            
            pinch = stage["pinch"]
            is_with_waves = stage.get("with_waves", False)
            scale = 1.0
            
            # Surface properties: white face color, black edges, opaque
            surf_kwargs = {
                'color': 'white', 
                'alpha': 1.0, 
                'linewidth': 0.8, # Increased linewidth for bold, crisp lines
                'edgecolors': 'black',
                'shade': False # This turns off the 3D lighting shadows completely!
            }
            
            # Normal Wormhole embedding (for all other panels)
            # Upper universe
            r_up = np.linspace(b0 * pinch, r_max, num_r)
            R_up, Theta = np.meshgrid(r_up, theta)
            X_up = R_up * np.cos(Theta)
            Y_up = R_up * np.sin(Theta)
            Z_up = wormhole_embedding(R_up, b0, scale=scale, throat_pinch=pinch)
            
            # Lower universe
            r_down = np.linspace(b0 * pinch, r_max, num_r)
            R_down, Theta = np.meshgrid(r_down, theta)
            X_down = R_down * np.cos(Theta)
            Y_down = R_down * np.sin(Theta)
            Z_down = -wormhole_embedding(R_down, b0, scale=scale, throat_pinch=pinch)
            
            # If it's the pinch-off stage, separate the two universes slightly
            if pinch <= 0.05:
                Z_up += 1.2
                Z_down -= 1.2
                
                # Draw event horizons (black spheres)
                u = np.linspace(0, 2 * np.pi, 20)
                v = np.linspace(0, np.pi, 20)
                hx = 0.3 * np.outer(np.cos(u), np.sin(v))
                hy = 0.3 * np.outer(np.sin(u), np.sin(v))
                hz_up = 0.3 * np.outer(np.ones(np.size(u)), np.cos(v)) + Z_up.min()
                hz_down = 0.3 * np.outer(np.ones(np.size(u)), np.cos(v)) + Z_down.max()
                
                ax.plot_surface(hx, hy, hz_up, color='black', alpha=1.0, shade=False)
                ax.plot_surface(hx, hy, hz_down, color='black', alpha=1.0, shade=False)
                
                if is_with_waves:
                    # Add some gravitational wave ripples for the case with waves
                    r_waves = np.linspace(1.5, r_max, 5) # Reduced wave radial points
                    Rw, Tw = np.meshgrid(r_waves, theta)
                    Xw = Rw * np.cos(Tw)
                    Yw = Rw * np.sin(Tw)
                    
                    # Quadrupolar wave pattern
                    wave_amp = 0.2 * np.exp(-(Rw-3.0)**2 / 1.0) * np.cos(2*Tw)
                    Zw_up = wormhole_embedding(Rw, b0, throat_pinch=pinch) + 1.2 + wave_amp
                    Zw_down = -wormhole_embedding(Rw, b0, throat_pinch=pinch) - 1.2 - wave_amp
                    
                    ax.plot_surface(Xw, Yw, Zw_up, color='white', alpha=1.0, 
                                   linewidth=0.8, edgecolors='black', antialiased=True, shade=False)
                    ax.plot_surface(Xw, Yw, Zw_down, color='white', alpha=1.0, 
                                   linewidth=0.8, edgecolors='black', antialiased=True, shade=False)
                    
                    # Plot the inner part without waves
                    ax.plot_surface(X_up, Y_up, Z_up, **surf_kwargs)
                    ax.plot_surface(X_down, Y_down, Z_down, **surf_kwargs)
                else:
                    # Spontaneous collapse (unperturbed) has no ripples
                    ax.plot_surface(X_up, Y_up, Z_up, **surf_kwargs)
                    ax.plot_surface(X_down, Y_down, Z_down, **surf_kwargs)
                
            else:
                # Plot continuous wormhole (both normal and expanding)
                ax.plot_surface(X_up, Y_up, Z_up, **surf_kwargs)
                ax.plot_surface(X_down, Y_down, Z_down, **surf_kwargs)
                
                # Draw throat indicator if it's not totally flattened
                throat_theta = np.linspace(0, 2*np.pi, 20) # Match azimuthal frequency
                throat_x = b0 * pinch * np.cos(throat_theta)
                throat_y = b0 * pinch * np.sin(throat_theta)
                throat_z = np.zeros_like(throat_theta)
                # Use a dark gray/black dashed line instead of red for B&W format
                ax.plot(throat_x, throat_y, throat_z, color='black', linewidth=2, linestyle='--')

            ax.set_title(stage["title"], pad=15) # Increased pad
            
            # Remove axes for cleaner look
            ax.set_axis_off()
            
            # Set consistent view
            ax.view_init(elev=20, azim=60) # Increased elevation slightly to show the hole better
            
            # Set consistent limits
            ax.set_xlim([-r_max, r_max])
            ax.set_ylim([-r_max, r_max])
            ax.set_zlim([-r_max, r_max])

    # Adjust layout to reduce vertical whitespace between rows and leave room for titles
    plt.tight_layout(h_pad=2.0)
    
    # Save the figure
    os.makedirs(output_dir, exist_ok=True)
    png_path = os.path.join(output_dir, "wormhole_collapse_stages.png")
    pdf_path = os.path.join(output_dir, "wormhole_collapse_stages.pdf")
    
    # Add facecolor='white' so it doesn't get messed up by transparent backgrounds in LaTeX
    plt.savefig(png_path, dpi=400, bbox_inches='tight', facecolor='white', transparent=False)
    plt.savefig(pdf_path, dpi=400, bbox_inches='tight', facecolor='white', transparent=False)
    
    print(f"Saved figures to:\n- {png_path}\n- {pdf_path}")

if __name__ == "__main__":
    output_dir = "/home/jovyan/nachevsky/test/simulation/GRTeclyn/src/visualisation/figures"
    plot_wormhole_stages(output_dir)
