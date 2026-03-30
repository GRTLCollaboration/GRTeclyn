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
    if throat_pinch < 1.0:
        r_safe = np.maximum(r, b0 * throat_pinch + 1e-5)
        z = b0 * np.arccosh(r_safe / (b0 * throat_pinch)) / scale
    else:
        B = b0 * throat_pinch
        r_safe = np.maximum(r, B + 1e-5)
        z = B * np.arccosh(r_safe / B) / scale / throat_pinch
        
    return z + z_shift

def plot_wormhole_stages(output_dir):
    """
    Creates a panel plot showing the stages of wormhole collapse
    """
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 18,
        "axes.titlesize": 22,
    })

    fig = plt.figure(figsize=(20, 12))
    
    b0 = 1.0
    r_max = 5.0
    theta = np.linspace(0, 2*np.pi, 16)
    num_r = 8
    
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
            ax = fig.add_subplot(2, 3, row_idx * 3 + col_idx + 1, projection='3d')
            
            pinch = stage["pinch"]
            is_with_waves = stage.get("with_waves", False)
            scale = 1.0
            
            surf_kwargs = {
                'color': 'white',
                'alpha': 1.0,
                'linewidth': 0.8,
                'edgecolors': 'black',
                'shade': False
            }
            
            r_up = np.linspace(b0 * pinch, r_max, num_r)
            R_up, Theta = np.meshgrid(r_up, theta)
            X_up = R_up * np.cos(Theta)
            Y_up = R_up * np.sin(Theta)
            Z_up = wormhole_embedding(R_up, b0, scale=scale, throat_pinch=pinch)
            
            r_down = np.linspace(b0 * pinch, r_max, num_r)
            R_down, Theta = np.meshgrid(r_down, theta)
            X_down = R_down * np.cos(Theta)
            Y_down = R_down * np.sin(Theta)
            Z_down = -wormhole_embedding(R_down, b0, scale=scale, throat_pinch=pinch)
            
            if pinch <= 0.05:
                Z_up += 1.2
                Z_down -= 1.2
                
                u = np.linspace(0, 2 * np.pi, 20)
                v = np.linspace(0, np.pi, 20)
                hx = 0.3 * np.outer(np.cos(u), np.sin(v))
                hy = 0.3 * np.outer(np.sin(u), np.sin(v))
                hz_up = 0.3 * np.outer(np.ones(np.size(u)), np.cos(v)) + Z_up.min()
                hz_down = 0.3 * np.outer(np.ones(np.size(u)), np.cos(v)) + Z_down.max()
                
                ax.plot_surface(hx, hy, hz_up, color='black', alpha=1.0, shade=False)
                ax.plot_surface(hx, hy, hz_down, color='black', alpha=1.0, shade=False)
                
                if is_with_waves:
                    r_waves = np.linspace(1.5, r_max, 5)
                    Rw, Tw = np.meshgrid(r_waves, theta)
                    Xw = Rw * np.cos(Tw)
                    Yw = Rw * np.sin(Tw)
                    
                    wave_amp = 0.2 * np.exp(-(Rw-3.0)**2 / 1.0) * np.cos(2*Tw)
                    Zw_up = wormhole_embedding(Rw, b0, throat_pinch=pinch) + 1.2 + wave_amp
                    Zw_down = -wormhole_embedding(Rw, b0, throat_pinch=pinch) - 1.2 - wave_amp
                    
                    ax.plot_surface(Xw, Yw, Zw_up, color='white', alpha=1.0, 
                                   linewidth=0.8, edgecolors='black', antialiased=True, shade=False)
                    ax.plot_surface(Xw, Yw, Zw_down, color='white', alpha=1.0, 
                                   linewidth=0.8, edgecolors='black', antialiased=True, shade=False)
                    
                    ax.plot_surface(X_up, Y_up, Z_up, **surf_kwargs)
                    ax.plot_surface(X_down, Y_down, Z_down, **surf_kwargs)
                else:
                    ax.plot_surface(X_up, Y_up, Z_up, **surf_kwargs)
                    ax.plot_surface(X_down, Y_down, Z_down, **surf_kwargs)
                
            else:
                ax.plot_surface(X_up, Y_up, Z_up, **surf_kwargs)
                ax.plot_surface(X_down, Y_down, Z_down, **surf_kwargs)
                
                throat_theta = np.linspace(0, 2*np.pi, 20)
                throat_x = b0 * pinch * np.cos(throat_theta)
                throat_y = b0 * pinch * np.sin(throat_theta)
                throat_z = np.zeros_like(throat_theta)
                ax.plot(throat_x, throat_y, throat_z, color='black', linewidth=2, linestyle='--')

            ax.set_title(stage["title"], pad=15)
            
            ax.set_axis_off()
            
            ax.view_init(elev=20, azim=60)
            
            ax.set_xlim([-r_max, r_max])
            ax.set_ylim([-r_max, r_max])
            ax.set_zlim([-r_max, r_max])

    plt.tight_layout(h_pad=2.0)
    
    os.makedirs(output_dir, exist_ok=True)
    png_path = os.path.join(output_dir, "wormhole_collapse_stages.png")
    pdf_path = os.path.join(output_dir, "wormhole_collapse_stages.pdf")
    
    plt.savefig(png_path, dpi=400, bbox_inches='tight', facecolor='white', transparent=False)
    plt.savefig(pdf_path, dpi=400, bbox_inches='tight', facecolor='white', transparent=False)
    
    print(f"Saved figures to:\n- {png_path}\n- {pdf_path}")

if __name__ == "__main__":
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
    plot_wormhole_stages(output_dir)
