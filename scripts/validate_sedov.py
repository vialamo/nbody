#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider
from scipy.integrate import solve_ivp

def get_sedov_exact_profile(t, E=1.0, rho_bg=1.0, gamma=5.0/3.0, n_points=500):
    """
    Computes the exact analytical profile for the 3D Sedov-Taylor blastwave 
    by numerically integrating the self-similar ODEs.
    """
    if abs(gamma - (5.0/3.0)) < 1e-3:
        alpha = 1.15167
    elif abs(gamma - 1.4) < 1e-3:
        alpha = 1.033
    else:
        alpha = 1.0  # Fallback approximation
        
    if t <= 1e-8:
        return np.array([0]), np.array([rho_bg]), np.array([0]), np.array([0]), np.array([0])
        
    Rs = alpha * (E * (t**2) / rho_bg)**0.2
    D = 0.4 * Rs / t
    
    # Boundary conditions at the shock front (xi = r/Rs = 1)
    U1 = 2.0 / (gamma + 1.0)
    G1 = (gamma + 1.0) / (gamma - 1.0)
    Pi1 = 2.0 / (gamma + 1.0)
    
    def sedov_odes(xi, y):
        U, G, Pi = y
        
        # Prevent division by zero near the origin
        if xi < 1e-6 or G < 1e-6:
            return [0, 0, 0]
            
        C2 = gamma * Pi / G
        
        num = 3.0 * Pi / G - 1.5 * U * (U - xi) - 2.0 * C2 * U / xi
        den = C2 - (U - xi)**2
        
        if abs(den) < 1e-10:
            dU = 0
        else:
            dU = num / den
            
        dG = -(G / (U - xi)) * (dU + 2.0 * U / xi)
        dPi = G * (1.5 * U - (U - xi) * dU)
        
        return [dU, dG, dPi]
        
    # Integrate backwards from the shock front (xi=1) towards the center
    # We stop at xi=1e-3 because the density drops to zero and creates a singularity
    res = solve_ivp(sedov_odes, [1.0, 1e-3], [U1, G1, Pi1], 
                    t_eval=np.linspace(1.0, 1e-3, n_points),
                    method='RK45', rtol=1e-6, atol=1e-8)
                    
    # Reverse arrays so they go from center (r=0) to shock front (r=Rs)
    xi = res.t[::-1]
    U = res.y[0][::-1]
    G = res.y[1][::-1]
    Pi = res.y[2][::-1]
    
    # Convert dimensionless similarity variables back to physical units
    r_exact = np.concatenate(([0.0], xi * Rs))
    v_exact = np.concatenate(([0.0], U * D))
    rho_exact = np.concatenate(([0.0], G * rho_bg))
    
    # Pressure tends to a constant non-zero value at the center
    P_exact = np.concatenate(([Pi[0] * rho_bg * D**2], Pi * rho_bg * D**2)) 
    
    # Internal energy (u = P / rho(gamma-1)) approaches infinity at the center as rho -> 0
    u_exact = np.zeros_like(P_exact)
    valid = rho_exact > 1e-8
    u_exact[valid] = P_exact[valid] / (rho_exact[valid] * (gamma - 1.0))
    u_exact[~valid] = np.nan # Mask the central singularity
    
    return r_exact, rho_exact, v_exact, P_exact, u_exact

def validate_sedov_interactive(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        setup_type = config.get('setup', b"").decode('utf-8')
        hydro_method = config.get('hydro_method', b"none").decode('utf-8')
        
        if setup_type != "sedov_blastwave":
            print(f"[ERROR] Expected config attribute 'setup'='sedov_blastwave'. Found: '{setup_type}'")
            sys.exit(1)
            
        gamma = config.get('gamma', 5.0/3.0)
        domain_size = config.get('domain_size', 1.0)
        center = domain_size / 2.0

        if abs(gamma - (5.0/3.0)) > 1e-5:
            print(f"[WARNING] Expected gamma = 5/3 for the standard Sedov test, but found {gamma:.5f}.")
            
        # Physical constraints for Sedov
        if bool(config.get('expanding_universe', 0)):
            print("[ERROR] Sedov test requires a static box, but 'expanding_universe' is TRUE.")
            sys.exit(1)
            
        if bool(config.get('enable_cooling', 0)):
            print("[ERROR] Sedov test requires adiabatic gas, but 'enable_cooling' is TRUE.")
            sys.exit(1)
            
        if bool(config.get('use_pm', 0)) or bool(config.get('use_pp', 0)):
            print(f"[ERROR] Sedov test must be pure hydrodynamics.")
            sys.exit(1)

    # Styling based on solver method
    label_prefix = "Eulerian Grid" if hydro_method == "eulerian" else "MFM Particles"
    color_prefix = "darkorange" if hydro_method == "eulerian" else "royalblue"
    
    # Plot Setup using GridSpec for the 2x3 layout
    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(2, 3, width_ratios=[1, 1, 1.2], figure=fig)
    plt.subplots_adjust(bottom=0.15, wspace=0.3)
    
    axs = np.empty((2, 2), dtype=object)
    axs[0, 0] = fig.add_subplot(gs[0, 0])
    axs[0, 1] = fig.add_subplot(gs[0, 1])
    axs[1, 0] = fig.add_subplot(gs[1, 0])
    axs[1, 1] = fig.add_subplot(gs[1, 1])
    
    # 2D Slice panel spanning both rows on the right
    ax_slice = fig.add_subplot(gs[:, 2])
    
    scatter_kwargs = {'s': 2, 'color': color_prefix, 'alpha': 0.3, 'label': f'{label_prefix} Sim'}
    exact_line_kwargs = {'color': 'black', 'lw': 2.0, 'linestyle': '--', 'label': 'Exact Solution'}
    
    # Density
    scat_rho = axs[0, 0].scatter([], [], **scatter_kwargs)
    line_rho, = axs[0, 0].plot([], [], **exact_line_kwargs)
    axs[0, 0].set_ylabel(r"Density ($\rho$)")
    axs[0, 0].set_title("Radial Density Profile")
    
    # Radial Velocity
    scat_v = axs[0, 1].scatter([], [], **scatter_kwargs)
    line_v, = axs[0, 1].plot([], [], **exact_line_kwargs)
    axs[0, 1].set_ylabel(r"Radial Velocity ($v_r$)")
    axs[0, 1].set_title("Radial Velocity Profile")
    
    # Pressure
    scat_P = axs[1, 0].scatter([], [], **scatter_kwargs)
    line_P, = axs[1, 0].plot([], [], **exact_line_kwargs)
    axs[1, 0].set_ylabel(r"Pressure ($P$)")
    axs[1, 0].set_xlabel("Radius (r)")
    axs[1, 0].set_title("Pressure Profile")
    
    # Internal Energy
    scat_u = axs[1, 1].scatter([], [], **scatter_kwargs)
    line_u, = axs[1, 1].plot([], [], **exact_line_kwargs)
    axs[1, 1].set_ylabel(r"Internal Energy ($u$)")
    axs[1, 1].set_xlabel("Radius (r)")
    axs[1, 1].set_title("Specific Internal Energy")

    for ax in axs.flat:
        ax.set_xlim(0, domain_size / 2.0)
        ax.grid(True, linestyle=':', alpha=0.6)
    axs[0, 0].legend(loc='upper right', fontsize=9)
    
    # Setup 2D Slice
    ax_slice.set_title("Specific Internal Energy (2D Slice z=0)")
    ax_slice.set_xlabel("x")
    ax_slice.set_ylabel("y")
    ax_slice.set_aspect('equal')
    ax_slice.set_facecolor('midnightblue') # Give empty space a dark background
    scat_slice = ax_slice.scatter([], [], c=[], cmap='jet', norm=mcolors.LogNorm(vmin=1, vmax=1000), s=15, marker='s', edgecolors='none')
    cbar = fig.colorbar(scat_slice, ax=ax_slice, fraction=0.046, pad=0.04)
    cbar.set_label(r"Internal Energy ($u$)")

    def update(val):
        idx = int(snap_slider.val)
        target_file = files[idx]
        
        with h5py.File(target_file, 'r') as f:
            time = f['Header'].attrs['simulation_time']
            cfg = f['Config'].attrs
            method = cfg.get('hydro_method', b'none').decode('utf-8')
            
            if method == "eulerian":
                N = cfg.get('mesh_size_1d', 32)
                dx = domain_size / float(N)
                
                coords = np.linspace(dx/2.0, domain_size - dx/2.0, N)
                X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')
                
                x_full, y_full, z_full = X.flatten(), Y.flatten(), Z.flatten()
                r_flat = np.sqrt((x_full - center)**2 + (y_full - center)**2 + (z_full - center)**2)
                
                rho_flat = f['Gas/density'][:].flatten()
                mom_x_flat = f['Gas/momentum_x'][:].flatten()
                mom_y_flat = f['Gas/momentum_y'][:].flatten()
                mom_z_flat = f['Gas/momentum_z'][:].flatten()
                P_flat = f['Gas/pressure'][:].flatten()
                
                v_x = mom_x_flat / rho_flat
                v_y = mom_y_flat / rho_flat
                v_z = mom_z_flat / rho_flat
                
                v_r_flat = (v_x*(x_full-center) + v_y*(y_full-center) + v_z*(z_full-center)) / np.maximum(r_flat, 1e-8)
                u_flat = P_flat / (rho_flat * (gamma - 1.0))
                
                # Eulerian slice extraction
                mid_idx = N // 2
                x_slice = X[:, :, mid_idx].flatten()
                y_slice = Y[:, :, mid_idx].flatten()
                u_slice_2d = u_flat.reshape((N, N, N))[:, :, mid_idx].flatten()
                
            elif method == "mfm":
                x_full = f['Gas/position_x'][:]
                y_full = f['Gas/position_y'][:]
                z_full = f['Gas/position_z'][:]
                rho_flat = f['Gas/density'][:]
                v_x = f['Gas/velocity_x'][:]
                v_y = f['Gas/velocity_y'][:]
                v_z = f['Gas/velocity_z'][:]
                u_flat = f['Gas/internal_energy'][:]
                P_flat = rho_flat * u_flat * (gamma - 1.0)
                
                r_flat = np.sqrt((x_full - center)**2 + (y_full - center)**2 + (z_full - center)**2)
                v_r_flat = (v_x*(x_full-center) + v_y*(y_full-center) + v_z*(z_full-center)) / np.maximum(r_flat, 1e-8)

                # MFM slice extraction (isolate particles near z=center)
                dz = domain_size * 0.03
                slice_mask = np.abs(z_full - center) < dz
                x_slice = x_full[slice_mask]
                y_slice = y_full[slice_mask]
                u_slice_2d = u_flat[slice_mask]
                
                # Sort slice to render high internal energy particles on top
                sort_slice = np.argsort(u_slice_2d)
                x_slice = x_slice[sort_slice]
                y_slice = y_slice[sort_slice]
                u_slice_2d = u_slice_2d[sort_slice]

        # Update 2D Slice panel
        if len(x_slice) > 0:
            scat_slice.set_offsets(np.c_[x_slice, y_slice])
            scat_slice.set_array(u_slice_2d)
            u_min = max(u_slice_2d.min(), 1e-1)
            u_max = max(u_slice_2d.max(), 10.0)
            scat_slice.set_clim(vmin=u_min, vmax=u_max)
            
        ax_slice.set_xlim(center - 0.5 * domain_size, center + 0.5 * domain_size)
        ax_slice.set_ylim(center - 0.5 * domain_size, center + 0.5 * domain_size)

        # Sub-sampling (For 1D Profiles only)
        max_points = 15000
        if len(r_flat) > max_points:
            sort_idx = np.argsort(r_flat)
            core_idx = sort_idx[:8000]
            bg_idx = sort_idx[8000::25] 
            sample_idx = np.concatenate([core_idx, bg_idx])
            
            r_plot = r_flat[sample_idx]
            rho_plot = rho_flat[sample_idx]
            v_r_plot = v_r_flat[sample_idx]
            P_plot = P_flat[sample_idx]
            u_plot = u_flat[sample_idx]
        else:
            r_plot, rho_plot, v_r_plot, P_plot, u_plot = r_flat, rho_flat, v_r_flat, P_flat, u_flat
            
        fig.suptitle(f"Sedov Blastwave Validation - Snapshot {idx} (t={time:.4f})", fontsize=16)
        
        # Update Scatter Points
        scat_rho.set_offsets(np.c_[r_plot, rho_plot])
        scat_v.set_offsets(np.c_[r_plot, v_r_plot])
        scat_P.set_offsets(np.c_[r_plot, P_plot])
        scat_u.set_offsets(np.c_[r_plot, u_plot])

        # Get Analytical Shock Profile
        r_ex, rho_ex, v_ex, P_ex, u_ex = get_sedov_exact_profile(time, E=1.0, rho_bg=1.0, gamma=gamma)
        
        # Update Analytical Lines
        line_rho.set_data(r_ex, rho_ex)
        line_v.set_data(r_ex, v_ex)
        line_P.set_data(r_ex, P_ex)
        line_u.set_data(r_ex, u_ex)

        # Dynamic Y-axis limits
        if len(r_ex) > 1 and r_ex[-1] > 0:
            axs[0, 0].set_ylim(-0.2, rho_ex[-1] * 1.3)
            axs[0, 1].set_ylim(-0.1, v_ex[-1] * 1.3)
            axs[1, 0].set_ylim(-0.1, P_ex[-1] * 1.3)
            u_shock_peak = u_ex[-1]
            axs[1, 1].set_ylim(-0.1, u_shock_peak * 5.0)
        else:
            axs[0, 0].set_ylim(-0.2, 5.0)
            axs[0, 1].set_ylim(-0.1, 1.0)
            axs[1, 0].set_ylim(-0.1, 1.0)
            axs[1, 1].set_ylim(-0.1, 1.0)

        fig.canvas.draw_idle()

    # Slider
    ax_slider = fig.add_axes([0.15, 0.05, 0.7, 0.03])
    snap_slider = Slider(
        ax=ax_slider,
        label='Snapshot ID',
        valmin=0,
        valmax=len(files) - 1,
        valinit=0,
        valfmt='%0.0f'
    )
    
    snap_slider.on_changed(update)
    update(0)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interactive MFM/Eulerian Sedov Blastwave validation.")
    parser.add_argument("path", type=str, help="Path to snapshot directory.")
    parser.add_argument("-l", "--latest", action="store_true", help="Load latest run_* directory")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    target_dir = args.path
    
    if args.latest:
        runs = sorted(glob.glob(os.path.join(target_dir, "run_*")))
        if runs: target_dir = runs[-1]
        
    validate_sedov_interactive(target_dir)