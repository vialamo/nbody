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
        return np.array([0]), np.array([rho_bg]), np.array([0]), np.array([0]), np.array([0]), np.array([0])
        
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

    # Entropy S = P / rho^gamma
    S_exact = np.zeros_like(P_exact)
    S_exact[valid] = P_exact[valid] / (rho_exact[valid]**gamma)
    S_exact[~valid] = np.nan
    
    return r_exact, rho_exact, v_exact, P_exact, u_exact, S_exact

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

    # Styling based on solver method
    label_prefix = "Eulerian Grid" if hydro_method == "eulerian" else "MFM Particles"
    color_prefix = "darkorange" if hydro_method == "eulerian" else "royalblue"
    
    # Plot Setup using GridSpec for the 3x4 layout
    fig = plt.figure(figsize=(24, 14))
    gs = GridSpec(3, 4, width_ratios=[1, 1, 1, 1], height_ratios=[1, 1, 1.2], figure=fig)
    plt.subplots_adjust(bottom=0.1, top=0.92, wspace=0.35, hspace=0.35)
    
    # Create the axes
    ax_rho   = fig.add_subplot(gs[0, 0])
    ax_v     = fig.add_subplot(gs[0, 1])
    ax_h     = fig.add_subplot(gs[0, 2])
    ax_S     = fig.add_subplot(gs[0, 3])
    
    ax_P     = fig.add_subplot(gs[1, 0])
    ax_u     = fig.add_subplot(gs[1, 1])
    ax_gradP = fig.add_subplot(gs[1, 2])
    ax_cond  = fig.add_subplot(gs[1, 3])
    ax_raw   = ax_cond.twinx()  # Twin axis for Raw Sum
    
    ax_slice = fig.add_subplot(gs[2, 1:3])  # Slice spans the two middle columns
    
    scatter_kwargs = {'s': 2, 'color': color_prefix, 'alpha': 0.3, 'label': f'{label_prefix} Sim'}
    exact_line_kwargs = {'color': 'black', 'lw': 2.0, 'linestyle': '--', 'label': 'Exact Solution'}
    
    # Density
    scat_rho = ax_rho.scatter([], [], **scatter_kwargs)
    line_rho, = ax_rho.plot([], [], **exact_line_kwargs)
    ax_rho.set_ylabel(r"Density ($\rho$)")
    ax_rho.set_title("Radial Density Profile")
    
    # Radial Velocity
    scat_v = ax_v.scatter([], [], **scatter_kwargs)
    line_v, = ax_v.plot([], [], **exact_line_kwargs)
    ax_v.set_ylabel(r"Radial Velocity ($v_r$)")
    ax_v.set_title("Radial Velocity Profile")
    
    # Smoothing Length
    scat_h = ax_h.scatter([], [], **scatter_kwargs)
    ax_h.set_ylabel(r"Smoothing Length ($h$)")
    ax_h.set_title("Smoothing Length Profile")

    # Pressure
    scat_P = ax_P.scatter([], [], **scatter_kwargs)
    line_P, = ax_P.plot([], [], **exact_line_kwargs)
    ax_P.set_ylabel(r"Pressure ($P$)")
    ax_P.set_title("Pressure Profile")
    
    # Internal Energy
    scat_u = ax_u.scatter([], [], **scatter_kwargs)
    line_u, = ax_u.plot([], [], **exact_line_kwargs)
    ax_u.set_ylabel(r"Internal Energy ($u$)")
    ax_u.set_title("Specific Internal Energy")

    # Entropy
    scat_S = ax_S.scatter([], [], **scatter_kwargs)
    line_S, = ax_S.plot([], [], **exact_line_kwargs)
    ax_S.set_ylabel(r"Entropy ($P / \rho^\gamma$)")
    ax_S.set_title("Entropy Profile")

    # Pressure Gradient Magnitude
    scat_gradP = ax_gradP.scatter([], [], **scatter_kwargs)
    ax_gradP.set_ylabel(r"Grad P Mag ($|\nabla P|$)")
    ax_gradP.set_xlabel("Radius (r)")
    ax_gradP.set_title("Pressure Gradient Magnitude")

    # Condition Number & Raw Sum
    scat_cond = ax_cond.scatter([], [], s=2, color='tab:red', alpha=0.5, label='Condition Num')
    scat_raw  = ax_raw.scatter([], [], s=2, color='tab:purple', alpha=0.5, label='|Raw Sum P|')
    ax_cond.set_ylabel("Condition Number", color='tab:red')
    ax_raw.set_ylabel(r"Raw Sum Magnitude ($|\Sigma \Delta P \cdot \mathbf{x} W V|$)", color='tab:purple')
    ax_cond.set_xlabel("Radius (r)")
    ax_cond.set_title("Matrix Condition & Raw Sum")
    ax_cond.tick_params(axis='y', labelcolor='tab:red')
    ax_raw.tick_params(axis='y', labelcolor='tab:purple')

    # Global axis limits styling
    plot_axes = [ax_rho, ax_v, ax_h, ax_P, ax_u, ax_S, ax_gradP, ax_cond]
    for ax in plot_axes:
        ax.set_xlim(0, domain_size / 2.0)
        ax.grid(True, linestyle=':', alpha=0.6)
    
    ax_rho.legend(loc='upper right', fontsize=9)
    
    # Setup 2D Slice
    ax_slice.set_title("Specific Internal Energy (2D Slice z=0)")
    ax_slice.set_xlabel("x")
    ax_slice.set_ylabel("y")
    ax_slice.set_aspect('equal')
    ax_slice.set_facecolor('midnightblue')
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
                S_flat = P_flat / (rho_flat**gamma)
                
                h_flat = np.full_like(r_flat, np.nan) 
                gradP_mag_flat = np.full_like(r_flat, np.nan)
                cond_flat = np.full_like(r_flat, np.nan)
                raw_mag_flat = np.full_like(r_flat, np.nan)
                
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

                h_flat = f['Gas/smoothing_length'][:] if 'smoothing_length' in f['Gas'] else np.full_like(rho_flat, np.nan)
                S_flat = f['Gas/entropy'][:] if 'entropy' in f['Gas'] else P_flat / (rho_flat**gamma)
                    
                if 'grad_p' in f['Gas']:
                    grad_p = f['Gas/grad_p'][:]
                    gradP_mag_flat = np.linalg.norm(grad_p, axis=1)
                else:
                    gradP_mag_flat = np.full_like(rho_flat, np.nan)
                    
                if 'condition_number' in f['Gas']:
                    cond_flat = f['Gas/condition_number'][:]
                else:
                    cond_flat = np.full_like(rho_flat, np.nan)

                if 'raw_sum_p' in f['Gas']:
                    raw_sum = f['Gas/raw_sum_p'][:]
                    raw_mag_flat = np.linalg.norm(raw_sum, axis=1)
                else:
                    raw_mag_flat = np.full_like(rho_flat, np.nan)

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
            h_plot = h_flat[sample_idx]
            S_plot = S_flat[sample_idx]
            gradP_plot = gradP_mag_flat[sample_idx]
            cond_plot = cond_flat[sample_idx]
            raw_plot = raw_mag_flat[sample_idx]
        else:
            r_plot, rho_plot, v_r_plot = r_flat, rho_flat, v_r_flat
            P_plot, u_plot, h_plot = P_flat, u_flat, h_flat
            S_plot, gradP_plot = S_flat, gradP_mag_flat
            cond_plot, raw_plot = cond_flat, raw_mag_flat
            
        fig.suptitle(f"Sedov Blastwave Validation - Snapshot {idx} (t={time:.4f})", fontsize=16)
        
        # Update Scatter Points
        scat_rho.set_offsets(np.c_[r_plot, rho_plot])
        scat_v.set_offsets(np.c_[r_plot, v_r_plot])
        scat_P.set_offsets(np.c_[r_plot, P_plot])
        scat_u.set_offsets(np.c_[r_plot, u_plot])
        scat_h.set_offsets(np.c_[r_plot, h_plot])
        scat_S.set_offsets(np.c_[r_plot, S_plot])
        scat_gradP.set_offsets(np.c_[r_plot, gradP_plot])
        scat_cond.set_offsets(np.c_[r_plot, cond_plot])
        scat_raw.set_offsets(np.c_[r_plot, raw_plot])

        # Get Analytical Shock Profile
        r_ex, rho_ex, v_ex, P_ex, u_ex, S_ex = get_sedov_exact_profile(time, E=1.0, rho_bg=1.0, gamma=gamma)
        
        # Update Analytical Lines
        line_rho.set_data(r_ex, rho_ex)
        line_v.set_data(r_ex, v_ex)
        line_P.set_data(r_ex, P_ex)
        line_u.set_data(r_ex, u_ex)
        line_S.set_data(r_ex, S_ex)

        # Dynamic Y-axis limits
        if len(r_ex) > 1 and r_ex[-1] > 0:
            ax_rho.set_ylim(-0.2, rho_ex[-1] * 1.3)
            ax_v.set_ylim(-0.1, v_ex[-1] * 1.3)
            ax_P.set_ylim(-0.1, P_ex[-1] * 1.3)
            u_shock_peak = u_ex[-1]
            ax_u.set_ylim(-0.1, u_shock_peak * 5.0)
        else:
            ax_rho.set_ylim(-0.2, 5.0)
            ax_v.set_ylim(-0.1, 1.0)
            ax_P.set_ylim(-0.1, 1.0)
            ax_u.set_ylim(-0.1, 1.0)
            
        # Optional attributes limit bounds (Percentiles avoid singularities)
        if np.count_nonzero(~np.isnan(h_plot)) > 0:
            ax_h.set_ylim(0, np.nanmax(h_plot) * 1.2)
        else:
            ax_h.set_ylim(0, 0.1)

        if np.count_nonzero(~np.isnan(S_plot)) > 0:
            s_max = np.nanpercentile(S_plot, 98)
            ax_S.set_ylim(-0.1, s_max * 1.5)
        else:
            ax_S.set_ylim(-0.1, 1.0)

        if np.count_nonzero(~np.isnan(gradP_plot)) > 0:
            ax_gradP.set_ylim(-0.1, np.nanmax(gradP_plot) * 1.1)
        else:
            ax_gradP.set_ylim(-0.1, 1.0)
            
        # Condition number and Raw Sum bounds
        if np.count_nonzero(~np.isnan(cond_plot)) > 0:
            ax_cond.set_ylim(-0.1, np.nanmax(cond_plot) * 1.1)
        else:
            ax_cond.set_ylim(-0.1, 10.0)

        if np.count_nonzero(~np.isnan(raw_plot)) > 0:
            ax_raw.set_ylim(-0.1, np.nanmax(raw_plot) * 1.1)
        else:
            ax_raw.set_ylim(-0.1, 10.0)

        fig.canvas.draw_idle()

    # Slider
    ax_slider = fig.add_axes([0.15, 0.03, 0.7, 0.03])
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