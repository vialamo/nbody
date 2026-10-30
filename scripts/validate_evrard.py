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

def validate_evrard_interactive(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5"))) 
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'") 
        return

    with h5py.File(files[0], 'r') as f: 
        config = f['Config'].attrs 
        hydro_method = config.get('hydro_method', b"none").decode('utf-8') 
        gamma = config.get('gamma', 5.0/3.0) 
        domain_size = config.get('domain_size', 2.0) 
        center = domain_size / 2.0 

    print("Pre-computing global time-series diagnostics...") 
    times = [] 
    e_errors = [] 
    switch_energies = [] 
    E_initial = None 

    for fname in files: 
        with h5py.File(fname, 'r') as f: 
            t = f['Header'].attrs['simulation_time'] 
            times.append(t) 
            
            method = f['Config'].attrs.get('hydro_method', b'none').decode('utf-8') 

            # Collect energy variables (Evrard involves gravity)
            ph_e = f['Gas'].attrs.get('accumulated_photoheating_energy', 0.0) 
            e_rad = f['Gas'].attrs.get('cumulative_radiated_energy', 0.0) 
            w_grav = f['Gas'].attrs.get('cumulative_gravitational_work', 0.0) 
            w_exp = f['Gas'].attrs.get('cumulative_expansion_work', 0.0) 

            if method == "mfm": 
                sw_e = f['Gas'].attrs.get('cumulative_entropy_switch_energy', 0.0) 
                mass = f['Gas/mass'][:] 
                
                if 'total_energy' in f['Gas']: 
                    u_tot_spec = f['Gas/total_energy'][:] 
                    E_tot = np.sum(mass * u_tot_spec) 
                else: 
                    u = f['Gas/internal_energy'][:] 
                    vx = f['Gas/velocity_x'][:] 
                    vy = f['Gas/velocity_y'][:] 
                    vz = f['Gas/velocity_z'][:] 
                    ke = 0.5 * (vx**2 + vy**2 + vz**2) 
                    E_tot = np.sum(mass * (u + ke)) 
            else: 
                sw_e = f['Gas'].attrs.get('cumulative_dual_energy_switch_energy', 0.0) 
                energy_dens = f['Gas/energy'][:] 
                N = f['Config'].attrs.get('mesh_size_1d', 32) 
                vol = (domain_size / N)**3 
                E_tot = np.sum(energy_dens) * vol 
                # Note: Exact gravitational potential integration varies by Eulerian code implementation
                
            switch_energies.append(sw_e) 
            
            if E_initial is None: 
                E_initial = E_tot 
                
            delta_e = E_tot - E_initial 
            absolute_error = delta_e - w_grav + w_exp + e_rad - ph_e 
            err = absolute_error / E_initial if E_initial != 0 else 0 
            e_errors.append(err) 

    times = np.array(times) 
    e_errors = np.array(e_errors) 
    switch_energies = np.array(switch_energies) 
    
    fig = plt.figure(figsize=(22, 14)) 
    gs = GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1.2], figure=fig) 
    plt.subplots_adjust(bottom=0.1, top=0.92, wspace=0.35, hspace=0.40) 
    
    ax_rho   = fig.add_subplot(gs[0, 0]) 
    ax_v     = fig.add_subplot(gs[0, 1]) 
    ax_S     = fig.add_subplot(gs[0, 2])
    
    ax_u     = fig.add_subplot(gs[1, 0]) 
    ax_gradP = fig.add_subplot(gs[1, 1]) 
    ax_cond  = fig.add_subplot(gs[1, 2]) 
    
    ax_ts    = fig.add_subplot(gs[2, 0]) 
    ax_slice = fig.add_subplot(gs[2, 1:3])  
    
    cmap_angle = 'turbo' 
    scatter_kwargs = {'s': 3, 'alpha': 0.6, 'color': 'tab:blue'} 
    
    scat_rho = ax_rho.scatter([], [], **scatter_kwargs) 
    ax_rho.set_ylabel(r"Density ($\rho$)") 
    ax_rho.set_yscale('log')
    ax_rho.set_title("Radial Density Profile") 
    
    scat_v = ax_v.scatter([], [], **scatter_kwargs) 
    ax_v.set_ylabel(r"Radial Velocity ($v_r$)") 
    ax_v.set_title("Radial Velocity Profile") 
    
    scat_S = ax_S.scatter([], [], **scatter_kwargs)
    ax_S.set_ylabel(r"Entropy ($P/\rho^\gamma$)")
    ax_S.set_title("Entropy Profile")
    
    scat_u = ax_u.scatter([], [], **scatter_kwargs) 
    ax_u.set_ylabel(r"Internal Energy ($u$)") 
    ax_u.set_yscale('log')
    ax_u.set_title("Specific Internal Energy") 

    scat_gradP = ax_gradP.scatter([], [], **scatter_kwargs) 
    ax_gradP.set_ylabel(r"Radial Grad P ($\nabla P \cdot \hat{r}$)") 
    ax_gradP.set_xlabel("Radius (r)") 
    ax_gradP.set_title("Pressure Gradient") 

    scat_cond = ax_cond.scatter([], [], **scatter_kwargs) 
    ax_cond.set_ylabel("Condition Number") 
    ax_cond.set_xlabel("Radius (r)") 
    ax_cond.set_title("Condition Num")

    plot_axes = [ax_rho, ax_v, ax_S, ax_u, ax_gradP, ax_cond]
    for ax in plot_axes: 
        ax.grid(True, linestyle=':', alpha=0.6) 
    
    ax_slice.set_title("Density (2D Slice z=0)")
    ax_slice.set_xlabel("x") 
    ax_slice.set_ylabel("y") 
    ax_slice.set_aspect('equal') 
    ax_slice.set_facecolor('midnightblue') 
    scat_slice = ax_slice.scatter([], [], c=[], cmap='jet', norm=mcolors.LogNorm(vmin=1e-3, vmax=1e3), s=15, marker='o', edgecolors='none') 
    cbar = fig.colorbar(scat_slice, ax=ax_slice, fraction=0.046, pad=0.04) 
    cbar.set_label(r"Density ($\rho$)")
    
    ax_ts.plot(times, e_errors, color='tab:red', lw=2, label='Unexplained Error') 
    ax_ts.plot(times, switch_energies, color='tab:blue', lw=2, linestyle='--', label='Switch Energy') 
    ax_ts.set_xlabel("Time (t)") 
    ax_ts.set_ylabel("Fractional Energy", color='black') 
    ax_ts.grid(True, linestyle=':', alpha=0.6) 
    ax_ts.legend(loc='upper left', fontsize=8) 
    ax_ts.set_title("Global Energy Diagnostics") 
    vl_ts = ax_ts.axvline(times[0], color='black', linestyle='-', lw=1.5) 

    def update(val): 
        idx = int(snap_slider.val) 
        target_file = files[idx] 
        
        with h5py.File(target_file, 'r') as f: 
            time = f['Header'].attrs['simulation_time'] 
            vl_ts.set_xdata([time, time]) 
            
            cfg = f['Config'].attrs 
            method = cfg.get('hydro_method', b'none').decode('utf-8') 
            
            if method == "eulerian": 
                N = cfg.get('mesh_size_1d', 32) 
                dx_cell = domain_size / float(N) 
                coords = np.linspace(dx_cell/2.0, domain_size - dx_cell/2.0, N) 
                X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij') 
                x_full, y_full, z_full = X.flatten(), Y.flatten(), Z.flatten() 
                
                dx = x_full - center 
                dy = y_full - center 
                dz = z_full - center 
                r_flat = np.sqrt(dx**2 + dy**2 + dz**2) 
                
                rho_flat = f['Gas/density'][:].flatten() 
                mom_x_flat = f['Gas/momentum_x'][:].flatten() 
                mom_y_flat = f['Gas/momentum_y'][:].flatten() 
                mom_z_flat = f['Gas/momentum_z'][:].flatten() 
                P_flat = f['Gas/pressure'][:].flatten() 
                
                v_x = mom_x_flat / rho_flat 
                v_y = mom_y_flat / rho_flat 
                v_z = mom_z_flat / rho_flat 
                v_r_flat = (v_x*dx + v_y*dy + v_z*dz) / np.maximum(r_flat, 1e-8) 
                u_flat = P_flat / (rho_flat * (gamma - 1.0)) 
                entropy_flat = P_flat / (rho_flat**gamma)
                
                gradP_mag_flat = np.full_like(r_flat, np.nan) 
                cond_flat = np.full_like(r_flat, np.nan) 
                
                mid_idx = N // 2 
                x_slice = X[:, :, mid_idx].flatten() 
                y_slice = Y[:, :, mid_idx].flatten() 
                rho_slice_2d = rho_flat.reshape((N, N, N))[:, :, mid_idx].flatten()
                
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
                entropy_flat = P_flat / (rho_flat**gamma)
                
                dx = x_full - center 
                dy = y_full - center 
                dz = z_full - center 
                r_flat = np.sqrt(dx**2 + dy**2 + dz**2) 
                v_r_flat = (v_x*dx + v_y*dy + v_z*dz) / np.maximum(r_flat, 1e-8) 
                    
                if 'grad_p' in f['Gas']: 
                    grad_p = f['Gas/grad_p'][:] 
                    gradP_mag_flat = (grad_p[:, 0]*dx + grad_p[:, 1]*dy + grad_p[:, 2]*dz) / np.maximum(r_flat, 1e-8) 
                else: 
                    gradP_mag_flat = np.full_like(rho_flat, np.nan) 
                    
                if 'condition_number' in f['Gas']: 
                    cond_flat = f['Gas/condition_number'][:] 
                else: 
                    cond_flat = np.full_like(rho_flat, np.nan) 

                dz_slice = domain_size * 0.03 
                slice_mask = np.abs(z_full - center) < dz_slice 
                x_slice = x_full[slice_mask] 
                y_slice = y_full[slice_mask] 
                rho_slice_2d = rho_flat[slice_mask]
                
                sort_slice = np.argsort(rho_slice_2d)
                x_slice = x_slice[sort_slice] 
                y_slice = y_slice[sort_slice] 
                rho_slice_2d = rho_slice_2d[sort_slice]

        if len(x_slice) > 0: 
            scat_slice.set_offsets(np.c_[x_slice, y_slice]) 
            scat_slice.set_array(rho_slice_2d)
            u_min = max(rho_slice_2d.min(), 1e-3)
            u_max = max(rho_slice_2d.max(), 1e3)
            scat_slice.set_clim(vmin=u_min, vmax=u_max) 
            
        ax_slice.set_xlim(center - 0.5 * domain_size, center + 0.5 * domain_size) 
        ax_slice.set_ylim(center - 0.5 * domain_size, center + 0.5 * domain_size) 

        max_points = 64**2 
        if len(r_flat) > max_points: 
            sort_idx = np.argsort(r_flat) 
            core_idx = sort_idx[:8000] 
            bg_idx = sort_idx[8000::25] 
            sample_idx = np.concatenate([core_idx, bg_idx]) 
            
            r_plot = r_flat[sample_idx] 
            rho_plot = rho_flat[sample_idx] 
            v_r_plot = v_r_flat[sample_idx] 
            S_plot = entropy_flat[sample_idx]
            u_plot = u_flat[sample_idx] 
            gradP_plot = gradP_mag_flat[sample_idx] 
            cond_plot = cond_flat[sample_idx]  
        else: 
            r_plot, rho_plot, v_r_plot = r_flat, rho_flat, v_r_flat 
            S_plot, u_plot = entropy_flat, u_flat
            gradP_plot = gradP_mag_flat 
            cond_plot = cond_flat 
            
        fig.canvas.manager.set_window_title(os.path.basename(os.path.normpath(snapshot_dir))) 
        fig.suptitle(f"Evrard Collapse Validation - Snapshot {idx} (t={time:.4f})", fontsize=16) 
        
        shuffle_idx = np.random.permutation(len(r_plot)) 

        scat_rho.set_offsets(np.c_[r_plot[shuffle_idx], rho_plot[shuffle_idx]]) 
        scat_v.set_offsets(np.c_[r_plot[shuffle_idx], v_r_plot[shuffle_idx]]) 
        scat_S.set_offsets(np.c_[r_plot[shuffle_idx], S_plot[shuffle_idx]])
        scat_u.set_offsets(np.c_[r_plot[shuffle_idx], u_plot[shuffle_idx]]) 
        scat_gradP.set_offsets(np.c_[r_plot[shuffle_idx], gradP_plot[shuffle_idx]]) 
        scat_cond.set_offsets(np.c_[r_plot[shuffle_idx], cond_plot[shuffle_idx]])

        if len(rho_plot) > 0:
            valid_r = r_plot[r_plot > 0]
            
            # Dynamically calculate the maximum and minimum radius
            r_max = np.nanmax(valid_r) if len(valid_r) > 0 else (domain_size / 2.0)
            r_min = np.nanmin(valid_r) if len(valid_r) > 0 else 1e-3
            
            # Apply dynamic x-axis limits to all radial profile plots
            for ax in [ax_rho, ax_v, ax_S, ax_u, ax_gradP, ax_cond]:
                ax.set_xscale('log')
                ax.set_xlim(r_min * 0.9, r_max * 1.1)
            
            # Dynamic Y-limits (Log plots require >0 filters)
            rho_min = np.nanmin(rho_plot[rho_plot > 0]) if np.any(rho_plot > 0) else 1e-2
            ax_rho.set_ylim(rho_min * 0.5, np.nanmax(rho_plot) * 1.5)
            
            # Linear scale padding
            v_pad = 0.1 * (np.nanmax(v_r_plot) - np.nanmin(v_r_plot))
            v_pad = v_pad if v_pad > 0 else 0.2
            ax_v.set_ylim(np.nanmin(v_r_plot) - v_pad, np.nanmax(v_r_plot) + v_pad)
            
            S_min = np.nanmin(S_plot[S_plot > 0]) if np.any(S_plot > 0) else 1e-3
            ax_S.set_ylim(S_min * 0.5, np.nanmax(S_plot) * 1.5)
            
            u_min = np.nanmin(u_plot[u_plot > 0]) if np.any(u_plot > 0) else 1e-2
            ax_u.set_ylim(u_min * 0.5, np.nanmax(u_plot) * 1.5)
            
        if np.count_nonzero(~np.isnan(gradP_plot)) > 0: 
            ax_gradP.set_ylim(np.nanmin(gradP_plot), np.nanmax(gradP_plot) * 1.1) 
        else: 
            ax_gradP.set_ylim(-0.1, 1.0) 
            
        if np.count_nonzero(~np.isnan(cond_plot)) > 0: 
            ax_cond.set_ylim(-0.1, np.nanmax(cond_plot) * 1.1) 
        else: 
            ax_cond.set_ylim(-0.1, 10.0) 

        fig.canvas.draw_idle() 

    ax_slider = fig.add_axes([0.15, 0.03, 0.7, 0.03]) 
    snap_slider = Slider( 
        ax=ax_slider, 
        label='Snapshot ID', 
        valmin=0, 
        valmax=len(files) - 1, 
        valinit=len(files) - 1, 
        valfmt='%0.0f' 
    ) 
    
    snap_slider.on_changed(update) 
    update(len(files) - 1) 
    plt.show() 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interactive Evrard Collapse validation.") 
    parser.add_argument("path", type=str, help="Path to snapshot directory.") 
    parser.add_argument("-l", "--latest", action="store_true", help="Load latest in directory") 
    
    if len(sys.argv) == 1: 
        parser.print_help() 
        sys.exit(1) 
        
    args = parser.parse_args() 
    target_dir = args.path 
    
    if args.latest: 
        runs = sorted(glob.glob(os.path.join(target_dir, "evrard_collapse_*"))) 
        if runs: target_dir = runs[-1] 
        
    validate_evrard_interactive(target_dir) 
