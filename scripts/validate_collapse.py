#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider

def validate_dm_collapse_interactive(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5"))) 
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'") 
        return

    with h5py.File(files[0], 'r') as f: 
        config = f['Config'].attrs 
        domain_size = config.get('domain_size_x', 2.0) 
        center = domain_size / 2.0 

    print("Pre-computing global time-series diagnostics...") 
    times = [] 
    e_errors = [] 
    K_initial = None 

    for fname in files: 
        with h5py.File(fname, 'r') as f: 
            t = f['Header'].attrs['simulation_time'] 
            times.append(t) 
            
            # Collect energy variables for DM
            w_grav = f['Particles'].attrs.get('cumulative_gravitational_work', 0.0) 
            w_exp = f['Particles'].attrs.get('cumulative_expansion_work', 0.0) 

            mass = f['Particles/mass'][:] 
            vx = f['Particles/velocity_x'][:] 
            vy = f['Particles/velocity_y'][:] 
            vz = f['Particles/velocity_z'][:] 
            
            K_current = np.sum(0.5 * mass * (vx**2 + vy**2 + vz**2))
            
            if K_initial is None: 
                K_initial = K_current 
                
            delta_K = K_current - K_initial 
            absolute_error = delta_K - w_grav + w_exp 
            
            # Match the C++ diagnostics logic for division by zero if starting at rest
            denominator = K_initial if K_initial > 1e-12 else K_current
            err = absolute_error / denominator if denominator > 1e-12 else 0.0
            e_errors.append(err) 

    times = np.array(times) 
    e_errors = np.array(e_errors) 
    
    fig = plt.figure(figsize=(16, 12)) 
    gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1], figure=fig) 
    plt.subplots_adjust(bottom=0.15, top=0.92, wspace=0.25, hspace=0.30) 
    
    ax_rho   = fig.add_subplot(gs[0, 0]) 
    ax_vr    = fig.add_subplot(gs[0, 1]) 
    ax_ts    = fig.add_subplot(gs[1, 0]) 
    ax_slice = fig.add_subplot(gs[1, 1])  
    
    # 1. Radial Density Profile Setup
    line_rho, = ax_rho.plot([], [], color='tab:blue', lw=2)
    ax_rho.set_ylabel(r"Density ($\rho$)") 
    ax_rho.set_xlabel("Radius (r)")
    ax_rho.set_yscale('log')
    ax_rho.set_xscale('log')
    ax_rho.set_xlim(1e-4, 10.0) # Add default positive bound
    ax_rho.set_ylim(1e-4, 10.0) # Add default positive bound
    ax_rho.set_title("Radial Density Profile (Binned)") 
    ax_rho.grid(True, linestyle=':', alpha=0.6)
    
    # 2. Phase Space Setup
    scatter_kwargs = {'s': 2, 'alpha': 0.5, 'color': 'tab:orange'} 
    scat_vr = ax_vr.scatter([], [], **scatter_kwargs) 
    ax_vr.set_ylabel(r"Radial Velocity ($v_r$)") 
    ax_vr.set_xlabel("Radius (r)")
    ax_vr.set_xscale('log')
    ax_vr.set_xlim(1e-4, 10.0) # Add default positive bound
    ax_vr.set_title("Phase Space: Radial Velocity") 
    ax_vr.grid(True, linestyle=':', alpha=0.6)
    
    # 3. Time Series Setup
    ax_ts.plot(times, e_errors, color='tab:red', lw=2, label='Energy Error') 
    ax_ts.set_xlabel("Time (t)") 
    ax_ts.set_ylabel("Fractional Energy Error", color='black') 
    ax_ts.grid(True, linestyle=':', alpha=0.6) 
    ax_ts.legend(loc='upper left', fontsize=10) 
    ax_ts.set_title("Global Energy Diagnostics") 
    vl_ts = ax_ts.axvline(times[0], color='black', linestyle='-', lw=1.5) 

    # 4. 2D Slice Setup
    ax_slice.set_title("Particle Distribution (2D Slice z≈0)")
    ax_slice.set_xlabel("x") 
    ax_slice.set_ylabel("y") 
    ax_slice.set_aspect('equal') 
    ax_slice.set_facecolor('black') 
    scat_slice = ax_slice.scatter([], [], s=1, color='cyan', alpha=0.3, edgecolors='none') 
    
    def update(val): 
        idx = int(snap_slider.val) 
        target_file = files[idx] 
        
        with h5py.File(target_file, 'r') as f: 
            time = f['Header'].attrs['simulation_time'] 
            vl_ts.set_xdata([time, time]) 
            
            x_full = f['Particles/position_x'][:] 
            y_full = f['Particles/position_y'][:] 
            z_full = f['Particles/position_z'][:] 
            v_x = f['Particles/velocity_x'][:] 
            v_y = f['Particles/velocity_y'][:] 
            v_z = f['Particles/velocity_z'][:] 
            mass = f['Particles/mass'][:]
            
        dx = x_full - center 
        dy = y_full - center 
        dz = z_full - center 
        r_flat = np.sqrt(dx**2 + dy**2 + dz**2) 
        v_r_flat = (v_x*dx + v_y*dy + v_z*dz) / np.maximum(r_flat, 1e-8) 

        # Slice data
        dz_slice = domain_size * 0.05 
        slice_mask = np.abs(z_full - center) < dz_slice 
        x_slice = x_full[slice_mask] 
        y_slice = y_full[slice_mask] 
            
        if len(x_slice) > 0: 
            scat_slice.set_offsets(np.c_[x_slice, y_slice]) 
            
        ax_slice.set_xlim(center - 0.5 * domain_size, center + 0.5 * domain_size) 
        ax_slice.set_ylim(center - 0.5 * domain_size, center + 0.5 * domain_size) 

        # Downsample phase space if there are too many particles
        max_points = 100000 
        if len(r_flat) > max_points: 
            shuffle_idx = np.random.choice(len(r_flat), max_points, replace=False)
        else:
            shuffle_idx = np.random.permutation(len(r_flat))

        r_plot = r_flat[shuffle_idx] 
        v_r_plot = v_r_flat[shuffle_idx] 
            
        fig.canvas.manager.set_window_title(os.path.basename(os.path.normpath(snapshot_dir))) 
        fig.suptitle(f"Collisionless Spherical Collapse - Snapshot {idx} (t={time:.4f})", fontsize=16) 
        
        # Update Phase Space (Filter out r=0 for the log scale)
        valid_phase = r_plot > 0
        scat_vr.set_offsets(np.c_[r_plot[valid_phase], v_r_plot[valid_phase]])
        
        valid_r = r_flat[r_flat > 0]
        if len(valid_r) > 0:
            r_max = np.nanmax(valid_r)
            r_min = np.nanmin(valid_r)
            
            # Adaptive limits for phase space
            ax_vr.set_xlim(max(r_min, 1e-4) * 0.9, r_max * 1.1)
            v_pad = 0.1 * (np.nanmax(v_r_plot) - np.nanmin(v_r_plot))
            v_pad = v_pad if v_pad > 0 else 0.2
            ax_vr.set_ylim(np.nanmin(v_r_plot) - v_pad, np.nanmax(v_r_plot) + v_pad)

            # Compute Radial Density Profile (Binning)
            bins = np.logspace(np.log10(max(r_min, 1e-3)), np.log10(r_max), 50)
            hist_mass, _ = np.histogram(r_flat, bins=bins, weights=mass)
            bin_volumes = (4.0 / 3.0) * np.pi * (bins[1:]**3 - bins[:-1]**3)
            rho_bins = hist_mass / bin_volumes
            bin_centers = np.sqrt(bins[1:] * bins[:-1])
            
            # Filter empty bins
            valid_bins = rho_bins > 0
            
            if np.any(valid_bins):
                line_rho.set_data(bin_centers[valid_bins], rho_bins[valid_bins])
                ax_rho.set_xlim(bins[0] * 0.9, bins[-1] * 1.1)
                ax_rho.set_ylim(np.min(rho_bins[valid_bins]) * 0.5, np.max(rho_bins[valid_bins]) * 1.5)

        fig.canvas.draw_idle() 

    ax_slider = fig.add_axes([0.15, 0.05, 0.7, 0.03]) 
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
    parser = argparse.ArgumentParser(description="Interactive DM Collapse validation.") 
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
        
    validate_dm_collapse_interactive(target_dir)
