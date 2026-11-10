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

def validate_kh_interactive(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5"))) 
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'") 
        return

    with h5py.File(files[0], 'r') as f: 
        config = f['Config'].attrs 
        gamma = config.get('gamma', 5.0/3.0) 
        domain_size = config.get('domain_size', 1.0) 
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

            ph_e = f['Gas'].attrs.get('accumulated_photoheating_energy', 0.0) 
            e_rad = f['Gas'].attrs.get('cumulative_radiated_energy', 0.0) 
            w_grav = f['Gas'].attrs.get('cumulative_gravitational_work', 0.0) 
            w_exp = f['Gas'].attrs.get('cumulative_expansion_work', 0.0) 

            
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
    
    # 2x2 Layout for 2D Slices + 1 Row for Time Series
    fig = plt.figure(figsize=(18, 14)) 
    gs = GridSpec(3, 2, height_ratios=[1, 1, 0.6], figure=fig) 
    plt.subplots_adjust(bottom=0.1, top=0.92, wspace=0.2, hspace=0.3) 
    
    ax_rho = fig.add_subplot(gs[0, 0])
    ax_u = fig.add_subplot(gs[0, 1])
    ax_vx = fig.add_subplot(gs[1, 0])
    ax_vy = fig.add_subplot(gs[1, 1])
    ax_ts = fig.add_subplot(gs[2, :]) 
    
    slice_axes = [ax_rho, ax_u, ax_vx, ax_vy]
    for ax in slice_axes:
        ax.set_aspect('equal') 
        ax.set_facecolor('midnightblue') 
        ax.set_xlim(0, domain_size)
        ax.set_ylim(0, domain_size)

    ax_rho.set_title("Density (2D Slice)")
    ax_u.set_title("Internal Energy (2D Slice)")
    ax_vx.set_title("Velocity X (2D Slice)")
    ax_vy.set_title("Velocity Y (2D Slice)")

    scatter_kwargs = {'s': 3, 'marker': 's', 'edgecolors': 'none'} 
    scat_rho = ax_rho.scatter([], [], c=[], cmap='viridis', **scatter_kwargs)
    scat_u = ax_u.scatter([], [], c=[], cmap='plasma', **scatter_kwargs)
    scat_vx = ax_vx.scatter([], [], c=[], cmap='coolwarm', **scatter_kwargs)
    scat_vy = ax_vy.scatter([], [], c=[], cmap='coolwarm', **scatter_kwargs)

    fig.colorbar(scat_rho, ax=ax_rho, fraction=0.046, pad=0.04) 
    fig.colorbar(scat_u, ax=ax_u, fraction=0.046, pad=0.04) 
    fig.colorbar(scat_vx, ax=ax_vx, fraction=0.046, pad=0.04) 
    fig.colorbar(scat_vy, ax=ax_vy, fraction=0.046, pad=0.04) 

    # Setup Time Series (Energy)[cite: 62]
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
            
            x_full = f['Gas/position_x'][:] 
            y_full = f['Gas/position_y'][:] 
            z_full = f['Gas/position_z'][:] 
            rho_full = f['Gas/density'][:] 
            vx_full = f['Gas/velocity_x'][:] 
            vy_full = f['Gas/velocity_y'][:] 
            u_full = f['Gas/internal_energy'][:] 

            dz_slice = domain_size * 0.05 
            slice_mask = np.abs(z_full - center) < dz_slice 
            
            x_slice = x_full[slice_mask] 
            y_slice = y_full[slice_mask] 
            rho_slice = rho_full[slice_mask]
            vx_slice = vx_full[slice_mask]
            vy_slice = vy_full[slice_mask]
            u_slice = u_full[slice_mask] 

        if len(x_slice) > 0:
            coords = np.c_[x_slice, y_slice] 
            
            scat_rho.set_offsets(coords)
            scat_rho.set_array(rho_slice)
            scat_rho.set_clim(vmin=rho_slice.min(), vmax=rho_slice.max())
            
            scat_u.set_offsets(coords)
            scat_u.set_array(u_slice)
            scat_u.set_clim(vmin=u_slice.min(), vmax=u_slice.max())

            scat_vx.set_offsets(coords)
            scat_vx.set_array(vx_slice)
            scat_vx.set_clim(vmin=vx_slice.min(), vmax=vx_slice.max())

            scat_vy.set_offsets(coords)
            scat_vy.set_array(vy_slice)
            scat_vy.set_clim(vmin=vy_slice.min(), vmax=vy_slice.max())

        fig.canvas.manager.set_window_title(os.path.basename(os.path.normpath(snapshot_dir)))
        fig.suptitle(f"Kelvin-Helmholtz Validation - Snapshot {idx} (t={time:.4f})", fontsize=16) 
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
    parser = argparse.ArgumentParser(description="Interactive KH Instability validation.") 
    parser.add_argument("path", type=str, help="Path to snapshot directory.") 
    parser.add_argument("-l", "--latest", action="store_true", help="Load latest run_* directory") 
    
    if len(sys.argv) == 1: 
        parser.print_help() 
        sys.exit(1) 
        
    args = parser.parse_args() 
    target_dir = args.path 
    
    if args.latest: 
        runs = sorted(glob.glob(os.path.join(target_dir, "kelvin_helmholtz_*"))) 
        if runs: target_dir = runs[-1] 
        
    validate_kh_interactive(target_dir) 
