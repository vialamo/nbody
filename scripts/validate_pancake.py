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

def get_analytical_zeldovich(z, config, T_i, z_i, unit_length_mpc):  
    Lx = config.get('domain_size_x', 1.0)
    q = np.linspace(-Lx / 2.0, Lx / 2.0, 2000)
    
    k_code = 2.0 * np.pi / Lx
    k_phys = 2.0 * np.pi / (Lx * unit_length_mpc)
    z_c = 1.0
    
    H_0_kms = 100.0 * config.get('hubble_h', 0.7)
    
    x = q - ((1.0 + z_c) / (1.0 + z)) * (np.sin(k_code * q) / k_code)
    rho_rel = np.abs(1.0 / (1.0 - ((1.0 + z_c) / (1.0 + z)) * np.cos(k_code * q)))
    v_x = -H_0_kms * ((1.0 + z_c) / np.sqrt(1.0 + z)) * (np.sin(k_code * q) / k_phys)
    
    # Temperature scales adiabatically with the density contrast
    T_theo = T_i * ((1.0 + z) / (1.0 + z_i))**2.0 * (rho_rel)**(2.0 / 3.0)
    
    # Filter out the multi-valued region (shell-crossing) after caustic formation
    if z < z_c:
        # 1. Find the Lagrangian coordinate (q) where the flow turns around (dx/dq = 0)
        q_turn = np.arccos((1.0 + z) / (1.0 + z_c)) / k_code
        
        # 2. Map this to the Eulerian boundary (x) of the folded region
        x_turn = q_turn - ((1.0 + z_c) / (1.0 + z)) * (np.sin(k_code * q_turn) / k_code)
        x_caustic = np.abs(x_turn)
        
        # 3. Replace the multi-stream region with NaNs to break the Matplotlib line
        invalid_mask = np.abs(x) < x_caustic
        
        rho_rel[invalid_mask] = np.nan
        v_x[invalid_mask] = np.nan
        T_theo[invalid_mask] = np.nan
        
    return x, rho_rel, v_x, T_theo

def validate_zeldovich_interactive(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5"))) 
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'") 
        return

    with h5py.File(files[0], 'r') as f: 
        config = f['Config'].attrs 
        units = f['Units'].attrs
        
        domain_size = config.get('domain_size_x', 1.0) 
        unit_length_mpc = units.get('unit_length_in_mpc', 91.428)
        unit_vel_kms = units.get('unit_velocity_in_kms', 1.0)
        
        a_i = f['Header'].attrs['scale_factor']
        z_i = (1.0 / a_i) - 1.0
        
        # Read the initial Temperature baseline directly
        T_i = np.mean(f['Gas/temperature'][:])
        rho_0 = np.mean(f['Gas/density'][:])

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
            w_grav = f['Gas'].attrs.get('cumulative_gravitational_work', 0.0) 
            w_exp = f['Gas'].attrs.get('cumulative_expansion_work', 0.0) 
            sw_e = f['Gas'].attrs.get('cumulative_entropy_switch_energy', 0.0)

            if method == "mfm": 
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
                energy_dens = f['Gas/energy'][:] 
                N = f['Config'].attrs.get('mesh_size_1d', 32) 
                vol = (domain_size / N)**3 
                E_tot = np.sum(energy_dens) * vol 
                
            if E_initial is None: 
                E_initial = E_tot 
                
            delta_e = E_tot - E_initial 
            absolute_error = delta_e - w_grav + w_exp
            err = absolute_error / E_initial if E_initial != 0 else 0 
            e_errors.append(err)
            
            sw_err = sw_e / E_initial if E_initial != 0 else 0
            switch_energies.append(sw_err)

    times = np.array(times) 
    e_errors = np.array(e_errors) 
    
    fig = plt.figure(figsize=(18, 12)) 
    gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1.2, 1], figure=fig) 
    plt.subplots_adjust(bottom=0.1, top=0.92, wspace=0.25, hspace=0.30) 
    
    ax_v     = fig.add_subplot(gs[0, 0]) 
    ax_rho   = fig.add_subplot(gs[0, 1]) 
    ax_T     = fig.add_subplot(gs[1, 0]) 
    ax_ts    = fig.add_subplot(gs[1, 1]) 
    
    scatter_kwargs = {'s': 3, 'alpha': 0.6, 'color': 'tab:blue', 'label': 'Simulation'} 
    line_kwargs = {'color': 'black', 'linestyle': '--', 'linewidth': 1.5, 'label': 'Analytic (Linear)'}
    
    scat_v = ax_v.scatter([], [], **scatter_kwargs) 
    line_v, = ax_v.plot([], [], **line_kwargs)
    ax_v.set_ylabel(r"Peculiar Velocity ($v_x$) [km/s]") 
    ax_v.set_xlabel(r"x [$h^{-1}$ Mpc]")
    ax_v.set_title("x-Velocity Profile") 
    ax_v.legend(loc='upper left')
    
    scat_rho = ax_rho.scatter([], [], **scatter_kwargs) 
    line_rho, = ax_rho.plot([], [], **line_kwargs)
    ax_rho.set_ylabel(r"Density Contrast ($\rho/\rho_0$)") 
    ax_rho.set_xlabel(r"x [$h^{-1}$ Mpc]")
    ax_rho.set_yscale('log')
    ax_rho.set_title("Density Profile") 
    ax_rho.legend(loc='upper right')
    
    scat_T = ax_T.scatter([], [], **scatter_kwargs) 
    line_T, = ax_T.plot([], [], **line_kwargs)
    ax_T.set_ylabel(r"Temperature ($K$)") 
    ax_T.set_xlabel(r"x [$h^{-1}$ Mpc]")
    ax_T.set_yscale('log')
    ax_T.set_title("Temperature Profile") 
    ax_T.legend(loc='upper right')

    for ax in [ax_rho, ax_v, ax_T]: 
        ax.grid(True, linestyle=':', alpha=0.6) 
        ax.set_xlim(-domain_size / 2.0, domain_size / 2.0)
    
    ax_ts.plot(times, e_errors, color='tab:red', lw=2, label='Energy Error') 
    ax_ts.plot(times, switch_energies, color='tab:blue', lw=2, linestyle='--', label='Switch Energy')
    ax_ts.set_xlabel("Time (t)") 
    ax_ts.set_ylabel("Fractional Energy Error", color='black') 
    ax_ts.grid(True, linestyle=':', alpha=0.6) 
    ax_ts.legend(loc='upper left', fontsize=10) 
    ax_ts.set_title("Global Energy Diagnostics") 
    vl_ts = ax_ts.axvline(times[0], color='black', linestyle='-', lw=1.5) 

    def update(val): 
        idx = int(snap_slider.val) 
        target_file = files[idx] 
        
        with h5py.File(target_file, 'r') as f: 
            time = f['Header'].attrs['simulation_time'] 
            a_current = f['Header'].attrs['scale_factor']
            z_current = (1.0 / a_current) - 1.0
            vl_ts.set_xdata([time, time]) 
            
            x_full = f['Gas/position_x'][:] 
            rho_flat = f['Gas/density'][:] 
            
            v_x = f['Gas/velocity_x'][:] * a_current * unit_vel_kms
            
            # Read Temperature directly
            T_flat = f['Gas/temperature'][:]
            
        x_plot = (x_full + domain_size / 2.0) % domain_size - (domain_size / 2.0)
        rho_rel = rho_flat / rho_0
        
        max_points = 64000 
        if len(x_plot) > max_points: 
            shuffle_idx = np.random.choice(len(x_plot), max_points, replace=False)
            x_plot = x_plot[shuffle_idx] 
            v_x = v_x[shuffle_idx] 
            rho_rel = rho_rel[shuffle_idx] 
            T_flat = T_flat[shuffle_idx] 
            
        fig.canvas.manager.set_window_title(os.path.basename(os.path.normpath(snapshot_dir))) 
        fig.suptitle(f"Zeldovich Pancake - Snapshot {idx} (z={z_current:.2f})", fontsize=16) 
        
        x_theo, rho_theo, v_theo, T_theo = get_analytical_zeldovich(
            z_current, config, T_i, z_i, unit_length_mpc
        )
        
        line_v.set_data(x_theo, v_theo)
        line_rho.set_data(x_theo, rho_theo)
        line_T.set_data(x_theo, T_theo)

        scat_v.set_offsets(np.c_[x_plot, v_x]) 
        scat_rho.set_offsets(np.c_[x_plot, rho_rel]) 
        scat_T.set_offsets(np.c_[x_plot, T_flat]) 

        v_pad = 0.1 * (np.nanmax(v_x) - np.nanmin(v_x)) if len(v_x) > 0 else 100
        ax_v.set_ylim(np.nanmin(v_x) - v_pad, np.nanmax(v_x) + v_pad)
        
        if len(rho_rel) > 0:
            rho_min = np.nanmin(rho_rel[rho_rel > 0]) if np.any(rho_rel > 0) else 1e-1
            ax_rho.set_ylim(rho_min * 0.5, np.nanmax(rho_rel) * 2.0)
            
            T_min = np.nanmin(T_flat[T_flat > 0]) if np.any(T_flat > 0) else 1e-1
            ax_T.set_ylim(T_min * 0.5, np.nanmax(T_flat) * 2.0)

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
    parser = argparse.ArgumentParser(description="Interactive Zeldovich Pancake validation.") 
    parser.add_argument("path", type=str, help="Path to snapshot directory.") 
    parser.add_argument("-l", "--latest", action="store_true", help="Load latest in directory") 
    
    if len(sys.argv) == 1: 
        parser.print_help() 
        sys.exit(1) 
        
    args = parser.parse_args() 
    target_dir = args.path 
    
    if args.latest: 
        runs = sorted(glob.glob(os.path.join(target_dir, "zeldovich_pancake_*"))) 
        if runs: target_dir = runs[-1] 
        
    validate_zeldovich_interactive(target_dir)
