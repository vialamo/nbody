#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def get_exact_soundwave_solution(x, t, domain_size_x, gamma, rho_0_emp, P_0_emp):
    """
    Computes the exact analytical solution for the Linear Soundwave test.
    Uses strict theoretical constants for wave speed/phase to prevent phase drift,
    but centers the wave amplitude on the empirical baseline to ignore static kernel bias.
    """
    # Theoretical constants dictate the physics of the wave propagation
    rho_theo = 1.0
    c_s_theo = 2.0 / 3.0
    
    delta_rho = 1e-4
    k = 2.0 * np.pi / domain_size_x
    
    # Left-propagating wave phase using the exact theoretical speed
    phase = k * (x + c_s_theo * t)
    sin_phase = np.sin(phase)
    
    # Wave centered on the empirical kernel baseline
    rho = rho_0_emp + delta_rho * sin_phase
    
    # Velocity perturbation is absolute, scaled by the theoretical density
    v_x = -c_s_theo * (delta_rho / rho_theo) * sin_phase
    
    P = P_0_emp + (c_s_theo**2) * delta_rho * sin_phase
    u = P / ((gamma - 1.0) * rho)
    
    return rho, v_x, P, u

def validate_soundwave_interactive(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    # Extract configuration and empirical baselines from the very first snapshot
    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        hydro_method = config.get('hydro_method', b"none").decode('utf-8')
        gamma = config.get('gamma', 5.0/3.0)
        
        # Extract rectangular domain sizes, falling back to cubic 'domain_size' if absent
        domain_base = config.get('domain_size', 1.0)
        domain_size_x = config.get('domain_size_x', domain_base)
        domain_size_y = config.get('domain_size_y', domain_base)
        domain_size_z = config.get('domain_size_z', domain_base)

        gas = f['Gas']
        if hydro_method == "mfm":
            rho_sim_init = gas['density'][:]
            P_sim_init = gas['pressure'][:]
        elif hydro_method == "eulerian":
            N_base = config.get('mesh_size_1d', 32)
            N_x = config.get('mesh_size_x', int(round(N_base * (domain_size_x / domain_base))))
            N_y = config.get('mesh_size_y', int(round(N_base * (domain_size_y / domain_base))))
            N_z = config.get('mesh_size_z', int(round(N_base * (domain_size_z / domain_base))))
            
            rho_3d = gas['density'][:].reshape((N_x, N_y, N_z))
            P_3d = gas['pressure'][:].reshape((N_x, N_y, N_z))
            
            mid_y, mid_z = N_y // 2, N_z // 2
            rho_sim_init = rho_3d[:, mid_y, mid_z]
            P_sim_init = P_3d[:, mid_y, mid_z]
        else:
            rho_sim_init, P_sim_init = np.array([1.0]), np.array([1.0])

        rho_0_emp = np.mean(rho_sim_init)
        P_0_emp = np.mean(P_sim_init)
        
    print(f"Empirical Baseline Detected: rho_0 = {rho_0_emp:.6f}, P_0 = {P_0_emp:.6f}")
    print(f"Domain Bounds: [{domain_size_x:.2f}, {domain_size_y:.2f}, {domain_size_z:.2f}]")
    print("Pre-calculating global L1 errors from all snapshots...")
    
    times = []
    l1_errors = []

    for f_name in files:
        with h5py.File(f_name, 'r') as f:
            t = f['Header'].attrs['simulation_time']
            times.append(t)
            
            if 'Gas' in f:
                gas = f['Gas']
                
                if hydro_method == "mfm":
                    x = gas['position_x'][:]
                    rho_sim = gas['density'][:]
                    mass = gas['mass'][:]
                elif hydro_method == "eulerian":
                    N_base = config.get('mesh_size_1d', 32)
                    N_x = config.get('mesh_size_x', int(round(N_base * (domain_size_x / domain_base))))
                    N_y = config.get('mesh_size_y', int(round(N_base * (domain_size_y / domain_base))))
                    N_z = config.get('mesh_size_z', int(round(N_base * (domain_size_z / domain_base))))
                    
                    dx = domain_size_x / float(N_x)
                    x = np.linspace(dx/2.0, domain_size_x - dx/2.0, N_x)
                    
                    rho_3d = gas['density'][:].reshape((N_x, N_y, N_z))
                    mid_y, mid_z = N_y // 2, N_z // 2
                    rho_sim = rho_3d[:, mid_y, mid_z]
                    mass = None
                else:
                    x, rho_sim, mass = [], [], []

                if len(x) > 0:
                    rho_exact = get_exact_soundwave_solution(x, t, domain_size_x, gamma, rho_0_emp, P_0_emp)[0]
                    if hydro_method == "mfm":
                        # Volume-weighted continuous L1 norm 
                        volume = mass / rho_sim
                        l1_err = np.sum(np.abs(rho_sim - rho_exact) * volume) / np.sum(volume)
                    else:
                        # Uniform cell volume L1 norm
                        l1_err = np.mean(np.abs(rho_sim - rho_exact))
                    l1_errors.append(l1_err)
                else:
                    l1_errors.append(0.0)
            else:
                l1_errors.append(0.0)

    label_prefix = "Eulerian Grid" if hydro_method == "eulerian" else "MFM Particles"
    color_prefix = "darkorange" if hydro_method == "eulerian" else "royalblue"
    
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.4])
    plt.subplots_adjust(bottom=0.10, hspace=0.3) 

    ax_rho = fig.add_subplot(gs[0, 0])
    ax_v   = fig.add_subplot(gs[0, 1])
    ax_P   = fig.add_subplot(gs[1, 0])
    ax_u   = fig.add_subplot(gs[1, 1])
    ax_err = fig.add_subplot(gs[2, :])
    
    spatial_axes = [ax_rho, ax_v, ax_P, ax_u]

    # Theoretical wave scaling for dynamic plot limits
    c_s_theo = 2.0 / 3.0
    u_0_emp = P_0_emp / ((gamma - 1.0) * rho_0_emp)
    
    delta_rho = 1e-4
    delta_v = c_s_theo * (delta_rho / 1.0)
    delta_P = (c_s_theo**2) * delta_rho
    delta_u = delta_P / ((gamma - 1.0) * rho_0_emp)

    ax_rho.set_ylim(rho_0_emp - 1.5*delta_rho, rho_0_emp + 1.5*delta_rho)
    ax_v.set_ylim(-1.5*delta_v, 1.5*delta_v)  
    ax_P.set_ylim(P_0_emp - 1.5*delta_P, P_0_emp + 1.5*delta_P) 
    ax_u.set_ylim(u_0_emp - 1.5*delta_u, u_0_emp + 1.5*delta_u)
    
    scatter_kwargs = {'s': 6 if hydro_method == "eulerian" else 4, 
                      'color': color_prefix, 'alpha': 0.8, 'label': f'{label_prefix} Sim'}
    exact_line_kwargs = {'color': 'black', 'lw': 1.5, 'linestyle': '--', 'label': 'Exact Solution'}
    
    scat_rho = ax_rho.scatter([], [], **scatter_kwargs)
    line_rho, = ax_rho.plot([], [], **exact_line_kwargs)
    ax_rho.set_ylabel(r"Density ($\rho$)")
    ax_rho.set_title("Density Profile")
    
    scat_v = ax_v.scatter([], [], **scatter_kwargs)
    line_v, = ax_v.plot([], [], **exact_line_kwargs)
    ax_v.set_ylabel(r"Velocity ($v_x$)")
    ax_v.set_title("Velocity Profile")
    
    scat_P = ax_P.scatter([], [], **scatter_kwargs)
    line_P, = ax_P.plot([], [], **exact_line_kwargs)
    ax_P.set_ylabel(r"Pressure ($P$)")
    ax_P.set_title("Pressure Profile")
    
    scat_u = ax_u.scatter([], [], **scatter_kwargs)
    line_u, = ax_u.plot([], [], **exact_line_kwargs)
    ax_u.set_ylabel(r"Internal Energy ($u$)")
    ax_u.set_title("Specific Internal Energy")

    ax_err.plot(times, l1_errors, label=r'$L_1$ Density Error', color='black', lw=1.5)
    vline_time = ax_err.axvline(0, color='gray', linestyle='--', alpha=0.7)
    if len(times) > 1:
        ax_err.set_xlim(min(times), max(times))
    ax_err.set_ylabel("L1 Error Norm")
    ax_err.set_xlabel("Time")
    ax_err.set_title("Density Convergence")
    ax_err.set_yscale('log')
    ax_err.grid(True, linestyle=':', alpha=0.6)
    ax_err.legend(loc='lower right', fontsize=9)

    # Use the final snapshot's error for the label
    final_l1_error = l1_errors[-1] if l1_errors else 0.0
    ax_err.text(0.05, 0.85, f"Final $L_1$ Error: {final_l1_error:.5e}", 
                transform=ax_err.transAxes,
                ha='left', va='center', fontsize=12, fontweight='bold',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5', alpha=0.8))

    for ax in spatial_axes:
        ax.set_xlim(0, domain_size_x)
        ax.grid(True, linestyle=':', alpha=0.6)
            
    ax_rho.legend(loc='upper right', fontsize=9)

    def update(val):
        idx = int(snap_slider.val)
        target_file = files[idx]
        
        with h5py.File(target_file, 'r') as f:
            time = f['Header'].attrs['simulation_time']
            cfg = f['Config'].attrs
            method = cfg.get('hydro_method', b'none').decode('utf-8')
            
            if method == "eulerian":
                N_base = cfg.get('mesh_size_1d', 32)
                N_x = cfg.get('mesh_size_x', int(round(N_base * (domain_size_x / domain_base))))
                N_y = cfg.get('mesh_size_y', int(round(N_base * (domain_size_y / domain_base))))
                N_z = cfg.get('mesh_size_z', int(round(N_base * (domain_size_z / domain_base))))
                
                dx = domain_size_x / float(N_x)
                x_v = np.linspace(dx/2.0, domain_size_x - dx/2.0, N_x)
                
                rho_3d = f['Gas/density'][:].reshape((N_x, N_y, N_z))
                mom_x_3d = f['Gas/momentum_x'][:].reshape((N_x, N_y, N_z))
                P_3d = f['Gas/pressure'][:].reshape((N_x, N_y, N_z))
                
                mid_y, mid_z = N_y // 2, N_z // 2
                rho_v = rho_3d[:, mid_y, mid_z]
                v_v = mom_x_3d[:, mid_y, mid_z] / rho_v
                P_v = P_3d[:, mid_y, mid_z]
                u_v = P_v / (rho_v * (gamma - 1.0))
                
            elif method == "mfm":
                x_v = f['Gas/position_x'][:]
                rho_v = f['Gas/density'][:]
                v_v = f['Gas/velocity_x'][:]
                u_v = f['Gas/internal_energy'][:]
                P_v = rho_v * u_v * (gamma - 1.0)
            
        fig.canvas.manager.set_window_title(os.path.basename(os.path.normpath(snapshot_dir)))
        fig.suptitle(f"Linear Soundwave Validation - Snapshot {idx} (t={time:.4f})", fontsize=16)
        
        sort_idx = np.argsort(x_v)
        x_v = x_v[sort_idx]
        rho_v = rho_v[sort_idx]
        v_v = v_v[sort_idx]
        P_v = P_v[sort_idx]
        u_v = u_v[sort_idx]
        
        scat_rho.set_offsets(np.c_[x_v, rho_v])
        scat_v.set_offsets(np.c_[x_v, v_v])
        scat_P.set_offsets(np.c_[x_v, P_v])
        scat_u.set_offsets(np.c_[x_v, u_v])

        x_exact = np.linspace(0, domain_size_x, 1000)
        rho_ex, v_ex, P_ex, u_ex = get_exact_soundwave_solution(x_exact, time, domain_size_x, gamma, rho_0_emp, P_0_emp)
        
        line_rho.set_data(x_exact, rho_ex)
        line_v.set_data(x_exact, v_ex)
        line_P.set_data(x_exact, P_ex)
        line_u.set_data(x_exact, u_ex)

        vline_time.set_xdata([time, time])
        fig.canvas.draw_idle()

    ax_slider = fig.add_axes([0.15, 0.02, 0.7, 0.03])
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
    parser = argparse.ArgumentParser(description="Interactive MFM/Eulerian Soundwave validation.")
    parser.add_argument("path", type=str, nargs='?', default=".", help="Path to snapshot directory.")
    parser.add_argument("-l", "--latest", action="store_true", help="Load latest run_* directory")
    
    args = parser.parse_args()
    target_dir = args.path
    
    if args.latest:
        runs = sorted(glob.glob(os.path.join(target_dir, "run_*")))
        if runs: target_dir = runs[-1]
        
    validate_soundwave_interactive(target_dir)
