#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def get_exact_soundwave_solution(x, t, domain_size, gamma, rho_0, P_0):
    """
    Computes the analytical solution for the Linear Soundwave test,
    centered around an empirical numerical baseline (rho_0, P_0).
    """
    delta_rho = 1e-4
    c_s = np.sqrt(gamma * P_0 / rho_0)
    k = 2.0 * np.pi / domain_size
    
    # Left-propagating wave phase
    phase = k * (x + c_s * t)
    sin_phase = np.sin(phase)
    
    rho = rho_0 + delta_rho * sin_phase
    
    # Velocity perturbation for a left-moving wave is negative
    v_x = -c_s * (delta_rho / rho_0) * sin_phase
    
    P = P_0 + (c_s**2) * delta_rho * sin_phase
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
        domain_size = config.get('domain_size', 1.0)

        gas = f['Gas']
        if hydro_method == "mfm":
            rho_sim_init = gas['density'][:]
            P_sim_init = gas['pressure'][:]
        elif hydro_method == "eulerian":
            N = config.get('mesh_size_1d', 32)
            rho_3d = gas['density'][:].reshape((N, N, N))
            P_3d = gas['pressure'][:].reshape((N, N, N))
            mid = N // 2
            rho_sim_init = rho_3d[:, mid, mid]
            P_sim_init = P_3d[:, mid, mid]
        else:
            rho_sim_init, P_sim_init = np.array([1.0]), np.array([1.0])

        rho_0_emp = np.mean(rho_sim_init)
        P_0_emp = np.mean(P_sim_init)
        
    print(f"Empirical Baseline Detected: rho_0 = {rho_0_emp:.6f}, P_0 = {P_0_emp:.6f}")
    print("Pre-calculating global conservation diagnostics from all snapshots...")
    
    times = []
    l1_errors = []
    vol_errors = []
    energy_errors = []
    mom_errors = []
    
    E0_total = None
    P0_total = None

    for f_name in files:
        with h5py.File(f_name, 'r') as f:
            t = f['Header'].attrs['simulation_time']
            times.append(t)
            
            if 'Gas' in f:
                gas = f['Gas']
                
                # Extract arrays
                if hydro_method == "mfm":
                    x = gas['position_x'][:]
                    rho_sim = gas['density'][:]
                    v_x = gas['velocity_x'][:]
                    u_ie = gas['internal_energy'][:]
                    mass = gas['mass'][:]
                elif hydro_method == "eulerian":
                    N = config.get('mesh_size_1d', 32)
                    dx = domain_size / float(N)
                    x = np.linspace(dx/2.0, domain_size - dx/2.0, N)
                    
                    rho_3d = gas['density'][:].reshape((N, N, N))
                    v_x_3d = (gas['momentum_x'][:] / gas['density'][:]).reshape((N, N, N))
                    P_3d = gas['pressure'][:].reshape((N, N, N))
                    
                    mid = N // 2
                    rho_sim = rho_3d[:, mid, mid]
                    v_x = v_x_3d[:, mid, mid]
                    u_ie = P_3d[:, mid, mid] / (rho_sim * (gamma - 1.0))
                    mass = gas['density'][:] * (dx**3)
                    
                    # For global conservation we need the full flat arrays
                    v_x_full = gas['momentum_x'][:] / gas['density'][:]
                    u_ie_full = gas['pressure'][:] / (gas['density'][:] * (gamma - 1.0))
                else:
                    x, rho_sim, v_x, u_ie, mass = [], [], [], [], []

                # 1. L1 Density Error
                if len(x) > 0:
                    rho_exact = get_exact_soundwave_solution(x, t, domain_size, gamma, rho_0_emp, P_0_emp)[0]
                    l1_err = np.mean(np.abs(rho_sim - rho_exact))
                    l1_errors.append(l1_err)
                else:
                    l1_errors.append(0.0)
                    
                # 2. Volume Partition Error
                box_volume = domain_size**3
                if hydro_method == "mfm":
                    total_sim_volume = np.sum(mass / gas['density'][:])
                else:
                    total_sim_volume = np.sum(mass / gas['density'][:])
                vol_errors.append((total_sim_volume / box_volume) - 1.0)
                
                # 3. Total Energy Error
                if hydro_method == "mfm":
                    total_energy = np.sum(mass * (gas['internal_energy'][:] + 0.5 * gas['velocity_x'][:]**2))
                else:
                    total_energy = np.sum(mass * (u_ie_full + 0.5 * v_x_full**2))
                    
                if E0_total is None:
                    E0_total = total_energy
                energy_errors.append((total_energy - E0_total) / E0_total if E0_total else 0.0)
                
                # 4. Total Momentum Error
                if hydro_method == "mfm":
                    total_momentum = np.sum(mass * gas['velocity_x'][:])
                else:
                    total_momentum = np.sum(mass * v_x_full)
                    
                if P0_total is None:
                    P0_total = total_momentum
                mom_errors.append(total_momentum - P0_total)

            else:
                l1_errors.append(0.0)
                vol_errors.append(0.0)
                energy_errors.append(0.0)
                mom_errors.append(0.0)

    label_prefix = "Eulerian Grid" if hydro_method == "eulerian" else "MFM Particles"
    color_prefix = "darkorange" if hydro_method == "eulerian" else "royalblue"
    
    fig, axs = plt.subplots(3, 2, figsize=(14, 14))
    plt.subplots_adjust(bottom=0.10, hspace=0.3) 

    # Dynamic limits based on empirical baseline and 1e-4 perturbation
    c_s_emp = np.sqrt(gamma * P_0_emp / rho_0_emp)
    u_0_emp = P_0_emp / ((gamma - 1.0) * rho_0_emp)
    
    delta_rho = 1e-4
    delta_v = c_s_emp * (delta_rho / rho_0_emp)
    delta_P = (c_s_emp**2) * delta_rho
    delta_u = delta_P / ((gamma - 1.0) * rho_0_emp)

    axs[0, 0].set_ylim(rho_0_emp - 1.5*delta_rho, rho_0_emp + 1.5*delta_rho)
    axs[0, 1].set_ylim(-1.5*delta_v, 1.5*delta_v)  
    axs[1, 0].set_ylim(P_0_emp - 1.5*delta_P, P_0_emp + 1.5*delta_P) 
    axs[1, 1].set_ylim(u_0_emp - 1.5*delta_u, u_0_emp + 1.5*delta_u)
    
    scatter_kwargs = {'s': 6 if hydro_method == "eulerian" else 4, 
                      'color': color_prefix, 'alpha': 0.8, 'label': f'{label_prefix} Sim'}
    exact_line_kwargs = {'color': 'black', 'lw': 1.5, 'linestyle': '--', 'label': 'Exact Solution'}
    
    scat_rho = axs[0, 0].scatter([], [], **scatter_kwargs)
    line_rho, = axs[0, 0].plot([], [], **exact_line_kwargs)
    axs[0, 0].set_ylabel(r"Density ($\rho$)")
    axs[0, 0].set_title("Density Profile")
    
    scat_v = axs[0, 1].scatter([], [], **scatter_kwargs)
    line_v, = axs[0, 1].plot([], [], **exact_line_kwargs)
    axs[0, 1].set_ylabel(r"Velocity ($v_x$)")
    axs[0, 1].set_title("Velocity Profile")
    
    scat_P = axs[1, 0].scatter([], [], **scatter_kwargs)
    line_P, = axs[1, 0].plot([], [], **exact_line_kwargs)
    axs[1, 0].set_ylabel(r"Pressure ($P$)")
    axs[1, 0].set_title("Pressure Profile")
    
    scat_u = axs[1, 1].scatter([], [], **scatter_kwargs)
    line_u, = axs[1, 1].plot([], [], **exact_line_kwargs)
    axs[1, 1].set_ylabel(r"Internal Energy ($u$)")
    axs[1, 1].set_title("Specific Internal Energy")

    # Time Series L1 Error
    axs[2, 0].plot(times, l1_errors, label=r'$L_1$ Density Error', color='black', lw=1.5)
    vline_time = axs[2, 0].axvline(0, color='gray', linestyle='--', alpha=0.7)
    if len(times) > 1:
        axs[2, 0].set_xlim(min(times), max(times))
    axs[2, 0].set_ylabel("L1 Error Norm")
    axs[2, 0].set_xlabel("Time")
    axs[2, 0].set_title("Density Convergence")
    axs[2, 0].set_yscale('log')
    axs[2, 0].legend(loc='lower right', fontsize=9)

    # Embed Maximum L1 Error directly in the L1 plot
    max_l1_error = np.max(l1_errors) if l1_errors else 0.0
    axs[2, 0].text(0.05, 0.85, f"Max $L_1$ Error: {max_l1_error:.5e}", 
                   transform=axs[2, 0].transAxes,
                   ha='left', va='center', fontsize=12, fontweight='bold',
                   bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5', alpha=0.8))

    # Macroscopic Conservation Tracking
    axs[2, 1].plot(times, vol_errors, label=r'Volume Error $(\sum V_i / L^3 - 1)$', color='blue', lw=1.5)
    axs[2, 1].plot(times, energy_errors, label=r'Energy Error $(\Delta E / E_0)$', color='red', lw=1.5)
    axs[2, 1].plot(times, mom_errors, label=r'Momentum Drift ($\Delta P_x$)', color='green', lw=1.5, linestyle=':')
    
    vline_time2 = axs[2, 1].axvline(0, color='gray', linestyle='--', alpha=0.7)
    if len(times) > 1:
        axs[2, 1].set_xlim(min(times), max(times))
    
    axs[2, 1].set_ylabel("Fractional Error")
    axs[2, 1].set_xlabel("Time")
    axs[2, 1].set_title("Macroscopic Conservation Diagnostics")
    axs[2, 1].grid(True, linestyle=':', alpha=0.6)
    axs[2, 1].legend(loc='best', fontsize=9)

    for ax in axs.flat:
        if ax != axs[2, 0] and ax != axs[2, 1]:
            ax.set_xlim(0, domain_size)
            ax.grid(True, linestyle=':', alpha=0.6)
    axs[0, 0].legend(loc='upper right', fontsize=9)

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
                x_v = np.linspace(dx/2.0, domain_size - dx/2.0, N)
                
                rho_3d = f['Gas/density'][:].reshape((N, N, N))
                mom_x_3d = f['Gas/momentum_x'][:].reshape((N, N, N))
                P_3d = f['Gas/pressure'][:].reshape((N, N, N))
                
                mid = N // 2
                rho_v = rho_3d[:, mid, mid]
                v_v = mom_x_3d[:, mid, mid] / rho_v
                P_v = P_3d[:, mid, mid]
                u_v = P_v / (rho_v * (gamma - 1.0))
                
            elif method == "mfm":
                x_v = f['Gas/position_x'][:]
                rho_v = f['Gas/density'][:]
                v_v = f['Gas/velocity_x'][:]
                u_v = f['Gas/internal_energy'][:]
                P_v = rho_v * u_v * (gamma - 1.0)
            
        fig.canvas.manager.set_window_title(os.path.basename(os.path.normpath(snapshot_dir)))
        fig.suptitle(f"Linear Soundwave Validation - Snapshot {idx} (t={time:.4f})", fontsize=16)
        
        scat_rho.set_offsets(np.c_[x_v, rho_v])
        scat_v.set_offsets(np.c_[x_v, v_v])
        scat_P.set_offsets(np.c_[x_v, P_v])
        scat_u.set_offsets(np.c_[x_v, u_v])

        x_exact = np.linspace(0, domain_size, 1000)
        rho_ex, v_ex, P_ex, u_ex = get_exact_soundwave_solution(x_exact, time, domain_size, gamma, rho_0_emp, P_0_emp)
        
        line_rho.set_data(x_exact, rho_ex)
        line_v.set_data(x_exact, v_ex)
        line_P.set_data(x_exact, P_ex)
        line_u.set_data(x_exact, u_ex)

        vline_time.set_xdata([time, time])
        vline_time2.set_xdata([time, time])
        fig.canvas.draw_idle()

    ax_slider = fig.add_axes([0.15, 0.02, 0.7, 0.03])
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
    parser = argparse.ArgumentParser(description="Interactive MFM/Eulerian Soundwave validation.")
    parser.add_argument("path", type=str, nargs='?', default=".", help="Path to snapshot directory.")
    parser.add_argument("-l", "--latest", action="store_true", help="Load latest run_* directory")
    
    args = parser.parse_args()
    target_dir = args.path
    
    if args.latest:
        runs = sorted(glob.glob(os.path.join(target_dir, "run_*")))
        if runs: target_dir = runs[-1]
        
    validate_soundwave_interactive(target_dir)