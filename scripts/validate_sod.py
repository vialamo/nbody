#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def get_exact_sod_solution(x, t, gamma=5.0/3.0, x0=0.5):
    """
    Computes the analytical solution for the standard Sod Shock Tube.
    Handles t=0 to return the initial step function.
    """
    P_L, rho_L, v_L = 1.0, 1.0, 0.0
    P_R, rho_R, v_R = 0.1, 0.125, 0.0
    
    if t <= 1e-8:
        rho = np.where(x <= x0, rho_L, rho_R)
        P = np.where(x <= x0, P_L, P_R)
        v = np.zeros_like(x)
        ie = P / (rho * (gamma - 1.0))
        return rho, v, P, ie

    mu2 = (gamma - 1.0) / (gamma + 1.0)
    c_L = np.sqrt(gamma * P_L / rho_L)
    c_R = np.sqrt(gamma * P_R / rho_R)
    
    P_star = 0.5 * (P_L + P_R) 
    for _ in range(100):
        term_L = (P_star / P_L)**((gamma - 1.0) / (2.0 * gamma))
        f_L = (2.0 * c_L / (gamma - 1.0)) * (term_L - 1.0)
        df_L_dP = (c_L / (gamma * P_L)) * (P_star / P_L)**(-(gamma + 1.0) / (2.0 * gamma))
        
        A_R = 2.0 / ((gamma + 1.0) * rho_R)
        B_R = mu2 * P_R
        f_R = (P_star - P_R) * np.sqrt(A_R / (P_star + B_R))
        df_R_dP = np.sqrt(A_R / (P_star + B_R)) * (1.0 - 0.5 * (P_star - P_R) / (P_star + B_R))
        
        f = f_L + f_R - (v_L - v_R)
        df_dP = df_L_dP + df_R_dP
        
        dP = -f / df_dP
        P_star += dP
        if abs(dP) < 1e-8:
            break

    P_3 = P_star
    u_3 = v_R + (P_star - P_R) * np.sqrt(A_R / (P_star + B_R))
    rho_3 = rho_L * (P_star / P_L)**(1.0 / gamma)
    rho_4 = rho_R * (P_star + mu2 * P_R) / (P_R + mu2 * P_star)
    
    c_3 = np.sqrt(gamma * P_3 / rho_3)
    V_shock = v_R + c_R * np.sqrt((gamma + 1.0) / (2.0 * gamma) * (P_star / P_R) + (gamma - 1.0) / (2.0 * gamma))
    
    x_head = x0 - c_L * t
    x_tail = x0 + (u_3 - c_3) * t
    x_contact = x0 + u_3 * t
    x_shock = x0 + V_shock * t

    rho = np.zeros_like(x)
    P = np.zeros_like(x)
    v = np.zeros_like(x)

    for i, pos in enumerate(x):
        if pos <= x_head:
            rho[i], P[i], v[i] = rho_L, P_L, v_L
        elif pos > x_head and pos <= x_tail:
            v_fan = (2.0 / (gamma + 1.0)) * (c_L + (pos - x0) / t)
            c_fan = c_L - 0.5 * (gamma - 1.0) * v_fan
            rho[i] = rho_L * (c_fan / c_L)**(2.0 / (gamma - 1.0))
            P[i] = P_L * (c_fan / c_L)**(2.0 * gamma / (gamma - 1.0))
            v[i] = v_fan
        elif pos > x_tail and pos <= x_contact:
            rho[i], P[i], v[i] = rho_3, P_3, u_3
        elif pos > x_contact and pos <= x_shock:
            rho[i], P[i], v[i] = rho_4, P_3, u_3
        else:
            rho[i], P[i], v[i] = rho_R, P_R, v_R
            
    ie = P / (rho * (gamma - 1.0))
    return rho, v, P, ie

def validate_sod_shock_interactive(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        setup_type = config.get('setup', b"").decode('utf-8')
        hydro_method = config.get('hydro_method', b"none").decode('utf-8')
        
        if setup_type != "sod_shock_tube":
            print(f"[ERROR] Expected config attribute 'setup'='sod_shock_tube'. Found: '{setup_type}'")
            sys.exit(1)
            
        gamma = config.get('gamma', 1.4)
        domain_size = config.get('domain_size', 1.0)
        domain_size_x = config.get('domain_size_x', 1.0)

        if abs(gamma - 1.4) > 1e-5:
            print(f"[WARNING] Expected gamma = 1.4 for this Shock test, but found {gamma:.5f}.")
        if bool(config.get('expanding_universe', 0)):
            print("[ERROR] Sod shock test requires a static box, but 'expanding_universe' is TRUE.")
            sys.exit(1)
        if bool(config.get('enable_cooling', 0)):
            print("[ERROR] Sod shock test requires adiabatic gas, but 'enable_cooling' is TRUE.")
            sys.exit(1)
        if bool(config.get('use_pm', 0)) or bool(config.get('use_pp', 0)):
            print(f"[ERROR] Sod shock test must be pure hydrodynamics.")
            sys.exit(1)

    print("Pre-calculating global energy tracking from all snapshots...")
    times = []
    err_tot = []
    err_sw = []
    E0 = None

    for f_name in files:
        with h5py.File(f_name, 'r') as f:
            t = f['Header'].attrs['simulation_time']
            times.append(t)
            
            if 'Gas' in f:
                gas = f['Gas']
                if hydro_method == "mfm":
                    if 'mass' in gas and 'internal_energy' in gas and 'velocity_x' in gas:
                        mass = gas['mass'][:]
                        u_ie = gas['internal_energy'][:]
                        v_x = gas['velocity_x'][:]
                        v_y = gas['velocity_y'][:] if 'velocity_y' in gas else 0.0
                        v_z = gas['velocity_z'][:] if 'velocity_z' in gas else 0.0
                        
                        if 'total_energy' in gas:
                            e_tot = gas['total_energy'][:]
                            E_current = np.sum(mass * e_tot)
                        else:
                            ke = 0.5 * (v_x**2 + v_y**2 + v_z**2)
                            E_current = np.sum(mass * (u_ie + ke))
                    else:
                        E_current = 0.0
                        
                elif hydro_method == "eulerian":
                    N = config.get('mesh_size_1d', 32)
                    dx = domain_size / N
                    vol = dx**3
                    rho_g = gas['density'][:]
                    p_g = gas['pressure'][:]
                    mom_x = gas['momentum_x'][:]
                    mom_y = gas['momentum_y'][:] if 'momentum_y' in gas else 0.0
                    mom_z = gas['momentum_z'][:] if 'momentum_z' in gas else 0.0
                    
                    kin = 0.5 * (mom_x**2 + mom_y**2 + mom_z**2) / rho_g
                    therm = p_g / (gamma - 1.0)
                    E_current = np.sum((therm + kin) * vol)
                else:
                    E_current = 0.0

                if E0 is None and E_current != 0:
                    E0 = E_current
                
                if E0 is not None and E0 != 0:
                    delta_e = E_current - E0
                    w_grav = gas.attrs.get('cumulative_gravitational_work', 0.0)
                    w_exp = gas.attrs.get('cumulative_expansion_work', 0.0)
                    e_rad = gas.attrs.get('cumulative_radiated_energy', 0.0)
                    e_heat = gas.attrs.get('cumulative_photoheating_energy', 0.0)
                    
                    absolute_error = delta_e - w_grav + w_exp + e_rad - e_heat
                    frac_err = absolute_error / abs(E0)
                else:
                    frac_err = 0.0
                err_tot.append(frac_err)
                
                sw_energy = gas.attrs.get('cumulative_entropy_switch_energy', 0.0)
                if E0 is not None and E0 != 0:
                    sw_err = sw_energy / abs(E0)
                else:
                    sw_err = 0.0
                err_sw.append(sw_err)
            else:
                err_tot.append(0.0)
                err_sw.append(0.0)

    label_prefix = "Eulerian Grid" if hydro_method == "eulerian" else "MFM Particles"
    color_prefix = "darkorange" if hydro_method == "eulerian" else "royalblue"
    
    # GridSpec setup: 3 rows, 2 columns. Bottom row spans both columns and is half the height.
    fig = plt.figure(figsize=(14, 14))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.5])
    plt.subplots_adjust(bottom=0.10, hspace=0.3) 

    ax_rho = fig.add_subplot(gs[0, 0])
    ax_v   = fig.add_subplot(gs[0, 1])
    ax_P   = fig.add_subplot(gs[1, 0])
    ax_u   = fig.add_subplot(gs[1, 1])
    ax_err = fig.add_subplot(gs[2, :])
    
    spatial_axes = [ax_rho, ax_v, ax_P, ax_u]

    ax_rho.set_ylim(-0.05, 1.15)
    ax_v.set_ylim(-0.1, 1.2)  
    ax_P.set_ylim(-0.05, 1.15) 
    
    scatter_kwargs = {'s': 6 if hydro_method == "eulerian" else 4, 'color': color_prefix, 'alpha': 0.8, 'label': f'{label_prefix} Sim'}
    exact_line_kwargs = {'color': 'black', 'lw': 1.5, 'linestyle': '--', 'label': 'Exact Solution'}
    
    # Density
    scat_rho = ax_rho.scatter([], [], **scatter_kwargs)
    line_rho, = ax_rho.plot([], [], **exact_line_kwargs)
    ax_rho.set_ylabel(r"Density ($\rho$)")
    ax_rho.set_title("Density Profile")
    
    # Velocity
    scat_v = ax_v.scatter([], [], **scatter_kwargs)
    line_v, = ax_v.plot([], [], **exact_line_kwargs)
    ax_v.set_ylabel(r"Velocity ($v_x$)")
    ax_v.set_title("Velocity Profile")
    
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

    # Global Fractional Errors
    ax_err.plot(times, err_tot, label=r'$\Delta E_{tot} / E_0$', color='black', lw=1.5)
    ax_err.plot(times, err_sw, label=r'$E_{switch} / E_0$', color='red', lw=1.5, linestyle='-.')
    vline_time = ax_err.axvline(0, color='gray', linestyle='--', alpha=0.7)
    if len(times) > 1:
        ax_err.set_xlim(min(times), max(times))
    ax_err.set_ylabel("Fractional Error")
    ax_err.set_xlabel("Time")
    ax_err.set_title("Global Energy Conservation")
    ax_err.legend(loc='upper left', fontsize=9)

    # Format bounds/grids
    for ax in spatial_axes:
        ax.set_xlim(0, domain_size_x)
        ax.grid(True, linestyle=':', alpha=0.6)
    ax_err.grid(True, linestyle=':', alpha=0.6)
    ax_rho.legend(loc='upper right', fontsize=9)

    def update(val):
        idx = int(snap_slider.val)
        target_file = files[idx]
        
        with h5py.File(target_file, 'r') as f:
            time = f['Header'].attrs['simulation_time']
            cfg = f['Config'].attrs
            method = cfg.get('hydro_method', b'none').decode('utf-8')
            
            left_bnd = 0.05 * domain_size_x
            right_bnd = 0.95 * domain_size_x
            
            if method == "eulerian":
                N = cfg.get('mesh_size_1d', 32)
                dx = domain_size / float(N)
                
                x_all = np.linspace(dx/2.0, domain_size - dx/2.0, N)
                rho_3d = f['Gas/density'][:].reshape((N, N, N))
                mom_x_3d = f['Gas/momentum_x'][:].reshape((N, N, N))
                P_3d = f['Gas/pressure'][:].reshape((N, N, N))
                
                mid = N // 2
                rho_line = rho_3d[:, mid, mid]
                v_line = mom_x_3d[:, mid, mid] / rho_line
                P_line = P_3d[:, mid, mid]
                u_line = P_line / (rho_line * (gamma - 1.0))
                
                valid_mask = (x_all >= left_bnd) & (x_all <= right_bnd)
                x_v = x_all[valid_mask]
                rho_v = rho_line[valid_mask]
                v_v = v_line[valid_mask]
                P_v = P_line[valid_mask]
                u_v = u_line[valid_mask]
                
            elif method == "mfm":
                x = f['Gas/position_x'][:]
                rho = f['Gas/density'][:]
                v_x = f['Gas/velocity_x'][:]
                u = f['Gas/internal_energy'][:]
                P = rho * u * (gamma - 1.0)
                
                valid_mask = (x >= left_bnd) & (x <= right_bnd)
                x_v = x[valid_mask]
                rho_v = rho[valid_mask]
                v_v = v_x[valid_mask]
                P_v = P[valid_mask]
                u_v = u[valid_mask]
            
        fig.canvas.manager.set_window_title(os.path.basename(os.path.normpath(snapshot_dir)))
        fig.suptitle(f"Sod Shock Validation - Snapshot {idx} (t={time:.4f})", fontsize=16)
        
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
        rho_ex, v_ex, P_ex, u_ex = get_exact_sod_solution(x_exact, time, gamma, x0=domain_size_x/2.0)
        
        line_rho.set_data(x_exact, rho_ex)
        line_v.set_data(x_exact, v_ex)
        line_P.set_data(x_exact, P_ex)
        line_u.set_data(x_exact, u_ex)

        vline_time.set_xdata([time, time])

        if len(u_v) > 0:
            current_u_min = min(np.min(u_v), np.min(u_ex))
            current_u_max = max(np.max(u_v), np.max(u_ex))
            u_range = current_u_max - current_u_min
            padding = u_range * 0.15 if u_range > 0 else 0.5
            ax_u.set_ylim(current_u_min - padding, current_u_max + padding)

        for ax in spatial_axes:
            [p.remove() for p in reversed(ax.patches)]
            ax.axvspan(0, min(left_bnd, domain_size_x), color='gray', alpha=0.2)
            if right_bnd < domain_size_x:
                ax.axvspan(max(right_bnd, 0), domain_size_x, color='gray', alpha=0.2)

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
    parser = argparse.ArgumentParser(description="Interactive MFM/Eulerian Sod Shock Tube validation.")
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
        
    validate_sod_shock_interactive(target_dir)
