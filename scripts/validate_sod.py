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
    # GADGET/GIZMO Initial States (Hernquist & Katz variant)
    P_L, rho_L, v_L = 1.0, 1.0, 0.0
    P_R, rho_R, v_R = 0.1795, 0.25, 0.0
    
    if t <= 1e-8:
        rho = np.where(x <= x0, rho_L, rho_R)
        P = np.where(x <= x0, P_L, P_R)
        v = np.zeros_like(x)
        ie = P / (rho * (gamma - 1.0))
        return rho, v, P, ie

    # Dynamic Riemann Solver for the Star State
    mu2 = (gamma - 1.0) / (gamma + 1.0)
    c_L = np.sqrt(gamma * P_L / rho_L)
    c_R = np.sqrt(gamma * P_R / rho_R)
    
    # Newton-Raphson iteration to find the intermediate pressure (P_star)
    P_star = 0.5 * (P_L + P_R) 
    for _ in range(100):
        # Left Rarefaction Wave
        term_L = (P_star / P_L)**((gamma - 1.0) / (2.0 * gamma))
        f_L = (2.0 * c_L / (gamma - 1.0)) * (term_L - 1.0)
        df_L_dP = (c_L / (gamma * P_L)) * (P_star / P_L)**(-(gamma + 1.0) / (2.0 * gamma))
        
        # Right Shock Wave
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

    # Calculate intermediate states from P_star
    P_3 = P_star
    u_3 = v_R + (P_star - P_R) * np.sqrt(A_R / (P_star + B_R))
    rho_3 = rho_L * (P_star / P_L)**(1.0 / gamma)
    rho_4 = rho_R * (P_star + mu2 * P_R) / (P_R + mu2 * P_star)
    
    c_3 = np.sqrt(gamma * P_3 / rho_3)
    V_shock = v_R + c_R * np.sqrt((gamma + 1.0) / (2.0 * gamma) * (P_star / P_R) + (gamma - 1.0) / (2.0 * gamma))
    
    # Wave Boundaries
    x_head = x0 - c_L * t
    x_tail = x0 + (u_3 - c_3) * t
    x_contact = x0 + u_3 * t
    x_shock = x0 + V_shock * t

    # Map regions
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
        
        if setup_type != "sod_shock_tube":
            print(f"[ERROR] Expected config attribute 'setup'='sod_shock_tube'. Found: '{setup_type}'")
            sys.exit(1)
            
        gamma = config.get('gamma', 1.4)
        domain_size = config.get('domain_size', 1.0)

        # Check Gamma
        if abs(gamma - (5.0/3.0)) > 1e-5:
            print(f"[WARNING] Expected gamma = 5/3 for this Shock test, but found {gamma:.5f}.")
            print("          The analytical solution may not visually align with your data!")
        
        # Check Cosmological Expansion
        if bool(config.get('expanding_universe', 0)):
            print("[ERROR] Sod shock test requires a static box, but 'expanding_universe' is TRUE.")
            sys.exit(1)
            
        # Check Radiative Cooling
        if bool(config.get('enable_cooling', 0)):
            print("[ERROR] Sod shock test requires adiabatic gas, but 'enable_cooling' is TRUE.")
            sys.exit(1)
            
        # Check Gravity (PM and PP)
        use_pm = bool(config.get('use_pm', 0))
        use_pp = bool(config.get('use_pp', 0))
        if use_pm or use_pp:
            print(f"[ERROR] Sod shock test must be pure hydrodynamics.")
            print(f"        Found gravity enabled: PM={use_pm}, PP={use_pp}")
            sys.exit(1)

    # Plot Setup
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    plt.subplots_adjust(bottom=0.15) 

    axs[0, 0].set_ylim(-0.05, 1.15)
    axs[0, 1].set_ylim(-0.1, 1.2)  
    axs[1, 0].set_ylim(-0.05, 1.15) 
    
    # Styling variables
    scatter_kwargs = {'s': 4, 'color': 'royalblue', 'alpha': 0.6, 'label': 'MFM Simulation'}
    sim_line_kwargs = {'color': 'royalblue', 'lw': 0.5, 'alpha': 0.8} # The thin connecting line
    exact_line_kwargs = {'color': 'black', 'lw': 1.5, 'linestyle': '--', 'label': 'Exact Solution'}
    
    # Density
    line_sim_rho, = axs[0, 0].plot([], [], **sim_line_kwargs)
    scat_rho = axs[0, 0].scatter([], [], **scatter_kwargs)
    line_rho, = axs[0, 0].plot([], [], **exact_line_kwargs)
    axs[0, 0].set_ylabel(r"Density ($\rho$)")
    axs[0, 0].set_title("Density Profile")
    
    # Velocity
    line_sim_v, = axs[0, 1].plot([], [], **sim_line_kwargs)
    scat_v = axs[0, 1].scatter([], [], **scatter_kwargs)
    line_v, = axs[0, 1].plot([], [], **exact_line_kwargs)
    axs[0, 1].set_ylabel(r"Velocity ($v_x$)")
    axs[0, 1].set_title("Velocity Profile")
    
    # Pressure
    line_sim_P, = axs[1, 0].plot([], [], **sim_line_kwargs)
    scat_P = axs[1, 0].scatter([], [], **scatter_kwargs)
    line_P, = axs[1, 0].plot([], [], **exact_line_kwargs)
    axs[1, 0].set_ylabel(r"Pressure ($P$)")
    axs[1, 0].set_xlabel("Position (x)")
    axs[1, 0].set_title("Pressure Profile")
    
    # Internal Energy
    line_sim_u, = axs[1, 1].plot([], [], **sim_line_kwargs)
    scat_u = axs[1, 1].scatter([], [], **scatter_kwargs)
    line_u, = axs[1, 1].plot([], [], **exact_line_kwargs)
    
    axs[1, 1].set_ylabel(r"Internal Energy ($u$)")
    axs[1, 1].set_xlabel("Position (x)")
    axs[1, 1].set_title("Specific Internal Energy")

    for ax in axs.flat:
        ax.set_xlim(0, domain_size)
        ax.grid(True, linestyle=':', alpha=0.6)
    axs[0, 0].legend(loc='upper right', fontsize=9)


    def update(val):
        idx = int(snap_slider.val)
        target_file = files[idx]
        
        with h5py.File(target_file, 'r') as f:
            time = f['Header'].attrs['simulation_time']
            x = f['Gas/position_x'][:]
            y = f['Gas/position_y'][:]
            z = f['Gas/position_z'][:]
            rho = f['Gas/density'][:]
            v_x = f['Gas/velocity_x'][:]
            u = f['Gas/internal_energy'][:]
            P = rho * u * (gamma - 1.0)
            
        fig.suptitle(f"Sod Shock Validation - Snapshot {idx} (t={time:.4f})", fontsize=16)

        # In the update(val) function:
        c_L = np.sqrt(gamma * 1.0 / 1.0)
        
        # Recalculate dynamic V_shock for the gray boundary boxes
        P_R_plot = 0.1795
        rho_R_plot = 0.25
        c_R_plot = np.sqrt(gamma * P_R_plot / rho_R_plot)
        _, _, P_star_plot, _ = get_exact_sod_solution(np.array([domain_size/2.0]), 1e-8, gamma)
        V_shock = c_R_plot * np.sqrt((gamma + 1.0) / (2.0 * gamma) * (P_star_plot[0] / P_R_plot) + (gamma - 1.0) / (2.0 * gamma))
        
        left_bnd = c_L * time
        right_bnd = domain_size - (V_shock * time)

        # We only keep the inner part of the Y-Z plane to ensure we are looking at pristine 1D flow.
        safe_y_min, safe_y_max = domain_size * 0.495, domain_size * 0.505
        safe_z_min, safe_z_max = domain_size * 0.495, domain_size * 0.505

        # Add a small buffer to account for the kernel smoothing length (h) "reaching" into the error
        h_buffer = 0.1 
        
        # Calculate dynamic intrusion depth
        transverse_intrusion = (c_L * time) + h_buffer
        
        # Ensure we don't accidentally invert the bounds if the waves cross the center
        #safe_y_min = min(transverse_intrusion, domain_size / 2.0)
        #safe_y_max = max(domain_size - transverse_intrusion, domain_size / 2.0)
        #safe_z_min = min(transverse_intrusion, domain_size / 2.0)
        #safe_z_max = max(domain_size - transverse_intrusion, domain_size / 2.0)
        
        valid_mask = (x >= left_bnd) & (x <= right_bnd) & \
                     (y >= safe_y_min) & (y <= safe_y_max) & \
                     (z >= safe_z_min) & (z <= safe_z_max)
        
        #valid_mask = (x >= left_bnd) & (x <= right_bnd)
        
        # Filter valid particles
        x_v = x[valid_mask]
        rho_v = rho[valid_mask]
        v_v = v_x[valid_mask]
        P_v = P[valid_mask]
        u_v = u[valid_mask]
        
        # Sort arrays strictly by x-coordinate to prevent zig-zag lines
        sort_idx = np.argsort(x_v)
        x_sorted = x_v[sort_idx]
        
        # Update Scatter Points
        scat_rho.set_offsets(np.c_[x_v, rho_v])
        scat_v.set_offsets(np.c_[x_v, v_v])
        scat_P.set_offsets(np.c_[x_v, P_v])
        scat_u.set_offsets(np.c_[x_v, u_v])

        # Thin Connecting Lines (using sorted data)
        line_sim_rho.set_data(x_sorted, rho_v[sort_idx])
        line_sim_v.set_data(x_sorted, v_v[sort_idx])
        line_sim_P.set_data(x_sorted, P_v[sort_idx])
        line_sim_u.set_data(x_sorted, u_v[sort_idx])

        # Update Analytical Lines
        x_exact = np.linspace(0, domain_size, 1000)
        rho_ex, v_ex, P_ex, u_ex = get_exact_sod_solution(x_exact, time, gamma, x0=domain_size/2.0)
        
        line_rho.set_data(x_exact, rho_ex)
        line_v.set_data(x_exact, v_ex)
        line_P.set_data(x_exact, P_ex)
        line_u.set_data(x_exact, u_ex)

        if len(u_v) > 0:
            # Find the min and max from both simulation and exact solution
            current_u_min = min(np.min(u_v), np.min(u_ex))
            current_u_max = max(np.max(u_v), np.max(u_ex))
            
            # Add a 15% padding
            u_range = current_u_max - current_u_min
            padding = u_range * 0.15 if u_range > 0 else 0.5
            
            axs[1, 1].set_ylim(current_u_min - padding, current_u_max + padding)

        for ax in axs.flat:
            [p.remove() for p in reversed(ax.patches)]
            ax.axvspan(0, min(left_bnd, domain_size), color='gray', alpha=0.2)
            if right_bnd < domain_size:
                ax.axvspan(max(right_bnd, 0), domain_size, color='gray', alpha=0.2)

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
    parser = argparse.ArgumentParser(description="Interactive MFM Sod Shock Tube validation.")
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