#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt

def audit_gravity_snapshots(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    print(f"==================================================================")
    print(f"               GRAVITY SOLVER AUDIT & VALIDATION                  ")
    print(f"==================================================================")
    print(f"Found {len(files)} snapshots. Extracting particle data...\n")

    # Read configuration from the first snapshot
    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        is_expanding = bool(config.get('expanding_universe', 1))
        G = config.get('G', 1.0)
        box_size = config.get('domain_size', 1.0)
        method = config.get('hydro_method', b"none").decode('utf-8')
        has_particle_hydro = bool(method == "mfm")
        has_gas = 'Gas' in f and has_particle_hydro

    print(f"Mode Auto-detected : {'EXPANDING UNIVERSE' if is_expanding else 'STATIC BOX (a=1.0)'}")
    print(f"Gas Particles Present : {has_gas}")
    print(f"Newtonian Constant G  : {G:.6e}")
    print("------------------------------------------------------------------\n")

    # Time series arrays
    times = []
    scale_factors = []
    
    # Kinetic Energies
    ke_dm = []
    ke_gas = []
    te_gas = []
    
    # Linear Momentum components
    momentum_x, momentum_y, momentum_z = [], [], []
    
    # Net Force (Newton's 3rd Law)
    f_net_mag = []
    
    # Gas Particle Cumulative Work (from HDF5 attributes)
    w_grav_gas_list = []
    w_exp_gas_list = []
    
    # DM Numerically Integrated Work
    w_grav_dm_list = []
    w_exp_dm_list = []
    
    for i, f_path in enumerate(files):
        sys.stdout.write(f"\rProcessing [{i+1}/{len(files)}] {os.path.basename(f_path)}")
        sys.stdout.flush()

        with h5py.File(f_path, 'r') as f:
            header = f['Header'].attrs
            a = header['scale_factor']
            t = header['simulation_time']
            dt = header['dt_macro']
            
            times.append(t)
            scale_factors.append(a)
            
            # Dark Matter Particles         
            dm_vx = f['Particles/velocity_x'][:]
            dm_vy = f['Particles/velocity_y'][:]
            dm_vz = f['Particles/velocity_z'][:]
            
            dm_ax = f['Particles/acceleration_x'][:]
            dm_ay = f['Particles/acceleration_y'][:]
            dm_az = f['Particles/acceleration_z'][:]
            
            dm_m  = f['Particles/mass'][:]

            w_grav_dm_list.append(f['Particles'].attrs.get('cumulative_gravitational_work', 0.0))
            w_exp_dm_list.append(f['Particles'].attrs.get('cumulative_expansion_work', 0.0))
            
            # Kinetic Energy & Momentum
            k_dm_step = 0.5 * np.sum(dm_m * (dm_vx**2 + dm_vy**2 + dm_vz**2))
            ke_dm.append(k_dm_step)
            
            px = np.sum(dm_m * dm_vx)
            py = np.sum(dm_m * dm_vy)
            pz = np.sum(dm_m * dm_vz)
            
            # Net Force on DM
            fnx = np.sum(dm_m * dm_ax)
            fny = np.sum(dm_m * dm_ay)
            fnz = np.sum(dm_m * dm_az)
            
            # Gas Particles (If enabled)
            if has_gas:
                gas_vx = f['Gas/velocity_x'][:]
                gas_vy = f['Gas/velocity_y'][:]
                gas_vz = f['Gas/velocity_z'][:]
                
                gas_ax = f['Gas/acceleration_x'][:]
                gas_ay = f['Gas/acceleration_y'][:]
                gas_az = f['Gas/acceleration_z'][:]
                
                gas_m  = f['Gas/mass'][:]
                gas_u  = f['Gas/internal_energy'][:]
                
                k_gas_step = 0.5 * np.sum(gas_m * (gas_vx**2 + gas_vy**2 + gas_vz**2))
                u_gas_step = np.sum(gas_m * gas_u)
                ke_gas.append(k_gas_step)
                te_gas.append(k_gas_step + u_gas_step)
                
                px += np.sum(gas_m * gas_vx)
                py += np.sum(gas_m * gas_vy)
                pz += np.sum(gas_m * gas_vz)
                
                fnx += np.sum(gas_m * gas_ax)
                fny += np.sum(gas_m * gas_ay)
                fnz += np.sum(gas_m * gas_az)
                
                # Gas Cumulative Work from HDF5 attributes
                w_grav_gas_list.append(f['Gas'].attrs.get('cumulative_gravitational_work', 0.0))
                w_exp_gas_list.append(f['Gas'].attrs.get('cumulative_expansion_work', 0.0))
            else:
                ke_gas.append(0.0)
                te_gas.append(0.0)
                w_grav_gas_list.append(0.0)
                w_exp_gas_list.append(0.0)

            momentum_x.append(px)
            momentum_y.append(py)
            momentum_z.append(pz)
            
            f_net_mag.append(np.sqrt(fnx**2 + fny**2 + fnz**2))

    print("\n\nData collection complete. Evaluating diagnostic metrics...")

    # Convert to Numpy Arrays
    times = np.array(times)
    scale_factors = np.array(scale_factors)
    ke_dm = np.array(ke_dm)
    ke_gas = np.array(ke_gas)
    te_gas = np.array(te_gas)
    
    p_mag = np.sqrt(np.array(momentum_x)**2 + np.array(momentum_y)**2 + np.array(momentum_z)**2)
    f_net_mag = np.array(f_net_mag)
    
    w_grav_gas = np.array(w_grav_gas_list)
    w_exp_gas = np.array(w_exp_gas_list)
    
    w_grav_dm = np.array(w_grav_dm_list)
    w_exp_dm = np.array(w_exp_dm_list)

    # Work-Energy Fractional Errors
    # Gas Error
    e0_gas = te_gas[0] if len(te_gas) > 0 and te_gas[0] > 0 else 1.0
    delta_e_gas = te_gas - te_gas[0]
    err_gas = (delta_e_gas - w_grav_gas + w_exp_gas) / e0_gas

    # DM Error
    k0_dm = ke_dm[0] if ke_dm[0] > 0 else 1.0
    delta_ke_dm = ke_dm - ke_dm[0]
    err_dm = (delta_ke_dm - w_grav_dm + w_exp_dm) / k0_dm

    # Print Report
    max_f_net = np.max(f_net_mag)
    max_err_gas = np.max(np.abs(err_gas)) if has_gas else 0.0
    max_err_dm = np.max(np.abs(err_dm))

    print(f"------------------------------------------------------------------")
    print(f"                       AUDIT SUMMARY REPORT                       ")
    print(f"------------------------------------------------------------------")
    print(f"  Max Net System Force |F_net| : {max_f_net:.6e}  ", end="")
    print("[PASS]" if max_f_net < 1e-10 else "[WARNING: Non-zero net force! Check PP/PM symmetry]")
    
    print(f"  Max Momentum Drift |P_total| : {np.max(p_mag):.6e}")
    
    if has_gas:
        print(f"  Max Gas Energy Error        : {max_err_gas * 100.0:.4f}% ", end="")
        print("[PASS]" if max_err_gas < 0.01 else "[FAIL: >1% drift. Check hydro force array leakage!]")
        
    print(f"  Max DM Energy Error         : {max_err_dm * 100.0:.4f}% ", end="")
    print("[PASS]" if max_err_dm < 0.01 else "[FAIL: >1% drift. Check gravity integration!]")
    print(f"------------------------------------------------------------------\n")

    # Plotting
    x_axis = scale_factors if is_expanding else times
    x_label = "Scale Factor (a)" if is_expanding else "Simulation Time [Code Units]"

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.canvas.manager.set_window_title(f"Gravity Validation Audit - {'Expanding' if is_expanding else 'Static'}")

    # Plot 1: Energy Evolution
    axs[0, 0].plot(x_axis, ke_dm, color='blue', lw=2, label='DM Kinetic Energy')
    if has_gas:
        axs[0, 0].plot(x_axis, ke_gas, color='green', lw=2, linestyle='--', label='Gas Kinetic Energy')
        axs[0, 0].plot(x_axis, te_gas, color='darkgreen', lw=2, label='Gas Total Energy (K+U)')
    axs[0, 0].set_xlabel(x_label)
    axs[0, 0].set_ylabel("Energy [Code Units]")
    axs[0, 0].set_title("Energy Evolution")
    axs[0, 0].grid(True, linestyle=':', alpha=0.6)
    axs[0, 0].legend()

    # Plot 2: Fractional Work-Energy Conservation Error
    axs[0, 1].plot(x_axis, err_dm, color='blue', lw=1.5, linestyle='-', label='DM Error (Numeric Work)')
    if has_gas:
        axs[0, 1].plot(x_axis, err_gas, color='green', lw=1.5, linestyle='--', label='Gas Error (Logged Work)')
    
    axs[0, 1].axhline(0, color='black', linestyle='-', linewidth=0.8)
    axs[0, 1].axhline(0.01, color='red', linestyle=':', alpha=0.5, label='±1% Threshold')
    axs[0, 1].axhline(-0.01, color='red', linestyle=':', alpha=0.5)
    
    axs[0, 1].set_xlabel(x_label)
    axs[0, 1].set_ylabel(r"Fractional Error $\frac{\Delta K - W_{\mathrm{grav}} + W_{\mathrm{exp}}}{K_0}$")
    axs[0, 1].set_title("Work-Energy Conservation Error")
    
    # Set symmetric y-limits
    max_plot_err = max(0.02, np.max(np.abs(err_dm)), np.max(np.abs(err_gas)) if has_gas else 0.0) * 1.2
    axs[0, 1].set_ylim(-max_plot_err, max_plot_err)
    axs[0, 1].grid(True, linestyle=':', alpha=0.6)
    axs[0, 1].legend()

    # Plot 3: Net Force (Newton's 3rd Law)
    axs[1, 0].plot(x_axis, f_net_mag, color='purple', lw=2, label=r'$|\sum m_i \vec{a}_i|$')
    axs[1, 0].set_yscale('log')
    axs[1, 0].set_xlabel(x_label)
    axs[1, 0].set_ylabel(r"Net Force Magnitude $|F_{\mathrm{net}}|$")
    axs[1, 0].set_title("Newton's 3rd Law Check (Net Force)")
    axs[1, 0].grid(True, linestyle=':', alpha=0.6)
    axs[1, 0].legend()

    # Plot 4: Total Momentum Magnitude
    axs[1, 1].plot(x_axis, p_mag, color='crimson', lw=2, label=r'$|\vec{P}_{\mathrm{total}}|$')
    axs[1, 1].set_xlabel(x_label)
    axs[1, 1].set_ylabel(r"Linear Momentum Magnitude $|P|$")
    axs[1, 1].set_title("System Momentum Conservation")
    axs[1, 1].grid(True, linestyle=':', alpha=0.6)
    axs[1, 1].legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Audit and validate gravity and work-energy conservation from simulation snapshots."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to the simulation snapshot directory."
    )
    
    parser.add_argument(
        "-l", "--latest", 
        action="store_true", 
        help="Automatically find and load the most recent 'run_*' directory inside the provided path"
    )

    # Show help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Parse the arguments
    args = parser.parse_args()
    target_dir = args.path

    # Handle the --latest flag
    if args.latest:
        search_pattern = os.path.join(target_dir, "run_*")
        runs = sorted(glob.glob(search_pattern))
        
        if runs:
            target_dir = runs[-1]
            print(f"Auto-selected latest run: {target_dir}")
        else:
            print(f"Error: No run directories found in '{args.path}'")
            sys.exit(1)

    audit_gravity_snapshots(target_dir)