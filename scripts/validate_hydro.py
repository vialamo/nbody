#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt

def audit_hydro_snapshots(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    print(f"==================================================================")
    print(f"               PURE HYDRODYNAMICS CONSERVATION AUDIT              ")
    print(f"==================================================================")
    print(f"Found {len(files)} snapshots. Extracting gas data...\n")

    # Read configuration from the first snapshot
    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        is_expanding = bool(config.get('expanding_universe', 0))
        use_pm = bool(config.get('use_pm', 0))
        use_pp = bool(config.get('use_pp', 0))
        method = config.get('hydro_method', b"none").decode('utf-8')
        gamma = config.get('gamma', 1.6666666)
        has_gas = 'Gas' in f

    if not has_gas:
        print("[ERROR] No 'Gas' group found in HDF5. Is hydro enabled?")
        sys.exit(1)

    print(f"Hydro Method  : {method.upper()}")
    print(f"Adiabatic Idx : {gamma:.3f}")
    print(f"Mode          : {'EXPANDING UNIVERSE' if is_expanding else 'STATIC BOX'}")
    print(f"Gravity       : {'ON (Warning: This audit assumes pure hydro)' if (use_pm or use_pp) else 'OFF (Pure Hydro)'}")
    print("------------------------------------------------------------------\n")

    # Time series arrays
    times = []
    scale_factors = []
    
    # Conservation Arrays
    total_mass = []
    ke_gas = []
    ie_gas = []
    te_gas = []
    
    momentum_x, momentum_y, momentum_z = [], [], []
    
    # Cumulative Work Arrays
    w_grav_gas_list = []
    w_exp_gas_list = []
    
    for i, f_path in enumerate(files):
        sys.stdout.write(f"\rProcessing [{i+1}/{len(files)}] {os.path.basename(f_path)}")
        sys.stdout.flush()

        with h5py.File(f_path, 'r') as f:
            header = f['Header'].attrs
            times.append(header['simulation_time'])
            scale_factors.append(header['scale_factor'])
            
            # Extract MFM arrays
            m  = f['Gas/mass'][:]
            vx = f['Gas/velocity_x'][:]
            vy = f['Gas/velocity_y'][:]
            vz = f['Gas/velocity_z'][:]
            u  = f['Gas/internal_energy'][:]
            
            # Extract cumulative work attributes
            w_grav_gas_list.append(f['Gas'].attrs.get('cumulative_gravitational_work', 0.0))
            w_exp_gas_list.append(f['Gas'].attrs.get('cumulative_expansion_work', 0.0))
            
            # Computations
            mass_tot = np.sum(m)
            px = np.sum(m * vx)
            py = np.sum(m * vy)
            pz = np.sum(m * vz)
            
            k_step = 0.5 * np.sum(m * (vx**2 + vy**2 + vz**2))
            u_step = np.sum(m * u)
            
            total_mass.append(mass_tot)
            momentum_x.append(px)
            momentum_y.append(py)
            momentum_z.append(pz)
            ke_gas.append(k_step)
            ie_gas.append(u_step)
            te_gas.append(k_step + u_step)

    print("\n\nData collection complete. Evaluating diagnostic metrics...")

    # Convert to Numpy Arrays
    times = np.array(times)
    scale_factors = np.array(scale_factors)
    total_mass = np.array(total_mass)
    ke_gas = np.array(ke_gas)
    ie_gas = np.array(ie_gas)
    te_gas = np.array(te_gas)
    
    w_grav_gas = np.array(w_grav_gas_list)
    w_exp_gas = np.array(w_exp_gas_list)
    
    p_mag = np.sqrt(np.array(momentum_x)**2 + np.array(momentum_y)**2 + np.array(momentum_z)**2)

    # Fractional Errors relative to Step 0
    err_mass = (total_mass - total_mass[0]) / total_mass[0]
    
    # Correct Energy Error using accumulated cosmological expansion and gravitational work
    e0_gas = te_gas[0] if te_gas[0] > 0 else 1.0
    delta_e_gas = te_gas - te_gas[0]
    err_energy = (delta_e_gas - w_grav_gas + w_exp_gas) / e0_gas
    
    # For momentum, use absolute magnitude drift if starting at 0
    p0 = p_mag[0] if p_mag[0] > 1e-10 else 1.0
    err_mom = (p_mag - p_mag[0]) / p0

    # Print Report
    max_err_mass = np.max(np.abs(err_mass))
    max_err_energy = np.max(np.abs(err_energy))
    max_err_mom = np.max(np.abs(err_mom))

    print(f"------------------------------------------------------------------")
    print(f"                       HYDRO AUDIT SUMMARY                        ")
    print(f"------------------------------------------------------------------")
    print(f"  Max Mass Drift       : {max_err_mass * 100.0:.4e}% ", end="")
    print("[PASS]" if max_err_mass < 1e-10 else "[FAIL]")
    
    print(f"  Max Momentum Drift   : {max_err_mom * 100.0:.4e}% ", end="")
    print("[PASS]" if max_err_mom < 1e-10 else "[FAIL]")

    print(f"  Max Total Energy Err : {max_err_energy * 100.0:.6f}% ", end="")
    print("[PASS]" if max_err_energy < 0.05 else "[FAIL: >0.05% drift. Check PdV work / fluxes!]")
    print(f"------------------------------------------------------------------\n")

    # Plotting 
    x_axis = scale_factors if is_expanding else times
    x_label = "Scale Factor (a)" if is_expanding else "Simulation Time [Code Units]"

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.canvas.manager.set_window_title(f"Hydrodynamics Audit - {'Expanding' if is_expanding else 'Static'}")

    # Plot 1: Energy Partitioning
    axs[0, 0].plot(x_axis, ke_gas, color='blue', lw=2, label='Kinetic Energy')
    axs[0, 0].plot(x_axis, ie_gas, color='red', lw=2, label='Internal Energy')
    axs[0, 0].plot(x_axis, te_gas, color='green', lw=2, linestyle='--', label='Total Energy')
    axs[0, 0].set_xlabel(x_label)
    axs[0, 0].set_ylabel("Energy")
    axs[0, 0].set_title("Energy Conversion (Shock Dissipation)")
    axs[0, 0].grid(True, linestyle=':', alpha=0.6)
    axs[0, 0].legend()

    # Plot 2: Total Energy Fractional Error
    axs[0, 1].plot(x_axis, err_energy * 100.0, color='darkgreen', lw=2, label=r'$\frac{\Delta E_{tot} - W_{grav} + W_{exp}}{E_0}$')
    axs[0, 1].axhline(0, color='black', linestyle='-', linewidth=0.8)
    axs[0, 1].axhline(0.05, color='red', linestyle=':', alpha=0.5, label='±0.05% Bound')
    axs[0, 1].axhline(-0.05, color='red', linestyle=':', alpha=0.5)
    axs[0, 1].set_xlabel(x_label)
    axs[0, 1].set_ylabel("Total Energy Error (%)")
    axs[0, 1].set_title("Corrected Total Energy Conservation")
    
    # Auto-scale Y-axis to fit the error
    y_limit = max(0.1, np.max(np.abs(err_energy*100))*1.2)
    axs[0, 1].set_ylim(-y_limit, y_limit)
    axs[0, 1].grid(True, linestyle=':', alpha=0.6)
    axs[0, 1].legend()

    # Plot 3: Momentum Drift
    axs[1, 0].plot(x_axis, p_mag, color='purple', lw=2, label=r'$|P_{tot}|$')
    axs[1, 0].set_xlabel(x_label)
    axs[1, 0].set_ylabel(r"Total Momentum Magnitude $|P|$")
    axs[1, 0].set_title("Linear Momentum Conservation")
    axs[1, 0].grid(True, linestyle=':', alpha=0.6)
    axs[1, 0].legend()

    # Plot 4: Kinetic / Internal Phase Space of Final Snapshot
    with h5py.File(files[-1], 'r') as f:
        rho = f['Gas/density'][:]
        u = f['Gas/internal_energy'][:]
        final_time = f['Header'].attrs['simulation_time']
        final_a = f['Header'].attrs['scale_factor']
        
    axs[1, 1].scatter(rho, u, s=1, color='crimson', alpha=0.5)
    axs[1, 1].set_xscale('log')
    axs[1, 1].set_yscale('log')
    axs[1, 1].set_xlabel("Density (rho)")
    axs[1, 1].set_ylabel("Internal Energy (u)")
    plot_title = f"Final Phase Plot (a={final_a:.3f})" if is_expanding else f"Final Phase Plot (t={final_time:.3f})"
    axs[1, 1].set_title(plot_title)
    axs[1, 1].grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Audit pure hydro conservation from HDF5 snapshots.")
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
    audit_hydro_snapshots(target_dir)