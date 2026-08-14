#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt

def audit_mfm_temperature(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    print(f"==================================================================")
    print(f"              MFM ADIABATIC THERMODYNAMICS VALIDATION             ")
    print(f"==================================================================")
    print(f"Found {len(files)} snapshots. Extracting Temperature data...\n")

    # Validate Config
    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        is_expanding = bool(config.get('expanding_universe', 1))
        
        # Hydro method is stored as a string, handle byte decoding safely
        method_raw = config.get('hydro_method', b"none")
        method = method_raw.decode('utf-8').lower() if isinstance(method_raw, bytes) else method_raw.lower()
        
        if method != "mfm":
            print(f"[ERROR] This script is exclusively for MFM validation. Detected hydro_method: '{method}'")
            return
            
        has_gas = 'Gas' in f
        gamma = config.get('gamma', 1.6666667)
        temp_floor_k = 1.0 #config.get('temp_floor_k', 1.0)

    if not has_gas:
        print("[ERROR] No 'Gas' group found in the HDF5 file.")
        return

    print(f"Mode Auto-detected : {'EXPANDING UNIVERSE' if is_expanding else 'STATIC BOX (a=1.0)'}")
    print(f"Adiabatic Index (γ): {gamma:.4f}")
    print(f"Temperature Floor  : {temp_floor_k:.2f} K")
    print("------------------------------------------------------------------\n")

    # Time series arrays
    times = []
    scale_factors = []
    
    max_temps = []
    mean_temps = []
    min_temps = []

    for i, f_path in enumerate(files):
        sys.stdout.write(f"\rProcessing [{i+1}/{len(files)}] {os.path.basename(f_path)}")
        sys.stdout.flush()

        with h5py.File(f_path, 'r') as f:
            header = f['Header'].attrs
            a = header['scale_factor']
            t = header['simulation_time']
            
            times.append(t)
            scale_factors.append(a)
            
            # Read Gas Temperature
            T_arr = f['Gas/temperature'][:]
            
            max_temps.append(np.max(T_arr))
            mean_temps.append(np.mean(T_arr))
            min_temps.append(np.min(T_arr))

    print("\n\nData collection complete. Evaluating diagnostic metrics...")

    # Convert to Numpy Arrays
    times = np.array(times)
    scale_factors = np.array(scale_factors)
    max_temps = np.array(max_temps)
    mean_temps = np.array(mean_temps)
    min_temps = np.array(min_temps)

    # Theoretical Curve Calculation
    a_start = scale_factors[0]
    T_start = max_temps[0]
    
    # Adiabatic expansion: T scales as a^(-3(gamma - 1))
    decay_power = -3.0 * (gamma - 1.0)
    
    if is_expanding:
        T_theory = T_start * (scale_factors / a_start)**decay_power
    else:
        # In a static box with no cooling/heating, temperature should remain constant
        T_theory = np.full_like(times, T_start)
        
    # Apply the thermodynamic floor
    T_theory = np.clip(T_theory, temp_floor_k, None)
    
    # Calculate Fractional Error
    # We use max_temps here assuming a uniform initial temperature field for the test
    frac_error = np.abs(max_temps - T_theory) / T_theory

    # Print Report
    max_err = np.max(frac_error)
    final_T = max_temps[-1]

    print(f"------------------------------------------------------------------")
    print(f"                       AUDIT SUMMARY REPORT                       ")
    print(f"------------------------------------------------------------------")
    
    print(f"  Max Temperature Deviation Error : {max_err * 100.0:.6f}% ", end="")
    print("[PASS]" if max_err < 1e-3 else "[WARNING: >0.1% drift from theoretical adiabatic curve!]")
    
    print(f"  Final Max Temperature           : {final_T:.4f} K")
    print(f"  Expected Theoretical Floor      : {temp_floor_k:.4f} K ", end="")
    
    floor_hit = np.isclose(final_T, temp_floor_k, rtol=1e-3)
    if is_expanding and T_theory[-1] == temp_floor_k:
        print("[PASS]" if floor_hit else "[FAIL: Floor not maintained]")
    else:
        print("[-]")
        
    print(f"------------------------------------------------------------------\n")

    # Plotting
    x_axis = scale_factors if is_expanding else times
    x_label = "Scale Factor (a)" if is_expanding else "Simulation Time [Code Units]"

    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    fig.canvas.manager.set_window_title(f"MFM Thermodynamics Validation Audit")

    # Plot 1: Temperature Evolution
    axs[0].plot(x_axis, max_temps, color='crimson', lw=3, label='Simulation (Max Temp)')
    axs[0].plot(x_axis, mean_temps, color='orange', lw=2, linestyle='-.', label='Simulation (Mean Temp)')
    axs[0].plot(x_axis, T_theory, color='black', lw=2, linestyle='--', label=f'Theory ($T \propto a^{{{decay_power:.2f}}}$)')
    
    axs[0].axhline(temp_floor_k, color='blue', linestyle=':', lw=2, label=f'Floor ({temp_floor_k} K)')
    
    axs[0].set_yscale('log')
    if is_expanding:
        axs[0].set_xscale('log')
        
    axs[0].set_xlabel(x_label, fontsize=12)
    axs[0].set_ylabel("Temperature [K]", fontsize=12)
    axs[0].set_title("Adiabatic Cooling Evolution", fontsize=14)
    axs[0].grid(True, which='both', linestyle=':', alpha=0.6)
    axs[0].legend(fontsize=11)

    # Plot 2: Fractional Error Evolution
    axs[1].plot(x_axis, frac_error, color='purple', lw=2, label='Relative Error')
    
    axs[1].axhline(0.0, color='black', linestyle='-', linewidth=0.8)
    axs[1].axhline(1e-3, color='red', linestyle=':', alpha=0.5, label='0.1% Threshold')
    
    if is_expanding:
        axs[1].set_xscale('log')
        
    # Scale Y-axis to show small errors cleanly, or expand if error is large
    axs[1].set_ylim(-1e-4, max(2e-3, np.max(frac_error) * 1.2))
    
    axs[1].set_xlabel(x_label, fontsize=12)
    axs[1].set_ylabel(r"Fractional Error $|T_{sim} - T_{theory}| / T_{theory}$", fontsize=12)
    axs[1].set_title("Deviation from Theoretical Expansion", fontsize=14)
    axs[1].grid(True, linestyle=':', alpha=0.6)
    axs[1].legend(fontsize=11)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Audit and validate MFM adiabatic temperature expansion from simulation snapshots."
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

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    target_dir = args.path

    if args.latest:
        search_pattern = os.path.join(target_dir, "run_*")
        runs = sorted(glob.glob(search_pattern))
        
        if runs:
            target_dir = runs[-1]
            print(f"Auto-selected latest run: {target_dir}")
        else:
            print(f"Error: No run directories found in '{args.path}'")
            sys.exit(1)

    audit_mfm_temperature(target_dir)