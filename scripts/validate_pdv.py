#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt

def validate_adiabatic_expansion_physical(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    # Extract configuration and units from the first snapshot
    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        header = f['Header'].attrs
        units = f['Units'].attrs
        
        gamma = config.get('gamma', 5.0/3.0)
        a_start = header['scale_factor']
        
        # Load CGS conversion factors
        unit_density_cgs = units['unit_density_in_cgs']
        unit_velocity_cgs = units['unit_velocity_in_cgs']
        
        # Validation Checks for Adiabatic Expansion
        if not bool(config.get('expanding_universe', 0)):
            print("[ERROR] Adiabatic test requires an expanding box, but 'expanding_universe' is FALSE.")
            sys.exit(1)
            
        if bool(config.get('enable_cooling', 0)):
            print("[ERROR] Adiabatic test requires adiabatic gas, but 'enable_cooling' is TRUE.")
            sys.exit(1)
            
        # Initial states (Code Units)
        rho_code_0 = np.mean(f['Gas/density'][:])
        P_code_0 = np.mean(f['Gas/pressure'][:])
        T_0 = np.mean(f['Gas/temperature'][:])  # Already converted to Kelvin by HDF5Writer
        
        # Convert to Initial Physical States (CGS)
        rho_phys_0 = (rho_code_0 / (a_start**3)) * unit_density_cgs
        P_phys_0 = (P_code_0 / a_start) * (unit_density_cgs * unit_velocity_cgs**2)

    # Pre-read the evolution history
    history_a = []
    history_rho = []
    history_T = []
    history_P = []
    
    for file in files:
        with h5py.File(file, 'r') as f:
            a_curr = f['Header'].attrs['scale_factor']
            
            # Read current code units]
            rho_code = np.mean(f['Gas/density'][:])
            P_code = np.mean(f['Gas/pressure'][:])
            T_curr = np.mean(f['Gas/temperature'][:]) 
            
            # Convert to physical CGS variables
            rho_phys = (rho_code / (a_curr**3)) * unit_density_cgs
            P_phys = (P_code / a_curr) * (unit_density_cgs * unit_velocity_cgs**2)
            
            history_a.append(a_curr)
            history_rho.append(rho_phys)
            history_T.append(T_curr)
            history_P.append(P_phys)
            
    history_a = np.array(history_a)
    history_rho = np.array(history_rho)
    history_T = np.array(history_T)
    history_P = np.array(history_P)
    
    # Analytical Solutions for Adiabatic Expansion
    # rho ~ a^-3
    # T ~ a^(-3(gamma-1))
    # P ~ a^(-3*gamma)
    a_exact = np.linspace(a_start, max(history_a), 500)
    rho_exact = rho_phys_0 * (a_exact / a_start)**(-3.0)
    T_exact = T_0 * (a_exact / a_start)**(-3.0 * (gamma - 1.0))
    P_exact = P_phys_0 * (a_exact / a_start)**(-3.0 * gamma)

    # Plot Setup
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    plt.subplots_adjust(bottom=0.15, wspace=0.25) 

    sim_kwargs = {'color': 'royalblue', 'lw': 3, 'label': 'MFM Simulation'}
    exact_kwargs = {'color': 'black', 'lw': 1.5, 'linestyle': '--', 'label': 'Exact Solution'}
    
    # 1. Density Panel
    axs[0].plot(a_exact, rho_exact, **exact_kwargs)
    axs[0].plot(history_a, history_rho, **sim_kwargs)
    axs[0].set_ylabel(r"Physical Density [g/cm$^3$]")
    axs[0].set_xlabel("Scale Factor (a)")
    axs[0].set_title(r"Density Evolution ($\rho \propto a^{-3}$)")
    
    # 2. Temperature Panel
    axs[1].plot(a_exact, T_exact, **exact_kwargs)
    axs[1].plot(history_a, history_T, **sim_kwargs)
    axs[1].set_ylabel(r"Temperature [K]")
    axs[1].set_xlabel("Scale Factor (a)")
    decay_exp = -3.0 * (gamma - 1.0)
    axs[1].set_title(f"Temperature Evolution ($T \\propto a^{{{decay_exp:.2f}}}$)")
    
    # 3. Pressure Panel
    axs[2].plot(a_exact, P_exact, **exact_kwargs)
    axs[2].plot(history_a, history_P, **sim_kwargs)
    axs[2].set_ylabel(r"Physical Pressure [dyne/cm$^2$]")
    axs[2].set_xlabel("Scale Factor (a)")
    axs[2].set_title(f"Pressure Evolution ($P \\propto a^{{-3\\gamma}}$)")
    
    for ax in axs:
        ax.set_yscale('log')
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='upper right', fontsize=10)

    fig.suptitle("Adiabatic Expansion Validation (Physical Units)", fontsize=16, fontweight='bold')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MFM Adiabatic Expansion physical validation.")
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
        
    validate_adiabatic_expansion_physical(target_dir)