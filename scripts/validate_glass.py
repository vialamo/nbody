#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt

def validate_glass_generation(snapshot_dir):
    # Load and sort files identically to the reference script
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    #last_file = files[-1]          # Save the absolute last snapshot
    #files = files[::8]             # Take every 8th snapshot
    #if files[-1] != last_file:     # Ensure the last snapshot is included
    #    files.append(last_file)

    # Validate config from the first file
    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        setup_type = config.get('setup', b"").decode('utf-8')
        domain_size = config.get('domain_size', 1.0)

        if setup_type != "glass":
            print(f"[WARNING] Expected config attribute 'setup'='glass'. Found: '{setup_type}'")

    print(f"Processing {len(files)} snapshots to evaluate glass convergence...")

    times = []
    max_velocities = []
    density_covs = [] # Coefficient of Variation (std / mean)
    density_max = []
    density_min = []

    # Process time evolution
    for f_name in files:
        with h5py.File(f_name, 'r') as f:
            t = f['Header'].attrs['simulation_time']
            times.append(t)

            gas = f['Gas']

            # 1. Velocity Metric
            v_x = gas['velocity_x'][:]
            v_y = gas['velocity_y'][:]
            v_z = gas['velocity_z'][:]
            v_mag = np.sqrt(v_x**2 + v_y**2 + v_z**2)
            max_velocities.append(np.max(v_mag))

            # 2. Density Uniformity Metric
            rho = gas['density'][:]
            mean_rho = np.mean(rho)
            std_rho = np.std(rho)

            density_covs.append(std_rho / mean_rho)
            density_max.append(np.max(rho) / mean_rho)
            density_min.append(np.min(rho) / mean_rho)

    # Process Final State (for spatial visualization)
    with h5py.File(files[-1], 'r') as f:
        final_x = f['Gas/position_x'][:]
        final_y = f['Gas/position_y'][:]
        final_z = f['Gas/position_z'][:]
        final_rho = f['Gas/density'][:]
        final_mean_rho = np.mean(final_rho)

    # --- Plotting ---
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("MFM Glass Generation Convergence", fontsize=16)
    plt.subplots_adjust(hspace=0.3, wspace=0.2)

    # Top Left: Maximum Velocity (Log Scale)
    axs[0, 0].plot(times, max_velocities, color='royalblue', lw=2)
    axs[0, 0].set_yscale('log')
    axs[0, 0].set_title("Kinetic Energy Damping (Max Velocity)")
    axs[0, 0].set_xlabel("Simulation Time")
    axs[0, 0].set_ylabel(r"Max Velocity ($|v|_{max}$)")
    axs[0, 0].grid(True, linestyle=':', alpha=0.6)

    # Top Right: Density Uniformity
    axs[0, 1].plot(times, density_max, color='darkorange', lw=1.5, linestyle='--', label='Max / Mean')
    axs[0, 1].plot(times, density_min, color='darkorange', lw=1.5, linestyle='-.', label='Min / Mean')
    axs[0, 1].plot(times, density_covs, color='forestgreen', lw=2, label=r'Std Dev ($\sigma_\rho / \mu_\rho$)')
    axs[0, 1].set_title("Density Convergence")
    axs[0, 1].set_xlabel("Simulation Time")
    axs[0, 1].set_ylabel("Density Ratio")
    axs[0, 1].set_ylim(0, 1.8)
    axs[0, 1].legend()
    axs[0, 1].grid(True, linestyle=':', alpha=0.6)

    # Bottom Left: 2D Spatial Slice (Checking for grid lines or voids)
    slice_thickness = domain_size * 0.1
    slice_mask = (final_z >= 0) & (final_z <= slice_thickness)
    axs[1, 0].scatter(final_x[slice_mask], final_y[slice_mask], s=2, color='black', alpha=0.5)
    axs[1, 0].set_title(f"Final Spatial Slice (Z < {slice_thickness:.2f})")
    axs[1, 0].set_xlabel("X Position")
    axs[1, 0].set_ylabel("Y Position")
    axs[1, 0].set_xlim(0, domain_size)
    axs[1, 0].set_ylim(0, domain_size)
    axs[1, 0].set_aspect('equal')

    # Bottom Right: Final Density Histogram
    axs[1, 1].hist(final_rho / final_mean_rho, bins=50, color='royalblue', alpha=0.7, edgecolor='black')
    axs[1, 1].set_title("Final Density Distribution")
    axs[1, 1].set_xlabel(r"Normalized Density ($\rho / \mu_\rho$)")
    axs[1, 1].set_ylabel("Particle Count")
    axs[1, 1].axvline(1.0, color='red', linestyle='dashed', linewidth=1)
    axs[1, 1].grid(True, linestyle=':', alpha=0.6)

    plt.show()

if __name__ == "__main__":
    # Argument parsing mirrors the reference script
    parser = argparse.ArgumentParser(description="Evaluate MFM Glass Generation convergence.")
    parser.add_argument("path", type=str, nargs='?', help="Path to snapshot directory.")
    parser.add_argument("-l", "--latest", action="store_true", help="Load latest run_* directory")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    target_dir = args.path
    if args.latest:
        # Resolve the latest run directory
        runs = sorted(glob.glob(os.path.join(target_dir if target_dir else ".", "run_*")))
        if runs: target_dir = runs[-1]

    if target_dir is None:
        print("[ERROR] No target directory specified and no run_* directories found.")
        sys.exit(1)

    validate_glass_generation(target_dir)
