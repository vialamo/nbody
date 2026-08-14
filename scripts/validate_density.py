#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

def cubic_spline_kernel(r, h):
    """
    Evaluates the 3D cubic spline kernel using the GIZMO/MFM convention 
    where the support radius is exactly 1h.
    """
    q = r / h
    norm = 8.0 / (np.pi * h**3)
    W = np.zeros_like(q)
    
    # 0 <= q < 0.5
    mask1 = q < 0.5
    W[mask1] = norm * (1.0 - 6.0 * q[mask1]**2 + 6.0 * q[mask1]**3)
    
    # 0.5 <= q < 1.0
    mask2 = (q >= 0.5) & (q < 1.0)
    W[mask2] = norm * 2.0 * (1.0 - q[mask2])**3
    
    return W

def audit_mfm_snapshots(snapshot_dir, sample_size=5000):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return

    print(f"==================================================================")
    print(f"                 MFM DENSITY SOLVER VALIDATION                    ")
    print(f"==================================================================")
    print(f"Found {len(files)} snapshots. Extracting MFM data...\n")

    # Validate Config
    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        is_expanding = bool(config.get('expanding_universe', 1))
        domain_size = config.get('domain_size', 1.0)
        method = config.get('hydro_method', b"none").decode('utf-8')
        
        if method != "mfm":
            print(f"[ERROR] This script is exclusively for MFM validation. Detected hydro_method: '{method}'")
            return
            
        has_gas = 'Gas' in f

    if not has_gas:
        print("[ERROR] No 'Gas' group found in the HDF5 file.")
        return

    print(f"Mode Auto-detected : {'EXPANDING UNIVERSE' if is_expanding else 'STATIC BOX (a=1.0)'}")
    print(f"Domain Size        : {domain_size}")
    print(f"Max Validation Size: {sample_size} particles per snapshot")
    print("------------------------------------------------------------------\n")

    # Time series arrays
    times = []
    scale_factors = []
    
    mean_frac_errors = []
    max_frac_errors = []
    mean_n_enc_list = []
    std_n_enc_list = []
    
    # To store the final snapshot data for the scatter plot
    final_rho_saved = None
    final_rho_calc = None
    final_n_enc = None

    for i, f_path in enumerate(files):
        sys.stdout.write(f"\rProcessing [{i+1}/{len(files)}] {os.path.basename(f_path)}")
        sys.stdout.flush()

        with h5py.File(f_path, 'r') as f:
            header = f['Header'].attrs
            a = header['scale_factor']
            t = header['simulation_time']
            
            times.append(t)
            scale_factors.append(a)
            
            # Read Gas data
            gas_x = f['Gas/position_x'][:]
            gas_y = f['Gas/position_y'][:]
            gas_z = f['Gas/position_z'][:]
            mass = f['Gas/mass'][:]
            h = f['Gas/smoothing_length'][:]
            rho_saved = f['Gas/density'][:]
            
        N_total = len(mass)
        
        # Build KDTree for periodic neighbor search
        pos = np.vstack((gas_x, gas_y, gas_z)).T
        tree = cKDTree(pos, boxsize=domain_size)
        
        # Subsample to keep processing fast
        if N_total > sample_size:
            np.random.seed(42)  # Deterministic seed
            indices = np.random.choice(N_total, sample_size, replace=False)
        else:
            indices = np.arange(N_total)
            
        rho_calc = np.zeros(len(indices))
        n_enc_calc = np.zeros(len(indices))
        rho_saved_sampled = rho_saved[indices]
        
        for idx_local, idx_global in enumerate(indices):
            h_i = h[idx_global]
            p_i = pos[idx_global]
            m_i = mass[idx_global]
            
            # Find neighbors within support radius h_i
            neighbors = tree.query_ball_point(p_i, r=h_i)
            
            # Calculate periodic distances manually
            delta = np.abs(pos[neighbors] - p_i)
            delta = np.where(delta > 0.5 * domain_size, domain_size - delta, delta)
            r = np.linalg.norm(delta, axis=1)
            
            # Evaluate kernel
            W = cubic_spline_kernel(r, h_i)
            
            # MFM Density: rho_i = m_i * sum(W)
            n_i = np.sum(W)
            rho_calc[idx_local] = m_i * n_i
            
            # Calculate Effective Enclosed Neighbors
            n_enc_calc[idx_local] = (4.0 / 3.0) * np.pi * (h_i**3) * n_i
            
        # Calculate Error Metrics
        frac_error = np.abs(rho_calc - rho_saved_sampled) / rho_saved_sampled
        mean_frac_errors.append(np.mean(frac_error))
        max_frac_errors.append(np.max(frac_error))
        
        mean_n_enc_list.append(np.mean(n_enc_calc))
        std_n_enc_list.append(np.std(n_enc_calc))
        
        # Cache the last snapshot for scatter plotting
        if i == len(files) - 1:
            final_rho_saved = rho_saved_sampled
            final_rho_calc = rho_calc
            final_n_enc = n_enc_calc

    print("\n\nData collection complete. Evaluating diagnostic metrics...")

    # Convert to Numpy Arrays
    times = np.array(times)
    scale_factors = np.array(scale_factors)
    mean_frac_errors = np.array(mean_frac_errors)
    max_frac_errors = np.array(max_frac_errors)
    mean_n_enc_list = np.array(mean_n_enc_list)
    std_n_enc_list = np.array(std_n_enc_list)

    # Print Terminal Report
    global_max_err = np.max(max_frac_errors)
    final_mean_n_enc = mean_n_enc_list[-1]
    final_std_n_enc = std_n_enc_list[-1]

    print(f"------------------------------------------------------------------")
    print(f"                       AUDIT SUMMARY REPORT                       ")
    print(f"------------------------------------------------------------------")
    
    print(f"  Max Density Fractional Error : {global_max_err * 100.0:.6f}% ", end="")
    print("[PASS]" if global_max_err < 1e-4 else "[WARNING: >0.01% drift. Check root-finder or kernel implementation!]")
    
    print(f"  Final Snapshot N_enc (Mean)  : {final_mean_n_enc:.4f}")
    print(f"  Final Snapshot N_enc (Std)   : {final_std_n_enc:.4f} ", end="")
    print("[PASS]" if final_std_n_enc < 0.5 else "[WARNING: High N_enc variance. Root-finder tolerance may be too loose!]")
    print(f"------------------------------------------------------------------\n")

    # Plotting
    x_axis = scale_factors if is_expanding else times
    x_label = "Scale Factor (a)" if is_expanding else "Simulation Time [Code Units]"

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.canvas.manager.set_window_title(f"MFM Density Validation Audit")

    # Plot 1: Fractional Error Evolution
    axs[0, 0].plot(x_axis, mean_frac_errors, color='blue', lw=2, label='Mean Fractional Error')
    axs[0, 0].plot(x_axis, max_frac_errors, color='crimson', lw=2, linestyle='--', label='Max Fractional Error')
    axs[0, 0].set_yscale('log')
    axs[0, 0].set_xlabel(x_label)
    axs[0, 0].set_ylabel(r"Fractional Error $|\rho_{calc} - \rho_{saved}| / \rho_{saved}$")
    axs[0, 0].set_title("Density Reconstruction Error Evolution")
    axs[0, 0].grid(True, linestyle=':', alpha=0.6)
    axs[0, 0].legend()

    # Plot 2: N_enc Evolution
    axs[0, 1].plot(x_axis, mean_n_enc_list, color='darkgreen', lw=2, label='Mean $N_{enc}$')
    axs[0, 1].fill_between(x_axis, mean_n_enc_list - std_n_enc_list, mean_n_enc_list + std_n_enc_list, color='green', alpha=0.2, label='±1 Std Dev')
    axs[0, 1].set_xlabel(x_label)
    axs[0, 1].set_ylabel(r"Effective Enclosed Neighbors ($N_{enc}$)")
    axs[0, 1].set_title(r"Root-Finder Stability: $N_{enc} = \frac{4}{3}\pi h_i^3 n_i$")
    axs[0, 1].grid(True, linestyle=':', alpha=0.6)
    axs[0, 1].legend()

    # Plot 3: Scatter Plot (Final Snapshot)
    axs[1, 0].scatter(final_rho_saved, final_rho_calc, alpha=0.5, s=10, color='purple')
    # Ideal 1:1 line
    min_val = min(np.min(final_rho_saved), np.min(final_rho_calc))
    max_val = max(np.max(final_rho_saved), np.max(final_rho_calc))
    axs[1, 0].plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', lw=1.5, label='Exact Match (1:1)')
    axs[1, 0].set_xscale('log')
    axs[1, 0].set_yscale('log')
    axs[1, 0].set_xlabel(r"Saved Density ($\rho_{saved}$)")
    axs[1, 0].set_ylabel(r"Reconstructed Density ($\rho_{calc}$)")
    axs[1, 0].set_title("Final Snapshot: Density Match")
    axs[1, 0].grid(True, linestyle=':', alpha=0.6)
    axs[1, 0].legend()

    # Plot 4: N_enc Histogram (Final Snapshot)
    axs[1, 1].hist(final_n_enc, bins=50, color='teal', edgecolor='black', alpha=0.7)
    axs[1, 1].axvline(final_mean_n_enc, color='red', linestyle='dashed', linewidth=2, label=f'Mean: {final_mean_n_enc:.2f}')
    axs[1, 1].set_xlabel(r"Effective Enclosed Neighbors ($N_{enc}$)")
    axs[1, 1].set_ylabel("Particle Count")
    axs[1, 1].set_title("Final Snapshot: $N_{enc}$ Distribution")
    axs[1, 1].grid(True, linestyle=':', alpha=0.6)
    axs[1, 1].legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Audit and validate MFM density and effective volume from simulation snapshots."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to the simulation snapshot directory."
    )
    parser.add_argument(
        "-s", "--sample",
        type=int,
        default=4096,
        help="Number of random particles to sample per snapshot for validation (default: 5000)."
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

    audit_mfm_snapshots(target_dir, args.sample)