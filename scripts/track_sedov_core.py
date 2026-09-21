#!/usr/bin/env python3
import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt

def find_target_particles(f, num_neighbors):
    """Finds the central particle and a few closest neighbors along the +x axis."""
    config = f['Config'].attrs
    domain_size = config.get('domain_size', 1.0)
    center = domain_size / 2.0
    
    x = f['Gas/position_x'][:]
    y = f['Gas/position_y'][:]
    z = f['Gas/position_z'][:]
    
    # Distances from center
    dx = x - center
    dy = y - center
    dz = z - center
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    
    # 1. Find Central Particle
    center_idx = np.argmin(r)
    
    # 2. Find +x Neighbors (Cone search to handle Glass ICs)
    # Avoid the exact center particle returning an invalid angle
    valid = r > 1e-6
    
    # Calculate angle relative to the +x axis (cos(theta) = dx / r)
    cos_theta = np.zeros_like(r)
    cos_theta[valid] = dx[valid] / r[valid]
    
    # Select particles tightly aligned with the +x axis (cos(theta) > 0.98 is ~11 degrees)
    cone_mask = valid & (cos_theta > 0.98)
    
    # Get indices of particles in the cone, sorted by distance from the center
    cone_indices = np.where(cone_mask)[0]
    cone_indices_sorted = cone_indices[np.argsort(r[cone_indices])]
    
    # Select the closest `num_neighbors`
    neighbor_indices = cone_indices_sorted[:num_neighbors]
    
    return [center_idx] + neighbor_indices.tolist()

def track_core_evolution(snapshot_dir, num_snaps, num_neighbors):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"[ERROR] No HDF5 snapshots found in directory: '{snapshot_dir}'")
        return
        
    # Limit to the requested number of early snapshots
    files = files[:num_snaps]
    
    times = []
    
    # We will track a list of particles. Index 0 is the center, 1..N are the +x neighbors.
    # Data structure: history[particle_id][property] = list_of_values
    num_tracked = 1 + num_neighbors
    history = {i: {'x': [], 'rho': [], 'P': [], 'u': [], 'vx': [], 'h': []} for i in range(num_tracked)}
    
    # To track particles across resorted arrays, we remember their last known positions
    last_positions = np.zeros((num_tracked, 3))

    for step, fname in enumerate(files):
        with h5py.File(fname, 'r') as f:
            t = f['Header'].attrs['simulation_time']
            times.append(t)
            
            x = f['Gas/position_x'][:]
            y = f['Gas/position_y'][:]
            z = f['Gas/position_z'][:]
            rho = f['Gas/density'][:]
            vx = f['Gas/velocity_x'][:]
            u = f['Gas/internal_energy'][:]
            h = f['Gas/smoothing_length'][:]
            gamma = f['Config'].attrs.get('gamma', 5.0/3.0)
            P = rho * u * (gamma - 1.0)
            
            current_positions = np.c_[x, y, z]
            
            if step == 0:
                # Initialization: Find targets based on geometry
                target_indices = find_target_particles(f, num_neighbors)
                for i, idx in enumerate(target_indices):
                    last_positions[i] = current_positions[idx]
            else:
                # Tracking: Find nearest particle to the last known position
                target_indices = []
                for i in range(num_tracked):
                    dist = np.linalg.norm(current_positions - last_positions[i], axis=1)
                    idx = np.argmin(dist)
                    target_indices.append(idx)
                    last_positions[i] = current_positions[idx] # Update expected position
            
            # Record data
            for i, idx in enumerate(target_indices):
                history[i]['x'].append(x[idx])
                history[i]['rho'].append(rho[idx])
                history[i]['P'].append(P[idx])
                history[i]['u'].append(u[idx])
                history[i]['vx'].append(vx[idx])
                history[i]['h'].append(h[idx])

    # Convert to numpy arrays for easier plotting
    for i in range(num_tracked):
        for key in history[i]:
            history[i][key] = np.array(history[i][key])
            
    plot_diagnostics(np.array(times), history, num_neighbors)

def plot_diagnostics(times, history, num_neighbors):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("Lagrangian Core Tracker (Early Startup Phase)", fontsize=16)
    
    # Styling
    center_style = {'color': 'red', 'linewidth': 2.5, 'label': 'Center Particle', 'zorder': 10}
    
    # Generate a colormap for the neighbors (getting lighter as they get further away)
    cmap = plt.get_cmap('Blues_r')
    colors = [cmap(0.2 + 0.6 * (i / max(1, num_neighbors-1))) for i in range(num_neighbors)]

    # Map variables to axes
    plots = [
        (axes[0,0], 'x', r'X Position ($x$)', 'linear'),
        (axes[0,1], 'vx', r'X Velocity ($v_x$)', 'linear'),
        (axes[0,2], 'h', r'Smoothing Length ($h$)', 'linear'),
        (axes[1,0], 'rho', r'Density ($\rho$)', 'log'),
        (axes[1,1], 'P', r'Pressure ($P$)', 'log'),
        (axes[1,2], 'u', r'Internal Energy ($u$)', 'log')
    ]
    
    for ax, var, ylabel, scale in plots:
        # Plot center particle
        ax.plot(times, history[0][var], **center_style)
        
        # Plot +x neighbor particles
        for i in range(1, num_neighbors + 1):
            ax.plot(times, history[i][var], color=colors[i-1], 
                    linewidth=1.5, label=f'Neighbor {i}' if var == 'x' else "")
            
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Time (t)")
        ax.set_yscale(scale)
        ax.grid(True, linestyle=':', alpha=0.6)
        
    # Only add legend to the first plot to save space
    axes[0,0].legend()
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Track Lagrangian core particles in early Sedov blastwave.")
    parser.add_argument("path", type=str, help="Path to snapshot directory.")
    parser.add_argument("-l", "--latest", action="store_true", help="Load latest run_* directory")
    parser.add_argument("-n", "--num_snaps", type=int, default=30, help="Number of early snapshots to read.")
    parser.add_argument("-k", "--neighbors", type=int, default=5, help="Number of +x neighbors to track.")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    target_dir = args.path
    
    if args.latest:
        runs = sorted(glob.glob(os.path.join(target_dir, "run_*")))
        if runs: 
            target_dir = runs[-1]
            print(f"Auto-selected latest run: {target_dir}")
            
    track_core_evolution(target_dir, args.num_snaps, args.neighbors)