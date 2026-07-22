import h5py
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt
from numba import njit, prange

# Terminal Color Codes
GREEN = '\033[92m'
RED = '\033[91m'
BLUE = '\033[94m'
RESET = '\033[0m'

def get_hubble_param(a, config):
    """Reconstructs the code-unit Hubble parameter at a given scale factor."""
    if not bool(config.get('expanding_universe', 1)):
        return 0.0
        
    om = config['omega_M']
    ol = config['omega_lambda']
    H0_code = 2.0 / (3.0 * np.sqrt(om))
    
    if ol == 0.0:
        return H0_code * np.sqrt(om / (a**3))
    else:
        return H0_code * np.sqrt(om / (a**3) + ol)

@njit(parallel=True, fastmath=True)
def compute_potential_energy(pos, mass, box_size, G, eps_sq):
    """
    Computes the EXACT O(N^2) Particle-Particle potential energy.
    Compiled to multithreaded C-code via Numba.
    """
    N = len(mass)
    W_sum = 0.0
    half_box = box_size / 2.0
    
    # prange acts like #pragma omp parallel for
    for i in prange(N - 1):
        w_local = 0.0
        p1x = pos[i, 0]
        p1y = pos[i, 1]
        p1z = pos[i, 2]
        m1 = mass[i]
        
        for j in range(i + 1, N):
            dx = abs(pos[j, 0] - p1x)
            dy = abs(pos[j, 1] - p1y)
            dz = abs(pos[j, 2] - p1z)
            
            # Minimum image convention (Fast branching)
            if dx > half_box: dx = box_size - dx
            if dy > half_box: dy = box_size - dy
            if dz > half_box: dz = box_size - dz
            
            dist_sq = dx*dx + dy*dy + dz*dz
            
            w_local += mass[j] / np.sqrt(dist_sq + eps_sq)
            
        W_sum += m1 * w_local
        
    return -G * W_sum

def verify_energy(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    
    if not files:
        print(f"{RED}Error: No HDF5 snapshots found in {snapshot_dir}!{RESET}")
        return

    print(f"\n{BLUE}=== Offline Exact N-Body Energy Verification ==={RESET}")
    print(f"Loading {len(files)} snapshots...")

    # Data arrays
    times = []
    scale_factors = []
    K_phys_list = []
    W_phys_list = []
    H_list = []
    
    # Read base config from first file
    with h5py.File(files[0], 'r') as f:
        config = dict(f['Config'].attrs)
        
        if bool(config.get('use_hydro', 0)):
            print(f"{RED}Warning: Simulation contains Gas. This script is designed for pure DM runs!{RESET}")
            
        G = config['G']
        L = config['domain_size']
        
        # Reconstruct base comoving softening
        num_parts = config['num_particles_1d']**3
        mean_spacing = L / (num_parts**(1/3.0))
        base_softening = mean_spacing * config['comoving_softening_factor']
        cap_a = config['softening_cap_scale_factor']
        expanding = bool(config.get('expanding_universe', 1))

    # Process all snapshots
    for i, fpath in enumerate(files):
        sys.stdout.write(f"\rProcessing Snapshot {i+1}/{len(files)}")
        sys.stdout.flush()
        
        with h5py.File(fpath, 'r') as f:
            a = f['Header'].attrs['scale_factor']
            t = f['Header'].attrs['simulation_time']
            
            pos = np.c_[f['Particles/position_x'][:], 
                        f['Particles/position_y'][:], 
                        f['Particles/position_z'][:]]
            vel = np.c_[f['Particles/velocity_x'][:], 
                        f['Particles/velocity_y'][:], 
                        f['Particles/velocity_z'][:]]
            mass = f['Particles/mass'][:]
            
            # Update Softening as integrator.cpp does
            if a > cap_a:
                eps_comoving = base_softening * (cap_a / a)
            else:
                eps_comoving = base_softening
            eps_sq = eps_comoving**2
            
            # Kinetic Energy
            v_sq = np.sum(vel**2, axis=1)
            K_code = 0.5 * np.sum(mass * v_sq)
            K_phys = K_code * (a**2)
            
            # Potential Energy (O(N^2))
            W_code = compute_potential_energy(pos, mass, L, G, eps_sq)
            W_phys = W_code / a
            
            # Store values
            times.append(t)
            scale_factors.append(a)
            K_phys_list.append(K_phys)
            W_phys_list.append(W_phys)
            H_list.append(get_hubble_param(a, config))

    print("\n\nIntegrating Expansion Work...")
    
    # Convert to numpy arrays for vectorized math
    times = np.array(times)
    scale_factors = np.array(scale_factors)
    K_phys = np.array(K_phys_list)
    W_phys = np.array(W_phys_list)
    E_phys = K_phys + W_phys
    H_vals = np.array(H_list)
    
    # Integrate Expansion Work (Trapezoidal Rule over the snapshots)
    # dW = H * (2K + W) * dt
    expansion_work = np.zeros_like(E_phys)
    expansion_power = H_vals * (2.0 * K_phys + W_phys)
    
    for i in range(1, len(times)):
        dt = times[i] - times[i-1]
        avg_power = 0.5 * (expansion_power[i] + expansion_power[i-1])
        expansion_work[i] = expansion_work[i-1] + avg_power * dt
        
    # The Invariant Layzer-Irvine Energy
    invariant_energy = E_phys + expansion_work
    
    # Calculate fractional error
    E0 = invariant_energy[0]
    fractional_error = (invariant_energy - E0) / abs(E0)
    
    print(f"\n{GREEN}=== Exact Energy Conservation Results ==={RESET}")
    print(f"Mode: {'Expanding Universe' if expanding else 'Static Universe'}")
    print(f"Initial Energy (E0): {E0:.5e}")
    print(f"Final Energy Error : {fractional_error[-1]*100:.5f} %")
    
    # =========================================================================
    # Plotting
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f"Offline Exact N-Body Energy Conservation\n{os.path.basename(os.path.normpath(snapshot_dir))}", fontsize=14)
    
    x_data = scale_factors if expanding else times
    x_label = "Scale Factor (a)" if expanding else "Simulation Time (Code Units)"
    
    # Plot Energy Inventory
    ax1.plot(x_data, K_phys, label="Kinetic Energy ($K$)", color="green", lw=2)
    ax1.plot(x_data, W_phys, label="Potential Energy ($W$)", color="red", lw=2)
    ax1.plot(x_data, E_phys, label="Physical Energy ($E = K+W$)", color="orange", lw=2, linestyle='--')
    ax1.plot(x_data, invariant_energy, label="Invariant Energy ($E + W_{exp}$)", color="black", lw=2)
    
    ax1.set_xlabel(x_label)
    ax1.set_ylabel("Energy (Code Units)")
    ax1.set_title("Energy Inventory")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot Fractional Error
    ax2.plot(x_data, fractional_error * 100.0, color="black", lw=2)
    ax2.axhline(0, color="gray", linestyle="--")
    ax2.set_xlabel(x_label)
    ax2.set_ylabel("Fractional Error [%]")
    ax2.set_title("Global Energy Error (Layzer-Irvine Invariant)")
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    snapshot_dir = sys.argv[1] if len(sys.argv) > 1 else "./output/"
    
    if not os.path.exists(snapshot_dir):
        print(f"{RED}Error: Directory '{snapshot_dir}' does not exist.{RESET}")
        sys.exit(1)
        
    verify_energy(snapshot_dir)