import h5py
import numpy as np
import glob
import os
import sys

# --- Terminal Color Codes ---
GREEN = '\033[92m'
RED = '\033[91m'
BLUE = '\033[94m'
RESET = '\033[0m'

def run_check(condition, pass_msg, fail_msg):
    """Evaluates a condition and prints formatted colored output."""
    if condition:
        print(f"    {GREEN}[PASS]{RESET} {pass_msg}")
    else:
        print(f"    {RED}[FAIL]{RESET} {fail_msg}")
        sys.exit(1)  # Gracefully stop the script without a traceback

def validate_snapshot(file_path, initial_mass=None):
    """Runs physical validations on a single HDF5 snapshot."""
    with h5py.File(file_path, 'r') as f:
        # --- 1. Read Root Attributes ---
        domain_size = f.attrs['domain_size']
        mesh_size = f.attrs['mesh_size']
        use_hydro = bool(f.attrs['use_hydro'])
        sim_time = f.attrs['simulation_time']
        scale_factor = f.attrs['scale_factor']
        
        print(f"{BLUE}>> Validating {os.path.basename(file_path)} {RESET}| t={sim_time:.4f} | a={scale_factor:.4f}")
        
        # --- 2. Load Particle Data ---
        p_x = f['particles/position_x'][:]
        p_y = f['particles/position_y'][:]
        p_z = f['particles/position_z'][:]
        p_mass = f['particles/mass'][:]
        
        total_mass = np.sum(p_mass)
        
        # Particle checks: NaNs
        run_check(not np.isnan(p_x).any(), "Particle X positions valid (No NaNs)", "NaNs detected in particle X positions!")
        run_check(not np.isnan(p_y).any(), "Particle Y positions valid (No NaNs)", "NaNs detected in particle Y positions!")
        run_check(not np.isnan(p_z).any(), "Particle Z positions valid (No NaNs)", "NaNs detected in particle Z positions!")
        
        # Particle checks: Periodic Boundaries
        out_of_bounds_x = (p_x < 0.0) | (p_x >= domain_size)
        out_of_bounds_y = (p_y < 0.0) | (p_y >= domain_size)
        out_of_bounds_z = (p_z < 0.0) | (p_z >= domain_size)
        
        total_escaped = np.sum(out_of_bounds_x | out_of_bounds_y | out_of_bounds_z)
        run_check(total_escaped == 0, "All particles strictly within periodic bounds", f"{total_escaped} particles escaped!")

        # --- 3. Conditionally Load and Check Gas Data ---
        if use_hydro:
            g_rho = f['gas/density'][:]
            g_eng = f['gas/energy'][:]
            
            cell_size = domain_size / mesh_size
            cell_volume = cell_size ** 3
            
            # Gas checks
            run_check(not np.isnan(g_rho).any(), "Gas density array valid (No NaNs)", "NaNs detected in gas density!")
            run_check(np.all(g_rho > 0), "Gas density strictly positive", "Negative gas density detected (Safety floor failed)!")
            run_check(np.all(g_eng > 0), "Gas energy strictly positive", "Negative gas energy detected (Safety floor failed)!")
            
            gas_mass = np.sum(g_rho) * cell_volume
            total_mass += gas_mass
            print(f"    {BLUE}[INFO]{RESET} Max Gas Density: {np.max(g_rho):.5e}")
        else:
            print(f"    {BLUE}[INFO]{RESET} Hydro disabled (Dark Matter only).")

        # --- 4. Conservation Checks ---
        if initial_mass is not None:
            mass_drift = abs(total_mass - initial_mass) / initial_mass
            run_check(mass_drift < 1e-10, f"Mass perfectly conserved (drift: {mass_drift:.2e})", f"Mass drifted by {mass_drift * 100}%")
            print(f"    {BLUE}[INFO]{RESET} Total Mass: {total_mass:.5e}")
        else:
            print(f"    {BLUE}[INFO]{RESET} Baseline Total Mass established: {total_mass:.5e}")
            
        print("") # Blank line for readability
        return total_mass

def run_validation_suite(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    
    if not files:
        print(f"{RED}Error: No HDF5 snapshots found in {snapshot_dir}!{RESET}")
        return
        
    print(f"\n{BLUE}=== Beginning Physics Validation Suite ({len(files)} snapshots) ==={RESET}\n")
    
    # Validate the first snapshot to establish the baseline total mass
    initial_mass = validate_snapshot(files[0])
    
    # Loop through the rest and ensure mass is strictly conserved
    for f in files[1:]:
        validate_snapshot(f, initial_mass=initial_mass)
        
    print(f"{GREEN}===================================================={RESET}")
    print(f"{GREEN}SUCCESS: All snapshots passed physical validation!{RESET}")
    print(f"{GREEN}===================================================={RESET}\n")

if __name__ == "__main__":
    snapshot_dir = sys.argv[1] if len(sys.argv) > 1 else "./output/"
    
    if not os.path.exists(snapshot_dir):
        print(f"{RED}Error: Directory '{snapshot_dir}' does not exist.{RESET}")
        sys.exit(1)
        
    run_validation_suite(snapshot_dir)