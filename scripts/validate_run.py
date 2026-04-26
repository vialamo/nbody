import h5py
import numpy as np
import glob
import os
import sys

# Terminal Color Codes
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
        sys.exit(1)  # Gracefully stop the script

def compute_cic_variance(p_x, p_y, p_z, p_mass, domain_size, mesh_size):
    """Calculates density variance using Cloud-In-Cell (CIC) mass assignment."""
    N = mesh_size
    grid = np.zeros((N, N, N), dtype=np.float64)
    
    # Scale continuous positions to grid indices
    x = (p_x / domain_size) * N
    y = (p_y / domain_size) * N
    z = (p_z / domain_size) * N
    
    # Integer cell indices
    ix = np.floor(x).astype(int) % N
    iy = np.floor(y).astype(int) % N
    iz = np.floor(z).astype(int) % N
    
    # Fractional distances
    fx = x - np.floor(x)
    fy = y - np.floor(y)
    fz = z - np.floor(z)
    
    # Neighboring cell indices (periodic)
    ix1 = (ix + 1) % N
    iy1 = (iy + 1) % N
    iz1 = (iz + 1) % N
    
    # Distribute mass to the 8 corners
    np.add.at(grid, (ix, iy, iz), (1-fx)*(1-fy)*(1-fz) * p_mass)
    np.add.at(grid, (ix1, iy, iz), fx*(1-fy)*(1-fz) * p_mass)
    np.add.at(grid, (ix, iy1, iz), (1-fx)*fy*(1-fz) * p_mass)
    np.add.at(grid, (ix1, iy1, iz), fx*fy*(1-fz) * p_mass)
    
    np.add.at(grid, (ix, iy, iz1), (1-fx)*(1-fy)*fz * p_mass)
    np.add.at(grid, (ix1, iy, iz1), fx*(1-fy)*fz * p_mass)
    np.add.at(grid, (ix, iy1, iz1), (1-fx)*fy*fz * p_mass)
    np.add.at(grid, (ix1, iy1, iz1), fx*fy*fz * p_mass)
    
    mean_mass = np.mean(grid)
    return np.var(grid / mean_mass)

def validate_snapshot(file_path, initial_mass=None):
    """Runs physical validations on a single HDF5 snapshot."""
    with h5py.File(file_path, 'r') as f:
        # Read Root Attributes
        domain_size = f.attrs['domain_size']
        mesh_size = f.attrs['mesh_size']
        use_hydro = bool(f.attrs['use_hydro'])
        sim_time = f.attrs['simulation_time']
        scale_factor = f.attrs['scale_factor']
        
        print(f"{BLUE}>> Validating {os.path.basename(file_path)} {RESET}| t={sim_time:.4f} | a={scale_factor:.4f}")
        
        # Load Particle Data
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

        # Compute Dark Matter Clumpiness (Density Variance)
        dm_variance = compute_cic_variance(p_x, p_y, p_z, p_mass, domain_size, mesh_size)
        print(f"    {BLUE}[INFO]{RESET} DM Density Variance (\u03c3\u00b2): {dm_variance:.5e}")

        # Conditionally Load and Check Gas Data
        if use_hydro:
            g_rho = f['gas/density'][:]
            g_eng = f['gas/energy'][:]
            
            cell_size = domain_size / mesh_size
            cell_volume = cell_size ** 3
            
            run_check(not np.isnan(g_rho).any(), "Gas density array valid (No NaNs)", "NaNs detected in gas density!")
            run_check(np.all(g_rho > 0), "Gas density strictly positive", "Negative gas density detected (Safety floor failed)!")
            run_check(np.all(g_eng > 0), "Gas energy strictly positive", "Negative gas energy detected (Safety floor failed)!")
            
            gas_mass = np.sum(g_rho) * cell_volume
            total_mass += gas_mass
            print(f"    {BLUE}[INFO]{RESET} Max Gas Density: {np.max(g_rho):.5e}")
        else:
            print(f"    {BLUE}[INFO]{RESET} Hydro disabled (Dark Matter only).")

        # Conservation Checks
        if initial_mass is not None:
            mass_drift = abs(total_mass - initial_mass) / initial_mass
            run_check(mass_drift < 1e-10, f"Mass perfectly conserved (drift: {mass_drift:.2e})", f"Mass drifted by {mass_drift * 100}%")
            print(f"    {BLUE}[INFO]{RESET} Total Mass: {total_mass:.5e}")
        else:
            print(f"    {BLUE}[INFO]{RESET} Baseline Total Mass established: {total_mass:.5e}")
        
        # Return state variables for timeline comparisons
        return {
            'total_mass': total_mass,
            'scale_factor': scale_factor,
            'dm_variance': dm_variance
        }

def run_validation_suite(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    
    if not files:
        print(f"{RED}Error: No HDF5 snapshots found in {snapshot_dir}!{RESET}")
        return
        
    print(f"\n{BLUE}=== Beginning Physics Validation Suite ({len(files)} snapshots) ==={RESET}\n")
    
    # Validate the first snapshot to establish the baseline
    initial_state = validate_snapshot(files[0])
    initial_mass = initial_state['total_mass']
    
    prev_state = initial_state
    
    # Loop through the rest to check conservation AND structure evolution
    for f in files[1:]:
        curr_state = validate_snapshot(f, initial_mass=initial_mass)
        
        prev_a = prev_state['scale_factor']
        curr_a = curr_state['scale_factor']
        prev_var = prev_state['dm_variance']
        curr_var = curr_state['dm_variance']
        
        # Clumpiness Evolution Check
        if curr_a > prev_a:
            # Did gravity work at all? (Variance must increase)
            run_check(curr_var > prev_var, 
                      f"Clumpiness increased (\u03c3\u00b2: {prev_var:.3e} -> {curr_var:.3e})", 
                      f"Universe got smoother! (\u03c3\u00b2 dropped to {curr_var:.3e})")
            
            # Did it grow at the right cosmological rate?
            linear_growth_sq = (curr_a / prev_a)**2
            min_expected_var = prev_var * linear_growth_sq * 0.7
            
            run_check(curr_var >= min_expected_var, 
                      "Structure growth rate consistent with gravity", 
                      f"Growth too slow! Expected \u03c3\u00b2 > {min_expected_var:.3e}, got {curr_var:.3e}")
        
        # Update state for the next iteration
        prev_state = curr_state
        print("") # Blank line for readability
        
    # Final Clumpiness Report
    final_state = prev_state
    print(f"{BLUE}=== Final Clumpiness Report ==={RESET}")
    print(f"    Initial Variance (a={initial_state['scale_factor']:.4f}): {initial_state['dm_variance']:.5e}")
    print(f"    Final Variance   (a={final_state['scale_factor']:.4f}): {final_state['dm_variance']:.5e}")
    
    if initial_state['dm_variance'] > 0:
        actual_growth = final_state['dm_variance'] / initial_state['dm_variance']
        theory_growth = (final_state['scale_factor'] / initial_state['scale_factor'])**2
        efficiency = (actual_growth / theory_growth) * 100
        
        print(f"    Total Growth Multiplier: {actual_growth:.2f}x")
        print(f"    Linear Theory Expectation: ~{theory_growth:.2f}x")
        print(f"    Simulation Efficiency vs Theory: {efficiency:.1f}%")
        
    print(f"\n{GREEN}===================================================={RESET}")
    print(f"{GREEN}SUCCESS: All snapshots passed physical validation!{RESET}")
    print(f"{GREEN}===================================================={RESET}\n")

if __name__ == "__main__":
    snapshot_dir = sys.argv[1] if len(sys.argv) > 1 else "./output/"
    
    if not os.path.exists(snapshot_dir):
        print(f"{RED}Error: Directory '{snapshot_dir}' does not exist.{RESET}")
        sys.exit(1)
        
    run_validation_suite(snapshot_dir)