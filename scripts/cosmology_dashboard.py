import h5py
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def compute_cic_variance(p_x, p_y, p_z, p_mass, domain_size, mesh_size):
    """Calculates density variance using Cloud-In-Cell (CIC) mass assignment."""
    N = mesh_size
    grid = np.zeros((N, N, N), dtype=np.float64)
    
    x = (p_x / domain_size) * N
    y = (p_y / domain_size) * N
    z = (p_z / domain_size) * N
    
    ix, iy, iz = np.floor(x).astype(int) % N, np.floor(y).astype(int) % N, np.floor(z).astype(int) % N
    fx, fy, fz = x - np.floor(x), y - np.floor(y), z - np.floor(z)
    ix1, iy1, iz1 = (ix + 1) % N, (iy + 1) % N, (iz + 1) % N
    
    # Distribute mass
    np.add.at(grid, (ix, iy, iz), (1-fx)*(1-fy)*(1-fz) * p_mass)
    np.add.at(grid, (ix1, iy, iz), fx*(1-fy)*(1-fz) * p_mass)
    np.add.at(grid, (ix, iy1, iz), (1-fx)*fy*(1-fz) * p_mass)
    np.add.at(grid, (ix1, iy1, iz), fx*fy*(1-fz) * p_mass)
    np.add.at(grid, (ix, iy, iz1), (1-fx)*(1-fy)*fz * p_mass)
    np.add.at(grid, (ix1, iy, iz1), fx*(1-fy)*fz * p_mass)
    np.add.at(grid, (ix, iy1, iz1), (1-fx)*fy*fz * p_mass)
    np.add.at(grid, (ix1, iy1, iz1), fx*fy*fz * p_mass)
    
    return np.var(grid / np.mean(grid))

def get_temperature(f):
    """Extracts or calculates temperature in Kelvin."""
    if 'Gas/temperature' in f:
        return f['Gas/temperature'][:]
    else:
        # Fallback for adiabatic runs
        pressure = f['Gas/pressure'][:]
        density = f['Gas/density'][:]
        
        gamma = f['Config'].attrs.get('gamma', 5.0/3.0)
        factor_u_to_t = f['Units'].attrs['factor_u_to_t']
        a = f['Header'].attrs['scale_factor']
        
        # Calculate specific internal energy: u = P / (rho * (gamma - 1))
        u_code = pressure / (density * (gamma - 1.0))
        
        return u_code * factor_u_to_t * (a**2)

def generate_dashboard(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"Error: No HDF5 snapshots found in {snapshot_dir}")
        return

    print(f"Processing {len(files)} snapshots for diagnostics...")

    # Time series data
    scale_factors = []
    times_gyr = []
    dm_variances = []
    
    # Hydro time series
    max_gas_densities = []
    max_temperatures = []
    kin_energies = []
    therm_energies = []
    rad_energies = []
    cold_gas_fracs = []
    
    # PDF Data (storing histograms for first, middle, and last snapshots)
    pdf_data = {}
    target_pdf_indices = [0, len(files)//2, len(files)-1]

    # Pre-calculate CGS Unit Multipliers
    with h5py.File(files[0], 'r') as f:
        domain_size = f['Config'].attrs.get('domain_size', 1.0)
        mesh_size = f['Config'].attrs.get('mesh_size', 32)


    # Volume of a single cell in code units (needed for integrating total energy)
    cell_vol_code = (domain_size / mesh_size)**3

    for i, f_path in enumerate(files):
        with h5py.File(f_path, 'r') as f:
            header = f['Header'].attrs
            config = f['Config'].attrs
            units  = f['Units'].attrs

            a = header['scale_factor']
            scale_factors.append(a)
            times_gyr.append(header.get('simulation_time', 0.0) * units.get('unit_time_in_gyr', 1.0))
            
            u_density_cgs = units['unit_density_in_cgs']
            u_energy_cgs = units['unit_velocity_in_cgs']**2 * (units['unit_mass_in_msun'] * 1.98847e33)
            
            # DM Variance
            p_x = f['Particles/position_x'][:]
            p_y = f['Particles/position_y'][:]
            p_z = f['Particles/position_z'][:]
            p_mass = f['Particles/mass'][:]
            dm_variances.append(compute_cic_variance(p_x, p_y, p_z, p_mass, config['domain_size'], config['mesh_size']))
            
            # Gas Stats
            if bool(config.get('use_hydro', 0)):
                rho = f['Gas/density'][:]
                px = f['Gas/momentum_x'][:]
                py = f['Gas/momentum_y'][:]
                pz = f['Gas/momentum_z'][:]
                e_tot = f['Gas/energy'][:]
                temp = get_temperature(f)
                
                # Convert Max Density to Hydrogen Number Density (n_H)
                # max_rho_cgs is the comoving density. Divide by a^3 to get physical density
                physical_rho_cgs = (np.max(rho) * u_density_cgs) / (a**3)
                M_P_CGS = 1.67262192e-24
                X_H = 0.76  # Primordial hydrogen mass fraction
                n_H = (physical_rho_cgs * X_H) / M_P_CGS
                max_gas_densities.append(n_H)
                max_temperatures.append(np.max(temp))
                
                # Calculate Total Energies in Ergs
                safe_rho = np.maximum(rho, 1e-10)
                ke_grid = (px**2 + py**2 + pz**2) / (2.0 * safe_rho)
                
                # Multiply the comoving energy by a^2 to get physical Ergs
                total_ke_ergs = np.sum(ke_grid) * cell_vol_code * u_energy_cgs * (a**2)
    
                # Do the same for thermal energy
                total_te_ergs = (np.sum(e_tot) - np.sum(ke_grid)) * cell_vol_code * u_energy_cgs * (a**2)
                
                kin_energies.append(total_ke_ergs)
                therm_energies.append(total_te_ergs)

                rad_energy_code = f['Gas'].attrs.get('cumulative_radiated_energy', 0.0)
                total_rad_ergs = rad_energy_code * u_energy_cgs
                rad_energies.append(total_rad_ergs)
                
                # Cold Dense Gas Fraction (Star formation precursor)
                mean_rho = np.mean(rho)
                overdensity = rho / mean_rho
                cold_dense_mask = (temp < 10000.0) & (overdensity > 100.0)
                mass_fraction = np.sum(rho[cold_dense_mask]) / np.sum(rho)
                cold_gas_fracs.append(mass_fraction)
                
                # Snapshot PDF Extraction
                if i in target_pdf_indices:
                    # Log10 of Overdensity
                    log_od = np.log10(np.maximum(overdensity, 1e-5))
                    hist, edges = np.histogram(log_od, bins=100, range=(-2, 5))
                    # Convert to volume fraction
                    pdf_data[a] = (hist / np.sum(hist), edges)

    # Info extraction (from the LAST snapshot)
    with h5py.File(files[-1], 'r') as f:
        config = f['Config'].attrs
        
        box_size = config.get('box_size_mpc', 'N/A')
        mesh_size = config.get('mesh_size', 'N/A')
        omega_m = config.get('omega_M', 'N/A')
        omega_b = config.get('omega_baryon', 'N/A')
        omega_l = config.get('omega_lambda', 'N/A')
        hubble = config.get('hubble_param', 'N/A')
        gamma = config.get('gamma', 'N/A')
        mu = config.get('primordial_mu', 'N/A')
        prim_index = config.get('spectral_index', 'N/A')
        
        num_particles = len(f['Particles/position_x']) if 'Particles/position_x' in f else 'N/A'
        
        if bool(config.get('use_hydro', 0)):
            final_rho = f['Gas/density'][:]
            final_overdensity = final_rho / np.mean(final_rho) 
            final_temp = get_temperature(f)
            x_data = final_overdensity.flatten()
            y_data = final_temp.flatten()

    # Plotting (3x3 Grid)
    fig, axs = plt.subplots(3, 3, figsize=(18, 14))
    fig.canvas.manager.set_window_title(f"Dashboard: {os.path.basename(os.path.normpath(snapshot_dir))}")
    fig.suptitle(f"Dashboard: {os.path.basename(os.path.normpath(snapshot_dir))}", fontsize=18)

    # Expansion History
    axs[0, 0].plot(times_gyr, scale_factors, color='purple', lw=2)
    axs[0, 0].set_xlabel('Simulation Time [Gyr]')
    axs[0, 0].set_ylabel('Scale Factor (a)')
    axs[0, 0].set_title('Cosmic Expansion History')

    # DM Variance Evolution (Structure Growth)
    axs[0, 1].plot(scale_factors, dm_variances, color='black', label='Simulated Variance')
    theory_growth = dm_variances[0] * (np.array(scale_factors) / scale_factors[0])**2
    axs[0, 1].plot(scale_factors, theory_growth, linestyle='--', color='red', label='Linear Theory')
    axs[0, 1].set_yscale('log')
    axs[0, 1].set_xlabel('Scale Factor (a)')
    axs[0, 1].set_ylabel(r'DM Variance ($\sigma^2$)')
    axs[0, 1].set_title('Structure Growth')
    axs[0, 1].legend()

    # Density PDF Evolution (Probability Density Function)
    if pdf_data:
        for a_val, (hist, edges) in pdf_data.items():
            centers = (edges[:-1] + edges[1:]) / 2
            axs[0, 2].plot(centers, hist, lw=2, label=f'a = {a_val:.2f}')
        axs[0, 2].set_yscale('log')
        axs[0, 2].set_xlabel(r'$\log_{10}$ Gas Overdensity ($\rho/\bar{\rho}$)')
        axs[0, 2].set_ylabel('Volume Fraction')
        axs[0, 2].set_title('Density PDF Evolution')
        axs[0, 2].legend()

    # Global Energy Inventory
    if kin_energies:
        axs[1, 0].plot(scale_factors, kin_energies, color='green', lw=2, label='Kinetic Energy')
        axs[1, 0].plot(scale_factors, therm_energies, color='orange', lw=2, label='Thermal Energy')
        axs[1, 0].plot(scale_factors, rad_energies, color='red', lw=2, linestyle=':', label='Radiated Energy (cumulative)')
        axs[1, 0].set_yscale('log')
        axs[1, 0].set_xlabel('Scale Factor (a)')
        axs[1, 0].set_ylabel('Total Energy [Ergs]')
        axs[1, 0].set_title('Global Gas Energy Inventory')
        axs[1, 0].legend()

    # Max Gas Density
    if max_gas_densities:
        axs[1, 1].plot(scale_factors, max_gas_densities, color='blue')
        axs[1, 1].set_yscale('log')
        axs[1, 1].set_xlabel('Scale Factor (a)')
        axs[1, 1].set_ylabel(r'Max Hydrogen Density $n_H$ [cm$^{-3}$]')
        axs[1, 1].set_title('Maximum Gas Density')

    # Max Temperature
    if max_temperatures:
        axs[1, 2].plot(scale_factors, max_temperatures, color='red')
        axs[1, 2].set_yscale('log')
        axs[1, 2].set_xlabel('Scale Factor (a)')
        axs[1, 2].set_ylabel('Max Temperature [K]')
        axs[1, 2].set_title('Maximum Temperature')

    # Phase Diagram
    if 'x_data' in locals():
        # Create logarithmically spaced bins
        x_bins = np.logspace(np.log10(1e-2), np.log10(np.max(x_data)), 100)
        y_bins = np.logspace(np.log10(10.0), np.log10(np.max(y_data)), 100)
        
        # Pass log bins into the hist2d function
        h = axs[2, 0].hist2d(x_data, y_data, bins=[x_bins, y_bins], norm=LogNorm(), cmap='plasma')
        axs[2, 0].set_xscale('log')
        axs[2, 0].set_yscale('log')
        axs[2, 0].set_xlabel(r'Gas Overdensity ($\rho / \bar{\rho}$)')
        axs[2, 0].set_ylabel('Temperature [K]')
        axs[2, 0].set_title('Phase Diagram (Final Snapshot)')
        fig.colorbar(h[3], ax=axs[2, 0], label='Number of Cells')
    else:
        axs[2, 0].text(0.5, 0.5, 'Hydro Disabled', ha='center', va='center')   

    # Cold Dense Gas Fraction 
    if cold_gas_fracs:
        axs[2, 1].plot(scale_factors, cold_gas_fracs, color='teal', lw=2)
        axs[2, 1].set_xlabel('Scale Factor (a)')
        axs[2, 1].set_ylabel('Mass Fraction')
        axs[2, 1].set_title(r'Cold Dense Gas ($T<10^4$ K, $\delta>100$)')

    # Info Card
    axs[2, 2].axis('off') # Hide the axes
    
    # Safely format numbers (so we don't crash if 'N/A' is returned)
    part_str = f"{num_particles:,}" if isinstance(num_particles, int) else str(num_particles)
    h_str = f"{hubble:.2f}" if isinstance(hubble, float) else str(hubble)
    g_str = f"{gamma:.3f}" if isinstance(gamma, float) else str(gamma)
    m_str = f"{mu:.2f}" if isinstance(mu, float) else str(mu)
    idx_str = f"{prim_index:.2f}" if isinstance(prim_index, float) else str(prim_index)

    # Use Python's < padding to perfectly align the two columns
    info_text = (
        f"SIMULATION DATA\n\n"
        f"{'Box Size':<9}: {str(box_size) + ' Mpc'}\n"
        f"{'Grid':<9}: {str(mesh_size) + '³'}\n"
        f"{'Particles':<9}: {part_str}\n\n"
        f"Cosmology & Physics:\n"
        f"  {'Ω_m':<3}: {str(omega_m):<8} | {'Hubble (h)':<12}: {h_str}\n"
        f"  {'Ω_b':<3}: {str(omega_b):<8} | {'Gamma (γ)':<12}: {g_str}\n"
        f"  {'Ω_Λ':<3}: {str(omega_l):<8} | {'Primord. μ':<12}: {m_str}\n"
        f"  {'':<3}  {'':<8} | {'Prim. Index':<12}: {idx_str}"
    )
    
    props = dict(boxstyle='round,pad=1.0', facecolor='#f8f9fa', alpha=1.0, edgecolor='#cccccc')
    
    # Place text in the center of the empty subplot
    axs[2, 2].text(0.5, 0.5, info_text, fontsize=12, color='black',
                   verticalalignment='center', horizontalalignment='center', 
                   multialignment='left',
                   bbox=props, family='monospace')

    plt.tight_layout(rect=[0, 0.01, 1, 0.97], h_pad=2.2, w_pad=2.0)
    
    try:
        plt.show()
    except Exception:
        pass

if __name__ == "__main__":
    target_dir = sys.argv[1] if len(sys.argv) > 1 else "./output/"
    generate_dashboard(target_dir)