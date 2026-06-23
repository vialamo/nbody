import h5py
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import camb

# Physics & Cosmological Constants
M_P_CGS = 1.67262192e-24
X_H = 0.76  # Primordial hydrogen mass fraction
G_CGS = 6.6743e-8
K_B_CGS = 1.3806e-16

def get_camb_pk(omega_m, omega_b, hubble, ns, target_sigma8, scale_factors, k_min=1.0, k_max=20.0):
    """Generates theoretical non-linear P(k) using CAMB and Halofit."""
    pars = camb.CAMBparams()
    
    if hubble > 10.0:
        hubble = hubble / 100.0
        
    H0 = hubble * 100.0
    ombh2 = omega_b * (hubble ** 2)
    omch2 = (omega_m - omega_b) * (hubble ** 2)
    
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(As=2e-9, ns=ns)
    pars.set_matter_power(redshifts=[0.0], kmax=2.0)
    default_results = camb.get_results(pars)
    default_sigma8 = default_results.get_sigma8()[-1]
    
    # Scale the primordial amplitude (As)
    target_As = 2e-9 * (target_sigma8 / default_sigma8)**2
    pars.InitPower.set_params(As=target_As, ns=ns)
    
    # Sanitize redshifts and sort for CAMB descending order
    redshifts_to_run = sorted(list({max(0.0, 1.0 / a - 1.0) for a in scale_factors}), reverse=True)
    
    pars.set_matter_power(redshifts=redshifts_to_run, kmax=k_max)
    pars.NonLinear = camb.model.NonLinear_both
    
    results = camb.get_results(pars)
    k_h, z_out, pk = results.get_matter_power_spectrum(minkh=k_min, maxkh=k_max, npoints=200)
    
    theory_pk = {}
    for i, z_val in enumerate(z_out):
        # Map the CAMB output back to the closest requested scale factor
        closest_a = min(scale_factors, key=lambda x: abs(max(0.0, 1.0/x - 1.0) - z_val))
        theory_pk[closest_a] = (k_h, pk[i, :])
        
    return theory_pk

def assign_cic_grid(p_x, p_y, p_z, p_mass, domain_size, mesh_size):
    """Helper function to perform 3D Cloud-In-Cell (CIC) mass assignment."""
    N = mesh_size
    grid = np.zeros((N, N, N), dtype=np.float64)
    
    x = (p_x / domain_size) * N
    y = (p_y / domain_size) * N
    z = (p_z / domain_size) * N
    
    ix, iy, iz = np.floor(x).astype(int) % N, np.floor(y).astype(int) % N, np.floor(z).astype(int) % N
    fx, fy, fz = x - np.floor(x), y - np.floor(y), z - np.floor(z)
    ix1, iy1, iz1 = (ix + 1) % N, (iy + 1) % N, (iz + 1) % N
    
    np.add.at(grid, (ix, iy, iz), (1-fx)*(1-fy)*(1-fz) * p_mass)
    np.add.at(grid, (ix1, iy, iz), fx*(1-fy)*(1-fz) * p_mass)
    np.add.at(grid, (ix, iy1, iz), (1-fx)*fy*(1-fz) * p_mass)
    np.add.at(grid, (ix1, iy1, iz), fx*fy*(1-fz) * p_mass)
    np.add.at(grid, (ix, iy, iz1), (1-fx)*(1-fy)*fz * p_mass)
    np.add.at(grid, (ix1, iy, iz1), fx*(1-fy)*fz * p_mass)
    np.add.at(grid, (ix, iy1, iz1), (1-fx)*fy*fz * p_mass)
    np.add.at(grid, (ix1, iy1, iz1), fx*fy*fz * p_mass)
    
    return grid

def compute_cic_variance(p_x, p_y, p_z, p_mass, domain_size, mesh_size):
    """Calculates density variance using Cloud-In-Cell (CIC) mass assignment."""
    grid = assign_cic_grid(p_x, p_y, p_z, p_mass, domain_size, mesh_size)
    return np.var(grid / np.mean(grid))

def compute_power_spectrum(p_x, p_y, p_z, p_mass, domain_size, mesh_size, gas_rho=None):
    """Calculates the 1D total matter power spectrum from particle and gas data."""
    N = mesh_size
    grid = assign_cic_grid(p_x, p_y, p_z, p_mass, domain_size, mesh_size)

    if gas_rho is not None:
        cell_vol = (domain_size / N)**3
        grid += gas_rho * cell_vol
    
    # Create overdensity field (delta)
    mean_rho = np.mean(grid)
    delta = (grid - mean_rho) / mean_rho if mean_rho > 0 else grid
    
    fft_delta = np.fft.fftn(delta)
    # Apply normalization
    box_vol = domain_size**3
    power_3d = (np.abs(fft_delta)**2) * (box_vol / (float(N)**6))
    
    # Calculate k modes
    k_1d = np.fft.fftfreq(N, d=domain_size/N) * 2.0 * np.pi
    kx, ky, kz = np.meshgrid(k_1d, k_1d, k_1d, indexing='ij')
    k_mag = np.sqrt(kx**2 + ky**2 + kz**2)
    
    # Flatten and remove k=0 mode
    k_mag_flat = k_mag.flatten()
    power_flat = power_3d.flatten()
    valid = k_mag_flat > 0
    k_mag_flat = k_mag_flat[valid]
    power_flat = power_flat[valid]
    
    # Binning
    k_min = 2.0 * np.pi / domain_size
    k_max = np.pi * N / domain_size
    bins = np.logspace(np.log10(k_min), np.log10(k_max), N//2)
    
    hist_power, bin_edges = np.histogram(k_mag_flat, bins=bins, weights=power_flat)
    hist_counts, _ = np.histogram(k_mag_flat, bins=bins)
    
    p_k = np.zeros_like(hist_power)
    mask = hist_counts > 0
    p_k[mask] = hist_power[mask] / hist_counts[mask]
    k_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    
    return k_centers[mask], p_k[mask]

def get_linear_growth(a, omega_m_0, omega_l_0):
    """Calculates the linear growth factor D(a) using the Carroll, Press & Turner (1992) approximation."""
    if a == 0: return 0.0
    a3 = a**3
    omega_k_0 = 1.0 - omega_m_0 - omega_l_0
    
    # E(a)^2 = (H(a)/H0)^2
    E2 = omega_m_0 / a3 + omega_k_0 / (a**2) + omega_l_0
    # Time-evolving density parameters at scale factor 'a'
    Om_a = (omega_m_0 / a3) / E2
    Ol_a = omega_l_0 / E2
    
    # Growth suppression factor g(a)
    g_a = (2.5 * Om_a) / (Om_a**(4.0/7.0) - Ol_a + (1.0 + Om_a / 2.0) * (1.0 + Ol_a / 70.0))

    # D(a) = a * g(a)
    return a * g_a

def get_temperature(f):
    """Extracts or calculates temperature in Kelvin."""
    if 'Gas/temperature' in f:
        return f['Gas/temperature'][:]
    
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

    # Pre-calculate Global Simulation Config & Units (from the first snapshot)
    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        units = f['Units'].attrs
        
        domain_size = config.get('domain_size', 1.0)
        mesh_size = config.get('mesh_size', 32)
        has_hydro = bool(config.get('use_hydro', 0))
        hubble = config.get('hubble_param', 70.0)
        h_val = hubble / 100.0 if hubble > 10.0 else hubble
        box_size_mpc = config.get('box_size_mpc', 1.0)
        box_h_mpc = box_size_mpc * h_val
        
        omega_m = config.get('omega_M', 1.0)
        omega_b = config.get('omega_baryon', 0.04)
        omega_l = config.get('omega_lambda', 0.0)
        gamma = config.get('gamma', 5.0/3.0)
        mu = config.get('primordial_mu', 0.6)
        prim_index = config.get('spectral_index', 0.96)
        sigma_8 = float(config.get('sigma_8', 0.81))
        T_floor = float(config.get('cooling_temp_floor', 1000.0))
        floor_k = config.get('temp_floor_k', 10.0)
        
        u_density_cgs = units.get('unit_density_in_cgs', 1.0)
        u_energy_cgs = units.get('unit_velocity_in_cgs', 1.0)**2 * (units.get('unit_mass_in_msun', 1.0) * 1.98847e33)
        u_time_gyr = units.get('unit_time_in_gyr', 1.0)
        
        # Volume of a single cell in code units
        cell_vol_code = (domain_size / mesh_size)**3

    # Time series lists
    scale_factors, times_gyr, dm_variances = [], [], []
    p999_gas_densities, max_gas_densities = [], []
    p999_temperatures, max_temperatures = [], []
    kin_energies, therm_energies, rad_energies, cold_gas_fracs, fractional_errors = [], [], [], [], []

    pdf_data, pk_data = {}, {}
    target_indices = [0, len(files)//2, len(files)-1]
    initial_e_code = None

    for i, f_path in enumerate(files):
        sys.stdout.write(f"\rProcessing snapshot {i+1}/{len(files)} [{(i + 1) / len(files) * 100:.1f}%]")
        sys.stdout.flush()

        with h5py.File(f_path, 'r') as f:
            header = f['Header'].attrs
            a = header['scale_factor']
            
            scale_factors.append(a)
            times_gyr.append(header.get('simulation_time', 0.0) * u_time_gyr)

            # Particles
            p_x = f['Particles/position_x'][:] * box_h_mpc
            p_y = f['Particles/position_y'][:] * box_h_mpc
            p_z = f['Particles/position_z'][:] * box_h_mpc
            p_mass = f['Particles/mass'][:]
            
            dm_variances.append(compute_cic_variance(p_x, p_y, p_z, p_mass, box_h_mpc, mesh_size))

            rho = f['Gas/density'][:] if has_hydro else None

            # Snapshot Extraction for PDF & Power Spectrum
            if i in target_indices:
                k_bins, pk = compute_power_spectrum(p_x, p_y, p_z, p_mass, box_h_mpc, mesh_size, gas_rho=rho)
                pk_data[a] = (k_bins, pk)
            
            # Gas Stats
            if has_hydro:
                px, py, pz = f['Gas/momentum_x'][:], f['Gas/momentum_y'][:], f['Gas/momentum_z'][:]
                e_tot = f['Gas/energy'][:]
                temp = get_temperature(f)
                
                # Convert Max Density to Hydrogen Number Density (n_H)
                # max_rho_cgs is the comoving density. Divide by a^3 to get physical density
                p999_rho, max_rho = np.percentile(rho, 99.9), np.max(rho)
                n_H_conv = (u_density_cgs / a**3 * X_H) / M_P_CGS
                p999_gas_densities.append(p999_rho * n_H_conv)
                max_gas_densities.append(max_rho * n_H_conv)

                p999_temperatures.append(np.percentile(temp, 99.9))
                max_temperatures.append(np.max(temp))
                
                # Calculate Total Energies in Ergs
                # Multiply the comoving energy by a^2 to get physical Ergs
                safe_rho = np.maximum(rho, 1e-10)
                ke_grid = (px**2 + py**2 + pz**2) / (2.0 * safe_rho)
                
                energy_conv = cell_vol_code * u_energy_cgs * (a**2)
                kin_energies.append(np.sum(ke_grid) * energy_conv)
                therm_energies.append((np.sum(e_tot) - np.sum(ke_grid)) * energy_conv)

                rad_energy_code = f['Gas'].attrs.get('cumulative_radiated_energy', 0.0)
                rad_energies.append(rad_energy_code * u_energy_cgs * (a**2))

                current_e_code = np.sum(e_tot) * cell_vol_code
                if initial_e_code is None:
                    initial_e_code = current_e_code

                w_grav_code = f['Gas'].attrs.get('cumulative_gravitational_work', 0.0)
                w_exp_code = f['Gas'].attrs.get('cumulative_expansion_work', 0.0)

                # The First Law of Thermodynamics (in Comoving Code Units)
                delta_e_code = current_e_code - initial_e_code
                absolute_error_code = delta_e_code - w_grav_code + w_exp_code + rad_energy_code
                fractional_errors.append(absolute_error_code / abs(initial_e_code) if initial_e_code != 0 else 0.0)
                
                # Cold Dense Gas Fraction
                mean_rho = np.mean(rho)
                overdensity = rho / mean_rho
                cold_dense_mask = (temp < 10000.0) & (overdensity > 100.0)
                cold_gas_fracs.append(np.sum(rho[cold_dense_mask]) / np.sum(rho))
                
                # Snapshot PDF Extraction
                if i in target_indices:
                    safe_od = np.maximum(overdensity, 1e-10)
                    hist, edges = np.histogram(np.log10(np.maximum(overdensity, 1e-5)), bins=100, range=(-2, 5))
                    pdf_data[a] = (hist / np.sum(hist), edges, np.var(safe_od))

    print("\nData extraction complete. Generating plots...")

    # Last snapshot phase data for Phase Diagram
    if has_hydro:
        with h5py.File(files[-1], 'r') as f:
            final_rho = f['Gas/density'][:]
            final_temp = get_temperature(f)
            x_data = (final_rho / np.mean(final_rho)).flatten()
            y_data = final_temp.flatten()
            num_particles = len(f['Particles/position_x']) if 'Particles/position_x' in f else 'N/A'

    # Plotting Grid Setup
    fig, axs = plt.subplots(3, 3, figsize=(18, 14))
    fig.canvas.manager.set_window_title(f"Dashboard: {os.path.basename(os.path.normpath(snapshot_dir))}")
    
    # Matter Power Spectrum (With CAMB Theory)
    if pk_data:
        try:
            ns = float(prim_index) if isinstance(prim_index, (float, int)) else 0.96
            all_k_bins = np.concatenate([k_bins for _, (k_bins, _) in pk_data.items()])
            camb_theory = get_camb_pk(omega_m, omega_b, hubble, ns, sigma_8, list(pk_data.keys()), 
                                      k_min=np.min(all_k_bins), k_max=np.max(all_k_bins))
        except Exception as e:
            camb_theory = None
            print(f"CAMB theory failed: {e}")

        colors = ['blue', 'green', 'red']
        for idx, (a_val, (k_bins, pk)) in enumerate(pk_data.items()):
            c = colors[idx % len(colors)]
            axs[0, 0].plot(k_bins, pk, lw=2, color=c, label=f'Sim a={a_val:.2f}')
            if camb_theory and a_val in camb_theory:
                theory_k, theory_pk = camb_theory[a_val]
                axs[0, 0].plot(theory_k, theory_pk, lw=1, linestyle='--', color=c, alpha=0.7, label=f'Theory a={a_val:.2f}')

        axs[0, 0].set(xscale='log', yscale='log', xlabel=r'k [$h$ Mpc$^{-1}$]', ylabel=r'$P(k)$ [($h^{-1}$ Mpc)$^3$]', title='Matter Power Spectrum')
        axs[0, 0].legend(fontsize=8)

    # 1-Point Volume-Weighted Density PDF
    if pdf_data:
        colors = ['blue', 'green', 'red']
        for idx, (a_val, (hist, edges, sigma2)) in enumerate(pdf_data.items()):
            c = colors[idx % len(colors)]
            centers = (edges[:-1] + edges[1:]) / 2
            dx = centers[1] - centers[0]
            
            # Plot the original simulation data
            # Filter by y > 1e-5
            # Keep points right before and after entering > 1e-5
            base_mask = hist > 1e-5
            mask_hist = base_mask.copy()
            mask_hist[:-1] |= base_mask[1:]
            mask_hist[1:] |= base_mask[:-1]
            axs[0, 1].plot(centers[mask_hist], hist[mask_hist], lw=2, color=c, label=f'Sim a={a_val:.2f}')

            # Determine the x-axis extent of the valid simulation data
            if np.any(mask_hist):
                x_min, x_max = centers[mask_hist][0], centers[mask_hist][-1]
            else:
                x_min, x_max = centers[0], centers[-1]      
            domain_mask = (centers >= x_min) & (centers <= x_max)

            Delta = 10**centers
            
            # Model A: Gaussian (early universe)
            if idx == 0:
                p_delta = (1.0 / np.sqrt(2.0 * np.pi * sigma2)) * np.exp(-0.5 * (Delta - 1.0)**2 / sigma2)
                y_gauss = p_delta * Delta * np.log(10) * dx
                mask_gauss = (y_gauss > 1e-5) & domain_mask
                axs[0, 1].plot(centers[mask_gauss], y_gauss[mask_gauss], color=c, linestyle=':', lw=1, label='Gaussian Model')
            
            # Model B: Lognormal (late universe)
            elif idx == len(pdf_data) - 1:
                sigma2_A, mu_A = np.log(1.0 + sigma2), -np.log(1.0 + sigma2) / 2.0
                p_A = (1.0 / np.sqrt(2.0 * np.pi * sigma2_A)) * np.exp(-0.5 * (centers * np.log(10) - mu_A)**2 / sigma2_A)
                y_lognorm = p_A * np.log(10) * dx
                mask_lognorm = (y_lognorm > 1e-5) & domain_mask
                axs[0, 1].plot(centers[mask_lognorm], y_lognorm[mask_lognorm], color=c, linestyle='--', lw=1, label='Lognormal Model')

        axs[0, 1].set(yscale='log', xlabel=r'$\log_{10}$ Gas Overdensity ($\rho/\bar{\rho}$)', ylabel='Volume Fraction', title='1-Point Volume-Weighted Density PDF')
        axs[0, 1].set_ylim(bottom=1e-5, top=1.5)
        axs[0, 1].legend(fontsize=8)

    # DM Variance Evolution (Structure Growth)
    axs[0, 2].plot(scale_factors, dm_variances, color='black', label='Simulated Variance')
    D_a = np.array([get_linear_growth(a, float(omega_m), float(omega_l)) for a in scale_factors])
    axs[0, 2].plot(scale_factors, dm_variances[0] * (D_a / D_a[0])**2, linestyle='--', color='red', label=r'$\Lambda$CDM Linear Theory')
    axs[0, 2].set(yscale='log', xlabel='Scale Factor (a)', ylabel=r'DM Variance ($\sigma^2$)', title='Structure Growth')
    axs[0, 2].legend()

    # Global Energy Inventory & Conservation Error
    if kin_energies:
        axs[1, 0].plot(scale_factors, kin_energies, color='green', lw=2, label='Kinetic')
        axs[1, 0].plot(scale_factors, therm_energies, color='orange', lw=2, label='Thermal')
        axs[1, 0].plot(scale_factors, rad_energies, color='red', lw=2, linestyle=':', label='Radiated')
        axs[1, 0].set(yscale='log', xlabel='Scale Factor (a)', ylabel='Energy Components [Ergs]', title='Energy Inventory & Conservation')
        
        ax_err = axs[1, 0].twinx()
        ax_err.plot(scale_factors, fractional_errors, color='black', lw=1.5, linestyle='--', label='Fractional Error')
        # Set symmetric limits around 0 for the error (e.g., +/- 1%)
        max_err = max(1e-4, np.max(np.abs(fractional_errors)) * 1.5)
        ax_err.set_ylim(-max_err, max_err)
        ax_err.set_ylabel(r'Fractional Error ($\Delta E / E_0$)', color='black')
        ax_err.axhline(0, color='gray', linestyle='-', linewidth=0.5)

        lines_1, labels_1 = axs[1, 0].get_legend_handles_labels()
        lines_2, labels_2 = ax_err.get_legend_handles_labels()
        axs[1, 0].legend(lines_1 + lines_2, labels_1 + labels_2, loc='center left', fontsize=8)

    # Extreme Gas States (99.9% to Max Envelope)
    if p999_gas_densities:
        # Plot Density on the primary left axis
        axs[1, 1].plot(scale_factors, p999_gas_densities, color='blue', lw=2, label='99.9% Density')
        axs[1, 1].fill_between(scale_factors, p999_gas_densities, max_gas_densities, color='blue', alpha=0.2)
        axs[1, 1].set(yscale='log', xlabel='Scale Factor (a)', ylabel=r'Hydrogen Density $n_H$ [cm$^{-3}$]')
        axs[1, 1].tick_params(axis='y', labelcolor='blue')
        
        # Plot Temperature on the secondary right axis
        ax_temp = axs[1, 1].twinx()
        ax_temp.plot(scale_factors, p999_temperatures, color='red', lw=2, label='99.9% Temp')
        ax_temp.fill_between(scale_factors, p999_temperatures, max_temperatures, color='red', alpha=0.2)
        ax_temp.set(yscale='log', ylabel='Temperature [K]')
        ax_temp.tick_params(axis='y', labelcolor='red')

        # Calculate the variance of the fundamental mode: sigma(L) = sigma_8 * (8/L)^0.9
        sigma_L = sigma_8 * (8.0 / box_size_mpc)**0.9
        a_limit = 0.3 / sigma_L
        if a_limit < 1.0:
            axs[1, 1].axvline(a_limit, color='black', linestyle='--', linewidth=1, label='Global Collapse Limit')
        
        lines_left, labels_left = axs[1, 1].get_legend_handles_labels()
        axs[1, 1].legend(lines_left, labels_left, loc='upper left', fontsize=8)
        axs[1, 1].set_title('Extreme States (99.9% - Max Envelope)')

    # Cosmic Expansion History
    axs[1, 2].plot(times_gyr, scale_factors, color='purple', lw=2, label='Simulation')
    axs[1, 2].set(xlabel='Simulation Time [Gyr]', ylabel='Scale Factor (a)', title='Cosmic Expansion History', ylim=(0.0, None))
    axs[1, 2].legend()

    # Phase Diagram
    if has_hydro and 'x_data' in locals():
        # Create logarithmically spaced bins
        x_bins, y_bins = np.logspace(-2, np.log10(np.max(x_data)), 100), np.logspace(1, np.log10(np.max(y_data)), 100)
        h = axs[2, 0].hist2d(x_data, y_data, bins=[x_bins, y_bins], norm=LogNorm(), cmap='plasma')
        
        # THEORETICAL BASELINE
        mask_mean = (x_data > 0.95) & (x_data < 1.05)
        if np.any(mask_mean):
            T_line = np.maximum(np.min(y_data[mask_mean]) * (x_bins ** (gamma - 1.0)), floor_k)
            axs[2, 0].plot(x_bins, T_line, color='cyan', linestyle='--', linewidth=2, label=r'Adiabatic Track')
        
        # Truelove Jeans Fragmentation Limit (Density Trust Barrier)
        # Read parameters
        H0_cgs = (hubble * 100.0) * 1e5 / 3.086e24
        rho_crit_cgs = (3.0 * H0_cgs**2) / (8.0 * np.pi * G_CGS)
        dx_comoving_cm = (box_size_mpc / mesh_size) * 3.086e24
        # Truelove math: lambda_J >= 4 * dx_phys
        # This isolates the maximum overdensity (Delta) before fragmentation occurs
        K = (np.pi * gamma * K_B_CGS) / (mu * M_P_CGS * G_CGS)
        max_safe_delta = (K * T_floor * scale_factors[-1]) / (16.0 * (dx_comoving_cm**2) * omega_m * rho_crit_cgs)
        
        axs[2, 0].axvline(x=max_safe_delta, color='red', linestyle=':', linewidth=2, label='Truelove Limit')
        axs[2, 0].axvspan(max_safe_delta, np.max(x_data) * 10, color='red', alpha=0.1)
        
        axs[2, 0].set(xscale='log', yscale='log', xlabel=r'Gas Overdensity ($\rho / \bar{\rho}$)', ylabel='Temperature [K]', title='Phase Diagram (Final Snapshot)')
        axs[2, 0].legend(loc='upper left', fontsize=8)
        fig.colorbar(h[3], ax=axs[2, 0], label='Number of Cells')
    else:
        axs[2, 0].text(0.5, 0.5, 'Hydro Disabled', ha='center', va='center')  

    # Cold Dense Gas Fraction 
    if cold_gas_fracs:
        axs[2, 1].plot(scale_factors, cold_gas_fracs, color='teal', lw=2)
        axs[2, 1].set(xlabel='Scale Factor (a)', ylabel='Mass Fraction', title=r'Cold Dense Gas ($T<10^4$ K, $\delta>100$)')

    # Info Card
    axs[2, 2].axis('off')
    info_text = (
        f"SIMULATION: {os.path.basename(os.path.normpath(snapshot_dir))}\n\n"
        f"{'Box Size':<9}: {str(box_size_mpc) + ' Mpc'}\n"
        f"{'Grid':<9}: {str(mesh_size) + '³'}\n"
        f"{'Particles':<9}: {num_particles if has_hydro else 'N/A'}\n\n"
        f"Cosmology & Physics:\n"
        f"  {'Ω_m':<3}: {omega_m:<8} | {'Hubble (h)':<12}: {h_val:.2f}\n"
        f"  {'Ω_b':<3}: {omega_b:<8} | {'Gamma (γ)':<12}: {gamma:.3f}\n"
        f"  {'Ω_Λ':<3}: {omega_l:<8} | {'Primord. μ':<12}: {mu:.2f}\n"
        f"  {'':<3}  {'':<8} | {'Prim. Index':<12}: {prim_index:.2f}"
    )
    
    axs[2, 2].text(0.5, 0.5, info_text, fontsize=12, color='black', va='center', ha='center', multialignment='left',
                   bbox=dict(boxstyle='round,pad=1.0', facecolor='#f8f9fa', alpha=1.0, edgecolor='#cccccc'), family='monospace')

    plt.tight_layout(rect=[0, 0.01, 1, 0.99], h_pad=2.2, w_pad=2.0)
    plt.show()

if __name__ == "__main__":
    generate_dashboard(sys.argv[1] if len(sys.argv) > 1 else "./output/")