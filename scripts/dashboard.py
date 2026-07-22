import h5py
import numpy as np
import glob
import os
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import camb
from scipy.spatial import cKDTree

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

def compute_power_spectrum(p_x, p_y, p_z, p_mass, domain_size, mesh_size, part_mesh_size, gas_rho=None):
    """Calculates the 1D total matter power spectrum from particle and gas data."""
    N = mesh_size
    grid = assign_cic_grid(p_x, p_y, p_z, p_mass, domain_size, mesh_size)

    # If the grid is oversampled, apply a Gaussian blur to the DM particles
    if mesh_size > part_mesh_size:
        from scipy.ndimage import gaussian_filter
        # Smear the discrete particles over the extra grid cells
        sigma = mesh_size / (2.0 * part_mesh_size)
        grid = gaussian_filter(grid, sigma=sigma, mode='wrap')

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

def compute_max_kdtree_density(p_x, p_y, p_z, p_mass, box_size, k=32):
    """Calculates the maximum comoving density of particles using a KD-Tree nearest-neighbor search."""
    points = np.vstack((p_x, p_y, p_z)).T
    
    # Build a periodic KD-Tree
    tree = cKDTree(points, boxsize=box_size)
    
    # Find the distance to the K-th nearest neighbor for every particle
    k_nn = min(k, len(p_x))
    dists, _ = tree.query(points, k=k_nn, workers=-1)
    
    # dists[:, -1] contains the distance to the K-th neighbor
    # Enforce a tiny minimum distance to prevent division by zero
    r_comoving = np.maximum(dists[:, -1], 1e-6)
    
    # Volume of the sphere enclosing the neighbors
    vol_comoving = (4.0 / 3.0) * np.pi * (r_comoving**3)
    
    # True local comoving density = Mass / Volume
    # (Assuming all DM particles have the same mass)
    rho_comoving = (k_nn * p_mass[0]) / vol_comoving
    
    return np.max(rho_comoving)

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

def generate_dashboard(snapshot_dir, pair_dir=None):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"Error: No HDF5 snapshots found in {snapshot_dir}")
        return

    print(f"Processing {len(files)} snapshots for diagnostics...")

    # Pre-calculate Global Simulation Config & Units (from the last snapshot)
    with h5py.File(files[-1], 'r') as f:
        config = f['Config'].attrs
        units = f['Units'].attrs
        
        domain_size = config.get('domain_size', 1.0)
        mesh_size = config.get('mesh_size_1d', 32)
        has_hydro = bool(config.get('use_hydro', 0))
        hubble = config.get('Hubble_h', 70.0)
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
        T_floor = float(config.get('cooling_cutoff_k', 1000.0))
        floor_k = config.get('temp_floor_k', 10.0)
        
        u_density_cgs = units.get('unit_density_in_cgs', 1.0)
        u_energy_cgs = units.get('unit_velocity_in_cgs', 1.0)**2 * (units.get('unit_mass_in_msun', 1.0) * 1.98847e33)
        u_time_gyr = units.get('unit_time_in_gyr', 1.0)
        
        # Volume of a single cell in code units
        cell_vol_code = (domain_size / mesh_size)**3

        rho = f['Gas/density'][:] if has_hydro else None
        flat_index = np.argmax(rho)
        target_coords = np.unravel_index(flat_index, rho.shape)
        Z, Y, X = np.indices(rho.shape)

    # Tracking the densest cell at z=0
    radius_cells = 0
    dz = np.abs(Z - target_coords[0])
    dy = np.abs(Y - target_coords[1])
    dx = np.abs(X - target_coords[2])

    dz = np.minimum(dz, mesh_size - dz)
    dy = np.minimum(dy, mesh_size - dy)
    dx = np.minimum(dx, mesh_size - dx)

    dist_sq = dx**2 + dy**2 + dz**2

    # Time series lists
    scale_factors, times_gyr, dm_variances, dm_scale_factors = [], [], [], []
    p999_gas_densities, max_gas_densities, max_dm_densities = [], [], []
    p999_temperatures, max_temperatures = [], []
    rho_densest_cell, gas_temp_densest_cell, thermal_timescale_densest = [], [], []
    dt_hydro = []
    kin_energies, therm_energies, rad_energies, heat_energies, cold_gas_fracs, fractional_errors = [], [], [], [], [], []
    max_metallicity = []

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
            dt_hydro.append(header['dt_hydro'])

            # Particles
            p_x = f['Particles/position_x'][:] * box_h_mpc
            p_y = f['Particles/position_y'][:] * box_h_mpc
            p_z = f['Particles/position_z'][:] * box_h_mpc
            p_mass = f['Particles/mass'][:]
            
            # Extract the native 1D particle resolution (e.g., 32 from 32768 particles)
            num_particles = len(p_mass)
            part_mesh_size = int(np.round(num_particles**(1.0/3.0)))
            
            # Evaluate DM variance exclusively at the native particle resolution
            dm_variances.append(compute_cic_variance(p_x, p_y, p_z, p_mass, box_h_mpc, part_mesh_size))

            # Particle density estimator (expensive, so don't process every snapshot)
            if i%8 == 0:
                max_rho_comoving = compute_max_kdtree_density(p_x, p_y, p_z, p_mass, box_h_mpc, k=32)

                # Convert to physical proton densities
                dm_phys_conv = (u_density_cgs / a**3) / M_P_CGS
                max_dm_densities.append(max_rho_comoving * dm_phys_conv)
                dm_scale_factors.append(a)

            rho = f['Gas/density'][:] if has_hydro else None
            metal_density = f['Gas/metal_density'][:] if has_hydro else None

            # Snapshot Extraction for PDF & Power Spectrum
            if i in target_indices:
                k_bins, pk = compute_power_spectrum(p_x, p_y, p_z, p_mass, box_h_mpc, mesh_size, part_mesh_size, gas_rho=rho)
                pk_data[a] = (k_bins, pk)
            
            # Gas Stats
            if has_hydro:
                px, py, pz = f['Gas/momentum_x'][:], f['Gas/momentum_y'][:], f['Gas/momentum_z'][:]
                e_tot = f['Gas/energy'][:]
                temp = get_temperature(f)
                
                # Convert Max Density to Hydrogen Number Density (n_H)
                # max_rho_cgs is the comoving density. Divide by a^3 to get physical density
                max_rho_flat_index = np.argmax(rho)
                max_rho_index = np.unravel_index(max_rho_flat_index, rho.shape)
                max_rho = rho[max_rho_index]
                p999_rho = np.percentile(rho, 99.9)
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
                heat_energy_code = f['Gas'].attrs.get('cumulative_photoheating_energy', 0.0)
                heat_energies.append(heat_energy_code * u_energy_cgs * (a**2))

                current_e_code = np.sum(e_tot) * cell_vol_code
                if initial_e_code is None:
                    initial_e_code = current_e_code

                w_grav_code = f['Gas'].attrs.get('cumulative_gravitational_work', 0.0)
                w_exp_code = f['Gas'].attrs.get('cumulative_expansion_work', 0.0)

                                # Temp of the max density cell
                local_rho = np.where(dist_sq <= radius_cells**2, rho, -1.0)
                local_max_flat = np.argmax(local_rho)
                local_max_coords = np.unravel_index(local_max_flat, rho.shape)
                rho_densest_cell.append(rho[local_max_coords]*n_H_conv)
                gas_temp_densest_cell.append(temp[local_max_coords])
                thermal_timescale = f['Gas/thermal_timescale'][:]
                thermal_timescale_densest.append(thermal_timescale[local_max_coords])

                # The First Law of Thermodynamics (in Comoving Code Units)
                delta_e_code = current_e_code - initial_e_code
                absolute_error_code = delta_e_code - w_grav_code + w_exp_code + rad_energy_code - heat_energy_code
                fractional_errors.append(absolute_error_code / abs(initial_e_code) if initial_e_code != 0 else 0.0)
                
                # Cold Dense Gas Fraction
                mean_rho = np.mean(rho)
                overdensity = rho / mean_rho
                cold_dense_mask = (temp < 10000.0) & (overdensity > 100.0)
                cold_gas_fracs.append(np.sum(rho[cold_dense_mask]) / np.sum(rho))

                # Metallicity
                metallicity = metal_density / rho
                max_metallicity.append(np.max(metallicity))
                
                # Snapshot PDF Extraction
                if i in target_indices:
                    safe_od = np.maximum(overdensity, 1e-10)
                    hist, edges = np.histogram(np.log10(np.maximum(overdensity, 1e-5)), bins=100, range=(-2, 5))
                    pdf_data[a] = (hist / np.sum(hist), edges, np.var(safe_od))
                    
    # Process Paired Simulation P(k) if provided
    if pair_dir:
        print(f"\nExtracting paired P(k) from {pair_dir} for variance suppression...")
        pair_files = sorted(glob.glob(os.path.join(pair_dir, "snapshot_*.hdf5")))
        if pair_files:
            pair_target_indices = [0, len(pair_files)//2, len(pair_files)-1]
            pk_data_pair = {}
            for i in pair_target_indices:
                with h5py.File(pair_files[i], 'r') as f:
                    a_pair = f['Header'].attrs['scale_factor']
                    p_x = f['Particles/position_x'][:] * box_h_mpc
                    p_y = f['Particles/position_y'][:] * box_h_mpc
                    p_z = f['Particles/position_z'][:] * box_h_mpc
                    p_mass = f['Particles/mass'][:]
                    rho = f['Gas/density'][:] if has_hydro else None
                    
                    k_bins_p, pk_p = compute_power_spectrum(p_x, p_y, p_z, p_mass, box_h_mpc, mesh_size, part_mesh_size, gas_rho=rho)
                    pk_data_pair[a_pair] = (k_bins_p, pk_p)
            
            # Average the power spectra
            for a_val in list(pk_data.keys()):
                # Find corresponding scale factor in pair
                closest_a = min(pk_data_pair.keys(), key=lambda x: abs(x - a_val))
                if abs(closest_a - a_val) < 0.05:  # Tolerance check
                    k_bins_primary, pk_primary = pk_data[a_val]
                    _, pk_paired = pk_data_pair[closest_a]
                    
                    # Apply Angulo & Pontzen paired averaging
                    pk_data[a_val] = (k_bins_primary, (pk_primary + pk_paired) / 2.0)

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
            label_prefix = "Paired Avg" if pair_dir else "Sim"
            axs[0, 0].plot(k_bins, pk, lw=2, color=c, label=f'{label_prefix} a={a_val:.2f}')
            if camb_theory and a_val in camb_theory:
                theory_k, theory_pk = camb_theory[a_val]
                axs[0, 0].plot(theory_k, theory_pk, lw=1, linestyle='--', color=c, alpha=0.7, label=f'Theory a={a_val:.2f}')

        title_suffix = " (Variance Suppressed)" if pair_dir else ""
        axs[0, 0].set(xscale='log', yscale='log', xlabel=r'k [$h$ Mpc$^{-1}$]', ylabel=r'$P(k)$ [($h^{-1}$ Mpc)$^3$]', title=f'Matter Power Spectrum{title_suffix}')
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

        axs[0, 1].set(yscale='log', xlabel=r'$\log_{10}$ Gas Overdensity ($\rho/\bar{\rho}$)', ylabel='Volume Fraction', title='1-Point Volume-Weighted Density PDF')
        axs[0, 1].set_ylim(bottom=1e-5, top=1.5)
        axs[0, 1].legend(fontsize=8)

    # DM Variance Evolution (Structure Growth)
    #axs[0, 2].plot(scale_factors, dm_variances, color='black', label='Simulated Variance')
    #D_a = np.array([get_linear_growth(a, float(omega_m), float(omega_l)) for a in scale_factors])
    #axs[0, 2].plot(scale_factors, dm_variances[0] * (D_a / D_a[0])**2, linestyle='--', color='red', label=r'$\Lambda$CDM Linear Theory')
    #axs[0, 2].set(yscale='log', xlabel='Scale Factor (a)', ylabel=r'DM Variance ($\sigma^2$)', title='Structure Growth')
    #axs[0, 2].legend()

    # Cosmic Expansion History and Metallicity
    line_a = axs[0, 2].plot(times_gyr, scale_factors, color='purple', lw=2, label='Scale factor')
    axs[0, 2].set(xlabel='Simulation Time [Gyr]', ylabel='Scale Factor (a)', title='Cosmic Expansion History', ylim=(0.0, None))
    ax_metallicity = axs[0, 2].twinx()
    line_z = ax_metallicity.plot(times_gyr, max_metallicity, color='red', lw=2, label='Max Metallicity')
    ax_metallicity.set(ylabel='Max Metallicity', ylim=(0.0, 1.0))
    lines = line_a + line_z
    labels = [line.get_label() for line in lines]
    axs[0, 2].legend(lines, labels, loc='best')

    # Global Energy Inventory & Conservation Error
    if kin_energies:
        axs[1, 0].plot(scale_factors, kin_energies, color='green', lw=2, label='Kinetic')
        axs[1, 0].plot(scale_factors, therm_energies, color='orange', lw=2, label='Thermal')
        axs[1, 0].plot(scale_factors, rad_energies, color='red', lw=2, linestyle=':', label='Radiated')
        axs[1, 0].plot(scale_factors, heat_energies, color='blue', lw=2, linestyle=':', label='Heated')
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
        axs[1, 0].legend(lines_1 + lines_2, labels_1 + labels_2, loc='center right', fontsize=8)

    # Extreme Gas States (99.9% to Max Envelope)
    if p999_gas_densities:
        # Primary Left Axis (Gas & DM Density)
        axs[1, 1].plot(scale_factors, max_gas_densities, color='blue', lw=2, label='Max Gas Dens')
        axs[1, 1].fill_between(scale_factors, p999_gas_densities, max_gas_densities, color='blue', alpha=0.2)
        
        # Plot Max DM Density on the SAME axis
        axs[1, 1].plot(dm_scale_factors, max_dm_densities, color='black', lw=1.5, linestyle='--', alpha=0.8, label='Max DM Dens')
        
        axs[1, 1].set(yscale='log', xlabel='Scale Factor (a)', ylabel=r'Physical Density [$m_p$ cm$^{-3}$]')
        axs[1, 1].tick_params(axis='y', labelcolor='black')
        
        # Right Axis (Gas Temperature)
        ax_temp = axs[1, 1].twinx()
        ax_temp.plot(scale_factors, max_temperatures, color='red', lw=2, label='Max Temp')
        ax_temp.fill_between(scale_factors, p999_temperatures, max_temperatures, color='red', alpha=0.2)
        ax_temp.set(yscale='log', ylabel='Temperature [K]')
        ax_temp.tick_params(axis='y', labelcolor='red')
      
        # Gather legends from both axes so they appear in one box
        lines_left, labels_left = axs[1, 1].get_legend_handles_labels()
        lines_temp, labels_temp = ax_temp.get_legend_handles_labels()
        
        axs[1, 1].legend(lines_left + lines_temp, 
                         labels_left + labels_temp, 
                         loc='upper left', fontsize=8)
        axs[1, 1].set_title('Extreme States (Gas vs DM)')

    # Densest cell tracking
    # Primary Left Axis (Density)
    axs[1, 2].plot(scale_factors, rho_densest_cell, color='blue', lw=2, label='Density')
    
    # Plot Max DM Density on the SAME axis
    #axs[1, 2].plot(dm_scale_factors, max_dm_densities, color='black', lw=1.5, linestyle='--', alpha=0.8, label='Max DM Dens')
    
    axs[1, 2].set(yscale='log', xlabel='Scale Factor (a)', ylabel=r'Physical Density [$m_p$ cm$^{-3}$]')
    axs[1, 2].tick_params(axis='y', labelcolor='black')
    
    # Right Axis 1 (Gas Temperature)
    ax_temp = axs[1, 2].twinx()
    ax_temp.plot(scale_factors, gas_temp_densest_cell, color='red', lw=2, label='Temp')
    ax_temp.set(yscale='log', ylabel='Temperature [K]')
    ax_temp.tick_params(axis='y', labelcolor='red')

    # Right Axis 2 (Thermal Timescale)
    ax_time = axs[1, 2].twinx()
    ax_time.spines['right'].set_position(('outward', 50))
    
    t_therm = np.array(thermal_timescale_densest)
    dt_hydro_arr = np.array(dt_hydro) # Ensure this is a numpy array
    
    t_heating = np.where(t_therm > 0, t_therm, np.nan)
    t_cooling = np.where(t_therm < 0, np.abs(t_therm), np.nan)
    
    # Plot the timescales
    ax_time.plot(scale_factors, t_heating, color='green', lw=2, linestyle='-', label='+t_therm (Heating)')
    ax_time.plot(scale_factors, t_cooling, color='green', lw=2, linestyle='--', label='-t_therm (Cooling)')
    
    # Plot dt_hydro as the baseline threshold
    #ax_time.plot(scale_factors, dt_hydro_arr, color='orange', lw=2, linestyle=':', label='dt_hydro (Courant Limit)')
    
    ax_time.set(yscale='log', ylabel='Time [Code Units]')
    ax_time.tick_params(axis='y', labelcolor='green')

    # Gather legends from all axes so they appear in one box
    lines_left, labels_left = axs[1, 2].get_legend_handles_labels()
    lines_temp, labels_temp = ax_temp.get_legend_handles_labels()
    lines_time, labels_time = ax_time.get_legend_handles_labels()
    
    axs[1, 2].legend(lines_left + lines_temp + lines_time, 
                     labels_left + labels_temp + labels_time, 
                     loc='upper left', fontsize=8)
                     
    axs[1, 2].set_title('Densest cell (z=0) evolution')

    # Phase Diagram
    if has_hydro and 'x_data' in locals():
        # Create logarithmically spaced bins
        x_bins, y_bins = np.logspace(-2, np.log10(np.max(x_data)), 100), np.logspace(1, np.log10(np.max(y_data)), 100)
        weights = final_rho.flatten() * cell_vol_code
        h = axs[2, 0].hist2d(x_data, y_data, bins=[x_bins, y_bins], weights=weights, norm=LogNorm(), cmap='plasma')
        
        axs[2, 0].set(xscale='log', yscale='log', xlabel=r'Gas Overdensity ($\rho / \bar{\rho}$)', ylabel='Temperature [K]', title='Phase Diagram (Final Snapshot)')
        fig.colorbar(h[3], ax=axs[2, 0], label='Total Gas Mass (Code Units)')
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
        f"{'Particles':<9}: {str(int(np.cbrt(num_particles))) + '³'} ({num_particles})\n\n"
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
    parser = argparse.ArgumentParser(
        description="Generate a cosmological simulation dashboard from HDF5 snapshots."
    )
    
    parser.add_argument(
        "path", 
        type=str, 
        help="Path to the simulation directory (or base directory if using --latest)"
    )
    
    parser.add_argument(
        "-p", "--pair", 
        type=str,
        default=None,
        help="Path to a paired simulation directory (invert_phases=true) to average the Power Spectrum."
    )
    
    parser.add_argument(
        "-l", "--latest", 
        action="store_true", 
        help="Automatically find and load the most recent 'run_*' directory inside the provided path"
    )

    # Show help if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Parse the arguments
    args = parser.parse_args()
    target_dir = args.path

    # Handle the --latest flag
    if args.latest:
        search_pattern = os.path.join(target_dir, "run_*")
        runs = sorted(glob.glob(search_pattern))
        
        if runs:
            target_dir = runs[-1]
            print(f"Auto-selected latest run: {target_dir}")
        else:
            print(f"Error: No run directories found in '{args.path}'")
            sys.exit(1)

    # Run the dashboard
    generate_dashboard(target_dir, pair_dir=args.pair)