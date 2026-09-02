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
    if len(p_x) == 0:
        return grid
        
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
    if len(p_x) == 0:
        return 0.0
    grid = assign_cic_grid(p_x, p_y, p_z, p_mass, domain_size, mesh_size)
    mean_grid = np.mean(grid)
    return np.var(grid / mean_grid) if mean_grid > 0 else 0.0

def compute_power_spectrum(p_x, p_y, p_z, p_mass, domain_size, mesh_size, part_mesh_size, gas_rho=None):
    """Calculates the 1D total matter power spectrum from particle and gas data."""
    N = mesh_size
    grid = assign_cic_grid(p_x, p_y, p_z, p_mass, domain_size, mesh_size)

    if part_mesh_size > 0 and mesh_size > part_mesh_size:
        from scipy.ndimage import gaussian_filter
        sigma = mesh_size / (2.0 * part_mesh_size)
        grid = gaussian_filter(grid, sigma=sigma, mode='wrap')

    if gas_rho is not None:
        cell_vol = (domain_size / N)**3
        grid += gas_rho * cell_vol
    
    mean_rho = np.mean(grid)
    delta = (grid - mean_rho) / mean_rho if mean_rho > 0 else grid
    
    fft_delta = np.fft.fftn(delta)
    box_vol = domain_size**3
    power_3d = (np.abs(fft_delta)**2) * (box_vol / (float(N)**6))
    
    k_1d = np.fft.fftfreq(N, d=domain_size/N) * 2.0 * np.pi
    kx, ky, kz = np.meshgrid(k_1d, k_1d, k_1d, indexing='ij')
    k_mag = np.sqrt(kx**2 + ky**2 + kz**2)
    
    k_mag_flat = k_mag.flatten()
    power_flat = power_3d.flatten()
    valid = k_mag_flat > 0
    k_mag_flat = k_mag_flat[valid]
    power_flat = power_flat[valid]
    
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
    if len(p_x) == 0:
        return 0.0
    points = np.vstack((p_x, p_y, p_z)).T
    tree = cKDTree(points, boxsize=box_size)
    k_nn = min(k, len(p_x))
    dists, _ = tree.query(points, k=k_nn, workers=-1)
    
    r_comoving = np.maximum(dists[:, -1], 1e-6)
    vol_comoving = (4.0 / 3.0) * np.pi * (r_comoving**3)
    rho_comoving = (k_nn * p_mass[0]) / vol_comoving
    
    return np.max(rho_comoving)

def compute_kdtree_density_array(p_x, p_y, p_z, p_mass, box_size, k=32):
    """Calculates the comoving density of all particles using a KD-Tree nearest-neighbor search."""
    if len(p_x) == 0:
        return np.array([])
    points = np.vstack((p_x, p_y, p_z)).T
    tree = cKDTree(points, boxsize=box_size)
    k_nn = min(k, len(p_x))
    dists, _ = tree.query(points, k=k_nn, workers=-1)
    
    r_comoving = np.maximum(dists[:, -1], 1e-6)
    vol_comoving = (4.0 / 3.0) * np.pi * (r_comoving**3)
    rho_comoving = (k_nn * p_mass[0]) / vol_comoving
    return rho_comoving

def get_temperature(f):
    """Extracts or calculates temperature in Kelvin."""
    if 'Gas/temperature' in f:
        return f['Gas/temperature'][:]
    
    pressure = f['Gas/pressure'][:]
    density = f['Gas/density'][:]
    gamma = f['Config'].attrs.get('gamma', 5.0/3.0)
    factor_u_to_t = f['Units'].attrs['factor_u_to_t']
    a = f['Header'].attrs['scale_factor']
    
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
        method = config.get('hydro_method', b"none").decode('utf-8')
        has_eulerian_hydro = bool(method == "eulerian")
        has_particle_hydro = bool(method == "mfm")
        has_hydro = has_eulerian_hydro or has_particle_hydro
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
        
        u_density_cgs = units.get('unit_density_in_cgs', 1.0)
        u_energy_cgs = units.get('unit_velocity_in_cgs', 1.0)**2 * (units.get('unit_mass_in_msun', 1.0) * 1.98847e33)
        u_time_gyr = units.get('unit_time_in_gyr', 1.0)
        cell_vol_code = (domain_size / mesh_size)**3

        if has_eulerian_hydro:
            rho = f['Gas/density'][:]
            flat_index = np.argmax(rho)
            target_coords = np.unravel_index(flat_index, rho.shape)
            Z, Y, X = np.indices(rho.shape)
            radius_cells = 0
            dz = np.minimum(np.abs(Z - target_coords[0]), mesh_size - np.abs(Z - target_coords[0]))
            dy = np.minimum(np.abs(Y - target_coords[1]), mesh_size - np.abs(Y - target_coords[1]))
            dx = np.minimum(np.abs(X - target_coords[2]), mesh_size - np.abs(X - target_coords[2]))
            dist_sq = dx**2 + dy**2 + dz**2
        elif has_particle_hydro:
            rho = f['Gas/density'][:]
            max_idx = np.argmax(rho)
            target_pos = np.array([f['Gas/position_x'][max_idx], 
                                   f['Gas/position_y'][max_idx], 
                                   f['Gas/position_z'][max_idx]])

    # Time series lists
    scale_factors, times_gyr, dm_variances, dm_scale_factors = [], [], [], []
    p999_gas_densities, max_gas_densities, max_dm_densities = [], [], []
    p999_temperatures, max_temperatures = [], []
    rho_densest_cell, gas_temp_densest_cell, thermal_timescale_densest = [], [], []
    dt_hydro = []
    
    # Energy Arrays
    kin_energies, therm_energies, rad_energies, heat_energies, switch_energies = [], [], [], [], []
    fractional_errors = []
    ke_dm_list, fractional_errors_dm = [], []
    cold_gas_fracs, max_metallicity = [], []

    pdf_data, pk_data = {}, {}
    target_indices = [0, len(files)//2, len(files)-1]
    initial_e_code = None
    initial_e_dm = None
    num_dm_particles_total = 0

    for i, f_path in enumerate(files):
        sys.stdout.write(f"\rProcessing snapshot {i+1}/{len(files)} [{(i + 1) / len(files) * 100:.1f}%]")
        sys.stdout.flush()

        with h5py.File(f_path, 'r') as f:
            header = f['Header'].attrs
            a = header['scale_factor']          
            scale_factors.append(a)
            times_gyr.append(header.get('simulation_time', 0.0) * u_time_gyr)
            dt_hydro.append(header['dt_hydro'])
            
            energy_conv = u_energy_cgs * (a**2)

            # DM Particles
            p_x = f['Particles/position_x'][:] * box_h_mpc
            p_y = f['Particles/position_y'][:] * box_h_mpc
            p_z = f['Particles/position_z'][:] * box_h_mpc
            p_mass = f['Particles/mass'][:]
            num_dm_particles = len(p_mass)
            num_dm_particles_total = num_dm_particles
            
            if num_dm_particles > 0:
                p_vx = f['Particles/velocity_x'][:]
                p_vy = f['Particles/velocity_y'][:]
                p_vz = f['Particles/velocity_z'][:]
                
                ke_dm_code = np.sum(0.5 * p_mass * (p_vx**2 + p_vy**2 + p_vz**2))
                ke_dm_list.append(ke_dm_code * energy_conv)

                w_grav_dm = f['Particles'].attrs.get('cumulative_gravitational_work', 0.0)
                w_exp_dm = f['Particles'].attrs.get('cumulative_expansion_work', 0.0)

                if initial_e_dm is None:
                    initial_e_dm = ke_dm_code
                
                delta_e_dm = ke_dm_code - initial_e_dm
                abs_err_dm = delta_e_dm - w_grav_dm + w_exp_dm
                fractional_errors_dm.append(abs_err_dm / abs(initial_e_dm) if initial_e_dm != 0 else 0.0)

                dm_part_mesh_size = int(np.round(num_dm_particles**(1.0/3.0)))
                dm_variances.append(compute_cic_variance(p_x, p_y, p_z, p_mass, box_h_mpc, dm_part_mesh_size))

                if i % 8 == 0:
                    max_rho_comoving = compute_max_kdtree_density(p_x, p_y, p_z, p_mass, box_h_mpc, k=32)
                    dm_phys_conv = (u_density_cgs / a**3) / M_P_CGS
                    max_dm_densities.append(max_rho_comoving * dm_phys_conv)
                    dm_scale_factors.append(a)

            # Total Power Spectrum (DM + Gas)
            tot_px, tot_py, tot_pz, tot_pmass = p_x, p_y, p_z, p_mass
            if has_particle_hydro:
                gx = f['Gas/position_x'][:] * box_h_mpc
                gy = f['Gas/position_y'][:] * box_h_mpc
                gz = f['Gas/position_z'][:] * box_h_mpc
                gm = f['Gas/mass'][:]
                tot_px = np.concatenate([p_x, gx]) if len(p_x) > 0 else gx
                tot_py = np.concatenate([p_y, gy]) if len(p_y) > 0 else gy
                tot_pz = np.concatenate([p_z, gz]) if len(p_z) > 0 else gz
                tot_pmass = np.concatenate([p_mass, gm]) if len(p_mass) > 0 else gm
                
            tot_particles = len(tot_pmass)
            tot_part_mesh_size = int(np.round(tot_particles**(1.0/3.0))) if tot_particles > 0 else 0
            rho_grid = f['Gas/density'][:] if has_eulerian_hydro else None

            if i in target_indices:
                k_bins, pk = compute_power_spectrum(tot_px, tot_py, tot_pz, tot_pmass, box_h_mpc, mesh_size, tot_part_mesh_size, gas_rho=rho_grid)
                pk_data[a] = (k_bins, pk)

                if not has_hydro and num_dm_particles > 0:
                    dm_densities = compute_kdtree_density_array(p_x, p_y, p_z, p_mass, box_h_mpc, k=32)
                    mean_dm_density = np.sum(p_mass) / (box_h_mpc**3)
                    if mean_dm_density > 0:
                        dm_overdensity = dm_densities / mean_dm_density
                        safe_od = np.maximum(dm_overdensity, 1e-10)
                        hist, edges = np.histogram(np.log10(np.maximum(dm_overdensity, 1e-5)), bins=100, range=(-2, 5))
                        pdf_data[a] = (hist / np.sum(hist), edges, np.var(safe_od))
            
            # Gas Stats
            num_gas_particles = 0
            if has_hydro:
                rho = f['Gas/density'][:]
                temp = get_temperature(f)
                n_H_conv = (u_density_cgs / a**3 * X_H) / M_P_CGS

                if has_eulerian_hydro:
                    px_g, py_g, pz_g = f['Gas/momentum_x'][:], f['Gas/momentum_y'][:], f['Gas/momentum_z'][:]
                    e_tot = f['Gas/energy'][:]
                    metal_density = f['Gas/metal_density'][:]
                    metallicity = metal_density / rho
                    gas_mass = rho * cell_vol_code
                    
                    safe_rho = np.maximum(rho, 1e-10)
                    ke_grid = (px_g**2 + py_g**2 + pz_g**2) / (2.0 * safe_rho)
                    
                    kin_energies.append(np.sum(ke_grid) * energy_conv * cell_vol_code)
                    therm_energies.append((np.sum(e_tot) - np.sum(ke_grid)) * energy_conv * cell_vol_code)
                    current_e_code = np.sum(e_tot) * cell_vol_code
                    
                    local_rho = np.where(dist_sq <= radius_cells**2, rho, -1.0)
                    local_max_coords = np.unravel_index(np.argmax(local_rho), rho.shape)
                    rho_densest_cell.append(rho[local_max_coords] * n_H_conv)
                    gas_temp_densest_cell.append(temp[local_max_coords])
                    thermal_timescale_densest.append(f['Gas/thermal_timescale'][local_max_coords])
                    mean_rho = np.mean(rho)

                    switch_energy_code = f['Gas'].attrs.get('cumulative_dual_energy_switch_energy', 0.0)
                    switch_energies.append(switch_energy_code)

                elif has_particle_hydro:
                    vx, vy, vz = f['Gas/velocity_x'][:], f['Gas/velocity_y'][:], f['Gas/velocity_z'][:]
                    u_int = f['Gas/internal_energy'][:]
                    gas_mass = f['Gas/mass'][:]
                    metallicity = f['Gas/metal_fraction'][:]
                    num_gas_particles = len(gas_mass)
                    
                    ke_array = 0.5 * gas_mass * (vx**2 + vy**2 + vz**2)
                    
                    kin_energies.append(np.sum(ke_array) * energy_conv)
                    therm_energies.append(np.sum(gas_mass * u_int) * energy_conv)
                    current_e_code = np.sum(ke_array) + np.sum(gas_mass * u_int)
                    
                    gx_g, gy_g, gz_g = f['Gas/position_x'][:], f['Gas/position_y'][:], f['Gas/position_z'][:]
                    d2 = (gx_g - target_pos[0])**2 + (gy_g - target_pos[1])**2 + (gz_g - target_pos[2])**2
                    closest_idx = np.argmin(d2)
                    rho_densest_cell.append(rho[closest_idx] * n_H_conv)
                    gas_temp_densest_cell.append(temp[closest_idx])
                    thermal_timescale_densest.append(np.nan) 
                    mean_rho = np.sum(gas_mass) / domain_size**3

                    switch_energy_code = f['Gas'].attrs.get('cumulative_entropy_switch_energy', 0.0)
                    switch_energies.append(switch_energy_code)

                p999_gas_densities.append(np.percentile(rho, 99.9) * n_H_conv)
                max_gas_densities.append(np.max(rho) * n_H_conv)
                p999_temperatures.append(np.percentile(temp, 99.9))
                max_temperatures.append(np.max(temp))
                max_metallicity.append(np.max(metallicity))

                rad_energy_code = f['Gas'].attrs.get('cumulative_radiated_energy', 0.0)
                rad_energies.append(rad_energy_code * energy_conv)
                heat_energy_code = f['Gas'].attrs.get('cumulative_photoheating_energy', 0.0)
                heat_energies.append(heat_energy_code * energy_conv)

                if initial_e_code is None:
                    initial_e_code = current_e_code

                w_grav_code = f['Gas'].attrs.get('cumulative_gravitational_work', 0.0)
                w_exp_code = f['Gas'].attrs.get('cumulative_expansion_work', 0.0)
                
                delta_e_code = current_e_code - initial_e_code
                absolute_error_code = delta_e_code - w_grav_code + w_exp_code + rad_energy_code - heat_energy_code
                fractional_errors.append(absolute_error_code / abs(initial_e_code) if initial_e_code != 0 else 0.0)

                overdensity = rho / mean_rho
                cold_dense_mask = (temp < 10000.0) & (overdensity > 100.0)
                cold_gas_fracs.append(np.sum(gas_mass[cold_dense_mask]) / np.sum(gas_mass))
                
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
                    
                    tot_px, tot_py, tot_pz, tot_pmass = p_x, p_y, p_z, p_mass
                    if has_particle_hydro:
                        gx = f['Gas/position_x'][:] * box_h_mpc
                        gy = f['Gas/position_y'][:] * box_h_mpc
                        gz = f['Gas/position_z'][:] * box_h_mpc
                        gm = f['Gas/mass'][:]
                        tot_px = np.concatenate([p_x, gx]) if len(p_x) > 0 else gx
                        tot_py = np.concatenate([p_y, gy]) if len(p_y) > 0 else gy
                        tot_pz = np.concatenate([p_z, gz]) if len(p_z) > 0 else gz
                        tot_pmass = np.concatenate([p_mass, gm]) if len(p_mass) > 0 else gm

                    tot_particles = len(tot_pmass)
                    tot_part_mesh_size = int(np.round(tot_particles**(1.0/3.0))) if tot_particles > 0 else 0
                    rho_grid = f['Gas/density'][:] if has_eulerian_hydro else None
                    
                    k_bins_p, pk_p = compute_power_spectrum(tot_px, tot_py, tot_pz, tot_pmass, box_h_mpc, mesh_size, tot_part_mesh_size, gas_rho=rho_grid)
                    pk_data_pair[a_pair] = (k_bins_p, pk_p)
            
            for a_val in list(pk_data.keys()):
                closest_a = min(pk_data_pair.keys(), key=lambda x: abs(x - a_val))
                if abs(closest_a - a_val) < 0.05:
                    k_bins_primary, pk_primary = pk_data[a_val]
                    _, pk_paired = pk_data_pair[closest_a]
                    pk_data[a_val] = (k_bins_primary, (pk_primary + pk_paired) / 2.0)

    print("\nData extraction complete. Generating plots...")

    # Last snapshot phase data
    if has_hydro:
        with h5py.File(files[-1], 'r') as f:
            final_rho = f['Gas/density'][:]
            final_temp = get_temperature(f)
            
            if has_eulerian_hydro:
                mean_rho = np.mean(final_rho)
                weights = (final_rho * cell_vol_code).flatten()
            else:
                final_mass = f['Gas/mass'][:]
                mean_rho = np.sum(final_mass) / domain_size**3
                weights = final_mass.flatten()
                
            x_data = (final_rho / mean_rho).flatten()
            y_data = final_temp.flatten()

    # Plotting Grid Setup
    fig, axs = plt.subplots(3, 3, figsize=(18, 14))
    fig.canvas.manager.set_window_title(f"Dashboard: {os.path.basename(os.path.normpath(snapshot_dir))}")
    
    # Matter Power Spectrum 
    if pk_data:
        try:
            ns = float(prim_index) if isinstance(prim_index, (float, int)) else 0.96
            all_k_bins = np.concatenate([k_bins for _, (k_bins, _) in pk_data.items()])
            camb_theory = get_camb_pk(omega_m, omega_b, hubble, ns, sigma_8, list(pk_data.keys()), 
                                      k_min=np.min(all_k_bins), k_max=np.max(all_k_bins))
        except Exception as e:
            camb_theory = None

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

    # Density PDF
    if pdf_data:
        colors = ['blue', 'green', 'red']
        for idx, (a_val, (hist, edges, sigma2)) in enumerate(pdf_data.items()):
            c = colors[idx % len(colors)]
            centers = (edges[:-1] + edges[1:]) / 2
            
            base_mask = hist > 1e-5
            mask_hist = base_mask.copy()
            mask_hist[:-1] |= base_mask[1:]
            mask_hist[1:] |= base_mask[:-1]
            axs[0, 1].plot(centers[mask_hist], hist[mask_hist], lw=2, color=c, label=f'Sim a={a_val:.2f}')

        x_label = r'$\log_{10}$ Gas Overdensity ($\rho/\bar{\rho}$)' if has_hydro else r'$\log_{10}$ DM Overdensity ($\rho/\bar{\rho}$)'
        axs[0, 1].set(yscale='log', xlabel=x_label, ylabel='Volume Fraction', title='1-Point Volume-Weighted Density PDF')
        axs[0, 1].set_ylim(bottom=1e-5, top=1.5)
        axs[0, 1].legend(fontsize=8)

    # Cosmic Expansion History
    line_a = axs[0, 2].plot(times_gyr, scale_factors, color='purple', lw=2, label='Scale factor')
    axs[0, 2].set(xlabel='Simulation Time [Gyr]', ylabel='Scale Factor (a)', title='Cosmic Expansion History', ylim=(0.0, None))
    lines = line_a
    if has_hydro and max_metallicity:
        ax_metallicity = axs[0, 2].twinx()
        line_z = ax_metallicity.plot(times_gyr, max_metallicity, color='red', lw=2, label='Max Metallicity')
        ax_metallicity.set(ylabel='Max Metallicity', ylim=(0.0, 1.0))
        lines = line_a + line_z
    labels = [line.get_label() for line in lines]
    axs[0, 2].legend(lines, labels, loc='best')

    # Global Energy Inventory 
    if kin_energies or ke_dm_list:
        if kin_energies:
            axs[1, 0].plot(scale_factors, kin_energies, color='green', lw=2, label='Gas Kinetic')
            axs[1, 0].plot(scale_factors, therm_energies, color='orange', lw=2, label='Gas Thermal')
            axs[1, 0].plot(scale_factors, rad_energies, color='blue', lw=2, linestyle=':', label='Gas Radiated')
            axs[1, 0].plot(scale_factors, heat_energies, color='red', lw=2, linestyle=':', label='Gas Heated')
        
        if ke_dm_list:
            axs[1, 0].plot(scale_factors, ke_dm_list, color='darkblue', lw=2, linestyle='-', label='DM Kinetic')

        axs[1, 0].set(yscale='log', xlabel='Scale Factor (a)', ylabel='Energy Components [Ergs]', title='Energy Inventory & Conservation')
        
        ax_err = axs[1, 0].twinx()
        all_errs = []
        if fractional_errors:
            ax_err.plot(scale_factors, fractional_errors, color='black', lw=1.5, linestyle='--', label='Gas Frac Error')
            all_errs.extend(fractional_errors)
        if fractional_errors_dm:
            ax_err.plot(scale_factors, fractional_errors_dm, color='purple', lw=1.5, linestyle='-.', label='DM Frac Error')
            all_errs.extend(fractional_errors_dm)

        if switch_energies:
            switch_arr = np.array(switch_energies)
            switch_frac = switch_arr / abs(initial_e_code) if initial_e_code != 0 else np.zeros_like(switch_arr)
            ax_err.plot(scale_factors, switch_frac, color='magenta', lw=1.5, linestyle=':', label='Energy Switch Drift')
            all_errs.extend(switch_frac)
            
        if all_errs:
            max_err = max(1e-4, np.max(np.abs(all_errs)) * 1.5)
            ax_err.set_ylim(-max_err, max_err)
        ax_err.set_ylabel(r'Fractional Error ($\Delta E / E_0$)', color='black')
        ax_err.axhline(0, color='gray', linestyle='-', linewidth=0.5)

        lines_1, labels_1 = axs[1, 0].get_legend_handles_labels()
        lines_2, labels_2 = ax_err.get_legend_handles_labels()
        axs[1, 0].legend(lines_1 + lines_2, labels_1 + labels_2, loc='center right', fontsize=8)
    else:
        axs[1, 0].text(0.5, 0.5, 'No Energy Data', ha='center', va='center')

    # Extreme States
    if p999_gas_densities or max_dm_densities:
        if p999_gas_densities:
            axs[1, 1].plot(scale_factors, max_gas_densities, color='blue', lw=2, label='Max Gas Dens')
            axs[1, 1].fill_between(scale_factors, p999_gas_densities, max_gas_densities, color='blue', alpha=0.2)
        
        if len(max_dm_densities) > 0:
            axs[1, 1].plot(dm_scale_factors, max_dm_densities, color='black', lw=1.5, linestyle='--', alpha=0.8, label='Max DM Dens')
        
        axs[1, 1].set(yscale='log', xlabel='Scale Factor (a)', ylabel=r'Physical Density [$m_p$ cm$^{-3}$]')
        axs[1, 1].tick_params(axis='y', labelcolor='black')
        
        lines_left, labels_left = axs[1, 1].get_legend_handles_labels()
        
        if p999_temperatures:
            ax_temp = axs[1, 1].twinx()
            ax_temp.plot(scale_factors, max_temperatures, color='red', lw=2, label='Max Temp')
            ax_temp.fill_between(scale_factors, p999_temperatures, max_temperatures, color='red', alpha=0.2)
            ax_temp.set(yscale='log', ylabel='Temperature [K]')
            ax_temp.tick_params(axis='y', labelcolor='red')
            lines_temp, labels_temp = ax_temp.get_legend_handles_labels()
            lines_left += lines_temp
            labels_left += labels_temp
        
        axs[1, 1].legend(lines_left, labels_left, loc='upper left', fontsize=8)
        title = 'Extreme States (Gas vs DM)' if (p999_gas_densities and max_dm_densities) else ('Extreme States (Gas)' if p999_gas_densities else 'Extreme States (DM)')
        axs[1, 1].set_title(title)
    else:
        axs[1, 1].text(0.5, 0.5, 'No Extreme States Data', ha='center', va='center')

    # Densest Cell Evolution
    if has_hydro and has_eulerian_hydro:
        axs[1, 2].plot(scale_factors, rho_densest_cell, color='blue', lw=2, label='Density')
        axs[1, 2].set(yscale='log', xlabel='Scale Factor (a)', ylabel=r'Physical Density [$m_p$ cm$^{-3}$]')
        axs[1, 2].tick_params(axis='y', labelcolor='black')
        
        ax_temp = axs[1, 2].twinx()
        ax_temp.plot(scale_factors, gas_temp_densest_cell, color='red', lw=2, label='Temp')
        ax_temp.set(yscale='log', ylabel='Temperature [K]')
        ax_temp.tick_params(axis='y', labelcolor='red')

        ax_time = axs[1, 2].twinx()
        ax_time.spines['right'].set_position(('outward', 50))
        
        t_therm = np.array(thermal_timescale_densest)
        t_heating = np.where(t_therm > 0, t_therm, np.nan)
        t_cooling = np.where(t_therm < 0, np.abs(t_therm), np.nan)
        
        ax_time.plot(scale_factors, t_heating, color='green', lw=2, linestyle='-', label='+t_therm (Heating)')
        ax_time.plot(scale_factors, t_cooling, color='green', lw=2, linestyle='--', label='-t_therm (Cooling)')
        
        ax_time.set(yscale='log', ylabel='Time [Code Units]')
        ax_time.tick_params(axis='y', labelcolor='green')

        lines_left, labels_left = axs[1, 2].get_legend_handles_labels()
        lines_temp, labels_temp = ax_temp.get_legend_handles_labels()
        lines_time, labels_time = ax_time.get_legend_handles_labels()
        
        axs[1, 2].legend(lines_left + lines_temp + lines_time, labels_left + labels_temp + labels_time, loc='upper left', fontsize=8)
        axs[1, 2].set_title('Densest cell (z=0) evolution')
    else:
        axs[1, 2].text(0.5, 0.5, 'Graph Disabled', ha='center', va='center')

    # Phase Diagram
    if has_hydro and 'x_data' in locals():
        max_x = max(np.max(x_data), 1.01e-2)
        max_y = max(np.max(y_data), 10.1)
        x_bins = np.logspace(-2, np.log10(max_x), 100)
        y_bins = np.logspace(1, np.log10(max_y), 100)
        
        counts, _, _ = np.histogram2d(x_data, y_data, bins=[x_bins, y_bins], weights=weights)
        valid_counts = counts[counts > 0]
        
        if len(valid_counts) > 0:
            c_min, c_max = np.min(valid_counts), np.max(valid_counts)
            if c_min == c_max:
                c_min = c_max * 0.1
        else:
            c_min, c_max = 1e-10, 1.0
            
        h = axs[2, 0].hist2d(x_data, y_data, bins=[x_bins, y_bins], weights=weights, 
                             norm=LogNorm(vmin=c_min, vmax=c_max), cmap='plasma')
        
        axs[2, 0].set(xscale='log', yscale='log', xlabel=r'Gas Overdensity ($\rho / \bar{\rho}$)', ylabel='Temperature [K]', title='Phase Diagram (Final Snapshot)')
        fig.colorbar(h[3], ax=axs[2, 0], label='Total Gas Mass (Code Units)')
    else:
        axs[2, 0].text(0.5, 0.5, 'Hydro Disabled', ha='center', va='center')

    # Cold Dense Gas Fraction
    if cold_gas_fracs:
        axs[2, 1].plot(scale_factors, cold_gas_fracs, color='teal', lw=2)
        axs[2, 1].set(xlabel='Scale Factor (a)', ylabel='Mass Fraction', title=r'Cold Dense Gas ($T<10^4$ K, $\delta>100$)')
    else:
        axs[2, 1].text(0.5, 0.5, 'Hydro Disabled', ha='center', va='center')

    # Info Card
    axs[2, 2].axis('off')
    n_dm_1d = str(int(np.cbrt(num_dm_particles_total))) + '³' if num_dm_particles_total > 0 else '0³'
    n_gas_1d = str(int(np.cbrt(num_gas_particles))) + '³' if num_gas_particles > 0 else '0³'

    info_text = (
        f"SIMULATION: {os.path.basename(os.path.normpath(snapshot_dir))}\n\n"
        f"{'Box Size':<9}: {str(box_size_mpc) + ' Mpc'}\n"
        f"{'Grid':<9}: {str(mesh_size) + '³'}\n"
        f"{'Particles':<9}: {n_dm_1d} ({num_dm_particles_total})\n"
        f"{'Gas particles':<9}: {n_gas_1d} ({num_gas_particles})\n"
        f"{'Hydro method':<9}: {method}\n\n"
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
    parser = argparse.ArgumentParser(description="Generate a cosmological simulation dashboard from HDF5 snapshots.")
    parser.add_argument("path", type=str, help="Path to the simulation directory")
    parser.add_argument("-p", "--pair", type=str, default=None, help="Path to a paired simulation directory")
    parser.add_argument("-l", "--latest", action="store_true", help="Automatically find the most recent 'run_*' directory")

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

    generate_dashboard(target_dir, pair_dir=args.pair)