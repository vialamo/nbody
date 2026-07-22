import h5py
import numpy as np
import glob
import os
import sys
import argparse
import camb
from scipy.spatial import cKDTree

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import tempfile
import webview

# Physics & Cosmological Constants
M_P_CGS = 1.67262192e-24
X_H = 0.76  # Primordial hydrogen mass fraction
G_CGS = 6.6743e-8
K_B_CGS = 1.3806e-16

USE_THEME = "plotly_dark"

def hex_to_rgba(hex_color, opacity=0.2):
    """Converts a theme's hex color into a transparent rgba string for fills."""
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    return f'rgba({r},{g},{b},{opacity})'

def maximize(window):
    window.maximize()

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
        sigma = mesh_size / (2.0 * part_mesh_size)
        grid = gaussian_filter(grid, sigma=sigma, mode='wrap')

    if gas_rho is not None:
        cell_vol = (domain_size / N)**3
        grid += gas_rho * cell_vol
    
    # Create overdensity field (delta)
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
    points = np.vstack((p_x, p_y, p_z)).T
    tree = cKDTree(points, boxsize=box_size)
    
    k_nn = min(k, len(p_x))
    dists, _ = tree.query(points, k=k_nn, workers=-1)
    
    r_comoving = np.maximum(dists[:, -1], 1e-6)
    vol_comoving = (4.0 / 3.0) * np.pi * (r_comoving**3)
    rho_comoving = (k_nn * p_mass[0]) / vol_comoving
    
    return np.max(rho_comoving)

def get_linear_growth(a, omega_m_0, omega_l_0):
    """Calculates the linear growth factor D(a) using the Carroll, Press & Turner (1992) approximation."""
    if a == 0: return 0.0
    a3 = a**3
    omega_k_0 = 1.0 - omega_m_0 - omega_l_0
    
    E2 = omega_m_0 / a3 + omega_k_0 / (a**2) + omega_l_0
    Om_a = (omega_m_0 / a3) / E2
    Ol_a = omega_l_0 / E2
    
    g_a = (2.5 * Om_a) / (Om_a**(4.0/7.0) - Ol_a + (1.0 + Om_a / 2.0) * (1.0 + Ol_a / 70.0))
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
    
    u_code = pressure / (density * (gamma - 1.0))
    return u_code * factor_u_to_t * (a**2)

def generate_dashboard(snapshot_dir):
    files = sorted(glob.glob(os.path.join(snapshot_dir, "snapshot_*.hdf5")))
    if not files:
        print(f"Error: No HDF5 snapshots found in {snapshot_dir}")
        return

    print(f"Processing {len(files)} snapshots for diagnostics...")

    with h5py.File(files[0], 'r') as f:
        config = f['Config'].attrs
        units = f['Units'].attrs
        
        domain_size = config.get('domain_size', 1.0)
        mesh_size = config.get('mesh_size_1d', 32)
        has_hydro = bool(config.get('use_hydro', 0))
        hubble = config.get('Hubble_h', 0.7)
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
        # T_floor = float(config.get('cooling_cutoff_k', 1000.0))
        floor_k = config.get('temp_floor_k', 10.0)
        
        u_density_cgs = units.get('unit_density_in_cgs', 1.0)
        u_energy_cgs = units.get('unit_velocity_in_cgs', 1.0)**2 * (units.get('unit_mass_in_msun', 1.0) * 1.98847e33)
        u_time_gyr = units.get('unit_time_in_gyr', 1.0)
        
        cell_vol_code = (domain_size / mesh_size)**3

    scale_factors, times_gyr, dm_variances, dm_scale_factors = [], [], [], []
    p999_gas_densities, max_gas_densities, max_dm_densities = [], [], []
    p999_temperatures, max_temperatures = [], []
    kin_energies, therm_energies, rad_energies, heat_energies, cold_gas_fracs, fractional_errors = [], [], [], [], [], []

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

            p_x = f['Particles/position_x'][:] * box_h_mpc
            p_y = f['Particles/position_y'][:] * box_h_mpc
            p_z = f['Particles/position_z'][:] * box_h_mpc
            p_mass = f['Particles/mass'][:]
            
            num_particles = len(p_mass)
            part_mesh_size = int(np.round(num_particles**(1.0/3.0)))
            
            dm_variances.append(compute_cic_variance(p_x, p_y, p_z, p_mass, box_h_mpc, part_mesh_size))

            if i%8 == 0:
                max_rho_comoving = compute_max_kdtree_density(p_x, p_y, p_z, p_mass, box_h_mpc, k=32)
                dm_phys_conv = (u_density_cgs / a**3) / M_P_CGS
                max_dm_densities.append(max_rho_comoving * dm_phys_conv)
                dm_scale_factors.append(a)

            rho = f['Gas/density'][:] if has_hydro else None

            if i in target_indices:
                k_bins, pk = compute_power_spectrum(p_x, p_y, p_z, p_mass, box_h_mpc, mesh_size, part_mesh_size, gas_rho=rho)
                pk_data[a] = (k_bins, pk)
            
            if has_hydro:
                px, py, pz = f['Gas/momentum_x'][:], f['Gas/momentum_y'][:], f['Gas/momentum_z'][:]
                e_tot = f['Gas/energy'][:]
                temp = get_temperature(f)
                
                p999_rho, max_rho = np.percentile(rho, 99.9), np.max(rho)
                n_H_conv = (u_density_cgs / a**3 * X_H) / M_P_CGS
                p999_gas_densities.append(p999_rho * n_H_conv)
                max_gas_densities.append(max_rho * n_H_conv)

                p999_temperatures.append(np.percentile(temp, 99.9))
                max_temperatures.append(np.max(temp))
                
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

                delta_e_code = current_e_code - initial_e_code
                absolute_error_code = delta_e_code - w_grav_code + w_exp_code + rad_energy_code - heat_energy_code
                fractional_errors.append(absolute_error_code / abs(initial_e_code) if initial_e_code != 0 else 0.0)
                
                mean_rho = np.mean(rho)
                overdensity = rho / mean_rho
                cold_dense_mask = (temp < 10000.0) & (overdensity > 100.0)
                cold_gas_fracs.append(np.sum(rho[cold_dense_mask]) / np.sum(rho))
                
                if i in target_indices:
                    safe_od = np.maximum(overdensity, 1e-10)
                    hist, edges = np.histogram(np.log10(np.maximum(overdensity, 1e-5)), bins=100, range=(-2, 5))
                    pdf_data[a] = (hist / np.sum(hist), edges, np.var(safe_od))

    print("\nData extraction complete. Generating Plotly dashboard...")

    if has_hydro:
        with h5py.File(files[-1], 'r') as f:
            final_rho = f['Gas/density'][:]
            final_temp = get_temperature(f)
            x_data = (final_rho / np.mean(final_rho)).flatten()
            y_data = final_temp.flatten()

    # Grid setup
    theme_colors = pio.templates[USE_THEME].layout.colorway
    bg_paper = pio.templates[USE_THEME].layout.paper_bgcolor
    bg_plot = pio.templates[USE_THEME].layout.plot_bgcolor
    
    fig = make_subplots(
        rows=3, cols=3,
        subplot_titles=(
            "Matter Power Spectrum",
            "1-Point Volume-Weighted Density PDF",
            "Structure Growth",
            "Energy Inventory & Conservation",
            "Extreme States (Gas vs DM)",
            "Cosmic Expansion History",
            "Phase Diagram (Final Snapshot)",
            "Cold Dense Gas Fraction",
            "" # Info card space
        ),
        specs=[
            [{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}],
            [{"secondary_y": True},  {"secondary_y": True},  {"secondary_y": False}],
            [{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}]
        ]
    )

    # Matter Power Spectrum
    if pk_data:
        try:
            ns = float(prim_index) if isinstance(prim_index, (float, int)) else 0.96
            all_k_bins = np.concatenate([k for _, (k, _) in pk_data.items()])
            camb_theory = get_camb_pk(omega_m, omega_b, hubble, ns, sigma_8, list(pk_data.keys()), 
                                      k_min=np.min(all_k_bins), k_max=np.max(all_k_bins))
        except Exception as e:
            camb_theory = None
            print(f"CAMB theory failed: {e}")

        for idx, (a_val, (k_bins, pk)) in enumerate(pk_data.items()):
            c = theme_colors[idx % len(theme_colors)]
            fig.add_trace(go.Scatter(x=k_bins, y=pk, mode='lines', name=f'Sim a={a_val:.2f}', line=dict(color=c, width=2), legend="legend"), row=1, col=1)
            if camb_theory and a_val in camb_theory:
                theory_k, theory_pk = camb_theory[a_val]
                fig.add_trace(go.Scatter(x=theory_k, y=theory_pk, mode='lines', name=f'Theory a={a_val:.2f}', line=dict(color=c, dash='dash', width=1), legend="legend"), row=1, col=1)

        fig.update_xaxes(type='log', title_text="k [h Mpc^-1]", row=1, col=1)
        fig.update_yaxes(type='log', title_text="P(k) [(h^-1 Mpc)^3]", row=1, col=1)

    # 1-Point Volume-Weighted Density PDF
    if pdf_data:
        for idx, (a_val, (hist, edges, sigma2)) in enumerate(pdf_data.items()):
            c = theme_colors[idx % len(theme_colors)]
            centers = (edges[:-1] + edges[1:]) / 2
            dx = centers[1] - centers[0]
            
            base_mask = hist > 1e-5
            mask_hist = base_mask.copy()
            mask_hist[:-1] |= base_mask[1:]
            mask_hist[1:] |= base_mask[:-1]
            fig.add_trace(go.Scatter(x=centers[mask_hist], y=hist[mask_hist], mode='lines', name=f'PDF a={a_val:.2f}', line=dict(color=c, width=2), legend="legend2"), row=1, col=2)

            if np.any(mask_hist):
                x_min, x_max = centers[mask_hist][0], centers[mask_hist][-1]
            else:
                x_min, x_max = centers[0], centers[-1]      
            domain_mask = (centers >= x_min) & (centers <= x_max)

            Delta = 10**centers
            if idx == 0:
                p_delta = (1.0 / np.sqrt(2.0 * np.pi * sigma2)) * np.exp(-0.5 * (Delta - 1.0)**2 / sigma2)
                y_gauss = p_delta * Delta * np.log(10) * dx
                mask_gauss = (y_gauss > 1e-5) & domain_mask
                fig.add_trace(go.Scatter(x=centers[mask_gauss], y=y_gauss[mask_gauss], mode='lines', name='Gaussian Model', line=dict(color=c, dash='dot', width=1), legend="legend2"), row=1, col=2)
            
            elif idx == len(pdf_data) - 1:
                sigma2_A, mu_A = np.log(1.0 + sigma2), -np.log(1.0 + sigma2) / 2.0
                p_A = (1.0 / np.sqrt(2.0 * np.pi * sigma2_A)) * np.exp(-0.5 * (centers * np.log(10) - mu_A)**2 / sigma2_A)
                y_lognorm = p_A * np.log(10) * dx
                mask_lognorm = (y_lognorm > 1e-5) & domain_mask
                fig.add_trace(go.Scatter(x=centers[mask_lognorm], y=y_lognorm[mask_lognorm], mode='lines', name='Lognormal Model', line=dict(color=c, dash='dash', width=1), legend="legend2"), row=1, col=2)

        fig.update_yaxes(type='log', range=[-5, 0.176], title_text="Volume Fraction", row=1, col=2) # 10^0.176 is ~1.5
        fig.update_xaxes(title_text="log10 Gas Overdensity (rho/bar_rho)", row=1, col=2)

    # Structure Growth
    fig.add_trace(go.Scatter(x=scale_factors, y=dm_variances, mode='lines', name='Simulated Variance', legend="legend3"), row=1, col=3)
    D_a = np.array([get_linear_growth(a, float(omega_m), float(omega_l)) for a in scale_factors])
    fig.add_trace(go.Scatter(x=scale_factors, y=dm_variances[0] * (D_a / D_a[0])**2, mode='lines', name='LCDM Linear Theory', line=dict(dash='dash'), legend="legend3"), row=1, col=3)
    fig.update_yaxes(type='log', title_text="DM Variance (sigma^2)", row=1, col=3)
    fig.update_xaxes(title_text="Scale Factor (a)", row=1, col=3)

    # Energy Inventory
    if kin_energies:
        fig.add_trace(go.Scatter(x=scale_factors, y=kin_energies, name='Kinetic', line=dict(width=2), legend="legend4"), row=2, col=1, secondary_y=False)
        fig.add_trace(go.Scatter(x=scale_factors, y=therm_energies, name='Thermal', line=dict(width=2), legend="legend4"), row=2, col=1, secondary_y=False)
        fig.add_trace(go.Scatter(x=scale_factors, y=rad_energies, name='Radiated', line=dict(dash='dot', width=2), legend="legend4"), row=2, col=1, secondary_y=False)
        fig.add_trace(go.Scatter(x=scale_factors, y=heat_energies, name='Heated', line=dict(dash='dot', width=2), legend="legend4"), row=2, col=1, secondary_y=False)
        
        fig.add_trace(go.Scatter(x=scale_factors, y=fractional_errors, name='Fractional Error', line=dict(dash='dash', width=1.5), legend="legend4"), row=2, col=1, secondary_y=True)
        
        fig.update_yaxes(type='log', title_text="Energy Components [Ergs]", row=2, col=1, secondary_y=False)
        max_err = max(1e-4, np.max(np.abs(fractional_errors)) * 1.5)
        fig.update_yaxes(range=[-max_err, max_err], title_text="Fractional Error", row=2, col=1, secondary_y=True)
        fig.update_xaxes(title_text="Scale Factor (a)", row=2, col=1)
        fig.add_hline(y=0, line_width=1, line_dash="solid", row=2, col=1, secondary_y=True)

    # Extreme States
    if p999_gas_densities:
        c_gas = theme_colors[0]
        c_dm = theme_colors[1]
        c_temp = theme_colors[2]
        
        # Fill envelope for Gas Density
        fig.add_trace(go.Scatter(x=scale_factors, y=max_gas_densities, line=dict(width=0), showlegend=False, hoverinfo='skip'), row=2, col=2, secondary_y=False)
        fig.add_trace(go.Scatter(x=scale_factors, y=p999_gas_densities, fill='tonexty', fillcolor=hex_to_rgba(c_gas, 0.2), line=dict(color=c_gas, width=2), legend="legend5", name='99.9% Gas Dens'), row=2, col=2, secondary_y=False)
        
        # Max DM
        fig.add_trace(go.Scatter(x=dm_scale_factors, y=max_dm_densities, line=dict(color=c_dm, dash='dash', width=1.5), legend="legend5", name='Max DM Dens'), row=2, col=2, secondary_y=False)
        
        # Fill envelope for Temp (Secondary Y)
        fig.add_trace(go.Scatter(x=scale_factors, y=max_temperatures, line=dict(width=0), showlegend=False, hoverinfo='skip'), row=2, col=2, secondary_y=True)
        fig.add_trace(go.Scatter(x=scale_factors, y=p999_temperatures, fill='tonexty', fillcolor=hex_to_rgba(c_temp, 0.2), line=dict(color=c_temp, width=2), legend="legend5", name='99.9% Temp'), row=2, col=2, secondary_y=True)

        fig.update_yaxes(type='log', title_text="Physical Density [m_p cm^-3]", row=2, col=2, secondary_y=False)
        fig.update_yaxes(type='log', title_text="Temperature [K]", row=2, col=2, secondary_y=True)
        fig.update_xaxes(title_text="Scale Factor (a)", row=2, col=2)

    # Cosmic Expansion
    fig.add_trace(go.Scatter(x=times_gyr, y=scale_factors, line=dict(width=2), name='Simulation', showlegend=False), row=2, col=3)
    fig.update_yaxes(range=[0, max(scale_factors)*1.1], title_text="Scale Factor (a)", row=2, col=3)
    fig.update_xaxes(title_text="Simulation Time [Gyr]", row=2, col=3)

    # Phase Diagram
    if has_hydro and 'x_data' in locals():
        x_bins = np.logspace(-2, np.log10(np.max(x_data)), 100)
        y_bins = np.logspace(1, np.log10(np.max(y_data)), 100)
        
        # Calculate 2D histogram manually to support LogNorm coloring
        H, xedges, yedges = np.histogram2d(x_data, y_data, bins=[x_bins, y_bins], weights=final_rho.flatten())
        H_log = np.log10(np.where(H > 0, H, 1e-10))
        
        fig.add_trace(go.Heatmap(
            z=H_log.T, 
            x=xedges[:-1] + np.diff(xedges)/2, 
            y=yedges[:-1] + np.diff(yedges)/2,
            colorscale='Plasma',
            colorbar=dict(title="Log10 Total Gas Mass", x=0.22, y=0.12, len=0.25),
            showscale=True,
            showlegend=False
        ), row=3, col=1)
        
        fig.update_xaxes(type='log', title_text="Gas Overdensity (rho / bar_rho)", row=3, col=1)
        fig.update_yaxes(type='log', title_text="Temperature [K]", row=3, col=1)
    else:
        fig.add_annotation(text="Hydro Disabled", x=0.5, y=0.5, showarrow=False, font=dict(size=20), row=3, col=1)

    # Cold Dense Gas Fraction
    if cold_gas_fracs:
        fig.add_trace(go.Scatter(x=scale_factors, y=cold_gas_fracs, line=dict(width=2), name='Cold Gas Fraction', showlegend=False), row=3, col=2)
        fig.update_xaxes(title_text="Scale Factor (a)", row=3, col=2)
        fig.update_yaxes(title_text="Mass Fraction", row=3, col=2)

    # Info Card
    info_html = (
        f"<b>SIMULATION:</b> {os.path.basename(os.path.normpath(snapshot_dir))}<br><br>"
        f"<b>Box Size:</b> {box_size_mpc} Mpc<br>"
        f"<b>Grid:</b> {mesh_size}³<br>"
        f"<b>Particles:</b> {f'{int(np.cbrt(num_particles))}³ ({num_particles})' if num_particles != 'N/A' else 'N/A'}<br><br>"
        f"<b>Cosmology & Physics:</b><br>"
        f"Ω_m: {omega_m:<6} | Hubble (h): {h_val:.2f}<br>"
        f"Ω_b: {omega_b:<6} | Gamma (γ): {gamma:.3f}<br>"
        f"Ω_Λ: {omega_l:<6} | Primord. μ: {mu:.2f}<br>"
        f"Prim. Index: {prim_index:.2f}"
    )
    
    fig.add_annotation(
        text=info_html,
        x=0.5, y=0.5, showarrow=False,
        align="left",
        font=dict(family="monospace", size=14),
        borderwidth=2, borderpad=15,
        row=3, col=3
    )
    fig.update_xaxes(visible=False, row=3, col=3)
    fig.update_yaxes(visible=False, row=3, col=3)
    
    legend_style = dict(
        bgcolor='rgba(0,0,0,0)', # Transparent background
        font=dict(size=10),      # Smaller font so it fits inside the subplot
        yanchor="top",
        xanchor="right"          # Anchoring to the right keeps them out of the y-axis
    )

    # Global Layout
    fig.update_layout(
        template=USE_THEME,
        hovermode="x unified",
        showlegend=True,
        
        # Manually placing a legend box inside the domain of each subplot
        legend=dict( x=0.3, y=0.99, **legend_style),   # (1,1) P(k)
        legend2=dict(x=0.63, y=0.99, **legend_style),   # (1,2) PDF
        legend3=dict(x=0.99, y=0.99, **legend_style),   # (1,3) Structure Growth
        legend4=dict(x=0.33, y=0.62, **legend_style),   # (2,1) Energy
        legend5=dict(x=0.63, y=0.62, **legend_style)   # (2,2) Extreme States
    )
    
    # Add subtle gridlines across all subplots
    #fig.update_xaxes(showgrid=True)
    #fig.update_yaxes(showgrid=True)
    
    raw_html = fig.to_html(full_html=True)

    # Inject CSS to kill the browser's default white margins
    custom_html = raw_html.replace(
        '<head>', 
        '<head><style>body { margin: 0 !important; padding: 0 !important; overflow: hidden; }</style>'
    )

    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.html') as f:
        f.write(custom_html)
        temp_url = 'file://' + f.name

    window = webview.create_window(
        title=f"Dashboard: {os.path.basename(os.path.normpath(snapshot_dir))}", 
        url=temp_url
    )
    webview.start(maximize, window)

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
        "-l", "--latest", 
        action="store_true", 
        help="Automatically find and load the most recent 'run_*' directory inside the provided path"
    )

    # Show help if no arguments are provided at all
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
    generate_dashboard(target_dir)
