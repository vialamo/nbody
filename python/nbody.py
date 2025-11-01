import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

import math
import time
import datetime
import numpy as np
import pygame
import matplotlib.cm
import scipy.ndimage  # Image scaling
import h5py # Snapshot saving
import configparser
import sys

# ------------------------
# Config Class
# ------------------------
class Config:
    """Loads all simulation parameters from a .ini file."""
    def __init__(self, config_filename):
        config = configparser.ConfigParser()
        if not os.path.exists(config_filename):
            raise FileNotFoundError(f"Config file '{config_filename}' not found.")
        
        config.read(config_filename)
        
        # --- Store all parameters ---
        self.params = {} # For saving to HDF5
        
        # [domain]
        self.domain_size = self._get(config, 'domain', 'domain_size', 'float')
        self.mesh_size = self._get(config, 'domain', 'mesh_size', 'int')
        
        # [matter]
        self.omega_baryon = self._get(config, 'matter', 'omega_baryon', 'float')
        self.num_dm_particles_side = self._get(config, 'matter', 'num_dm_particles_side', 'int')
        self.enable_hydro = self._get(config, 'matter', 'enable_hydro', 'bool')
        
        # [physics]
        self.expansion_start_t = self._get(config, 'physics', 'expansion_start_t', 'float')
        self.expanding_universe = self._get(config, 'physics', 'expanding_universe', 'bool')
        self.g_const = self._get(config, 'physics', 'g_const', 'float')
        self.gamma = self._get(config, 'physics', 'gamma', 'float')
        self.initial_power_spectrum_index = self._get(config, 'physics', 'initial_power_spectrum_index', 'float')

        # [p3m]
        self.cutoff_radius_cells = self._get(config, 'p3m', 'cutoff_radius_cells', 'float')
        self.cutoff_transition_width_factor = self._get(config, 'p3m', 'cutoff_transition_width_factor', 'float')

        # [time]
        self.dt_factor = self._get(config, 'time', 'dt_factor', 'float')
        
        # [output]
        self.debug_info_every_cycles = self._get(config, 'output', 'debug_info_every_cycles', 'int')
        self.save_snapshot_every_cycles = self._get(config, 'output', 'save_snapshot_every_cycles', 'int')
        self.seed = self._get(config, 'output', 'seed', 'int')

        # [visualization]
        self.render_size = self._get(config, 'visualization', 'render_size', 'int')

        # --- Calculate derived parameters ---
        self.cell_size = self.domain_size / self.mesh_size
        self.omega_dm = 1.0 - self.omega_baryon
        
        self.cutoff_radius = self.cutoff_radius_cells * self.cell_size
        self.cutoff_radius_squared = self.cutoff_radius**2
        self.cutoff_transition_width = self.cutoff_transition_width_factor * self.cutoff_radius
        self.r_switch_start = self.cutoff_radius - self.cutoff_transition_width
        self.r_switch_start_sq = self.r_switch_start**2
        
        self.num_dm_particles = self.num_dm_particles_side ** 2
        self.dm_particle_mass = self.omega_dm / self.num_dm_particles
        self.gas_total_mass = self.omega_baryon
        
        self.mean_interparticle_spacing = self.domain_size / math.sqrt(self.num_dm_particles)
        self.softening_squared = (self.mean_interparticle_spacing / 50.0)**2
        
        self.initial_hubble_param = (2.0/3.0) / self.expansion_start_t if self.expanding_universe else 0.0
        self.initial_scale_factor = self.expansion_start_t**(2.0/3.0) if self.expanding_universe else 0.0
        
        self.dynamical_time = 1.0 / math.sqrt(self.g_const)
        self.dt = self.dt_factor * self.dynamical_time
        
        self.render_scale = self.render_size / self.domain_size

        # Set seed
        np.random.seed(self.seed)

    def _get(self, config, section, option, dtype='str'):
        """Helper to get config value and store it for HDF5 saving."""
        if dtype == 'float':
            val = config.getfloat(section, option)
        elif dtype == 'int':
            val = config.getint(section, option)
        elif dtype == 'bool':
            val = config.getboolean(section, option)
        else:
            val = config.get(section, option)
        
        self.params[f"{section}_{option}"] = val
        return val

    def get_all_params(self):
        """Get all *base* parameters for saving to HDF5."""
        return self.params

# ------------------------
# Simulation Classes
# ------------------------

class Particle:
    def __init__(self, x, y, mass, vx=0, vy=0):
        self.x = x
        self.y = y
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.ax = 0
        self.ay = 0

class GasGrid:
    """Represents the Eulerian grid for the baryonic gas."""
    def __init__(self, config):
        self.config = config
        mesh_size = config.mesh_size
        
        # Conservative variables
        self.density = np.zeros((mesh_size, mesh_size))
        self.momentum_x = np.zeros((mesh_size, mesh_size))
        self.momentum_y = np.zeros((mesh_size, mesh_size))
        self.energy = np.zeros((mesh_size, mesh_size)) # Total energy density E = rho*e + 0.5*rho*v^2
        
        # Primitive variables (derived)
        self.pressure = np.zeros((mesh_size, mesh_size))
        self.velocity_x = np.zeros((mesh_size, mesh_size))
        self.velocity_y = np.zeros((mesh_size, mesh_size))

        self.accel_x = np.zeros((mesh_size, mesh_size))
        self.accel_y = np.zeros((mesh_size, mesh_size))

    def update_primitive_variables(self):
        """Calculates pressure and velocity from the conservative variables."""
        # Avoid division by zero in empty cells
        non_zero_density = self.density > 1e-9
        
        self.velocity_x.fill(0.0)
        self.velocity_y.fill(0.0)
        self.velocity_x[non_zero_density] = self.momentum_x[non_zero_density] / self.density[non_zero_density]
        self.velocity_y[non_zero_density] = self.momentum_y[non_zero_density] / self.density[non_zero_density]

        kinetic_energy_density = np.zeros_like(self.density)
        kinetic_energy_density[non_zero_density] = 0.5 * (self.momentum_x[non_zero_density]**2 + self.momentum_y[non_zero_density]**2) / self.density[non_zero_density]
        
        internal_energy_density = self.energy - kinetic_energy_density
        
        self.pressure = (self.config.gamma - 1.0) * internal_energy_density
        self.pressure[self.pressure < 1e-9] = 1e-9 # Pressure floor

    def calculate_fluxes(self, axis):
        """Internal helper to calculate HLL fluxes along a given axis (0 for x, 1 for y)."""
        
        if axis == 1: # Permute axes for y-direction
            density, mom_n, mom_t, energy, v_n, v_t, pressure = \
                self.density.T, self.momentum_y.T, self.momentum_x.T, self.energy.T, \
                self.velocity_y.T, self.velocity_x.T, self.pressure.T
        else: # x-direction
            density, mom_n, mom_t, energy, v_n, v_t, pressure = \
                self.density, self.momentum_x, self.momentum_y, self.energy, \
                self.velocity_x, self.velocity_y, self.pressure

        # Left and Right states for all cells at once using roll
        rho_L, p_L, vn_L, vt_L, E_L, mom_n_L, mom_t_L = density, pressure, v_n, v_t, energy, mom_n, mom_t
        rho_R, p_R, vn_R, vt_R, E_R, mom_n_R, mom_t_R = \
            np.roll(rho_L, -1, axis=0), np.roll(p_L, -1, axis=0), np.roll(vn_L, -1, axis=0), \
            np.roll(vt_L, -1, axis=0), np.roll(E_L, -1, axis=0), np.roll(mom_n_L, -1, axis=0), \
            np.roll(mom_t_L, -1, axis=0)

        # Sound speed
        gamma = self.config.gamma
        cs_L = np.sqrt(gamma * p_L / rho_L, where=rho_L > 1e-12, out=np.zeros_like(rho_L))
        cs_R = np.sqrt(gamma * p_R / rho_R, where=rho_R > 1e-12, out=np.zeros_like(rho_R))

        # Wave speeds (vectorized)
        S_L = np.minimum(vn_L - cs_L, vn_R - cs_R)
        S_R = np.maximum(vn_L + cs_L, vn_R + cs_R)

        # Fluxes in L and R states (normal direction)
        F_dens_L, F_dens_R = rho_L * vn_L, rho_R * vn_R
        F_momn_L, F_momn_R = rho_L * vn_L**2 + p_L, rho_R * vn_R**2 + p_R
        F_momt_L, F_momt_R = rho_L * vn_L * vt_L, rho_R * vn_R * vt_R
        F_en_L, F_en_R = (E_L + p_L) * vn_L, (E_R + p_R) * vn_R

        # HLL flux calculation (vectorized)
        S_R_minus_S_L = S_R - S_L
        S_R_minus_S_L[np.abs(S_R_minus_S_L) < 1e-9] = 1e-9 # Avoid division by zero

        flux_density = np.where(S_L >= 0, F_dens_L, 
                    np.where(S_R <= 0, F_dens_R, 
                                (S_R*F_dens_L - S_L*F_dens_R + S_L*S_R*(rho_R-rho_L))/S_R_minus_S_L))
        flux_mom_n = np.where(S_L >= 0, F_momn_L,
                    np.where(S_R <= 0, F_momn_R,
                            (S_R*F_momn_L - S_L*F_momn_R + S_L*S_R*(mom_n_R-mom_n_L))/S_R_minus_S_L))
        flux_mom_t = np.where(S_L >= 0, F_momt_L,
                    np.where(S_R <= 0, F_momt_R,
                            (S_R*F_momt_L - S_L*F_momt_R + S_L*S_R*(mom_t_R-mom_t_L))/S_R_minus_S_L))
        flux_energy = np.where(S_L >= 0, F_en_L,
                    np.where(S_R <= 0, F_en_R,
                            (S_R*F_en_L - S_L*F_en_R + S_L*S_R*(E_R-E_L))/S_R_minus_S_L))

        # Transpose back if we were working on y-axis
        if axis == 1:
            return flux_density.T, flux_mom_t.T, flux_mom_n.T, flux_energy.T
        else:
            return flux_density, flux_mom_n, flux_mom_t, flux_energy

    def hydro_step(self, dt):
        """Performs one first-order finite-volume hydrodynamics step using operator splitting."""
        self.update_primitive_variables()

        factor = dt / self.config.cell_size
        
        # --- X-direction sweep ---
        flux_density_x, flux_mom_x_x, flux_mom_y_x, flux_energy_x = self.calculate_fluxes(axis=0)
        self.density    -= factor * (flux_density_x - np.roll(flux_density_x, 1, axis=0))
        self.momentum_x -= factor * (flux_mom_x_x - np.roll(flux_mom_x_x, 1, axis=0))
        self.momentum_y -= factor * (flux_mom_y_x - np.roll(flux_mom_y_x, 1, axis=0))
        self.energy     -= factor * (flux_energy_x - np.roll(flux_energy_x, 1, axis=0))
        self.update_primitive_variables()

        # --- Y-direction sweep ---
        flux_density_y, flux_mom_x_y, flux_mom_y_y, flux_energy_y = self.calculate_fluxes(axis=1)
        self.density    -= factor * (flux_density_y - np.roll(flux_density_y, 1, axis=1))
        self.momentum_x -= factor * (flux_mom_x_y - np.roll(flux_mom_x_y, 1, axis=1))
        self.momentum_y -= factor * (flux_mom_y_y - np.roll(flux_mom_y_y, 1, axis=1))
        self.energy     -= factor * (flux_energy_y - np.roll(flux_energy_y, 1, axis=1))
        self.update_primitive_variables()

# ------------------------
# Simulation Functions
# ------------------------

# minimal-image displacement helper (periodic)
def displacement(dx, config):
    L = config.domain_size
    dx = (dx + 0.5*L) % L - 0.5*L
    return dx

def cic_mass_assignment(particles, config):
    mesh_size = config.mesh_size
    cell_size = config.cell_size
    
    # Mass per cell
    rho = np.zeros((mesh_size, mesh_size))
    cic_data = []
    for p in particles:
        ix = int(p.x / cell_size)
        iy = int(p.y / cell_size)

        # Find the particle's fractional position within that cell (0 to 1)
        frac_x = (p.x / cell_size) - ix
        frac_y = (p.y / cell_size) - iy

        # Calculate weights for the 4 corners
        w1 = (1 - frac_x) * (1 - frac_y) # (ix, iy)
        w2 = frac_x * (1 - frac_y)       # (ix+1, iy)
        w3 = (1 - frac_x) * frac_y       # (ix, iy+1)
        w4 = frac_x * frac_y             # (ix+1, iy+1)

        cic_data.append((ix, iy, w1, w2, w3, w4))

        # Distribute mass to the 4 grid points (with periodic boundaries)
        rho[ix % mesh_size, iy % mesh_size] += p.mass * w1
        rho[(ix + 1) % mesh_size, iy % mesh_size] += p.mass * w2
        rho[ix % mesh_size, (iy + 1) % mesh_size] += p.mass * w3
        rho[(ix + 1) % mesh_size, (iy + 1) % mesh_size] += p.mass * w4
    
    rho /= (cell_size * cell_size) # Convert mass to density
    return rho, cic_data

def cic_force_interpolation(particles, ax_grid, ay_grid, cic_data, config):
    mesh_size = config.mesh_size
    mesh_forces = []
    for i, p in enumerate(particles):
        ix, iy, w1, w2, w3, w4 = cic_data[i]

        fx = (ax_grid[ix % mesh_size, iy % mesh_size] * w1 +
            ax_grid[(ix + 1) % mesh_size, iy % mesh_size] * w2 +
            ax_grid[ix % mesh_size, (iy + 1) % mesh_size] * w3 +
            ax_grid[(ix + 1) % mesh_size, (iy + 1) % mesh_size] * w4)
            
        fy = (ay_grid[ix % mesh_size, iy % mesh_size] * w1 +
            ay_grid[(ix + 1) % mesh_size, iy % mesh_size] * w2 +
            ay_grid[ix % mesh_size, (iy + 1) % mesh_size] * w3 +
            ay_grid[(ix + 1) % mesh_size, (iy + 1) % mesh_size] * w4)
        
        fx *= p.mass
        fy *= p.mass
        mesh_forces.append((fx, fy))
    
    return mesh_forces

def compute_gravitational_acceleration(particles, gas, config):
    mesh_size = config.mesh_size
    cell_size = config.cell_size
    
    # Cloud-In-Cell interpolation
    rho, cic_data = cic_mass_assignment(particles, config)
    if config.enable_hydro:
        rho += gas.density

    # FFT of density
    rho_k = np.fft.fftn(rho)

    # Poisson solver in k-space
    kx = np.fft.fftfreq(mesh_size, d=cell_size) * 2 * np.pi
    ky = np.fft.fftfreq(mesh_size, d=cell_size) * 2 * np.pi
    kx, ky = np.meshgrid(kx, ky, indexing='ij')
    k = np.sqrt(kx*kx + ky*ky)
    k[0,0] = 1.0  # avoid div by zero

    phi_k = -2 * np.pi * config.g_const * rho_k / k
    phi_k[0,0] = 0.0 # set the average value of the potential field to zero

    # Inverse FFT the POTENTIAL back to the real-space grid
    phi = np.fft.ifftn(phi_k).real

    # Gradient to get field (finite differences)
    ax_grid = (np.roll(phi, 1, axis=0) - np.roll(phi, -1, axis=0)) / (2*cell_size)
    ay_grid = (np.roll(phi, 1, axis=1) - np.roll(phi, -1, axis=1)) / (2*cell_size)

    return ax_grid, ay_grid, cic_data

def compute_PM_short_range_approx(dist_sq, p1_mass, p2_mass, dx, dy, config):
    softening_epsilon_squared = (0.5 * config.cell_size)**2 
    soft_dist_sq = dist_sq + softening_epsilon_squared
    soft_dist = math.sqrt(soft_dist_sq)
    f_pm_short = config.g_const * p1_mass * p2_mass / soft_dist_sq
    f_pm_short_x = f_pm_short * dx / soft_dist
    f_pm_short_y = f_pm_short * dy / soft_dist
    return f_pm_short_x, f_pm_short_y

def compute_switching_parameter(dist_sq, config):
    if dist_sq < config.r_switch_start_sq:
        S = 1.0
    else:
        # In the transition zone, smoothly fade the correction to zero
        dist = math.sqrt(dist_sq)
        x = (dist - config.r_switch_start) / config.cutoff_transition_width
        S = 2*x**3 - 3*x**2 + 1
    return S

def compute_PP_forces(particles, config):
    cell_size = config.cell_size
    mesh_size = config.mesh_size

    # Bin particles into a grid (cell list)
    cells = {} # Key: (ix, iy), Value: list of particle indices
    particle_cells = [] # Store (ix, iy) for each particle i
    for i, p in enumerate(particles):
        ix = int(p.x / cell_size)
        iy = int(p.y / cell_size)
        particle_cells.append((ix, iy))
        key = (ix, iy)
        if key not in cells:
            cells[key] = []
        cells[key].append(i)

    forces = [[0.0, 0.0] for _ in particles]
    
    # Base search radius in cells on cutoff
    search_radius_cells = int(math.ceil(config.cutoff_radius / cell_size))

    # Iterate over particles and their neighbors
    for i, p1 in enumerate(particles):
        ix1, iy1 = particle_cells[i]
        
        # Check blocks of cells around the particle
        for dx_cell in range(-search_radius_cells, search_radius_cells + 1):
            for dy_cell in range(-search_radius_cells, search_radius_cells + 1):
                
                neighbor_ix = (ix1 + dx_cell) % mesh_size
                neighbor_iy = (iy1 + dy_cell) % mesh_size
                cell_key = (neighbor_ix, neighbor_iy)
                
                if cell_key not in cells:
                    continue # Empty neighbor cells

                # Compute interactions with particles in the neighbor cell
                for j in cells[cell_key]:
                    # Avoid double counting (j > i) and self-interaction (j != i)
                    if i >= j:
                        continue
                            
                    p2 = particles[j]
                    
                    dx = displacement(p2.x - p1.x, config)
                    dy = displacement(p2.y - p1.y, config)
                    dist_sq = dx*dx + dy*dy

                    if dist_sq > config.cutoff_radius_squared:
                        continue
                    
                    # Subtractive scheme
                    fx_sub, fy_sub = compute_PM_short_range_approx(dist_sq, p1.mass, p2.mass, dx, dy, config)
                    S = compute_switching_parameter(dist_sq, config)
                    dist_sq_soft = dist_sq + config.softening_squared
                    dist = math.sqrt(dist_sq_soft)
                    f = config.g_const * p1.mass * p2.mass / dist_sq_soft
                    fx = S * (f * dx / dist - fx_sub)
                    fy = S * (f * dy / dist - fy_sub)
            
                    # Apply forces to both particles
                    forces[i][0] += fx
                    forces[i][1] += fy

                    forces[j][0] -= fx
                    forces[j][1] -= fy
            
    # Convert back to a list of tuples for compatibility
    return [tuple(f) for f in forces]

def update_cosmology(total_time, config):
    if config.expanding_universe:
        expansion_time = config.expansion_start_t + total_time
        scale_factor = expansion_time**(2.0/3.0)
        hubble_param = (2.0/3.0) / expansion_time
    else:
        scale_factor = 1.0
        hubble_param = 0.0
    return scale_factor, hubble_param

def KDK_step(total_time, particles, gas, config):
    scale_factor, hubble_param = update_cosmology(total_time, config)
    dt = config.dt
    
    # 1. KICK (First half-step for velocity)
    for p in particles:
        total_ax = (p.ax / scale_factor**3) - (2 * hubble_param * p.vx)
        total_ay = (p.ay / scale_factor**3) - (2 * hubble_param * p.vy)
        p.vx += total_ax * dt / 2.0
        p.vy += total_ay * dt / 2.0

    # Apply Kick to Gas Grid
    if config.enable_hydro:
        gas.update_primitive_variables()
        gas_accel_x = (gas.accel_x / scale_factor**3) - (2 * hubble_param * gas.velocity_x)
        gas_accel_y = (gas.accel_y / scale_factor**3) - (2 * hubble_param * gas.velocity_y)
        
        # Calculate power term P = F · v = (ρa) · v = ρ(a · v)
        power_density = gas.density * (gas_accel_x * gas.velocity_x + gas_accel_y * gas.velocity_y)
        
        gas.momentum_x += gas.density * gas_accel_x * dt / 2.0
        gas.momentum_y += gas.density * gas_accel_y * dt / 2.0
        gas.energy += power_density * dt / 2.0
    
    # 2. DRIFT (Full step for position)
    for p in particles:
        p.x = ( p.x + p.vx * dt ) % config.domain_size
        p.y = ( p.y + p.vy * dt ) % config.domain_size

    # 3. UPDATE COSMOLOGY
    scale_factor, hubble_param = update_cosmology(total_time + dt, config)

    # 4. COMPUTE FORCES at new positions
    pp_forces = compute_PP_forces(particles, config)
    ax_grid, ay_grid, cic_data = compute_gravitational_acceleration(particles, gas, config)
    pm_forces = cic_force_interpolation(particles, ax_grid, ay_grid, cic_data, config)
    forces = [(pp[0]+pm[0], pp[1]+pm[1]) for pp,pm in zip(pp_forces, pm_forces)]
    
    if config.enable_hydro:
        gas.hydro_step(dt)

    # 5. KICK (Second half-step for velocity)
    for i, p in enumerate(particles):
        fx, fy = forces[i]
        p.ax = fx / p.mass
        p.ay = fy / p.mass
        total_ax = (p.ax / scale_factor**3) - (2 * hubble_param * p.vx)
        total_ay = (p.ay / scale_factor**3) - (2 * hubble_param * p.vy)
        p.vx += total_ax * dt / 2.0
        p.vy += total_ay * dt / 2.0

    # Apply Kick to Gas Grid
    if config.enable_hydro:
        gas.update_primitive_variables()
        gas.accel_x = ax_grid
        gas.accel_y = ay_grid
        gas_accel_x = (gas.accel_x / scale_factor**3) - (2 * hubble_param * gas.velocity_x)
        gas_accel_y = (gas.accel_y / scale_factor**3) - (2 * hubble_param * gas.velocity_y)
        
        power_density = gas.density * (gas_accel_x * gas.velocity_x + gas_accel_y * gas.velocity_y)

        gas.momentum_x += gas.density * gas_accel_x * dt / 2.0
        gas.momentum_y += gas.density * gas_accel_y * dt / 2.0
        gas.energy += power_density * dt / 2.0

# Create initial conditions using the simplified Zel'dovich Approximation
def create_zeldovich_ics(config):
    n_per_side = config.num_dm_particles_side
    ic_mesh_size = n_per_side
    cell_size = config.domain_size / ic_mesh_size

    # Build k-grid
    kx = np.fft.fftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi
    ky = np.fft.fftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi
    kx_grid, ky_grid = np.meshgrid(kx, ky, indexing='ij')
    k2 = kx_grid**2 + ky_grid**2
    k2[0, 0] = 1.0  # avoid division by zero

    # Generate random density field δ(k) with desired P(k) ∝ k^n
    real_space_random_field = np.random.randn(ic_mesh_size, ic_mesh_size)
    random_k = np.fft.fft2(real_space_random_field)
    
    # Power spectrum P(k) ~ k^n
    # Amplitude delta(k) ~ sqrt(P(k)) ~ k^(n/2)
    power_spectrum_index = config.initial_power_spectrum_index
    delta_k = random_k * (k2 ** (power_spectrum_index / 4.0)) # n/2 / 2 = n/4
    delta_k[0, 0] = 0.0 # No DC offset

    # Compute potential Φ(k) = -δ(k)/k² (in 2D, Poisson is ∇²Φ = δ, so k²Φ = δ)
    # G factor is absorbed into normalization
    phi_k = -delta_k / k2
    phi_k[0, 0] = 0.0

    # Displacement field s(k) = -i k Φ(k)
    disp_x_k = -1j * kx_grid * phi_k
    disp_y_k = -1j * ky_grid * phi_k

    disp_x = np.fft.ifft2(disp_x_k).real
    disp_y = np.fft.ifft2(disp_y_k).real

    # Normalize & apply amplitude
    disp_x /= np.std(disp_x)
    disp_y /= np.std(disp_y)
    disp_x *= config.initial_scale_factor * cell_size
    disp_y *= config.initial_scale_factor * cell_size

    # Create displaced particles
    particles = []
    spacing = config.domain_size / n_per_side
    for i in range(n_per_side):
        for j in range(n_per_side):
            qx = (i + 0.5) * spacing
            qy = (j + 0.5) * spacing
            x = (qx + disp_x[i, j]) % config.domain_size
            y = (qy + disp_y[i, j]) % config.domain_size
            vx = config.initial_hubble_param * disp_x[i, j]
            vy = config.initial_hubble_param * disp_y[i, j]
            particles.append(Particle(x, y, config.dm_particle_mass, vx, vy))

    return particles

def calculate_total_energy(particles, scale_factor, config):
    kinetic_energy = 0.0
    potential_energy = 0.0
    g_const = config.g_const
    softening_squared = config.softening_squared

    for p in particles:
        vel_sq = (scale_factor*p.vx)**2 + (scale_factor*p.vy)**2
        kinetic_energy += 0.5 * p.mass * vel_sq

    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            p1 = particles[i]
            p2 = particles[j]

            dx = displacement(p2.x - p1.x, config)
            dy = displacement(p2.y - p1.y, config)
            
            dist_sq = (scale_factor**2) * (dx * dx + dy * dy + softening_squared)
            dist = math.sqrt(dist_sq)
            if dist > 0:
                potential_energy -= g_const * p1.mass * p2.mass / dist

    return kinetic_energy + potential_energy

# ------------------------
# Visualization helpers
# ------------------------
def init_window(render_size, caption):
    pygame.init()
    screen = pygame.display.set_mode((render_size, render_size))
    pygame.display.set_caption(caption)
    return screen

def process_events():
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            return False
    return True

def render(screen, particles, gas, config):
    render_size = config.render_size
    render_scale = config.render_scale
    mesh_size = config.mesh_size
    
    surface_array = np.zeros((render_size, render_size, 3), dtype=np.uint8)

    if config.enable_hydro:
        # Render Gas Pressure
        log_field = np.log10(gas.pressure + 1e-12) 
        min_log, max_log = np.min(log_field), np.max(log_field)
        
        if max_log > min_log:
            norm_field = (log_field - min_log) / (max_log - min_log)
        else:
            norm_field = np.zeros_like(log_field)

        colormap = matplotlib.cm.plasma
        colored_field = colormap(norm_field)
        
        # Manually apply brightness scaling
        colored_field[:, :, 0:3] *= norm_field[..., np.newaxis]
        
        zoom_factor = render_size / mesh_size
        zoomed_field = scipy.ndimage.zoom(colored_field[:, :, 0:3], 
                                        (zoom_factor, zoom_factor, 1), 
                                        order=1) # order=1 is bilinear

        surface_array = (zoomed_field * 255).astype(np.uint8)

    # Draw DM particles on top
    coords = np.array([(int(p.x*render_scale), int(p.y*render_scale)) for p in particles])
    
    # Clip coordinates to be within the render size
    coords = np.clip(coords, 0, render_size - 1)
    
    surface_array[coords[:, 0], coords[:, 1]] = [255, 255, 255] # White
    
    pygame.surfarray.blit_array(screen, surface_array)
    pygame.display.flip()

def save_snapshot(h5_file, snapshot_name, particles, gas, sim_time, scale_factor, config):
    """Saves the current simulation state to a group in the HDF5 file."""
    try:
        grp = h5_file.create_group(snapshot_name)
        
        # Store Metadata
        grp.attrs['time'] = sim_time
        grp.attrs['scale_factor'] = scale_factor
        
        # Store Particle Data
        num_p = len(particles)
        positions = np.zeros((num_p, 2))
        velocities = np.zeros((num_p, 2))
        masses = np.zeros(num_p)
        
        for i, p in enumerate(particles):
            positions[i] = [p.x, p.y]
            velocities[i] = [p.vx, p.vy]
            masses[i] = p.mass
            
        p_grp = grp.create_group("particles")
        p_grp.create_dataset("positions", data=positions, compression="gzip")
        p_grp.create_dataset("velocities", data=velocities, compression="gzip")
        p_grp.create_dataset("masses", data=masses, compression="gzip")
        
        # Store Gas Data
        if config.enable_hydro:
            g_grp = grp.create_group("gas")
            g_grp.create_dataset("density", data=gas.density, compression="gzip")
            g_grp.create_dataset("momentum_x", data=gas.momentum_x, compression="gzip")
            g_grp.create_dataset("momentum_y", data=gas.momentum_y, compression="gzip")
            g_grp.create_dataset("energy", data=gas.energy, compression="gzip")
            g_grp.create_dataset("pressure", data=gas.pressure, compression="gzip") 
            g_grp.create_dataset("velocity_x", data=gas.velocity_x, compression="gzip")
            g_grp.create_dataset("velocity_y", data=gas.velocity_y, compression="gzip")
            
    except Exception as e:
        print(f"Error saving snapshot {snapshot_name}: {e}")

# ------------------------
# Main entry point
# ------------------------
def main():
    # --- Load Configuration ---
    if len(sys.argv) > 1:
        config_filename = sys.argv[1]
    else:
        config_filename = "simulation.ini"
        
    try:
        config = Config(config_filename)
    except Exception as e:
        print(f"Error loading config file: {e}")
        return

    print(f"Loaded parameters from: {config_filename}")
    # --- End Configuration ---

    screen = init_window(config.render_size, 'N-Body + Hydrodynamics Simulation')

    particles = create_zeldovich_ics(config)
    
    # Uniform initialization for gas
    gas = GasGrid(config)
    gas.density.fill(config.gas_total_mass / (config.domain_size**2))
    initial_internal_energy = 1e-6
    gas.energy = initial_internal_energy * gas.density
    gas.update_primitive_variables()
    
    if config.enable_hydro:
        # Calculate initial gravitational field for the gas
        ax_grid, ay_grid, _ = compute_gravitational_acceleration(particles, gas, config)
        gas.accel_x = ax_grid
        gas.accel_y = ay_grid

    initial_energy = calculate_total_energy(particles, config.initial_scale_factor, config)

    cycle_count = 0
    simulation_start_time = time.time()
    total_simulation_time = 0.0

    # Create HDF5 file
    run_timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    snapshot_filename = f"sim_{run_timestamp}.hdf5"
    print(f"Saving snapshots to: {snapshot_filename}")

    h5_file = h5py.File(snapshot_filename, 'w')
    param_grp = h5_file.create_group("parameters")
    
    # Save all base parameters to the HDF5 file
    for key, val in config.get_all_params().items():
        param_grp.attrs[key] = val
        
    snapshot_count = 0

    running = True
    while running:
        running = process_events()

        total_simulation_time = config.dt * cycle_count
        KDK_step(total_simulation_time, particles, gas, config)
        render(screen, particles, gas, config)

        cycle_count += 1

        if config.save_snapshot_every_cycles > 0 and cycle_count % config.save_snapshot_every_cycles == 0:
            snapshot_name = f"snapshot_{snapshot_count:04d}"
            scale_factor, _ = update_cosmology(total_simulation_time, config)
            save_snapshot(h5_file, snapshot_name, particles, gas, total_simulation_time, scale_factor, config)
            print(f"Saved snapshot: {snapshot_name}")
            snapshot_count += 1

        if config.debug_info_every_cycles > 0 and cycle_count % config.debug_info_every_cycles == 0:
            total_elapsed_time = time.time() - simulation_start_time
            average_time = total_elapsed_time / cycle_count if cycle_count > 0 else 0
            scale_factor, _ = update_cosmology(total_simulation_time, config)
            current_energy = calculate_total_energy(particles, scale_factor, config)
            energy_error = abs(current_energy - initial_energy) / abs(initial_energy) if abs(initial_energy) > 1e-12 else 0.0
            print(f"Time:{total_elapsed_time:.0f}s Cycles:{cycle_count} Avg. dt:{average_time:.3f}s Scale:{scale_factor:.2f} Energy Error:{energy_error:.1%}")

    h5_file.close()
    pygame.quit()
    print("Simulation finished. HDF5 file closed.")

if __name__ == "__main__":
    main()
