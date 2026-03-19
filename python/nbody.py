import os
import math
import time
import datetime
import numpy as np
import h5py
import configparser
import sys
import csv

# VisPy imports
from vispy import app, scene
from vispy.color import get_colormap

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
        self.cfl_safety_factor = self._get(config, 'time', 'cfl_safety_factor', 'float')
        self.gravity_dt_factor = self._get(config, 'time', 'gravity_dt_factor', 'float')
        self.use_adaptive_dt = self._get(config, 'time', 'use_adaptive_dt', 'bool')

        # [output]
        self.debug_info_every_cycles = self._get(config, 'output', 'debug_info_every_cycles', 'int')
        self.save_snapshot_every_cycles = self._get(config, 'output', 'save_snapshot_every_cycles', 'int')
        self.seed = self._get(config, 'output', 'seed', 'int')

        # [visualization]
        self.render_size = self._get(config, 'visualization', 'render_size', 'int')

        # --- Calculate derived parameters ---
        self.cell_size = self.domain_size / self.mesh_size
        self.cell_volume = self.cell_size**3
        self.omega_dm = 1.0 - self.omega_baryon
        
        self.cutoff_radius = self.cutoff_radius_cells * self.cell_size
        self.cutoff_radius_squared = self.cutoff_radius**2
        self.cutoff_transition_width = self.cutoff_transition_width_factor * self.cutoff_radius
        self.r_switch_start = self.cutoff_radius - self.cutoff_transition_width
        self.r_switch_start_sq = self.r_switch_start**2
        
        self.num_dm_particles = self.num_dm_particles_side**3
        self.dm_particle_mass = self.omega_dm / self.num_dm_particles
        self.gas_total_mass = self.omega_baryon
        
        self.mean_interparticle_spacing = self.domain_size / (self.num_dm_particles**(1/3))
        self.softening_squared = (self.mean_interparticle_spacing / 50.0)**2
        
        self.initial_hubble_param = (2.0/3.0) / self.expansion_start_t if self.expanding_universe else 0.0
        self.initial_scale_factor = self.expansion_start_t**(2.0/3.0) if self.expanding_universe else 0.0
        
        self.dynamical_time = 1.0 / math.sqrt(4.0 * math.pi * self.g_const) # Assuming rho_crit=1
        self.fixed_dt = self.dt_factor * self.dynamical_time
        
        self.render_scale = self.render_size / self.domain_size

        # Set seed
        np.random.seed(self.seed)

    def _get(self, config, section, option, dtype='str', default=None):
        """Helper to get config value and store it for HDF5 saving."""
        try:
            if dtype == 'float':
                val = config.getfloat(section, option)
            elif dtype == 'int':
                val = config.getint(section, option)
            elif dtype == 'bool':
                val = config.getboolean(section, option)
            else:
                val = config.get(section, option)
        except (configparser.NoOptionError, configparser.NoSectionError):
            if default is not None:
                val = default
            else:
                raise
        
        self.params[f"{section}_{option}"] = val
        return val

    def get_all_params(self):
        """Get all *base* parameters for saving to HDF5."""
        return self.params

# ------------------------
# Diagnostics & Logging
# ------------------------
class Diagnostics:
    """A data bucket to hold all diagnostics for a given step."""
    def __init__(self):
        # Simulation state
        self.cycle = 0
        self.sim_time = 0.0
        self.scale_factor = 1.0
        
        # Conservation
        self.total_mass_gas = 0.0
        self.total_mass_dm = 0.0
        self.total_momentum_gas = (0.0, 0.0, 0.0)
        self.total_momentum_dm = (0.0, 0.0, 0.0)
        self.ke_gas = 0.0
        self.ke_dm = 0.0
        self.pe_dm = 0.0
        self.ie_gas = 0.0
        
        # Stability
        self.dt_cfl = 0.0
        self.dt_gravity = 0.0
        self.dt_final = 0.0
        self.max_gas_density = 0.0
        self.max_gas_pressure = 0.0
        self.max_gas_velocity = 0.0
        
        # Performance
        self.wall_time_total = 0.0
        self.wall_time_pm = 0.0
        self.wall_time_pp = 0.0
        self.wall_time_hydro = 0.0
        self.wall_time_io = 0.0
        
    @property
    def total_mass(self):
        return self.total_mass_gas + self.total_mass_dm

    @property
    def total_momentum(self):
        return (self.total_momentum_gas[0] + self.total_momentum_dm[0],
                self.total_momentum_gas[1] + self.total_momentum_dm[1],
                self.total_momentum_gas[2] + self.total_momentum_dm[2])
    
    @property
    def total_energy(self):
        # Note: PE for gas and gas-DM interaction is not included
        return self.ke_gas + self.ke_dm + self.pe_dm + self.ie_gas

class Logger:
    """Handles writing diagnostics to stdout and a CSV log file."""
    def __init__(self, log_filename="diagnostics.csv"):
        self.log_filename = log_filename
        self.log_file = open(log_filename, 'w', newline='')
        self.writer = None
        self.start_time = time.time()

    def write_header(self, diag):
        """Writes the header row to the CSV file."""
        # Manually define header order
        headers = [
            'cycle', 'sim_time', 'scale_factor',
            'total_mass_gas', 'total_mass_dm',
            'total_momentum_gas_x', 'total_momentum_gas_y', 'total_momentum_gas_z',
            'total_momentum_dm_x', 'total_momentum_dm_y', 'total_momentum_dm_z',
            'ke_gas', 'ke_dm', 'pe_dm', 'ie_gas',
            'dt_cfl', 'dt_gravity', 'dt_final',
            'max_gas_density', 'max_gas_pressure', 'max_gas_velocity',
            'wall_time_total', 'wall_time_pm', 'wall_time_pp', 'wall_time_hydro', 'wall_time_io'
        ]
        
        # Create a flat dictionary for writing
        flat_diag = self._flatten_diag(diag)
        
        # Filter headers to only include keys present in the flat_diag
        self.fieldnames = [h for h in headers if h in flat_diag]
        
        self.writer = csv.DictWriter(self.log_file, fieldnames=self.fieldnames)
        self.writer.writeheader()
        self.log_file.flush()

    def _flatten_diag(self, diag):
        """Flattens the diagnostics object for CSV writing."""
        flat = diag.__dict__.copy()
        
        # Flatten tuples
        mom_g = flat.pop('total_momentum_gas')
        flat['total_momentum_gas_x'] = mom_g[0]
        flat['total_momentum_gas_y'] = mom_g[1]
        flat['total_momentum_gas_z'] = mom_g[2]
        
        mom_dm = flat.pop('total_momentum_dm')
        flat['total_momentum_dm_x'] = mom_dm[0]
        flat['total_momentum_dm_y'] = mom_dm[1]
        flat['total_momentum_dm_z'] = mom_dm[2]
        
        return flat

    def log(self, diag):
        """Logs a diagnostic step to both the CSV file and stdout."""
        if self.writer is None:
            self.write_header(diag)
            
        # Write to CSV
        flat_diag = self._flatten_diag(diag)
        # Filter out keys not in header
        csv_row = {k: v for k, v in flat_diag.items() if k in self.fieldnames}
        self.writer.writerow(csv_row)
        self.log_file.flush()
        
        # Write to stdout
        wall_time_str = f"Wall: {time.time() - self.start_time:.1f}s"
        sim_time_str = f"SimTime: {diag.sim_time:.3f}"
        scale_str = f"a: {diag.scale_factor:.3f}"
        
        print(f"[Cycle {diag.cycle}] {wall_time_str} | {sim_time_str} | {scale_str}")
        
        mass_err = (diag.total_mass - 1.0) # Assuming total mass is 1.0
        print(f"  [Physics]")
        print(f"    - Mass (P/G/T):   {diag.total_mass_dm:.4f} | {diag.total_mass_gas:.4f} | {diag.total_mass:.4f} (Err: {mass_err:+.1e})")
        print(f"    - Momentum (P/G): ({diag.total_momentum_dm[0]:.1e}, {diag.total_momentum_dm[1]:.1e}, {diag.total_momentum_dm[2]:.1e}) | ({diag.total_momentum_gas[0]:.1e}, {diag.total_momentum_gas[1]:.1e}, {diag.total_momentum_gas[2]:.1e})")
        print(f"    - Energy (KE/PE/IE): {diag.ke_dm + diag.ke_gas:.3e} | {diag.pe_dm:.3e} | {diag.ie_gas:.3e} (Total: {diag.total_energy:.3e})")
        
        print(f"  [Stability]")
        print(f"    - Timestep (CFL): {diag.dt_cfl:.2e} | (Grav): {diag.dt_gravity:.2e} | (Final): {diag.dt_final:.2e}")
        print(f"    - Max(rho): {diag.max_gas_density:.2e} | Max(press): {diag.max_gas_pressure:.2e} | Max(vel): {diag.max_gas_velocity:.2e}")

        print(f"  [Performance (ms)]")
        print(f"    - PM: {diag.wall_time_pm*1000:.1f} | PP: {diag.wall_time_pp*1000:.1f} | Hydro: {diag.wall_time_hydro*1000:.1f} | I/O: {diag.wall_time_io*1000:.1f} | Total: {diag.wall_time_total*1000:.1f}")
        print("-" * 70)

    def close(self):
        self.log_file.close()

# ------------------------
# Simulation Classes
# ------------------------

class Particle:
    def __init__(self, x, y, z, mass, vx=0, vy=0, vz=0):
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.ax = 0
        self.ay = 0
        self.az = 0

class GasGrid:
    """Represents the Eulerian grid for the baryonic gas."""
    def __init__(self, config):
        self.config = config
        mesh_size = config.mesh_size
        grid_shape = (mesh_size, mesh_size, mesh_size)
        
        # Conservative variables
        self.density = np.zeros(grid_shape)
        self.momentum_x = np.zeros(grid_shape)
        self.momentum_y = np.zeros(grid_shape)
        self.momentum_z = np.zeros(grid_shape)
        self.energy = np.zeros(grid_shape)
        
        # Primitive variables (derived)
        self.pressure = np.zeros(grid_shape)
        self.velocity_x = np.zeros(grid_shape)
        self.velocity_y = np.zeros(grid_shape)
        self.velocity_z = np.zeros(grid_shape)

        # Acceleration field
        self.accel_x = np.zeros(grid_shape)
        self.accel_y = np.zeros(grid_shape)
        self.accel_z = np.zeros(grid_shape)

    def update_primitive_variables(self):
        """Calculates pressure and velocity from the conservative variables."""
        non_zero_density = self.density > 1e-12
        
        self.velocity_x.fill(0.0)
        self.velocity_y.fill(0.0)
        self.velocity_z.fill(0.0)
        
        self.velocity_x[non_zero_density] = self.momentum_x[non_zero_density] / self.density[non_zero_density]
        self.velocity_y[non_zero_density] = self.momentum_y[non_zero_density] / self.density[non_zero_density]
        self.velocity_z[non_zero_density] = self.momentum_z[non_zero_density] / self.density[non_zero_density]

        kinetic_energy_density = np.zeros_like(self.density)
        kinetic_energy_density[non_zero_density] = 0.5 * (
            self.momentum_x[non_zero_density]**2 + 
            self.momentum_y[non_zero_density]**2 +
            self.momentum_z[non_zero_density]**2
        ) / self.density[non_zero_density]
        
        internal_energy_density = self.energy - kinetic_energy_density
        
        self.pressure = (self.config.gamma - 1.0) * internal_energy_density
        self.pressure[self.pressure < 1e-12] = 1e-12 # Pressure floor

    def calculate_fluxes(self, axis):
        """
        Internal helper to calculate HLL fluxes along a given axis (0=x, 1=y, 2=z).
        Handles 3D by considering 1 normal (n) and 2 tangential (t1, t2) components.
        """
        
        # Permute axes to make the calculation axis=0
        if axis == 0: # x-sweep
            density, mom_n, mom_t1, mom_t2, energy, v_n, v_t1, v_t2, pressure = \
                self.density, self.momentum_x, self.momentum_y, self.momentum_z, self.energy, \
                self.velocity_x, self.velocity_y, self.velocity_z, self.pressure
        elif axis == 1: # y-sweep
            density = np.transpose(self.density, (1, 0, 2))
            mom_n = np.transpose(self.momentum_y, (1, 0, 2))
            mom_t1 = np.transpose(self.momentum_x, (1, 0, 2))
            mom_t2 = np.transpose(self.momentum_z, (1, 0, 2))
            energy = np.transpose(self.energy, (1, 0, 2))
            v_n = np.transpose(self.velocity_y, (1, 0, 2))
            v_t1 = np.transpose(self.velocity_x, (1, 0, 2))
            v_t2 = np.transpose(self.velocity_z, (1, 0, 2))
            pressure = np.transpose(self.pressure, (1, 0, 2))
        elif axis == 2: # z-sweep
            density = np.transpose(self.density, (2, 0, 1))
            mom_n = np.transpose(self.momentum_z, (2, 0, 1))
            mom_t1 = np.transpose(self.momentum_x, (2, 0, 1))
            mom_t2 = np.transpose(self.momentum_y, (2, 0, 1))
            energy = np.transpose(self.energy, (2, 0, 1))
            v_n = np.transpose(self.velocity_z, (2, 0, 1))
            v_t1 = np.transpose(self.velocity_x, (2, 0, 1))
            v_t2 = np.transpose(self.velocity_y, (2, 0, 1))
            pressure = np.transpose(self.pressure, (2, 0, 1))


        # Left and Right states for all cells at once using roll (on axis 0)
        rho_L, p_L, vn_L, vt1_L, vt2_L, E_L, mom_n_L, mom_t1_L, mom_t2_L = \
            density, pressure, v_n, v_t1, v_t2, energy, mom_n, mom_t1, mom_t2
        
        rho_R, p_R, vn_R, vt1_R, vt2_R, E_R, mom_n_R, mom_t1_R, mom_t2_R = \
            np.roll(rho_L, -1, axis=0), np.roll(p_L, -1, axis=0), np.roll(vn_L, -1, axis=0), \
            np.roll(vt1_L, -1, axis=0), np.roll(vt2_L, -1, axis=0), np.roll(E_L, -1, axis=0), \
            np.roll(mom_n_L, -1, axis=0), np.roll(mom_t1_L, -1, axis=0), np.roll(mom_t2_L, -1, axis=0)

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
        F_momt1_L, F_momt1_R = rho_L * vn_L * vt1_L, rho_R * vn_R * vt1_R
        F_momt2_L, F_momt2_R = rho_L * vn_L * vt2_L, rho_R * vn_R * vt2_R
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
        flux_mom_t1 = np.where(S_L >= 0, F_momt1_L,
                    np.where(S_R <= 0, F_momt1_R,
                            (S_R*F_momt1_L - S_L*F_momt1_R + S_L*S_R*(mom_t1_R-mom_t1_L))/S_R_minus_S_L))
        flux_mom_t2 = np.where(S_L >= 0, F_momt2_L,
                    np.where(S_R <= 0, F_momt2_R,
                            (S_R*F_momt2_L - S_L*F_momt2_R + S_L*S_R*(mom_t2_R-mom_t2_L))/S_R_minus_S_L))
        flux_energy = np.where(S_L >= 0, F_en_L,
                    np.where(S_R <= 0, F_en_R,
                            (S_R*F_en_L - S_L*F_en_R + S_L*S_R*(E_R-E_L))/S_R_minus_S_L))

        # Transpose back if necessary
        if axis == 0: # x-sweep
            return flux_density, flux_mom_n, flux_mom_t1, flux_mom_t2, flux_energy
        elif axis == 1: # y-sweep
            return flux_density.transpose(1, 0, 2), flux_mom_t1.transpose(1, 0, 2), \
                   flux_mom_n.transpose(1, 0, 2), flux_mom_t2.transpose(1, 0, 2), flux_energy.transpose(1, 0, 2)
        elif axis == 2: # z-sweep
            return flux_density.transpose(1, 2, 0), flux_mom_t1.transpose(1, 2, 0), \
                   flux_mom_t2.transpose(1, 2, 0), flux_mom_n.transpose(1, 2, 0), flux_energy.transpose(1, 2, 0)


    def hydro_step(self, dt):
        """
        Performs a 1st-order Godunov hydrodynamics step using an HLL Riemann solver and directional operator splitting.
        """
        self.update_primitive_variables()
        factor = dt / self.config.cell_size
        
        # --- X-direction sweep ---
        flux_density, flux_mom_x, flux_mom_y, flux_mom_z, flux_energy = self.calculate_fluxes(axis=0)
        self.density    -= factor * (flux_density - np.roll(flux_density, 1, axis=0))
        self.momentum_x -= factor * (flux_mom_x - np.roll(flux_mom_x, 1, axis=0))
        self.momentum_y -= factor * (flux_mom_y - np.roll(flux_mom_y, 1, axis=0))
        self.momentum_z -= factor * (flux_mom_z - np.roll(flux_mom_z, 1, axis=0))
        self.energy     -= factor * (flux_energy - np.roll(flux_energy, 1, axis=0))
        self.update_primitive_variables()

        # --- Y-direction sweep ---
        flux_density, flux_mom_x, flux_mom_y, flux_mom_z, flux_energy = self.calculate_fluxes(axis=1)
        self.density    -= factor * (flux_density - np.roll(flux_density, 1, axis=1))
        self.momentum_x -= factor * (flux_mom_x - np.roll(flux_mom_x, 1, axis=1))
        self.momentum_y -= factor * (flux_mom_y - np.roll(flux_mom_y, 1, axis=1))
        self.momentum_z -= factor * (flux_mom_z - np.roll(flux_mom_z, 1, axis=1))
        self.energy     -= factor * (flux_energy - np.roll(flux_energy, 1, axis=1))
        self.update_primitive_variables()
        
        # --- Z-direction sweep ---
        flux_density, flux_mom_x, flux_mom_y, flux_mom_z, flux_energy = self.calculate_fluxes(axis=2)
        self.density    -= factor * (flux_density - np.roll(flux_density, 1, axis=2))
        self.momentum_x -= factor * (flux_mom_x - np.roll(flux_mom_x, 1, axis=2))
        self.momentum_y -= factor * (flux_mom_y - np.roll(flux_mom_y, 1, axis=2))
        self.momentum_z -= factor * (flux_mom_z - np.roll(flux_mom_z, 1, axis=2))
        self.energy     -= factor * (flux_energy - np.roll(flux_energy, 1, axis=2))
        self.update_primitive_variables()

# ------------------------
# Simulation Functions
# ------------------------

# minimal-image displacement helper (periodic)
def displacement(dx, config):
    L = config.domain_size
    dx = (dx + 0.5*L) % L - 0.5*L
    return dx

def particle_binning_and_mass_assignment(particles, config):
    """
    Performs CIC mass assignment and particle binning in one pass.
    """
    mesh_size = config.mesh_size
    cell_size = config.cell_size
    
    # PM Data
    dm_rho = np.zeros((mesh_size, mesh_size, mesh_size))
    cic_data = [] # List of (ix, iy, iz, w1..w8)
    
    # PP Data
    cells = {} # Key: (ix, iy, iz), Value: list of particle indices
    particle_cells = [] # List of (ix, iy, iz) for each particle i

    for i, p in enumerate(particles):
        ix = int(p.x / cell_size)
        iy = int(p.y / cell_size)
        iz = int(p.z / cell_size)

        frac_x = (p.x / cell_size) - ix
        frac_y = (p.y / cell_size) - iy
        frac_z = (p.z / cell_size) - iz
        
        # CIC weights (8 corners)
        w1 = (1-frac_x) * (1-frac_y) * (1-frac_z)
        w2 = (  frac_x) * (1-frac_y) * (1-frac_z)
        w3 = (1-frac_x) * (  frac_y) * (1-frac_z)
        w4 = (1-frac_x) * (1-frac_y) * (  frac_z)
        w5 = (  frac_x) * (  frac_y) * (1-frac_z)
        w6 = (  frac_x) * (1-frac_y) * (  frac_z)
        w7 = (1-frac_x) * (  frac_y) * (  frac_z)
        w8 = (  frac_x) * (  frac_y) * (  frac_z)

        # 1. For PM (CIC weights)
        cic_data.append((ix, iy, iz, w1, w2, w3, w4, w5, w6, w7, w8))

        # 2. For PP (Cell list)
        particle_cells.append((ix, iy, iz))
        key = (ix % mesh_size, iy % mesh_size, iz % mesh_size)
        if key not in cells:
            cells[key] = []
        cells[key].append(i)

        # 3. For PM (Mass assignment)
        dm_rho[ix % mesh_size,     iy % mesh_size,     iz % mesh_size    ] += p.mass * w1
        dm_rho[(ix+1) % mesh_size, iy % mesh_size,     iz % mesh_size    ] += p.mass * w2
        dm_rho[ix % mesh_size,     (iy+1) % mesh_size, iz % mesh_size    ] += p.mass * w3
        dm_rho[ix % mesh_size,     iy % mesh_size,     (iz+1) % mesh_size] += p.mass * w4
        dm_rho[(ix+1) % mesh_size, (iy+1) % mesh_size, iz % mesh_size    ] += p.mass * w5
        dm_rho[(ix+1) % mesh_size, iy % mesh_size,     (iz+1) % mesh_size] += p.mass * w6
        dm_rho[ix % mesh_size,     (iy+1) % mesh_size, (iz+1) % mesh_size] += p.mass * w7
        dm_rho[(ix+1) % mesh_size, (iy+1) % mesh_size, (iz+1) % mesh_size] += p.mass * w8
    
    dm_rho /= (cell_size**3) # Convert mass to density
    
    return dm_rho, cic_data, cells, particle_cells

def cic_force_interpolation(particles, ax_grid, ay_grid, az_grid, cic_data, config):
    mesh_size = config.mesh_size
    mesh_forces = []
    for i, p in enumerate(particles):
        ix, iy, iz, w1, w2, w3, w4, w5, w6, w7, w8 = cic_data[i]

        ax = (ax_grid[ix % mesh_size,     iy % mesh_size,     iz % mesh_size    ] * w1 +
              ax_grid[(ix+1) % mesh_size, iy % mesh_size,     iz % mesh_size    ] * w2 +
              ax_grid[ix % mesh_size,     (iy+1) % mesh_size, iz % mesh_size    ] * w3 +
              ax_grid[ix % mesh_size,     iy % mesh_size,     (iz+1) % mesh_size] * w4 +
              ax_grid[(ix+1) % mesh_size, (iy+1) % mesh_size, iz % mesh_size    ] * w5 +
              ax_grid[(ix+1) % mesh_size, iy % mesh_size,     (iz+1) % mesh_size] * w6 +
              ax_grid[ix % mesh_size,     (iy+1) % mesh_size, (iz+1) % mesh_size] * w7 +
              ax_grid[(ix+1) % mesh_size, (iy+1) % mesh_size, (iz+1) % mesh_size] * w8)

        ay = (ay_grid[ix % mesh_size,     iy % mesh_size,     iz % mesh_size    ] * w1 +
              ay_grid[(ix+1) % mesh_size, iy % mesh_size,     iz % mesh_size    ] * w2 +
              ay_grid[ix % mesh_size,     (iy+1) % mesh_size, iz % mesh_size    ] * w3 +
              ay_grid[ix % mesh_size,     iy % mesh_size,     (iz+1) % mesh_size] * w4 +
              ay_grid[(ix+1) % mesh_size, (iy+1) % mesh_size, iz % mesh_size    ] * w5 +
              ay_grid[(ix+1) % mesh_size, iy % mesh_size,     (iz+1) % mesh_size] * w6 +
              ay_grid[ix % mesh_size,     (iy+1) % mesh_size, (iz+1) % mesh_size] * w7 +
              ay_grid[(ix+1) % mesh_size, (iy+1) % mesh_size, (iz+1) % mesh_size] * w8)
              
        az = (az_grid[ix % mesh_size,     iy % mesh_size,     iz % mesh_size    ] * w1 +
              az_grid[(ix+1) % mesh_size, iy % mesh_size,     iz % mesh_size    ] * w2 +
              az_grid[ix % mesh_size,     (iy+1) % mesh_size, iz % mesh_size    ] * w3 +
              az_grid[ix % mesh_size,     iy % mesh_size,     (iz+1) % mesh_size] * w4 +
              az_grid[(ix+1) % mesh_size, (iy+1) % mesh_size, iz % mesh_size    ] * w5 +
              az_grid[(ix+1) % mesh_size, iy % mesh_size,     (iz+1) % mesh_size] * w6 +
              az_grid[ix % mesh_size,     (iy+1) % mesh_size, (iz+1) % mesh_size] * w7 +
              az_grid[(ix+1) % mesh_size, (iy+1) % mesh_size, (iz+1) % mesh_size] * w8)
        
        fx = ax * p.mass
        fy = ay * p.mass
        fz = az * p.mass
        mesh_forces.append((fx, fy, fz))
    
    return mesh_forces

def compute_gravitational_acceleration(gas, config, dm_rho):
    """
    Solves the Poisson equation for the gravitational potential using Fast Fourier Transforms (FFT) in k-space.
    """
    mesh_size = config.mesh_size
    cell_size = config.cell_size
    
    if config.enable_hydro:
        total_rho = dm_rho + gas.density
    else:
        total_rho = dm_rho

    # FFT of density
    rho_k = np.fft.rfftn(total_rho)

    # Poisson solver in k-space
    kx = np.fft.fftfreq(mesh_size, d=cell_size) * 2 * np.pi
    ky = np.fft.fftfreq(mesh_size, d=cell_size) * 2 * np.pi
    kz = np.fft.rfftfreq(mesh_size, d=cell_size) * 2 * np.pi # Note: rfftfreq
    
    kx_grid, ky_grid, kz_grid = np.meshgrid(kx, ky, kz, indexing='ij')
    k2 = kx_grid**2 + ky_grid**2 + kz_grid**2
    k2[0,0,0] = 1.0  # avoid div by zero
    
    # Green's function: phi_k = -4*pi*G * rho_k / k^2
    phi_k = -4 * math.pi * config.g_const * rho_k / k2 
    phi_k[0,0,0] = 0.0

    # Inverse FFT the POTENTIAL
    phi = np.fft.irfftn(phi_k, s=total_rho.shape)

    # Gradient to get field (finite differences)
    ax_grid = (np.roll(phi, 1, axis=0) - np.roll(phi, -1, axis=0)) / (2*cell_size)
    ay_grid = (np.roll(phi, 1, axis=1) - np.roll(phi, -1, axis=1)) / (2*cell_size)
    az_grid = (np.roll(phi, 1, axis=2) - np.roll(phi, -1, axis=2)) / (2*cell_size)

    return ax_grid, ay_grid, az_grid, total_rho

def compute_PM_short_range_approx(dist_sq, p1_mass, p2_mass, dx, dy, dz, config):
    softening_epsilon_squared = (0.5 * config.cell_size)**2 
    soft_dist_sq = dist_sq + softening_epsilon_squared
    soft_dist = math.sqrt(soft_dist_sq)
    f_pm_short = config.g_const * p1_mass * p2_mass / soft_dist_sq
    f_pm_short_x = f_pm_short * dx / soft_dist
    f_pm_short_y = f_pm_short * dy / soft_dist
    f_pm_short_z = f_pm_short * dz / soft_dist
    return f_pm_short_x, f_pm_short_y, f_pm_short_z

def compute_switching_parameter(dist_sq, config):
    if dist_sq < config.r_switch_start_sq:
        S = 1.0
    else:
        dist = math.sqrt(dist_sq)
        x = (dist - config.r_switch_start) / config.cutoff_transition_width
        S = 2*x**3 - 3*x**2 + 1
    return S

def compute_PP_forces(particles, config, cells, particle_cells):
    """
    Calculates PP forces.
    """
    mesh_size = config.mesh_size
    forces = [[0.0, 0.0, 0.0] for _ in particles]
    
    search_radius_cells = int(math.ceil(config.cutoff_radius / config.cell_size))

    for i, p1 in enumerate(particles):
        ix1, iy1, iz1 = particle_cells[i]
        
        # Check block of cells around the particle
        for dx_cell in range(-search_radius_cells, search_radius_cells + 1):
            for dy_cell in range(-search_radius_cells, search_radius_cells + 1):
                for dz_cell in range(-search_radius_cells, search_radius_cells + 1):
                
                    neighbor_ix = (ix1 + dx_cell) % mesh_size
                    neighbor_iy = (iy1 + dy_cell) % mesh_size
                    neighbor_iz = (iz1 + dz_cell) % mesh_size
                    cell_key = (neighbor_ix, neighbor_iy, neighbor_iz)
                    
                    if cell_key not in cells:
                        continue # Empty neighbor cells

                    for j in cells[cell_key]:
                        if i >= j:
                            continue
                                
                        p2 = particles[j]
                        
                        dx = displacement(p2.x - p1.x, config)
                        dy = displacement(p2.y - p1.y, config)
                        dz = displacement(p2.z - p1.z, config)
                        dist_sq = dx*dx + dy*dy + dz*dz

                        if dist_sq > config.cutoff_radius_squared:
                            continue
                        
                        # Subtractive scheme
                        fx_sub, fy_sub, fz_sub = compute_PM_short_range_approx(dist_sq, p1.mass, p2.mass, dx, dy, dz, config)
                        S = compute_switching_parameter(dist_sq, config)
                        dist_sq_soft = dist_sq + config.softening_squared
                        dist = math.sqrt(dist_sq_soft)
                        f = config.g_const * p1.mass * p2.mass / dist_sq_soft
                        
                        fx = S * (f * dx / dist - fx_sub)
                        fy = S * (f * dy / dist - fy_sub)
                        fz = S * (f * dz / dist - fz_sub)
                
                        forces[i][0] += fx
                        forces[i][1] += fy
                        forces[i][2] += fz

                        forces[j][0] -= fx
                        forces[j][1] -= fy
                        forces[j][2] -= fz
            
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

def KDK_step(total_time: float, particles: list, gas: GasGrid, dt: float, config: Config) -> tuple:
    
    timings = {} # To store performance breakdown
    
    scale_factor, hubble_param = update_cosmology(total_time, config)
    
    # 1. KICK (First half-step for velocity)
    for p in particles:
        total_ax = (p.ax / scale_factor**3) - (2 * hubble_param * p.vx)
        total_ay = (p.ay / scale_factor**3) - (2 * hubble_param * p.vy)
        total_az = (p.az / scale_factor**3) - (2 * hubble_param * p.vz)
        p.vx += total_ax * dt / 2.0
        p.vy += total_ay * dt / 2.0
        p.vz += total_az * dt / 2.0

    if config.enable_hydro:
        gas.update_primitive_variables()
        gas_accel_x = (gas.accel_x / scale_factor**3) - (2 * hubble_param * gas.velocity_x)
        gas_accel_y = (gas.accel_y / scale_factor**3) - (2 * hubble_param * gas.velocity_y)
        gas_accel_z = (gas.accel_z / scale_factor**3) - (2 * hubble_param * gas.velocity_z)
        
        power_density = gas.density * (gas_accel_x * gas.velocity_x + 
                                       gas_accel_y * gas.velocity_y + 
                                       gas_accel_z * gas.velocity_z)
        
        gas.momentum_x += gas.density * gas_accel_x * dt / 2.0
        gas.momentum_y += gas.density * gas_accel_y * dt / 2.0
        gas.momentum_z += gas.density * gas_accel_z * dt / 2.0
        gas.energy += power_density * dt / 2.0
    
    # 2. DRIFT (Full step for position)
    for p in particles:
        p.x = ( p.x + p.vx * dt ) % config.domain_size
        p.y = ( p.y + p.vy * dt ) % config.domain_size
        p.z = ( p.z + p.vz * dt ) % config.domain_size

    # 3. UPDATE COSMOLOGY
    scale_factor_new, hubble_param_new = update_cosmology(total_time + dt, config)

    # 4. COMPUTE FORCES at new positions
    t_pm_start = time.time()
    dm_rho, cic_data, cells, particle_cells = particle_binning_and_mass_assignment(particles, config)
    ax_grid, ay_grid, az_grid, total_rho = compute_gravitational_acceleration(gas, config, dm_rho)
    pm_forces = cic_force_interpolation(particles, ax_grid, ay_grid, az_grid, cic_data, config)
    timings['pm'] = time.time() - t_pm_start
    
    t_pp_start = time.time()
    pp_forces = compute_PP_forces(particles, config, cells, particle_cells)
    timings['pp'] = time.time() - t_pp_start

    forces = [(pp[0]+pm[0], pp[1]+pm[1], pp[2]+pm[2]) for pp,pm in zip(pp_forces, pm_forces)]
    
    t_hydro_start = time.time()
    if config.enable_hydro:
        gas.hydro_step(dt)
    timings['hydro'] = time.time() - t_hydro_start

    # 5. KICK (Second half-step for velocity)
    for i, p in enumerate(particles):
        fx, fy, fz = forces[i]
        p.ax = fx / p.mass
        p.ay = fy / p.mass
        p.az = fz / p.mass
        total_ax = (p.ax / scale_factor_new**3) - (2 * hubble_param_new * p.vx)
        total_ay = (p.ay / scale_factor_new**3) - (2 * hubble_param_new * p.vy)
        total_az = (p.az / scale_factor_new**3) - (2 * hubble_param_new * p.vz)
        p.vx += total_ax * dt / 2.0
        p.vy += total_ay * dt / 2.0
        p.vz += total_az * dt / 2.0

    if config.enable_hydro:
        gas.update_primitive_variables()
        gas.accel_x = ax_grid
        gas.accel_y = ay_grid
        gas.accel_z = az_grid
        gas_accel_x = (gas.accel_x / scale_factor_new**3) - (2 * hubble_param_new * gas.velocity_x)
        gas_accel_y = (gas.accel_y / scale_factor_new**3) - (2 * hubble_param_new * gas.velocity_y)
        gas_accel_z = (gas.accel_z / scale_factor_new**3) - (2 * hubble_param_new * gas.velocity_z)
        
        power_density = gas.density * (gas_accel_x * gas.velocity_x + 
                                       gas_accel_y * gas.velocity_y +
                                       gas_accel_z * gas.velocity_z)

        gas.momentum_x += gas.density * gas_accel_x * dt / 2.0
        gas.momentum_y += gas.density * gas_accel_y * dt / 2.0
        gas.momentum_z += gas.density * gas_accel_z * dt / 2.0
        gas.energy += power_density * dt / 2.0
        
    return timings, total_rho

# Create initial conditions
def create_zeldovich_ics(config):
    n_per_side = config.num_dm_particles_side
    ic_mesh_size = n_per_side
    cell_size = config.domain_size / ic_mesh_size
    grid_shape = (ic_mesh_size, ic_mesh_size, ic_mesh_size)

    # Build k-grid
    kx = np.fft.fftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi
    ky = np.fft.fftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi
    kz = np.fft.rfftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi # Real FFT
    
    kx_grid, ky_grid, kz_grid = np.meshgrid(kx, ky, kz, indexing='ij')
    k2 = kx_grid**2 + ky_grid**2 + kz_grid**2
    k2[0, 0, 0] = 1.0

    # Generate random density field
    real_space_random_field = np.random.randn(*grid_shape)
    random_k = np.fft.rfftn(real_space_random_field)
    
    power_spectrum_index = config.initial_power_spectrum_index
    delta_k = random_k * (k2 ** (power_spectrum_index / 4.0)) # delta_k ~ P(k)^1/2 ~ k^(n/2)
    delta_k[0, 0, 0] = 0.0

    # Poisson solver for the Zel'dovich potential: phi_k = -4*pi*G * delta_k / k^2
    phi_k = -4.0 * math.pi * config.g_const * delta_k / k2
    phi_k[0, 0, 0] = 0.0

    # Displacement field s_k = i * k * (delta_k / k^2)
    disp_x_k = 1j * kx_grid * phi_k
    disp_y_k = 1j * ky_grid * phi_k
    disp_z_k = 1j * kz_grid * phi_k

    disp_x = np.fft.irfftn(disp_x_k, s=grid_shape)
    disp_y = np.fft.irfftn(disp_y_k, s=grid_shape)
    disp_z = np.fft.irfftn(disp_z_k, s=grid_shape)

    # Normalize & apply amplitude
    disp_x /= np.std(disp_x)
    disp_y /= np.std(disp_y)
    disp_z /= np.std(disp_z)
    
    disp_x *= config.initial_scale_factor * cell_size
    disp_y *= config.initial_scale_factor * cell_size
    disp_z *= config.initial_scale_factor * cell_size

    # Create displaced particles
    particles = []
    spacing = config.domain_size / n_per_side
    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                qx = (i + 0.5) * spacing
                qy = (j + 0.5) * spacing
                qz = (k + 0.5) * spacing
                
                x = (qx + disp_x[i, j, k]) % config.domain_size
                y = (qy + disp_y[i, j, k]) % config.domain_size
                z = (qz + disp_z[i, j, k]) % config.domain_size
                
                vx = config.initial_hubble_param * disp_x[i, j, k]
                vy = config.initial_hubble_param * disp_y[i, j, k]
                vz = config.initial_hubble_param * disp_z[i, j, k]
                
                particles.append(Particle(x, y, z, config.dm_particle_mass, vx, vy, vz))

    return particles

def calculate_particle_energies(particles, scale_factor, config):
    kinetic_energy = 0.0
    potential_energy = 0.0
    g_const = config.g_const
    softening_squared = config.softening_squared

    for p in particles:
        vel_sq = (scale_factor*p.vx)**2 + (scale_factor*p.vy)**2 + (scale_factor*p.vz)**2
        kinetic_energy += 0.5 * p.mass * vel_sq

    # Note: This O(N^2) PE calculation is slow! Only for diagnostics.
    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            p1 = particles[i]
            p2 = particles[j]

            dx = displacement(p2.x - p1.x, config)
            dy = displacement(p2.y - p1.y, config)
            dz = displacement(p2.z - p1.z, config)
            
            dist_sq = (scale_factor**2) * (dx*dx + dy*dy + dz*dz + softening_squared)
            dist = math.sqrt(dist_sq)
            if dist > 0:
                potential_energy -= g_const * p1.mass * p2.mass / dist

    return kinetic_energy, potential_energy

def get_cfl_timestep(gas, config):
    """Calculates the CFL timestep limit."""
    if not config.enable_hydro:
        return np.inf
        
    non_zero_density = gas.density > 1e-12
    if not np.any(non_zero_density):
        return np.inf
        
    cs = np.sqrt(config.gamma * gas.pressure[non_zero_density] / gas.density[non_zero_density])
    v_mag = np.sqrt(gas.velocity_x[non_zero_density]**2 + 
                    gas.velocity_y[non_zero_density]**2 +
                    gas.velocity_z[non_zero_density]**2)
    
    max_signal_vel = np.max(v_mag + cs)
    if max_signal_vel < 1e-9:
        return np.inf
        
    dt_cfl = config.cell_size / max_signal_vel
    return dt_cfl * config.cfl_safety_factor
    
def get_gravity_timestep(particles, config):
    """Calculates the gravity acceleration timestep limit."""
    if not particles:
        return np.inf
        
    max_accel_sq = 1e-9
    for p in particles:
        accel_sq = p.ax**2 + p.ay**2 + p.az**2
        if accel_sq > max_accel_sq:
            max_accel_sq = accel_sq
            
    dt_grav = np.sqrt(config.softening_squared) / np.sqrt(max_accel_sq)
    return dt_grav * config.gravity_dt_factor

def calculate_diagnostics(particles, gas, timings, dt, cycle, sim_time, scale_factor, config):
    """Fills a Diagnostics object with data from the simulation state."""
    
    diag = Diagnostics()
    diag.cycle = cycle
    diag.sim_time = sim_time
    diag.scale_factor = scale_factor

    # Performance
    diag.wall_time_pm = timings.get('pm', 0.0)
    diag.wall_time_pp = timings.get('pp', 0.0)
    diag.wall_time_hydro = timings.get('hydro', 0.0)
    diag.wall_time_io = timings.get('io', 0.0)
    diag.wall_time_total = diag.wall_time_pm + diag.wall_time_pp + diag.wall_time_hydro + diag.wall_time_io

    # Stability
    diag.dt_cfl = get_cfl_timestep(gas, config)
    diag.dt_gravity = get_gravity_timestep(particles, config)
    diag.dt_final = dt
    
    if config.enable_hydro:
        diag.max_gas_density = np.max(gas.density)
        diag.max_gas_pressure = np.max(gas.pressure)
        v_mag = np.sqrt(gas.velocity_x**2 + gas.velocity_y**2 + gas.velocity_z**2)
        diag.max_gas_velocity = np.max(v_mag)

    # Conservation
    diag.total_mass_dm = sum(p.mass for p in particles)
    diag.total_momentum_dm = (sum(p.mass * p.vx for p in particles),
                             sum(p.mass * p.vy for p in particles),
                             sum(p.mass * p.vz for p in particles))
    diag.ke_dm, diag.pe_dm = calculate_particle_energies(particles, scale_factor, config)

    if config.enable_hydro:
        diag.total_mass_gas = np.sum(gas.density) * config.cell_volume
        diag.total_momentum_gas = (np.sum(gas.momentum_x), 
                                  np.sum(gas.momentum_y), 
                                  np.sum(gas.momentum_z))
        
        ke_gas_density = 0.5 * (gas.momentum_x**2 + gas.momentum_y**2 + gas.momentum_z**2) / (gas.density + 1e-12)
        diag.ke_gas = np.sum(ke_gas_density) * config.cell_volume
        
        internal_energy_density = gas.pressure / (config.gamma - 1.0)
        diag.ie_gas = np.sum(internal_energy_density) * config.cell_volume
        
    return diag

# ------------------------
# Visualization (VisPy)
# ------------------------

class SimulationApp:
    def __init__(self, config, particles, gas, h5_file, logger):
        self.config = config
        self.particles = particles
        self.gas = gas
        self.h5_file = h5_file
        self.logger = logger
        
        self.cycle_count = 0
        self.total_simulation_time = 0.0
        self.snapshot_count = 0
        
        # --- VisPy Scene Setup ---
        self.canvas = scene.SceneCanvas(keys='interactive', 
                                        size=(config.render_size, config.render_size), 
                                        show=True, 
                                        title='N-Body + Hydrodynamics Simulation')
        
        self.view = self.canvas.central_widget.add_view()
        
        # --- Camera Setup ---
        self.view.camera = 'turntable' # orbital camera
        self.view.camera.fov = 60
        # Set camera center and distance
        self.view.camera.center = (config.domain_size / 2, 
                                   config.domain_size / 2, 
                                   config.domain_size / 2)
        self.view.camera.distance = config.domain_size * 2

        # --- Create Visuals ---
        self.colormap = get_colormap('plasma')
        
        # Gas: Render the full Volume
        if config.enable_hydro:
            # We must pass data in (Z, Y, X) order
            initial_volume_data = np.zeros((config.mesh_size, config.mesh_size, config.mesh_size))
            
            self.gas_volume = scene.visuals.Volume(
                initial_volume_data,
                parent=self.view.scene,
                cmap=self.colormap,
                method='mip', # Maximum Intensity Projection
                threshold=0.1 # Start by only showing brighter parts
            )
            # Set transform to scale the volume to the domain size
            self.gas_volume.transform = scene.transforms.STTransform(
                scale=(config.cell_size, config.cell_size, config.cell_size),
                translate=(0, 0, 0)
            )
        
        # Particles: Render as markers
        particle_pos = np.zeros((len(particles), 3)) 
        self.particles_visual = scene.visuals.Markers(parent=self.view.scene)
        self.particles_visual.set_data(particle_pos, 
                                       face_color='white', 
                                       edge_color=None, 
                                       size=2)
        
        # Add a box outline
        box = scene.visuals.Box(width=config.domain_size, 
                          height=config.domain_size, 
                          depth=config.domain_size, 
                          color=(1, 1, 1, 0), # Transparent faces
                          edge_color=(1, 1, 1, 0.8), # Dim white edges
                          parent=self.view.scene)
        
        box.transform = scene.transforms.STTransform(
            translate=(config.domain_size / 2, 
                       config.domain_size / 2, 
                       config.domain_size / 2)
        )

        # --- Timestep & Initial State ---
        self.dt_cfl = get_cfl_timestep(self.gas, self.config)
        self.dt_grav = get_gravity_timestep(self.particles, self.config)
        self.current_dt = min(self.dt_cfl, self.dt_grav, self.config.fixed_dt) if self.config.use_adaptive_dt else self.config.fixed_dt

        # Log initial state (Cycle 0)
        scale_factor, _ = update_cosmology(self.total_simulation_time, self.config)
        diag = calculate_diagnostics(self.particles, self.gas, {}, self.current_dt, 0, 0.0, scale_factor, self.config)
        self.logger.log(diag)
        
        self.update_visuals() # Initial render

        # --- Setup Simulation Timer ---
        self.timer = app.Timer(interval='auto', connect=self.update_simulation, start=True)
        self.canvas.events.close.connect(self.on_close)

    def update_visuals(self):
        """Update the visuals on the GPU with the new simulation data."""
        
        # 1. Update Gas Volume
        if self.config.enable_hydro:
            # Normalize the pressure field for visualization
            log_field = np.log10(self.gas.pressure + 1e-12)
            min_log, max_log = np.min(log_field), np.max(log_field)
            
            if max_log > min_log:
                norm_field = (log_field - min_log) / (max_log - min_log)
            else:
                norm_field = np.zeros_like(log_field)
            
            # VisPy's Volume visual expects data in (Z, Y, X) order.
            # Our simulation data is in (X, Y, Z) order. We must transpose.
            volume_data = norm_field.transpose(2, 1, 0).astype(np.float32)
            
            self.gas_volume.set_data(volume_data)
            self.gas_volume.clim = (0.0, 1.0) # We've already normalized

        # 2. Update Particle Positions
        particle_pos = np.array([[p.x, p.y, p.z] for p in self.particles])
        self.particles_visual.set_data(particle_pos, 
                                       face_color='white', 
                                       edge_color='white', 
                                       size=1)
        
        self.canvas.update()

    def update_simulation(self, event):
        """This is the main simulation loop, called by the VisPy timer."""
        
        timings, total_rho = KDK_step(self.total_simulation_time, self.particles, self.gas, self.current_dt, self.config)
        
        self.total_simulation_time += self.current_dt
        self.cycle_count += 1
        scale_factor, _ = update_cosmology(self.total_simulation_time, self.config)
        
        if self.config.use_adaptive_dt:
            self.dt_cfl = get_cfl_timestep(self.gas, self.config)
            self.dt_grav = get_gravity_timestep(self.particles, self.config)
            self.current_dt = min(self.dt_cfl, self.dt_grav, self.config.fixed_dt)
        else:
            self.current_dt = self.config.fixed_dt

        # Update visuals (every 5 steps to run faster)
        if self.cycle_count % 5 == 0:
            self.update_visuals()

        # I/O and Logging
        t_io_start = time.time()
        if self.config.save_snapshot_every_cycles > 0 and self.cycle_count % self.config.save_snapshot_every_cycles == 0:
            snapshot_name = f"snapshot_{self.snapshot_count:04d}"
            save_snapshot(self.h5_file, snapshot_name, self.particles, self.gas, self.total_simulation_time, scale_factor, self.config)
            print(f"Saved snapshot: {snapshot_name}")
            self.snapshot_count += 1
        timings['io'] = time.time() - t_io_start

        if self.config.debug_info_every_cycles > 0 and self.cycle_count % self.config.debug_info_every_cycles == 0:
            diag = calculate_diagnostics(self.particles, self.gas, timings, self.current_dt, self.cycle_count, self.total_simulation_time, scale_factor, self.config)
            self.logger.log(diag)

    def on_close(self, event):
        """Called when the user closes the VisPy window."""
        print("Closing simulation...")
        self.timer.stop()
        self.h5_file.close()
        self.logger.close()
        print("Simulation finished. HDF5 and log files closed.")
        self.canvas.app.quit()

    def run(self):
        """Starts the VisPy application event loop."""
        self.canvas.app.run()


def save_snapshot(h5_file, snapshot_name, particles, gas, sim_time, scale_factor, config):
    """Saves the simulation state to HDF5."""
    try:
        grp = h5_file.create_group(snapshot_name)
        
        grp.attrs['time'] = sim_time
        grp.attrs['scale_factor'] = scale_factor
        
        # Store Particle Data
        num_p = len(particles)
        positions = np.zeros((num_p, 3))
        velocities = np.zeros((num_p, 3))
        masses = np.zeros(num_p)
        
        for i, p in enumerate(particles):
            positions[i] = [p.x, p.y, p.z]
            velocities[i] = [p.vx, p.vy, p.vz]
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
            g_grp.create_dataset("momentum_z", data=gas.momentum_z, compression="gzip")
            g_grp.create_dataset("energy", data=gas.energy, compression="gzip")
            g_grp.create_dataset("pressure", data=gas.pressure, compression="gzip") 
            g_grp.create_dataset("velocity_x", data=gas.velocity_x, compression="gzip")
            g_grp.create_dataset("velocity_y", data=gas.velocity_y, compression="gzip")
            g_grp.create_dataset("velocity_z", data=gas.velocity_z, compression="gzip")
            
    except Exception as e:
        print(f"Error saving snapshot {snapshot_name}: {e}")

# ------------------------
# Main entry point
# ------------------------
def main():
    # Load Configuration
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
    
    run_timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = Logger(log_filename=f"sim_{run_timestamp}_diagnostics.csv")
    print(f"Logging diagnostics to: {logger.log_filename}")

    # --- Initialize Simulation State ---
    particles = create_zeldovich_ics(config)
    gas = GasGrid(config)
    if config.enable_hydro:
        gas.density.fill(config.gas_total_mass / (config.domain_size**3))
        initial_internal_energy = 1e-6
        gas.energy = initial_internal_energy * gas.density
        gas.update_primitive_variables()
    
    # --- Calculate initial forces ---
    dm_rho, cic_data, cells, particle_cells = particle_binning_and_mass_assignment(particles, config)
    ax_grid, ay_grid, az_grid, total_rho = compute_gravitational_acceleration(gas, config, dm_rho)
    pm_forces = cic_force_interpolation(particles, ax_grid, ay_grid, az_grid, cic_data, config)
    pp_forces = compute_PP_forces(particles, config, cells, particle_cells)
    
    forces = [(pp[0]+pm[0], pp[1]+pm[1], pp[2]+pm[2]) for pp,pm in zip(pp_forces, pm_forces)]
    for i, p in enumerate(particles):
        p.ax = forces[i][0] / p.mass
        p.ay = forces[i][1] / p.mass
        p.az = forces[i][2] / p.mass
    if config.enable_hydro:
        gas.accel_x = ax_grid
        gas.accel_y = ay_grid
        gas.accel_z = az_grid

    # --- Create HDF5 file ---
    snapshot_filename = f"sim_{run_timestamp}.hdf5"
    print(f"Saving snapshots to: {snapshot_filename}")
    h5_file = h5py.File(snapshot_filename, 'w')
    param_grp = h5_file.create_group("parameters")
    for key, val in config.get_all_params().items():
        param_grp.attrs[key] = val
        
    # --- Create and run the VisPy Application ---
    sim_app = SimulationApp(config, particles, gas, h5_file, logger)
    
    print("Starting simulation...")
    if (sys.flags.interactive != 1):
        sim_app.run()

if __name__ == "__main__":
    main()