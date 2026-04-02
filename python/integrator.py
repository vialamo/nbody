import math
import time
import numpy as np

from gas import get_cfl_timestep
from particles import (particle_binning_and_mass_assignment, cic_force_interpolation, 
                       compute_PP_forces, calculate_particle_energies, get_gravity_timestep)
from utils import Diagnostics

def compute_gravitational_acceleration(gas, config, dm_rho):
    mesh_size, cell_size = config.mesh_size, config.cell_size
    total_rho = dm_rho + gas.density if config.enable_hydro else dm_rho

    rho_k = np.fft.rfftn(total_rho)
    kx = np.fft.fftfreq(mesh_size, d=cell_size) * 2 * np.pi
    ky = np.fft.fftfreq(mesh_size, d=cell_size) * 2 * np.pi
    kz = np.fft.rfftfreq(mesh_size, d=cell_size) * 2 * np.pi
    
    kx_grid, ky_grid, kz_grid = np.meshgrid(kx, ky, kz, indexing='ij')
    k2 = kx_grid**2 + ky_grid**2 + kz_grid**2
    k2[0,0,0] = 1.0  
    
    phi_k = -4 * math.pi * config.g_const * rho_k / k2 
    phi_k[0,0,0] = 0.0
    phi = np.fft.irfftn(phi_k, s=total_rho.shape)

    ax_grid = (np.roll(phi, 1, axis=0) - np.roll(phi, -1, axis=0)) / (2*cell_size)
    ay_grid = (np.roll(phi, 1, axis=1) - np.roll(phi, -1, axis=1)) / (2*cell_size)
    az_grid = (np.roll(phi, 1, axis=2) - np.roll(phi, -1, axis=2)) / (2*cell_size)
    return ax_grid, ay_grid, az_grid, total_rho

def update_cosmology(total_time, config):
    if not config.expanding_universe: return 1.0, 0.0
    expansion_time = config.expansion_start_t + total_time
    return expansion_time**(2.0/3.0), (2.0/3.0) / expansion_time

def KDK_step(state, dt, config) -> tuple:
    """
    Performs a single Kick-Drift-Kick timestep for both particles and gas.
    Updates the simulation state (including time and cosmology) in-place.
    """
    timings = {}
    
    # Read current state cosmology for the first kick
    scale_factor, hubble_param = update_cosmology(state.total_time, config)
    
    # --- 1. KICK (First half-step for velocity) ---
    for p in state.particles:
        total_ax = (p.ax / scale_factor**3) - (2 * hubble_param * p.vx)
        total_ay = (p.ay / scale_factor**3) - (2 * hubble_param * p.vy)
        total_az = (p.az / scale_factor**3) - (2 * hubble_param * p.vz)
        p.vx += total_ax * dt / 2.0
        p.vy += total_ay * dt / 2.0
        p.vz += total_az * dt / 2.0

    if config.enable_hydro:
        state.gas.update_primitive_variables()
        gas_accel_x = (state.gas.accel_x / scale_factor**3) - (2 * hubble_param * state.gas.velocity_x)
        gas_accel_y = (state.gas.accel_y / scale_factor**3) - (2 * hubble_param * state.gas.velocity_y)
        gas_accel_z = (state.gas.accel_z / scale_factor**3) - (2 * hubble_param * state.gas.velocity_z)
        
        power_density = state.gas.density * (gas_accel_x * state.gas.velocity_x + 
                                             gas_accel_y * state.gas.velocity_y + 
                                             gas_accel_z * state.gas.velocity_z)
        
        state.gas.momentum_x += state.gas.density * gas_accel_x * dt / 2.0
        state.gas.momentum_y += state.gas.density * gas_accel_y * dt / 2.0
        state.gas.momentum_z += state.gas.density * gas_accel_z * dt / 2.0
        state.gas.energy += power_density * dt / 2.0
    
    # --- 2. DRIFT (Full step for position) ---
    for p in state.particles:
        p.x = ( p.x + p.vx * dt ) % config.domain_size
        p.y = ( p.y + p.vy * dt ) % config.domain_size
        p.z = ( p.z + p.vz * dt ) % config.domain_size

    # --- 3. UPDATE TIME AND COSMOLOGY ---
    # The integrator owns the progression of time and expansion
    state.total_time += dt
    state.scale_factor, hubble_param_new = update_cosmology(state.total_time, config)

    # --- 4. COMPUTE FORCES at new positions ---
    t_pm_start = time.time()
    dm_rho, cic_data, cells, particle_cells = particle_binning_and_mass_assignment(state.particles, config)
    ax_grid, ay_grid, az_grid, total_rho = compute_gravitational_acceleration(state.gas, config, dm_rho)
    pm_forces = cic_force_interpolation(state.particles, ax_grid, ay_grid, az_grid, cic_data, config)
    timings['pm'] = time.time() - t_pm_start
    
    t_pp_start = time.time()
    pp_forces = compute_PP_forces(state.particles, config, cells, particle_cells)
    timings['pp'] = time.time() - t_pp_start

    # Combine long-range (PM) and short-range (PP) forces
    forces = [(pp[0]+pm[0], pp[1]+pm[1], pp[2]+pm[2]) for pp,pm in zip(pp_forces, pm_forces)]
    
    t_hydro_start = time.time()
    if config.enable_hydro:
        state.gas.hydro_step(dt)
    timings['hydro'] = time.time() - t_hydro_start

    # --- 5. KICK (Second half-step for velocity) ---
    for i, p in enumerate(state.particles):
        p.ax = forces[i][0] / p.mass
        p.ay = forces[i][1] / p.mass
        p.az = forces[i][2] / p.mass
        
        total_ax = (p.ax / state.scale_factor**3) - (2 * hubble_param_new * p.vx)
        total_ay = (p.ay / state.scale_factor**3) - (2 * hubble_param_new * p.vy)
        total_az = (p.az / state.scale_factor**3) - (2 * hubble_param_new * p.vz)
        
        p.vx += total_ax * dt / 2.0
        p.vy += total_ay * dt / 2.0
        p.vz += total_az * dt / 2.0

    if config.enable_hydro:
        state.gas.update_primitive_variables()
        state.gas.accel_x = ax_grid
        state.gas.accel_y = ay_grid
        state.gas.accel_z = az_grid
        
        gas_accel_x = (state.gas.accel_x / state.scale_factor**3) - (2 * hubble_param_new * state.gas.velocity_x)
        gas_accel_y = (state.gas.accel_y / state.scale_factor**3) - (2 * hubble_param_new * state.gas.velocity_y)
        gas_accel_z = (state.gas.accel_z / state.scale_factor**3) - (2 * hubble_param_new * state.gas.velocity_z)
        
        power_density = state.gas.density * (gas_accel_x * state.gas.velocity_x + 
                                             gas_accel_y * state.gas.velocity_y +
                                             gas_accel_z * state.gas.velocity_z)

        state.gas.momentum_x += state.gas.density * gas_accel_x * dt / 2.0
        state.gas.momentum_y += state.gas.density * gas_accel_y * dt / 2.0
        state.gas.momentum_z += state.gas.density * gas_accel_z * dt / 2.0
        state.gas.energy += power_density * dt / 2.0
        
    return timings

def calculate_diagnostics(particles, gas, timings, dt, cycle, sim_time, scale_factor, config):
    diag = Diagnostics()
    diag.cycle, diag.sim_time, diag.scale_factor = cycle, sim_time, scale_factor

    diag.wall_time_pm = timings.get('pm', 0.0)
    diag.wall_time_pp = timings.get('pp', 0.0)
    diag.wall_time_hydro = timings.get('hydro', 0.0)
    diag.wall_time_io = timings.get('io', 0.0)
    diag.wall_time_total = diag.wall_time_pm + diag.wall_time_pp + diag.wall_time_hydro + diag.wall_time_io

    diag.dt_cfl, diag.dt_gravity, diag.dt_final = get_cfl_timestep(gas, config), get_gravity_timestep(particles, config), dt
    
    if config.enable_hydro:
        diag.max_gas_density = np.max(gas.density)
        diag.max_gas_pressure = np.max(gas.pressure)
        diag.max_gas_velocity = np.max(np.sqrt(gas.velocity_x**2 + gas.velocity_y**2 + gas.velocity_z**2))

    diag.total_mass_dm = sum(p.mass for p in particles)
    diag.total_momentum_dm = (sum(p.mass * p.vx for p in particles), sum(p.mass * p.vy for p in particles), sum(p.mass * p.vz for p in particles))
    diag.ke_dm, diag.pe_dm = calculate_particle_energies(particles, scale_factor, config)

    if config.enable_hydro:
        diag.total_mass_gas = np.sum(gas.density) * config.cell_volume
        diag.total_momentum_gas = (np.sum(gas.momentum_x), np.sum(gas.momentum_y), np.sum(gas.momentum_z))
        ke_gas_density = 0.5 * (gas.momentum_x**2 + gas.momentum_y**2 + gas.momentum_z**2) / (gas.density + 1e-12)
        diag.ke_gas = np.sum(ke_gas_density) * config.cell_volume
        diag.ie_gas = np.sum(gas.pressure / (config.gamma - 1.0)) * config.cell_volume
        
    return diag