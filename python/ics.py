import math
import numpy as np

from gas import GasGrid
from particles import (Particle, particle_binning_and_mass_assignment, 
                       cic_force_interpolation, compute_PP_forces)
from integrator import compute_gravitational_acceleration
from state import SimState
from integrator import update_cosmology

def create_zeldovich_ics(config):
    """Generates particles and gas velocity grids using the Zel'dovich approximation."""
    # Calculate the primordial fields at the Eulerian grid resolution
    ic_mesh_size = config.mesh_size
    cell_size = config.domain_size / ic_mesh_size
    grid_shape = (ic_mesh_size, ic_mesh_size, ic_mesh_size)

    kx = np.fft.fftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi
    ky = np.fft.fftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi
    kz = np.fft.rfftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi 
    
    kx_grid, ky_grid, kz_grid = np.meshgrid(kx, ky, kz, indexing='ij')
    k2 = kx_grid**2 + ky_grid**2 + kz_grid**2
    k2[0, 0, 0] = 1.0

    random_k = np.fft.rfftn(np.random.randn(*grid_shape))
    delta_k = random_k * (k2 ** (config.initial_power_spectrum_index / 4.0)) 
    delta_k[0, 0, 0] = 0.0

    phi_k = -4.0 * math.pi * config.g_const * delta_k / k2
    phi_k[0, 0, 0] = 0.0

    disp_x = np.fft.irfftn(1j * kx_grid * phi_k, s=grid_shape)
    disp_y = np.fft.irfftn(1j * ky_grid * phi_k, s=grid_shape)
    disp_z = np.fft.irfftn(1j * kz_grid * phi_k, s=grid_shape)
    
    disp_x = disp_x / np.std(disp_x)
    disp_y = disp_y / np.std(disp_y)
    disp_z = disp_z / np.std(disp_z)
    
    disp_x *= config.initial_scale_factor * cell_size
    disp_y *= config.initial_scale_factor * cell_size
    disp_z *= config.initial_scale_factor * cell_size

    gas_vx = config.initial_hubble_param * disp_x
    gas_vy = config.initial_hubble_param * disp_y
    gas_vz = config.initial_hubble_param * disp_z

    # Sample the high-resolution field to create the particles
    particles = []
    n_per_side = config.num_dm_particles_side
    spacing = config.domain_size / n_per_side

    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                # Particle's perfectly uniform starting coordinate
                qx, qy, qz = (i + 0.5) * spacing, (j + 0.5) * spacing, (k + 0.5) * spacing
                
                # Map physical coordinate to the high-res grid index
                ix = int(qx / cell_size) % ic_mesh_size
                iy = int(qy / cell_size) % ic_mesh_size
                iz = int(qz / cell_size) % ic_mesh_size

                # Displace the particle using that cell's vector
                x = (qx + disp_x[ix, iy, iz]) % config.domain_size
                y = (qy + disp_y[ix, iy, iz]) % config.domain_size
                z = (qz + disp_z[ix, iy, iz]) % config.domain_size
                
                # Assign the particle velocity
                vx = gas_vx[ix, iy, iz]
                vy = gas_vy[ix, iy, iz]
                vz = gas_vz[ix, iy, iz]
                
                particles.append(Particle(x, y, z, config.dm_particle_mass, vx, vy, vz))

    return particles, gas_vx, gas_vy, gas_vz

def initialize_state(config):
    """Initializes particles, gas grid, and computes the 'Step 0' forces."""
    particles, gas_vx, gas_vy, gas_vz = create_zeldovich_ics(config)
    gas = GasGrid(config)
    
    # Bin Dark Matter to get the primordial density structure FIRST
    dm_rho, cic_data, cells, particle_cells = particle_binning_and_mass_assignment(particles, config)
    
    # Initialize Gas to perfectly trace the Dark Matter
    if config.enable_hydro:
        # Gas density is directly proportional to DM density
        total_dm_mass = len(particles) * config.dm_particle_mass
        mass_ratio = config.gas_total_mass / total_dm_mass
        
        gas.density = dm_rho * mass_ratio
        
        # Prevent hydro crashes in deep voids by setting a density floor
        gas.density = np.maximum(gas.density, 1e-12)

        gas.momentum_x = gas.density * gas_vx
        gas.momentum_y = gas.density * gas_vy
        gas.momentum_z = gas.density * gas_vz

        initial_internal_energy = 1e-6
        kin_energy = 0.5 * (gas.momentum_x**2 + gas.momentum_y**2 + gas.momentum_z**2) / gas.density
        
        gas.energy = (gas.density * initial_internal_energy) + kin_energy
        gas.update_primitive_variables()
        
    # Compute Step 0 Gravity Forces
    ax_grid, ay_grid, az_grid, total_rho = compute_gravitational_acceleration(gas, config, dm_rho)
    pm_forces = cic_force_interpolation(particles, ax_grid, ay_grid, az_grid, cic_data, config)
    pp_forces = compute_PP_forces(particles, config, cells, particle_cells)
    
    forces = [(pp[0]+pm[0], pp[1]+pm[1], pp[2]+pm[2]) for pp,pm in zip(pp_forces, pm_forces)]
    for i, p in enumerate(particles):
        p.ax, p.ay, p.az = forces[i][0] / p.mass, forces[i][1] / p.mass, forces[i][2] / p.mass
        
    if config.enable_hydro:
        gas.accel_x, gas.accel_y, gas.accel_z = ax_grid, ay_grid, az_grid

    initial_scale_factor, _ = update_cosmology(0.0, config)
    
    return SimState(particles=particles, gas=gas, total_time=0.0, scale_factor=initial_scale_factor)