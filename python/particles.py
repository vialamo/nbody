import math
import numpy as np

class Particle:
    def __init__(self, x, y, z, mass, vx=0, vy=0, vz=0):
        self.x, self.y, self.z = x, y, z
        self.mass = mass
        self.vx, self.vy, self.vz = vx, vy, vz
        self.ax, self.ay, self.az = 0, 0, 0

def displacement(dx, config):
    L = config.domain_size
    return (dx + 0.5*L) % L - 0.5*L

def particle_binning_and_mass_assignment(particles, config):
    mesh_size, cell_size = config.mesh_size, config.cell_size
    dm_rho = np.zeros((mesh_size, mesh_size, mesh_size))
    cic_data = [] 
    cells = {} 
    particle_cells = [] 

    for i, p in enumerate(particles):
        ix, iy, iz = int(p.x / cell_size), int(p.y / cell_size), int(p.z / cell_size)
        frac_x, frac_y, frac_z = (p.x / cell_size) - ix, (p.y / cell_size) - iy, (p.z / cell_size) - iz
        
        w1 = (1-frac_x) * (1-frac_y) * (1-frac_z)
        w2 = (  frac_x) * (1-frac_y) * (1-frac_z)
        w3 = (1-frac_x) * (  frac_y) * (1-frac_z)
        w4 = (1-frac_x) * (1-frac_y) * (  frac_z)
        w5 = (  frac_x) * (  frac_y) * (1-frac_z)
        w6 = (  frac_x) * (1-frac_y) * (  frac_z)
        w7 = (1-frac_x) * (  frac_y) * (  frac_z)
        w8 = (  frac_x) * (  frac_y) * (  frac_z)

        cic_data.append((ix, iy, iz, w1, w2, w3, w4, w5, w6, w7, w8))
        particle_cells.append((ix, iy, iz))
        key = (ix % mesh_size, iy % mesh_size, iz % mesh_size)
        if key not in cells: cells[key] = []
        cells[key].append(i)

        dm_rho[ix % mesh_size,     iy % mesh_size,     iz % mesh_size    ] += p.mass * w1
        dm_rho[(ix+1) % mesh_size, iy % mesh_size,     iz % mesh_size    ] += p.mass * w2
        dm_rho[ix % mesh_size,     (iy+1) % mesh_size, iz % mesh_size    ] += p.mass * w3
        dm_rho[ix % mesh_size,     iy % mesh_size,     (iz+1) % mesh_size] += p.mass * w4
        dm_rho[(ix+1) % mesh_size, (iy+1) % mesh_size, iz % mesh_size    ] += p.mass * w5
        dm_rho[(ix+1) % mesh_size, iy % mesh_size,     (iz+1) % mesh_size] += p.mass * w6
        dm_rho[ix % mesh_size,     (iy+1) % mesh_size, (iz+1) % mesh_size] += p.mass * w7
        dm_rho[(ix+1) % mesh_size, (iy+1) % mesh_size, (iz+1) % mesh_size] += p.mass * w8
    
    dm_rho /= (cell_size**3) 
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
        
        mesh_forces.append((ax * p.mass, ay * p.mass, az * p.mass))
    return mesh_forces

def compute_PM_short_range_approx(dist_sq, p1_mass, p2_mass, dx, dy, dz, config):
    soft_dist_sq = dist_sq + (0.5 * config.cell_size)**2 
    soft_dist = math.sqrt(soft_dist_sq)
    f_pm_short = config.g_const * p1_mass * p2_mass / soft_dist_sq
    return f_pm_short * dx / soft_dist, f_pm_short * dy / soft_dist, f_pm_short * dz / soft_dist

def compute_switching_parameter(dist_sq, config):
    if dist_sq < config.r_switch_start_sq: return 1.0
    x = (math.sqrt(dist_sq) - config.r_switch_start) / config.cutoff_transition_width
    return 2*x**3 - 3*x**2 + 1

def compute_PP_forces(particles, config, cells, particle_cells):
    mesh_size = config.mesh_size
    forces = [[0.0, 0.0, 0.0] for _ in particles]
    search_radius_cells = int(math.ceil(config.cutoff_radius / config.cell_size))

    for i, p1 in enumerate(particles):
        ix1, iy1, iz1 = particle_cells[i]
        for dx_cell in range(-search_radius_cells, search_radius_cells + 1):
            for dy_cell in range(-search_radius_cells, search_radius_cells + 1):
                for dz_cell in range(-search_radius_cells, search_radius_cells + 1):
                    cell_key = ((ix1 + dx_cell) % mesh_size, (iy1 + dy_cell) % mesh_size, (iz1 + dz_cell) % mesh_size)
                    if cell_key not in cells: continue

                    for j in cells[cell_key]:
                        if i >= j: continue
                        p2 = particles[j]
                        dx, dy, dz = displacement(p2.x - p1.x, config), displacement(p2.y - p1.y, config), displacement(p2.z - p1.z, config)
                        dist_sq = dx*dx + dy*dy + dz*dz

                        if dist_sq > config.cutoff_radius_squared: continue
                        
                        fx_sub, fy_sub, fz_sub = compute_PM_short_range_approx(dist_sq, p1.mass, p2.mass, dx, dy, dz, config)
                        S = compute_switching_parameter(dist_sq, config)
                        dist_sq_soft = dist_sq + config.softening_squared
                        dist = math.sqrt(dist_sq_soft)
                        f = config.g_const * p1.mass * p2.mass / dist_sq_soft
                        
                        fx, fy, fz = S * (f * dx / dist - fx_sub), S * (f * dy / dist - fy_sub), S * (f * dz / dist - fz_sub)
                        forces[i][0] += fx; forces[i][1] += fy; forces[i][2] += fz
                        forces[j][0] -= fx; forces[j][1] -= fy; forces[j][2] -= fz
    return [tuple(f) for f in forces]

def calculate_particle_energies(particles, scale_factor, config):
    kinetic_energy, potential_energy = 0.0, 0.0
    for p in particles:
        kinetic_energy += 0.5 * p.mass * ((scale_factor*p.vx)**2 + (scale_factor*p.vy)**2 + (scale_factor*p.vz)**2)

    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            p1, p2 = particles[i], particles[j]
            dx, dy, dz = displacement(p2.x - p1.x, config), displacement(p2.y - p1.y, config), displacement(p2.z - p1.z, config)
            dist = math.sqrt((scale_factor**2) * (dx*dx + dy*dy + dz*dz + config.softening_squared))
            if dist > 0: potential_energy -= config.g_const * p1.mass * p2.mass / dist
    return kinetic_energy, potential_energy

def get_gravity_timestep(particles, config):
    if not particles: return np.inf
    max_accel_sq = max([p.ax**2 + p.ay**2 + p.az**2 for p in particles] + [1e-9])
    return np.sqrt(config.softening_squared) / np.sqrt(max_accel_sq) * config.gravity_dt_factor