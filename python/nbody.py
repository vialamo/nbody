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

# ------------------------
# Simulation parameters
# ------------------------
# Domain and mesh
DOMAIN_SIZE = 1.0
MESH_SIZE = 64
CELL_SIZE = DOMAIN_SIZE / MESH_SIZE

# Matter components
OMEGA_BARYON = 0.15 # Baryon fraction
OMEGA_DM = 1.0 - OMEGA_BARYON

# P3M parameters
CUTOFF_RADIUS = 2.5 * CELL_SIZE
CUTOFF_RADIUS_SQUARED = CUTOFF_RADIUS**2
CUTOFF_TRANSITION_WIDTH = 0.2 * CUTOFF_RADIUS 
R_SWITCH_START = CUTOFF_RADIUS - CUTOFF_TRANSITION_WIDTH
R_SWITCH_START_SQ = R_SWITCH_START**2

# Particles and gas
NUM_DM_PARTICLES = 25 ** 2 # should be a perfect square for simple IC grid
DM_PARTICLE_MASS = OMEGA_DM / NUM_DM_PARTICLES
GAS_TOTAL_MASS = OMEGA_BARYON
ENABLE_HYDRO = True

# Physics
MEAN_INTERPARTICLE_SPACING = DOMAIN_SIZE / math.sqrt(NUM_DM_PARTICLES)
SOFTENING_SQUARED = (MEAN_INTERPARTICLE_SPACING / 50.0)**2
EXPANSION_START_T = 0.5
EXPANDING_UNIVERSE = True
INITIAL_HUBBLE_PARAM = (2.0/3.0) / EXPANSION_START_T if EXPANDING_UNIVERSE else 0.0
INITIAL_SCALE_FACTOR = EXPANSION_START_T**(2.0/3.0) if EXPANDING_UNIVERSE else 0.0
G = 1.0 / (6 * math.pi)
GAMMA = 5.0/3.0 # Adiabatic index for monatomic gas

# Time integration
DYNAMICAL_TIME = 1.0 / math.sqrt(G)
DT = 2e-4 * DYNAMICAL_TIME

# Visualization
RENDER_SIZE = 800
RENDER_SCALE = RENDER_SIZE / DOMAIN_SIZE

# Misc
DEBUG_INFO_EVERY_CYCLES = 40
SAVE_SNAPSHOT_EVERY_CYCLES = 100
SEED = 42
np.random.seed(SEED)

class Particle:
    def __init__(self, x, y, mass = DM_PARTICLE_MASS, vx=0, vy=0):
        self.x = x
        self.y = y
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.ax = 0
        self.ay = 0

class GasGrid:
    """Represents the Eulerian grid for the baryonic gas."""
    def __init__(self):
        # Conservative variables
        self.density = np.zeros((MESH_SIZE, MESH_SIZE))
        self.momentum_x = np.zeros((MESH_SIZE, MESH_SIZE))
        self.momentum_y = np.zeros((MESH_SIZE, MESH_SIZE))
        self.energy = np.zeros((MESH_SIZE, MESH_SIZE)) # Total energy density E = rho*e + 0.5*rho*v^2
        
        # Primitive variables (derived)
        self.pressure = np.zeros((MESH_SIZE, MESH_SIZE))
        self.velocity_x = np.zeros((MESH_SIZE, MESH_SIZE))
        self.velocity_y = np.zeros((MESH_SIZE, MESH_SIZE))

        self.accel_x = np.zeros((MESH_SIZE, MESH_SIZE))
        self.accel_y = np.zeros((MESH_SIZE, MESH_SIZE))

    def update_primitive_variables(self):
        """Calculates pressure and velocity from the conservative variables."""
        # Avoid division by zero in empty cells
        non_zero_density = self.density > 1e-9
        
        self.velocity_x[non_zero_density] = self.momentum_x[non_zero_density] / self.density[non_zero_density]
        self.velocity_y[non_zero_density] = self.momentum_y[non_zero_density] / self.density[non_zero_density]

        kinetic_energy_density = 0.5 * (self.momentum_x**2 + self.momentum_y**2) / self.density
        internal_energy_density = self.energy - kinetic_energy_density
        
        self.pressure = (GAMMA - 1.0) * internal_energy_density
        self.pressure[self.pressure < 0] = 1e-9 # Pressure floor

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
        cs_L = np.sqrt(GAMMA * p_L / rho_L, where=rho_L > 0, out=np.zeros_like(rho_L))
        cs_R = np.sqrt(GAMMA * p_R / rho_R, where=rho_R > 0, out=np.zeros_like(rho_R))

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

        # --- X-direction sweep ---
        flux_density_x, flux_mom_x_x, flux_mom_y_x, flux_energy_x = self.calculate_fluxes(axis=0)
        factor = dt / CELL_SIZE
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

# minimal-image displacement helper (periodic)
def displacement(dx):
    L = DOMAIN_SIZE
    dx = (dx + 0.5*L) % L - 0.5*L
    return dx

def cic_mass_assignment(particles):
    # Mass per cell
    rho = np.zeros((MESH_SIZE, MESH_SIZE))
    cic_data = []
    for p in particles:
        ix = int(p.x / CELL_SIZE)
        iy = int(p.y / CELL_SIZE)

        # Find the particle's fractional position within that cell (0 to 1)
        frac_x = (p.x / CELL_SIZE) - ix
        frac_y = (p.y / CELL_SIZE) - iy

        # Calculate weights for the 4 corners
        w1 = (1 - frac_x) * (1 - frac_y) # (ix, iy)
        w2 = frac_x * (1 - frac_y)       # (ix+1, iy)
        w3 = (1 - frac_x) * frac_y       # (ix, iy+1)
        w4 = frac_x * frac_y             # (ix+1, iy+1)

        cic_data.append((ix, iy, w1, w2, w3, w4))

        # Distribute mass to the 4 grid points (with periodic boundaries)
        rho[ix % MESH_SIZE, iy % MESH_SIZE] += p.mass * w1
        rho[(ix + 1) % MESH_SIZE, iy % MESH_SIZE] += p.mass * w2
        rho[ix % MESH_SIZE, (iy + 1) % MESH_SIZE] += p.mass * w3
        rho[(ix + 1) % MESH_SIZE, (iy + 1) % MESH_SIZE] += p.mass * w4
    
    rho /= (CELL_SIZE * CELL_SIZE) # Convert mass to density
    return rho, cic_data

def cic_force_interpolation(particles, ax_grid, ay_grid, cic_data):
    mesh_forces = []
    for i, p in enumerate(particles):
        ix, iy, w1, w2, w3, w4 = cic_data[i]

        fx = (ax_grid[ix % MESH_SIZE, iy % MESH_SIZE] * w1 +
            ax_grid[(ix + 1) % MESH_SIZE, iy % MESH_SIZE] * w2 +
            ax_grid[ix % MESH_SIZE, (iy + 1) % MESH_SIZE] * w3 +
            ax_grid[(ix + 1) % MESH_SIZE, (iy + 1) % MESH_SIZE] * w4)
            
        fy = (ay_grid[ix % MESH_SIZE, iy % MESH_SIZE] * w1 +
            ay_grid[(ix + 1) % MESH_SIZE, iy % MESH_SIZE] * w2 +
            ay_grid[ix % MESH_SIZE, (iy + 1) % MESH_SIZE] * w3 +
            ay_grid[(ix + 1) % MESH_SIZE, (iy + 1) % MESH_SIZE] * w4)
        
        fx *= p.mass
        fy *= p.mass
        mesh_forces.append((fx, fy))
    
    return mesh_forces

def compute_gravitational_acceleration(particles, gas):
    # Cloud-In-Cell interpolation instead of NGP (Nearest Grid Point)
    # Mass per cell
    rho, cic_data = cic_mass_assignment(particles)
    if ENABLE_HYDRO:
        rho += gas.density

    # FFT of density
    rho_k = np.fft.fftn(rho)

    # Poisson solver in k-space
    kx = np.fft.fftfreq(MESH_SIZE, d=CELL_SIZE) * 2 * np.pi
    ky = np.fft.fftfreq(MESH_SIZE, d=CELL_SIZE) * 2 * np.pi
    kx, ky = np.meshgrid(kx, ky, indexing='ij')
    k = np.sqrt(kx*kx + ky*ky)
    k[0,0] = 1.0  # avoid div by zero, k[0,0] represents a wave with zero frequency

    phi_k = -2 * np.pi * G * rho_k / k
    phi_k[0,0] = 0.0 # set the average value of the potential field to zero

    # Inverse FFT the POTENTIAL back to the real-space grid
    phi = np.fft.ifftn(phi_k).real

    # Gradient to get field (finite differences)
    ax_grid = (np.roll(phi, 1, axis=0) - np.roll(phi, -1, axis=0)) / (2*CELL_SIZE)
    ay_grid = (np.roll(phi, 1, axis=1) - np.roll(phi, -1, axis=1)) / (2*CELL_SIZE)

    return ax_grid, ay_grid, cic_data

def compute_PM_short_range_approx(dist_sq, p1_mass, p2_mass, dx, dy):
    softening_epsilon_squared = (0.5 * CELL_SIZE)**2 
    soft_dist_sq = dist_sq + softening_epsilon_squared
    soft_dist = math.sqrt(soft_dist_sq)
    f_pm_short = G * p1_mass * p2_mass / soft_dist_sq
    f_pm_short_x = f_pm_short * dx / soft_dist
    f_pm_short_y = f_pm_short * dy / soft_dist
    return f_pm_short_x, f_pm_short_y

def compute_switching_parameter(dist_sq):
    if dist_sq < R_SWITCH_START_SQ:
        S = 1.0
    else:
        # In the transition zone, smoothly fade the correction to zero
        # Normalize distance within the transition zone (x goes from 0 to 1)
        dist = math.sqrt(dist_sq)
        x = (dist - R_SWITCH_START) / CUTOFF_TRANSITION_WIDTH
        # Use a smoothstep polynomial (2x³ - 3x² + 1) for a smooth transition
        S = 2*x**3 - 3*x**2 + 1
    return S

def compute_PP_forces(particles):
    # Bin particles into a grid (cell list)
    cells = {} # Key: (ix, iy), Value: list of particle indices
    particle_cells = [] # Store (ix, iy) for each particle i
    for i, p in enumerate(particles):
        ix = int(p.x / CELL_SIZE)
        iy = int(p.y / CELL_SIZE)
        particle_cells.append((ix, iy))
        key = (ix, iy)
        if key not in cells:
            cells[key] = []
        cells[key].append(i)

    forces = [[0.0, 0.0] for _ in particles]
    
    # Base search radius in cells on cutoff
    search_radius_cells = int(math.ceil(CUTOFF_RADIUS / CELL_SIZE))

    # Iterate over particles and their neighbors
    for i, p1 in enumerate(particles):
        ix1, iy1 = particle_cells[i]
        
        # Check blocks of cells around the particle
        for dx_cell in range(-search_radius_cells, search_radius_cells + 1):
            for dy_cell in range(-search_radius_cells, search_radius_cells + 1):
                
                neighbor_ix = (ix1 + dx_cell) % MESH_SIZE
                neighbor_iy = (iy1 + dy_cell) % MESH_SIZE
                cell_key = (neighbor_ix, neighbor_iy)
                
                if cell_key not in cells:
                    continue # Empty neighbor cells

                # Compute interactions with particles in the neighbor cell
                for j in cells[cell_key]:
                    # Avoid double counting (j > i) and self-interaction (j != i)
                    if i >= j:
                        continue
                            
                    p2 = particles[j]
                    
                    dx = displacement(p2.x - p1.x)
                    dy = displacement(p2.y - p1.y)
                    dist_sq = dx*dx + dy*dy

                    if dist_sq > CUTOFF_RADIUS_SQUARED:
                        continue
                    
                    # Subtractive scheme
                    fx_sub, fy_sub = compute_PM_short_range_approx(dist_sq, p1.mass, p2.mass, dx, dy)
                    S = compute_switching_parameter(dist_sq)
                    dist_sq_soft = dist_sq + SOFTENING_SQUARED
                    dist = math.sqrt(dist_sq_soft)
                    f = G * p1.mass * p2.mass / dist_sq_soft
                    fx = S * (f * dx / dist - fx_sub)
                    fy = S * (f * dy / dist - fy_sub)
            
                    # Apply forces to both particles
                    forces[i][0] += fx
                    forces[i][1] += fy

                    forces[j][0] -= fx
                    forces[j][1] -= fy
            
    # Convert back to a list of tuples for compatibility
    return [tuple(f) for f in forces]

def update_cosmology(total_time):
    if EXPANDING_UNIVERSE:
        expansion_time = EXPANSION_START_T + total_time
        scale_factor = expansion_time**(2.0/3.0)
        hubble_param = (2.0/3.0) / expansion_time
    else:
        scale_factor = 1.0
        hubble_param = 0.0
    return scale_factor, hubble_param

def KDK_step(total_time, particles, gas):
    scale_factor, hubble_param = update_cosmology(total_time)

    # 1. KICK (First half-step for velocity)
    for p in particles:
        total_ax = (p.ax / scale_factor**3) - (2 * hubble_param * p.vx)
        total_ay = (p.ay / scale_factor**3) - (2 * hubble_param * p.vy)
        p.vx += total_ax * DT / 2.0
        p.vy += total_ay * DT / 2.0

    # Apply Kick to Gas Grid
    if ENABLE_HYDRO:
        gas.update_primitive_variables()
        gas_accel_x = (gas.accel_x / scale_factor**3) - (2 * hubble_param * gas.velocity_x)
        gas_accel_y = (gas.accel_y / scale_factor**3) - (2 * hubble_param * gas.velocity_y)
        gas.momentum_x += gas.density * gas_accel_x * DT / 2.0
        gas.momentum_y += gas.density * gas_accel_y * DT / 2.0
        gas.energy += (gas.momentum_x * gas_accel_x + gas.momentum_y * gas_accel_y) * DT / 2.0
    
    # 2. DRIFT (Full step for position)
    for p in particles:
        p.x = ( p.x + p.vx * DT ) % DOMAIN_SIZE
        p.y = ( p.y + p.vy * DT ) % DOMAIN_SIZE

    # 3. UPDATE COSMOLOGY
    scale_factor, hubble_param = update_cosmology(total_time + DT)

    # 4. COMPUTE FORCES at new positions
    # P³M (PP+PM) Gravity Solver
    pp_forces = compute_PP_forces(particles)
    ax_grid, ay_grid, cic_data = compute_gravitational_acceleration(particles, gas)
    pm_forces = cic_force_interpolation(particles, ax_grid, ay_grid, cic_data)
    forces = [(pp[0]+pm[0], pp[1]+pm[1]) for pp,pm in zip(pp_forces, pm_forces)]
    if ENABLE_HYDRO:
        gas.hydro_step(DT)

    # 5. KICK (Second half-step for velocity)
    for i, p in enumerate(particles):
        fx, fy = forces[i]
        p.ax = fx / p.mass
        p.ay = fy / p.mass
        total_ax = (p.ax / scale_factor**3) - (2 * hubble_param * p.vx)
        total_ay = (p.ay / scale_factor**3) - (2 * hubble_param * p.vy)
        p.vx += total_ax * DT / 2.0
        p.vy += total_ay * DT / 2.0

    # Apply Kick to Gas Grid
    if ENABLE_HYDRO:
        gas.update_primitive_variables()
        gas.accel_x = ax_grid
        gas.accel_y = ay_grid
        gas_accel_x = (gas.accel_x / scale_factor**3) - (2 * hubble_param * gas.velocity_x)
        gas_accel_y = (gas.accel_y / scale_factor**3) - (2 * hubble_param * gas.velocity_y)
        gas.momentum_x += gas.density * gas_accel_x * DT / 2.0
        gas.momentum_y += gas.density * gas_accel_y * DT / 2.0
        gas.energy += (gas.momentum_x * gas_accel_x + gas.momentum_y * gas_accel_y) * DT / 2.0

# Create initial conditions using the simplified Zel'dovich Approximation
def create_zeldovich_ics(power_spectrum_index):
    n_per_side = int(round(NUM_DM_PARTICLES ** 0.5))
    if n_per_side**2 != NUM_DM_PARTICLES:
        print(f"Warning: {NUM_DM_PARTICLES} is not a perfect square. Using {n_per_side**2} particles.")

    ic_mesh_size = n_per_side
    cell_size = DOMAIN_SIZE / ic_mesh_size

    # Build k-grid
    kx = np.fft.fftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi
    ky = np.fft.fftfreq(ic_mesh_size, d=cell_size) * 2 * np.pi
    kx_grid, ky_grid = np.meshgrid(kx, ky, indexing='ij')
    k2 = kx_grid**2 + ky_grid**2
    k2[0, 0] = 1.0  # avoid division by zero

    # Generate random density field δ(k) with desired P(k) ∝ k^n
    real_space_random_field = np.random.randn(ic_mesh_size, ic_mesh_size)
    random_k = np.fft.fft2(real_space_random_field)
    delta_k = random_k * (k2 ** (power_spectrum_index / 4.0))

    # Compute potential Φ(k) = -δ(k)/k²
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
    disp_x *= INITIAL_SCALE_FACTOR * cell_size
    disp_y *= INITIAL_SCALE_FACTOR * cell_size

    # Create displaced particles
    particles = []
    spacing = DOMAIN_SIZE / n_per_side
    for i in range(n_per_side):
        for j in range(n_per_side):
            qx = (i + 0.5) * spacing
            qy = (j + 0.5) * spacing
            x = (qx + disp_x[i, j]) % DOMAIN_SIZE
            y = (qy + disp_y[i, j]) % DOMAIN_SIZE
            vx = INITIAL_HUBBLE_PARAM * disp_x[i, j]
            vy = INITIAL_HUBBLE_PARAM * disp_y[i, j]
            particles.append(Particle(x, y, DM_PARTICLE_MASS, vx, vy))

    return particles

def calculate_total_energy(particles, scale_factor=1.0):
    kinetic_energy = 0.0
    potential_energy = 0.0

    for p in particles:
        vel_sq = (scale_factor*p.vx)**2 + (scale_factor*p.vy)**2
        kinetic_energy += 0.5 * p.mass * vel_sq

    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            p1 = particles[i]
            p2 = particles[j]

            dx = displacement(p2.x - p1.x)
            dy = displacement(p2.y - p1.y)
            
            dist_sq = (scale_factor**2) * (dx * dx + dy * dy + SOFTENING_SQUARED)
            dist = math.sqrt(dist_sq)
            potential_energy -= G * p1.mass * p2.mass / dist

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

def render(screen, particles, gas):
    surface_array = np.zeros((RENDER_SIZE, RENDER_SIZE, 3), dtype=np.uint8)

    if ENABLE_HYDRO:
        # Render Gas Pressure
        # Normalize pressure for visualization using a logarithmic scale
        log_field = np.log10(gas.pressure + 1e-12) # Add a small epsilon to avoid log(0)
        min_log, max_log = np.min(log_field), np.max(log_field)
        
        if max_log > min_log:
            norm_field = (log_field - min_log) / (max_log - min_log)
        else:
            norm_field = np.zeros_like(log_field)

        # Apply a colormap
        colormap = matplotlib.cm.plasma
        colored_field = colormap(norm_field)
        colored_field[:, :, 0:3] *= norm_field[..., np.newaxis]
        
        # Upscale to render size using bilinear interpolation
        zoom_factor = RENDER_SIZE / MESH_SIZE
        zoomed_field = scipy.ndimage.zoom(colored_field[:, :, 0:3], 
                                        (zoom_factor, zoom_factor, 1), 
                                        order=1)

        surface_array = (zoomed_field * 255).astype(np.uint8)
        pygame.surfarray.blit_array(screen, surface_array)

    coords = np.array([(int(p.x*RENDER_SCALE), int(p.y*RENDER_SCALE)) for p in particles])
    surface_array[coords[:,0], coords[:,1]] = [255, 255, 255]
    pygame.surfarray.blit_array(screen, surface_array)
    pygame.display.flip()

def save_snapshot(h5_file, snapshot_name, particles, gas, sim_time, scale_factor):
    """Saves the current simulation state to a group in the HDF5 file."""
    try:
        grp = h5_file.create_group(snapshot_name)
        
        # Store Metadata
        grp.attrs['time'] = sim_time
        grp.attrs['scale_factor'] = scale_factor
        
        # Store Particle Data
        # Convert list of objects to NumPy arrays for efficient storage
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
        if ENABLE_HYDRO:
            g_grp = grp.create_group("gas")
            g_grp.create_dataset("density", data=gas.density, compression="gzip")
            g_grp.create_dataset("momentum_x", data=gas.momentum_x, compression="gzip")
            g_grp.create_dataset("momentum_y", data=gas.momentum_y, compression="gzip")
            g_grp.create_dataset("energy", data=gas.energy, compression="gzip")
            g_grp.create_dataset("pressure", data=gas.pressure, compression="gzip") # Also save primitive
            g_grp.create_dataset("velocity_x", data=gas.velocity_x, compression="gzip")
            g_grp.create_dataset("velocity_y", data=gas.velocity_y, compression="gzip")
            
    except Exception as e:
        print(f"Error saving snapshot {snapshot_name}: {e}")

# ------------------------
# Main entry point
# ------------------------
def main():
    screen = init_window(RENDER_SIZE, 'N-Body + Hydrodynamics Simulation')

    particles = create_zeldovich_ics(power_spectrum_index=-2.0)
    
    # Uniform initialization for gas
    gas = GasGrid()
    gas.density.fill(GAS_TOTAL_MASS / (DOMAIN_SIZE**2))
    initial_internal_energy = 1e-6
    gas.energy = initial_internal_energy * gas.density
    
    if ENABLE_HYDRO:
        ax_grid, ay_grid, _ = compute_gravitational_acceleration(particles, gas)
        gas.accel_x = ax_grid
        gas.accel_y = ay_grid

    initial_energy = calculate_total_energy(particles)

    cycle_count = 0
    simulation_start_time = time.time()
    total_simulation_time = 0.0

    # Create HDF5 file
    run_timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    snapshot_filename = f"sim_{run_timestamp}.hdf5"
    print(f"Saving snapshots to: {snapshot_filename}")

    h5_file = h5py.File(snapshot_filename, 'w')
    param_grp = h5_file.create_group("parameters")
    param_grp.attrs['DOMAIN_SIZE'] = DOMAIN_SIZE
    param_grp.attrs['MESH_SIZE'] = MESH_SIZE
    param_grp.attrs['NUM_DM_PARTICLES'] = NUM_DM_PARTICLES
    param_grp.attrs['OMEGA_BARYON'] = OMEGA_BARYON
    param_grp.attrs['OMEGA_DM'] = OMEGA_DM
    param_grp.attrs['ENABLE_HYDRO'] = ENABLE_HYDRO
    param_grp.attrs['G'] = G
    param_grp.attrs['GAMMA'] = GAMMA
    param_grp.attrs['DT'] = DT
    param_grp.attrs['EXPANDING_UNIVERSE'] = EXPANDING_UNIVERSE
    param_grp.attrs['INITIAL_HUBBLE_PARAM'] = INITIAL_HUBBLE_PARAM
    param_grp.attrs['SOFTENING_SQUARED'] = SOFTENING_SQUARED
    snapshot_count = 0

    while process_events():
        total_simulation_time = DT * cycle_count
        KDK_step(total_simulation_time, particles, gas)
        render(screen, particles, gas)

        cycle_count += 1

        if SAVE_SNAPSHOT_EVERY_CYCLES > 0 and cycle_count % SAVE_SNAPSHOT_EVERY_CYCLES == 0:
            snapshot_name = f"snapshot_{snapshot_count:04d}"
            scale_factor, _ = update_cosmology(total_simulation_time)
            save_snapshot(h5_file, snapshot_name, particles, gas, total_simulation_time, scale_factor)
            print(f"Savied snapshot: {snapshot_name}")
            snapshot_count += 1

        if cycle_count % DEBUG_INFO_EVERY_CYCLES == 0:
            total_elapsed_time = time.time() - simulation_start_time
            average_time = total_elapsed_time / cycle_count
            scale_factor, _ = update_cosmology(total_simulation_time)
            current_energy = calculate_total_energy(particles, scale_factor)
            energy_error = abs(current_energy - initial_energy) / abs(initial_energy)
            print(f"Time:{total_elapsed_time:.0f}s Cycles:{cycle_count} Avg. dt:{average_time:.3f}s Scale:{scale_factor:.2f} Energy Error:{energy_error:.1%}")

    h5_file.close()
    pygame.quit()

if __name__ == "__main__":
    main()