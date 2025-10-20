import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

import math
import time
import numpy as np
import pygame

# ------------------------
# Simulation parameters
# ------------------------
# Domain and mesh
DOMAIN_SIZE = 1
MESH_SIZE = 16
CELL_SIZE = DOMAIN_SIZE / MESH_SIZE

# P3M parameters
CUTOFF_RADIUS = 2.5 * CELL_SIZE
CUTOFF_RADIUS_SQUARED = CUTOFF_RADIUS**2
CUTOFF_TRANSITION_WIDTH = 0.2 * CUTOFF_RADIUS 
R_SWITCH_START = CUTOFF_RADIUS - CUTOFF_TRANSITION_WIDTH
R_SWITCH_START_SQ = R_SWITCH_START**2

# Particles
NUM_PARTICLES = 25 ** 2 # should be a perfect square for simple IC grid
PARTICLE_MASS = 1.0 / NUM_PARTICLES

# Physics
MEAN_INTERPARTICLE_SPACING = DOMAIN_SIZE / math.sqrt(NUM_PARTICLES)
SOFTENING_SQUARED = (MEAN_INTERPARTICLE_SPACING / 50.0)**2
EXPANSION_START_T = 0.5
EXPANDING_UNIVERSE = True
INITIAL_HUBBLE_PARAM = (2.0/3.0) / EXPANSION_START_T if EXPANDING_UNIVERSE else 0.0
INITIAL_SCALE_FACTOR = EXPANSION_START_T**(2.0/3.0) if EXPANDING_UNIVERSE else 0.0
G = 1.0 / (6 * math.pi)

# Time integration
DYNAMICAL_TIME = 1.0 / math.sqrt(G)
DT = 5e-4 * DYNAMICAL_TIME

# Visualization
RENDER_SIZE = 800
RENDER_SCALE = RENDER_SIZE / DOMAIN_SIZE

# Misc
DEBUG_INFO_EVERY_CYCLES = 40
SEED = 42
np.random.seed(SEED)

class Particle:
    def __init__(self, x, y, mass = PARTICLE_MASS, vx=0, vy=0):
        self.x = x
        self.y = y
        self.mass = mass
        self.vx = vx
        self.vy = vy
        self.ax = 0
        self.ay = 0

# minimal-image displacement helper (periodic)
def displacement(dx):
    L = DOMAIN_SIZE
    dx = (dx + 0.5*L) % L - 0.5*L
    return dx

def cic_mass_assignment(particles):
    # Cloud-In-Cell interpolation
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

def cic_force_interpolation(particles, fx_grid, fy_grid, cic_data):
    mesh_forces = []
    for i, p in enumerate(particles):
        ix, iy, w1, w2, w3, w4 = cic_data[i]

        fx = (fx_grid[ix % MESH_SIZE, iy % MESH_SIZE] * w1 +
            fx_grid[(ix + 1) % MESH_SIZE, iy % MESH_SIZE] * w2 +
            fx_grid[ix % MESH_SIZE, (iy + 1) % MESH_SIZE] * w3 +
            fx_grid[(ix + 1) % MESH_SIZE, (iy + 1) % MESH_SIZE] * w4)
            
        fy = (fy_grid[ix % MESH_SIZE, iy % MESH_SIZE] * w1 +
            fy_grid[(ix + 1) % MESH_SIZE, iy % MESH_SIZE] * w2 +
            fy_grid[ix % MESH_SIZE, (iy + 1) % MESH_SIZE] * w3 +
            fy_grid[(ix + 1) % MESH_SIZE, (iy + 1) % MESH_SIZE] * w4)
        
        fx *= p.mass
        fy *= p.mass
        mesh_forces.append((fx, fy))
    
    return mesh_forces

def compute_mesh_forces(particles):
    # Cloud-In-Cell interpolation instead of NGP (Nearest Grid Point)
    # Mass per cell
    rho, cic_data = cic_mass_assignment(particles)

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
    fx_grid = (np.roll(phi, 1, axis=0) - np.roll(phi, -1, axis=0)) / (2*CELL_SIZE)
    fy_grid = (np.roll(phi, 1, axis=1) - np.roll(phi, -1, axis=1)) / (2*CELL_SIZE)

    # Interpolate back to particles
    mesh_forces = cic_force_interpolation(particles, fx_grid, fy_grid, cic_data)
    return mesh_forces

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
    # Initialize forces as a list of mutable lists [fx, fy]
    forces = [[0.0, 0.0] for _ in particles]

    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            p1 = particles[i]
            p2 = particles[j]

            dx = displacement(p2.x - p1.x)
            dy = displacement(p2.y - p1.y)
            dist_sq = dx*dx + dy*dy

            if dist_sq > CUTOFF_RADIUS_SQUARED:
                continue
            
            # Substractive scheme
            fx_sub, fy_sub = compute_PM_short_range_approx(dist_sq, p1.mass, p2.mass, dx, dy)

            S = compute_switching_parameter(dist_sq)
            dist_sq += SOFTENING_SQUARED
            dist = math.sqrt(dist_sq)
            f = G * p1.mass * p2.mass / dist_sq
            fx = S * (f * dx / dist - fx_sub)
            fy = S * (f * dy / dist - fy_sub)
    
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

def KDK_step(total_time, particles):
    # Initial cosmology values
    scale_factor, hubble_param = update_cosmology(total_time)

    # 1. KICK (First half-step for velocity)
    for p in particles:
        total_ax = (p.ax / scale_factor**3) - (2 * hubble_param * p.vx)
        total_ay = (p.ay / scale_factor**3) - (2 * hubble_param * p.vy)
        p.vx += total_ax * DT / 2.0
        p.vy += total_ay * DT / 2.0
    
    # 2. DRIFT (Full step for position)
    for p in particles:
        p.x = ( p.x + p.vx * DT ) % DOMAIN_SIZE
        p.y = ( p.y + p.vy * DT ) % DOMAIN_SIZE

    # 3. UPDATE COSMOLOGY
    scale_factor, hubble_param = update_cosmology(total_time + DT)

    # 4. COMPUTE FORCES at new positions
    # P³M (PP+PM) Gravity Solver
    pp_forces = compute_PP_forces(particles)
    pm_forces = compute_mesh_forces(particles)
    forces = [(pp[0]+pm[0], pp[1]+pm[1]) for pp,pm in zip(pp_forces, pm_forces)]

    # 5. KICK (Second half-step for velocity)
    for i, p in enumerate(particles):
        fx, fy = forces[i]
        p.ax = fx / p.mass
        p.ay = fy / p.mass
        total_ax = (p.ax / scale_factor**3) - (2 * hubble_param * p.vx)
        total_ay = (p.ay / scale_factor**3) - (2 * hubble_param * p.vy)
        p.vx += total_ax * DT / 2.0
        p.vy += total_ay * DT / 2.0

# Create initial conditions using the simplified Zel'dovich Approximation
def create_zeldovich_ics(power_spectrum_index):
    n_per_side = int(round(NUM_PARTICLES ** 0.5))
    if n_per_side**2 != NUM_PARTICLES:
        print(f"Warning: {NUM_PARTICLES} is not a perfect square. Using {n_per_side**2} particles.")

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
            particles.append(Particle(x, y, PARTICLE_MASS, vx, vy))

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

def render(screen, particles):
    pixels = np.zeros((RENDER_SIZE, RENDER_SIZE, 3), dtype=np.uint8)
    coords = np.array([(int(p.x*RENDER_SCALE), int(p.y*RENDER_SCALE)) for p in particles])
    pixels[coords[:,0], coords[:,1]] = [255, 255, 255]
    pygame.surfarray.blit_array(screen, pixels)    
    pygame.display.flip()

# ------------------------
# Main entry point
# ------------------------
def main():
    screen = init_window(RENDER_SIZE, 'N-Body Simulation')

    particles = create_zeldovich_ics(power_spectrum_index=-2.0)
    initial_energy = calculate_total_energy(particles)

    cycle_count = 0
    simulation_start_time = time.time()
    total_simulation_time = 0.0

    while process_events():
        total_simulation_time = DT * cycle_count
        KDK_step(total_simulation_time, particles)
        render(screen, particles)

        cycle_count += 1

        if cycle_count % DEBUG_INFO_EVERY_CYCLES == 0:
            total_elapsed_time = time.time() - simulation_start_time
            average_time = total_elapsed_time / cycle_count
            scale_factor, _ = update_cosmology(total_simulation_time)
            current_energy = calculate_total_energy(particles, scale_factor)
            energy_error = abs(current_energy - initial_energy) / abs(initial_energy)
            print(f"Time:{total_elapsed_time:.0f}s Cycles:{cycle_count} Avg. dt:{average_time:.3f}s Energy Error:{energy_error:.1%}")

    pygame.quit()

if __name__ == "__main__":
    main()