import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from numpy.linalg import eigh

class MFMFluidSandbox:
    def __init__(self, num_particles=32, domain_size=1.0):
        self.N = num_particles
        self.L = domain_size
        self.gamma = 5.0 / 3.0
        self.target_neighbors = 12.0
        
        # Geometry
        grid_1d = int(np.ceil(np.sqrt(self.N)))
        x, y = np.meshgrid(np.linspace(0.2, 0.8, grid_1d), np.linspace(0.2, 0.8, grid_1d))
        self.pos = np.vstack([x.ravel(), y.ravel()]).T[:self.N]
        self.pos += np.random.normal(0, 0.01, self.pos.shape)
        
        # Fluid States
        self.mass = np.ones(self.N) * (1.0 / self.N)
        self.h = np.ones(self.N) * 0.15
        self.rho = np.zeros(self.N)
        self.vel = np.zeros((self.N, 2))
        
        # Setup a Sedov-like high energy center
        self.u = np.ones(self.N) * 0.1
        center = np.array([self.L/2, self.L/2])
        dists = np.linalg.norm(self.pos - center, axis=1)
        self.u[dists < 0.15] = 5.0 
        self.P = (self.gamma - 1.0) * self.rho * self.u
        
        # Geometry matrices & Gradients
        self.B = np.zeros((self.N, 2, 2))
        self.grad_u = np.zeros((self.N, 2))
        self.grad_v = np.zeros((self.N, 2, 2)) # Shape: (N, 2_vel_components, 2_spatial_dims)
        
        self.selected_idx = self.N // 2
        self.dragging_idx = None
        
        self.compute_density_and_h()
        self.compute_geometry_and_gradients()
        self.setup_plot()

    def cubic_spline_2d(self, r, h):
        q = r / h
        sigma = 10.0 / (7.0 * np.pi * h**2)
        if q <= 1.0:
            return sigma * (1.0 - 1.5 * q**2 + 0.75 * q**3)
        elif q <= 2.0:
            return sigma * 0.25 * (2.0 - q)**3
        return 0.0

    def compute_b_matrix_tikhonov(self, E):
        eigenvalues, eigenvectors = eigh(E)
        max_eigenvalue = np.max(eigenvalues)
        
        # Hard floor + relative threshold
        epsilon = max(1e-12, 1e-6 * max_eigenvalue)
        eps_sq = epsilon**2
        
        # Smooth Tikhonov regularization
        inv_eigenvalues = eigenvalues / (eigenvalues**2 + eps_sq)
        B = eigenvectors @ np.diag(inv_eigenvalues) @ eigenvectors.T
        return B

    def compute_density_and_h(self):
        # Iteratively solve for h to maintain target neighbors, and compute density
        for _ in range(4): # 4 Newton-Raphson-ish iterations
            for i in range(self.N):
                p_i = self.pos[i]
                rho_i = 0.0
                N_neighbors = 0.0
                for j in range(self.N):
                    r = np.linalg.norm(self.pos[j] - p_i)
                    if r < 2.0 * self.h[i]:
                        W = self.cubic_spline_2d(r, self.h[i])
                        rho_i += self.mass[j] * W
                        # Approximate neighbor count enclosed
                        N_neighbors += W * (np.pi * (2*self.h[i])**2)
                
                self.rho[i] = rho_i
                
                # Relax h towards target
                ratio = self.target_neighbors / max(N_neighbors, 1e-4)
                # Dampen the change to prevent oscillation
                self.h[i] *= ratio ** 0.25 
                
        self.P = (self.gamma - 1.0) * self.rho * self.u

    def compute_geometry_and_gradients(self):
        for i in range(self.N):
            p_i = self.pos[i]
            E = np.zeros((2, 2))
            V_i = self.mass[i] / self.rho[i]
            
            # Build Shape Matrix E
            for j in range(self.N):
                if i == j: continue
                dx = self.pos[j] - p_i
                r = np.linalg.norm(dx)
                if r < 2.0 * self.h[i]:
                    W = self.cubic_spline_2d(r, self.h[i])
                    V_j = self.mass[j] / self.rho[j]
                    weight = W * V_j
                    E += weight * np.outer(dx, dx)
                    
            self.B[i] = self.compute_b_matrix_tikhonov(E)
            
            # Reconstruct Gradients using B
            g_u = np.zeros(2)
            g_vx, g_vy = np.zeros(2), np.zeros(2)
            
            for j in range(self.N):
                if i == j: continue
                dx = self.pos[j] - p_i
                r = np.linalg.norm(dx)
                if r < 2.0 * self.h[i]:
                    W = self.cubic_spline_2d(r, self.h[i])
                    V_j = self.mass[j] / self.rho[j]
                    # Base weight vector: W_ij * B_i * dx_ij
                    w_vec = V_j * W * (self.B[i] @ dx)
                    
                    g_u += (self.u[j] - self.u[i]) * w_vec
                    g_vx += (self.vel[j, 0] - self.vel[i, 0]) * w_vec
                    g_vy += (self.vel[j, 1] - self.vel[i, 1]) * w_vec
                    
            self.grad_u[i] = g_u
            self.grad_v[i, 0] = g_vx
            self.grad_v[i, 1] = g_vy

    def advance_cycle(self, event):
        self.compute_density_and_h()
        self.compute_geometry_and_gradients()
        
        dt = 0.005 # Fixed for sandbox stability
        du_dt = np.zeros(self.N)
        dv_dt = np.zeros((self.N, 2))
        
        # Loop over all pairs to compute fluxes
        for i in range(self.N):
            for j in range(i+1, self.N):
                dx = self.pos[j] - self.pos[i]
                r = np.linalg.norm(dx)
                
                # Check interaction radius (symmetric max)
                if r > 2.0 * max(self.h[i], self.h[j]):
                    continue
                    
                # 1. Effective Face Area Vector (A_ij)
                W_i = self.cubic_spline_2d(r, self.h[i])
                W_j = self.cubic_spline_2d(r, self.h[j])
                W_avg = 0.5 * (W_i + W_j)
                
                V_i = self.mass[i] / self.rho[i]
                V_j = self.mass[j] / self.rho[j]
                
                B_avg = 0.5 * (self.B[i] + self.B[j])
                A_ij = V_i * V_j * W_avg * (B_avg @ dx)
                area_mag = np.linalg.norm(A_ij)
                if area_mag < 1e-12:
                    continue
                n_ij = A_ij / area_mag
                
                # 2. Reconstruct states at the face (Second-order)
                # Extrapolate halfway
                dx_face_i = 0.5 * dx
                dx_face_j = -0.5 * dx
                
                u_L = self.u[i] + np.dot(self.grad_u[i], dx_face_i)
                u_R = self.u[j] + np.dot(self.grad_u[j], dx_face_j)
                
                v_L = self.vel[i] + self.grad_v[i] @ dx_face_i
                v_R = self.vel[j] + self.grad_v[j] @ dx_face_j
                
                # Limit internal energy to physical bounds
                u_L, u_R = max(u_L, 0.01), max(u_R, 0.01)
                
                P_L = (self.gamma - 1.0) * self.rho[i] * u_L
                P_R = (self.gamma - 1.0) * self.rho[j] * u_R
                
                c_L = np.sqrt(self.gamma * P_L / self.rho[i])
                c_R = np.sqrt(self.gamma * P_R / self.rho[j])
                c_max = max(c_L, c_R)
                rho_avg = 0.5 * (self.rho[i] + self.rho[j])
                
                # 3. Acoustic Riemann Solver (Rusanov/Godunov approximation)
                v_n_L = np.dot(v_L, n_ij)
                v_n_R = np.dot(v_R, n_ij)
                
                # Interface Pressure and Normal Velocity
                P_star = 0.5 * (P_L + P_R) + 0.5 * c_max * rho_avg * (v_n_L - v_n_R)
                P_star = max(P_star, 0.0)
                
                v_star_n = 0.5 * (v_n_L + v_n_R) + 0.5 * (P_L - P_R) / (rho_avg * c_max + 1e-8)
                v_star = v_star_n * n_ij # Face velocity vector
                
                # 4. Lagrangian Fluxes
                # dv/dt = - P * A
                momentum_flux = P_star * A_ij
                dv_dt[i] -= momentum_flux / self.mass[i]
                dv_dt[j] += momentum_flux / self.mass[j]
                
                # du/dt = - P * v * A
                energy_flux = P_star * np.dot(v_star, A_ij)
                du_dt[i] -= energy_flux / self.mass[i]
                du_dt[j] += energy_flux / self.mass[j]

        # Update States
        self.vel += dv_dt * dt
        self.u += du_dt * dt
        self.pos += self.vel * dt
        
        # Enforce reflective boundaries for sandbox
        for d in range(2):
            mask_low = self.pos[:, d] < 0.05
            self.pos[mask_low, d] = 0.05
            self.vel[mask_low, d] *= -0.5
            
            mask_high = self.pos[:, d] > self.L - 0.05
            self.pos[mask_high, d] = self.L - 0.05
            self.vel[mask_high, d] *= -0.5

        self.update_physics_and_plot()

    def update_physics_and_plot(self):
        self.ax.clear()
        self.ax.set_xlim(0, self.L)
        self.ax.set_ylim(0, self.L)
        self.ax.set_title("MFM Sandbox: Pressure & Face Vectors")
        
        # Plot particles colored by Pressure
        sc = self.ax.scatter(self.pos[:, 0], self.pos[:, 1], c=self.P, cmap='jet', s=60, zorder=2, vmin=0)
        
        # Highlight focal particle
        p_i = self.pos[self.selected_idx]
        self.ax.scatter(p_i[0], p_i[1], c='none', s=150, zorder=3, edgecolors='black', linewidths=2)
        
        # Draw dynamic smoothing circle for focal particle
        h_i = self.h[self.selected_idx]
        circle = plt.Circle((p_i[0], p_i[1]), h_i * 2.0, color='black', fill=False, linestyle='--', alpha=0.5)
        self.ax.add_patch(circle)

        # Draw Face Vectors for focal particle
        V_i = self.mass[self.selected_idx] / self.rho[self.selected_idx]
        origins_x, origins_y, vectors_u, vectors_v = [], [], [], []

        for j in range(self.N):
            if self.selected_idx == j: continue
            dx = self.pos[j] - p_i
            r = np.linalg.norm(dx)
            
            if r < 2.0 * max(h_i, self.h[j]):
                W_avg = 0.5 * (self.cubic_spline_2d(r, h_i) + self.cubic_spline_2d(r, self.h[j]))
                V_j = self.mass[j] / self.rho[j]
                B_avg = 0.5 * (self.B[self.selected_idx] + self.B[j])
                
                A_ij = V_i * V_j * W_avg * (B_avg @ dx)
                
                midpoint = p_i + 0.5 * dx
                origins_x.append(midpoint[0])
                origins_y.append(midpoint[1])
                # Scale up arrow for visibility
                vectors_u.append(A_ij[0] * 5.0) 
                vectors_v.append(A_ij[1] * 5.0)

        if origins_x:
            self.ax.quiver(origins_x, origins_y, vectors_u, vectors_v, 
                           color='green', angles='xy', scale_units='xy', scale=1.0, zorder=4)

        # HUD Info
        hud_text = (f"Focal Particle {self.selected_idx}\n"
                    f"Rho: {self.rho[self.selected_idx]:.2f}\n"
                    f"h: {self.h[self.selected_idx]:.3f}\n"
                    f"P: {self.P[self.selected_idx]:.2f}")
        self.ax.text(0.02, 0.98, hud_text, transform=self.ax.transAxes, 
                     verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        self.fig.canvas.draw_idle()

    def setup_plot(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        plt.subplots_adjust(bottom=0.15)
        
        # Advance Button
        ax_btn = plt.axes([0.4, 0.05, 0.2, 0.05])
        self.btn = Button(ax_btn, 'Advance Cycle')
        self.btn.on_clicked(self.advance_cycle)
        
        self.update_physics_and_plot()
        
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        plt.show()

    def on_press(self, event):
        if event.inaxes != self.ax: return
        click_pos = np.array([event.xdata, event.ydata])
        distances = np.linalg.norm(self.pos - click_pos, axis=1)
        closest_idx = np.argmin(distances)
        
        if distances[closest_idx] < 0.05:
            if event.button == 1: # Left click to drag
                self.dragging_idx = closest_idx
            elif event.button == 3: # Right click to select focal particle
                self.selected_idx = closest_idx
                self.update_physics_and_plot()

    def on_motion(self, event):
        if self.dragging_idx is not None and event.inaxes == self.ax:
            self.pos[self.dragging_idx] = [event.xdata, event.ydata]
            # Real-time recalculation of density and matrices while dragging
            self.compute_density_and_h()
            self.compute_geometry_and_gradients()
            self.update_physics_and_plot()

    def on_release(self, event):
        self.dragging_idx = None

if __name__ == "__main__":
    app = MFMFluidSandbox(num_particles=36)
