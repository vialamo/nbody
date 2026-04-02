import numpy as np

class GasGrid:
    def __init__(self, config):
        self.config = config
        mesh_size = config.mesh_size
        grid_shape = (mesh_size, mesh_size, mesh_size)
        
        self.density = np.zeros(grid_shape)
        self.momentum_x = np.zeros(grid_shape)
        self.momentum_y = np.zeros(grid_shape)
        self.momentum_z = np.zeros(grid_shape)
        self.energy = np.zeros(grid_shape)
        
        self.pressure = np.zeros(grid_shape)
        self.velocity_x = np.zeros(grid_shape)
        self.velocity_y = np.zeros(grid_shape)
        self.velocity_z = np.zeros(grid_shape)

        self.accel_x = np.zeros(grid_shape)
        self.accel_y = np.zeros(grid_shape)
        self.accel_z = np.zeros(grid_shape)

    def update_primitive_variables(self):
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
        self.pressure[self.pressure < 1e-12] = 1e-12

    def calculate_fluxes(self, axis):
        if axis == 0: 
            density, mom_n, mom_t1, mom_t2, energy, v_n, v_t1, v_t2, pressure = \
                self.density, self.momentum_x, self.momentum_y, self.momentum_z, self.energy, \
                self.velocity_x, self.velocity_y, self.velocity_z, self.pressure
        elif axis == 1: 
            density = np.transpose(self.density, (1, 0, 2))
            mom_n = np.transpose(self.momentum_y, (1, 0, 2))
            mom_t1 = np.transpose(self.momentum_x, (1, 0, 2))
            mom_t2 = np.transpose(self.momentum_z, (1, 0, 2))
            energy = np.transpose(self.energy, (1, 0, 2))
            v_n = np.transpose(self.velocity_y, (1, 0, 2))
            v_t1 = np.transpose(self.velocity_x, (1, 0, 2))
            v_t2 = np.transpose(self.velocity_z, (1, 0, 2))
            pressure = np.transpose(self.pressure, (1, 0, 2))
        elif axis == 2: 
            density = np.transpose(self.density, (2, 0, 1))
            mom_n = np.transpose(self.momentum_z, (2, 0, 1))
            mom_t1 = np.transpose(self.momentum_x, (2, 0, 1))
            mom_t2 = np.transpose(self.momentum_y, (2, 0, 1))
            energy = np.transpose(self.energy, (2, 0, 1))
            v_n = np.transpose(self.velocity_z, (2, 0, 1))
            v_t1 = np.transpose(self.velocity_x, (2, 0, 1))
            v_t2 = np.transpose(self.velocity_y, (2, 0, 1))
            pressure = np.transpose(self.pressure, (2, 0, 1))

        rho_L, p_L, vn_L, vt1_L, vt2_L, E_L, mom_n_L, mom_t1_L, mom_t2_L = \
            density, pressure, v_n, v_t1, v_t2, energy, mom_n, mom_t1, mom_t2
        
        rho_R, p_R, vn_R, vt1_R, vt2_R, E_R, mom_n_R, mom_t1_R, mom_t2_R = \
            np.roll(rho_L, -1, axis=0), np.roll(p_L, -1, axis=0), np.roll(vn_L, -1, axis=0), \
            np.roll(vt1_L, -1, axis=0), np.roll(vt2_L, -1, axis=0), np.roll(E_L, -1, axis=0), \
            np.roll(mom_n_L, -1, axis=0), np.roll(mom_t1_L, -1, axis=0), np.roll(mom_t2_L, -1, axis=0)

        gamma = self.config.gamma
        cs_L = np.sqrt(gamma * p_L / rho_L, where=rho_L > 1e-12, out=np.zeros_like(rho_L))
        cs_R = np.sqrt(gamma * p_R / rho_R, where=rho_R > 1e-12, out=np.zeros_like(rho_R))

        S_L = np.minimum(vn_L - cs_L, vn_R - cs_R)
        S_R = np.maximum(vn_L + cs_L, vn_R + cs_R)

        F_dens_L, F_dens_R = rho_L * vn_L, rho_R * vn_R
        F_momn_L, F_momn_R = rho_L * vn_L**2 + p_L, rho_R * vn_R**2 + p_R
        F_momt1_L, F_momt1_R = rho_L * vn_L * vt1_L, rho_R * vn_R * vt1_R
        F_momt2_L, F_momt2_R = rho_L * vn_L * vt2_L, rho_R * vn_R * vt2_R
        F_en_L, F_en_R = (E_L + p_L) * vn_L, (E_R + p_R) * vn_R

        S_R_minus_S_L = S_R - S_L
        S_R_minus_S_L[np.abs(S_R_minus_S_L) < 1e-9] = 1e-9 

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

        if axis == 0: 
            return flux_density, flux_mom_n, flux_mom_t1, flux_mom_t2, flux_energy
        elif axis == 1: 
            return flux_density.transpose(1, 0, 2), flux_mom_t1.transpose(1, 0, 2), \
                   flux_mom_n.transpose(1, 0, 2), flux_mom_t2.transpose(1, 0, 2), flux_energy.transpose(1, 0, 2)
        elif axis == 2: 
            return flux_density.transpose(1, 2, 0), flux_mom_t1.transpose(1, 2, 0), \
                   flux_mom_t2.transpose(1, 2, 0), flux_mom_n.transpose(1, 2, 0), flux_energy.transpose(1, 2, 0)

    def hydro_step(self, dt):
        self.update_primitive_variables()
        factor = dt / self.config.cell_size
        
        flux_density, flux_mom_x, flux_mom_y, flux_mom_z, flux_energy = self.calculate_fluxes(axis=0)
        self.density    -= factor * (flux_density - np.roll(flux_density, 1, axis=0))
        self.momentum_x -= factor * (flux_mom_x - np.roll(flux_mom_x, 1, axis=0))
        self.momentum_y -= factor * (flux_mom_y - np.roll(flux_mom_y, 1, axis=0))
        self.momentum_z -= factor * (flux_mom_z - np.roll(flux_mom_z, 1, axis=0))
        self.energy     -= factor * (flux_energy - np.roll(flux_energy, 1, axis=0))
        self.update_primitive_variables()

        flux_density, flux_mom_x, flux_mom_y, flux_mom_z, flux_energy = self.calculate_fluxes(axis=1)
        self.density    -= factor * (flux_density - np.roll(flux_density, 1, axis=1))
        self.momentum_x -= factor * (flux_mom_x - np.roll(flux_mom_x, 1, axis=1))
        self.momentum_y -= factor * (flux_mom_y - np.roll(flux_mom_y, 1, axis=1))
        self.momentum_z -= factor * (flux_mom_z - np.roll(flux_mom_z, 1, axis=1))
        self.energy     -= factor * (flux_energy - np.roll(flux_energy, 1, axis=1))
        self.update_primitive_variables()
        
        flux_density, flux_mom_x, flux_mom_y, flux_mom_z, flux_energy = self.calculate_fluxes(axis=2)
        self.density    -= factor * (flux_density - np.roll(flux_density, 1, axis=2))
        self.momentum_x -= factor * (flux_mom_x - np.roll(flux_mom_x, 1, axis=2))
        self.momentum_y -= factor * (flux_mom_y - np.roll(flux_mom_y, 1, axis=2))
        self.momentum_z -= factor * (flux_mom_z - np.roll(flux_mom_z, 1, axis=2))
        self.energy     -= factor * (flux_energy - np.roll(flux_energy, 1, axis=2))
        self.update_primitive_variables()

def get_cfl_timestep(gas, config):
    if not config.enable_hydro or not np.any(gas.density > 1e-12): return np.inf
    non_zero = gas.density > 1e-12
    cs = np.sqrt(config.gamma * gas.pressure[non_zero] / gas.density[non_zero])
    v_mag = np.sqrt(gas.velocity_x[non_zero]**2 + gas.velocity_y[non_zero]**2 + gas.velocity_z[non_zero]**2)
    max_signal = np.max(v_mag + cs)
    return config.cell_size / max_signal * config.cfl_safety_factor if max_signal >= 1e-9 else np.inf