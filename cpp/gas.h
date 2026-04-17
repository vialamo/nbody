#pragma once
#include <limits>

#include "config.h"
#include "types.h"

class GasGrid;
struct SimState;
struct ZeldovichField;

class RiemannSolver {
   private:
    Grid3D density, mom_n, mom_t1, mom_t2, energy, v_n, v_t1, v_t2, pressure;
    Grid3D rho_L, p_L, vn_L, vt1_L, vt2_L, E_L, mom_n_L, mom_t1_L, mom_t2_L;
    Grid3D rho_R, p_R, vn_R, vt1_R, vt2_R, E_R, mom_n_R, mom_t1_R, mom_t2_R;
    Grid3D cs_L, cs_R, S_L, S_R, S_R_minus_S_L;
    Grid3D F_dens_L, F_dens_R, F_momn_L, F_momn_R, F_momt1_L, F_momt1_R,
        F_momt2_L, F_momt2_R, F_en_L, F_en_R;
    Grid3D flux_density, flux_mom_n, flux_mom_t1, flux_mom_t2, flux_energy;
    Grid3D flux_density_sh, flux_mom_n_sh, flux_mom_t1_sh, flux_mom_t2_sh,
        flux_energy_sh;

   public:
    RiemannSolver(int mesh_size);
    Grid3D solve_hll(const Grid3D& FL, const Grid3D& FR, const Grid3D& UL,
                     const Grid3D& UR);
    void compute_fluxes(const GasGrid& grid, int axis, double gamma);

    const Grid3D& get_flux_density() const { return flux_density; }
    const Grid3D& get_flux_density_sh() const { return flux_density_sh; }
    const Grid3D& get_flux_mom_n() const { return flux_mom_n; }
    const Grid3D& get_flux_mom_n_sh() const { return flux_mom_n_sh; }
    const Grid3D& get_flux_mom_t1() const { return flux_mom_t1; }
    const Grid3D& get_flux_mom_t1_sh() const { return flux_mom_t1_sh; }
    const Grid3D& get_flux_mom_t2() const { return flux_mom_t2; }
    const Grid3D& get_flux_mom_t2_sh() const { return flux_mom_t2_sh; }
    const Grid3D& get_flux_energy() const { return flux_energy; }
    const Grid3D& get_flux_energy_sh() const { return flux_energy_sh; }
};

class GasGrid {
    Grid3D density, momentum_x, momentum_y, momentum_z, energy;
    Grid3D pressure, velocity_x, velocity_y, velocity_z;

    RiemannSolver solver;
    const Config& config;

    friend void initialize_gas(SimState& state, const Config& config,
                               const ZeldovichField& zf);
    friend void apply_gas_kick(GasGrid& gas, const Grid3D& grav_x,
                               const Grid3D& grav_y, const Grid3D& grav_z,
                               double dt, double a, double H,
                               const Config& config);

   public:
    GasGrid(const Config& conf);

    void hydro_step(double dt);
    double get_cfl_timestep() const;

    const Grid3D& get_density() const { return density; }
    const Grid3D& get_momentum_x() const { return momentum_x; }
    const Grid3D& get_momentum_y() const { return momentum_y; }
    const Grid3D& get_momentum_z() const { return momentum_z; }
    const Grid3D& get_energy() const { return energy; }

    const Grid3D& get_pressure() const { return pressure; }
    const Grid3D& get_velocity_x() const { return velocity_x; }
    const Grid3D& get_velocity_y() const { return velocity_y; }
    const Grid3D& get_velocity_z() const { return velocity_z; }

   private:
    void update_primitive_variables();
};