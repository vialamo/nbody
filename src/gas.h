#pragma once

#include "config.h"
#include "cooling.h"
#include "types.h"

class GasGrid;
struct SimState;
struct ZeldovichField;

class RiemannSolver {
   private:
    Grid3D rho_L, p_L, vn_L, vt1_L, vt2_L, E_L, mom_n_L, mom_t1_L, mom_t2_L;
    Grid3D rho_R, p_R, vn_R, vt1_R, vt2_R, E_R, mom_n_R, mom_t1_R, mom_t2_R;
    Grid3D cs_L, cs_R, S_L, S_R;
    Grid3D F_dens_L, F_dens_R, F_momn_L, F_momn_R, F_momt1_L, F_momt1_R,
        F_momt2_L, F_momt2_R, F_en_L, F_en_R;
    Grid3D flux_density, flux_mom_n, flux_mom_t1, flux_mom_t2, flux_energy;
    Grid3D flux_density_sh, flux_mom_n_sh, flux_mom_t1_sh, flux_mom_t2_sh,
        flux_energy_sh;
    Grid3D ie_L, ie_R, F_ie_L, F_ie_R, flux_ie, flux_ie_sh;
    Grid3D q_minus, q_plus, dq_L, dq_R, slope;
    Grid3D flux_metal, flux_metal_sh;
    Grid3D Z, Z_L, Z_R;

    Eigen::ArrayXd den_star, S_star, omega_L, omega_R, denom_L, denom_R;
    Eigen::ArrayXd mom_n_star_L, mom_n_star_R, mom_t1_star_L, mom_t1_star_R;
    Eigen::ArrayXd mom_t2_star_L, mom_t2_star_R, E_star_L, E_star_R, ie_star_L,
        ie_star_R;

    friend struct RiemannSolverTestAccess;

    inline void reconstruct_muscl(const Grid3D& q, Grid3D& q_L, Grid3D& q_R,
                                  int axis);

    inline void apply_hllc_flux(const Grid3D& F_L, const Grid3D& F_R,
                                const Grid3D& U_L, const Grid3D& U_R,
                                const Eigen::ArrayXd& U_star_L,
                                const Eigen::ArrayXd& U_star_R,
                                const Eigen::ArrayXd& S_star, Grid3D& F_out);

   public:
    RiemannSolver(int mesh_size);
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
    const Grid3D& get_flux_ie() const { return flux_ie; }
    const Grid3D& get_flux_ie_sh() const { return flux_ie_sh; }
    const Grid3D& get_flux_metal() const { return flux_metal; }
    const Grid3D& get_flux_metal_sh() const { return flux_metal_sh; }
};

class GasGrid {
   private:
    Grid3D density, momentum_x, momentum_y, momentum_z, energy;
    Grid3D pressure, velocity_x, velocity_y, velocity_z;
    Grid3D internal_energy;
    Grid3D metal_density;
    size_t cooling_failed_cells;
    size_t cooling_total_cycles;
    double accumulated_radiated_energy;
    double accumulated_photoheating_energy;
    double accumulated_gravitational_work;
    double accumulated_expansion_work;
    double pressure_floor;

    RiemannSolver solver;
    Cooling cooling_module;
    const Config& config;

    friend struct GasGridTestAccess;
    friend void initialize_gas(SimState& state, const Config& config,
                               const ZeldovichField& zf);
    friend void apply_gas_kick(GasGrid& gas, const Grid3D& grav_x,
                               const Grid3D& grav_y, const Grid3D& grav_z,
                               double dt, double a, double H,
                               const Config& config);

    void update_primitive_variables();
    void compute_and_apply_fluxes(double dt);

   public:
    GasGrid(const Config& conf);

    void hydro_step(double dt);
    void apply_cooling(double dt, double a);
    double get_cfl_timestep() const;
    double get_cooling_timestep(double a) const;

    const Grid3D& get_density() const { return density; }
    const Grid3D& get_momentum_x() const { return momentum_x; }
    const Grid3D& get_momentum_y() const { return momentum_y; }
    const Grid3D& get_momentum_z() const { return momentum_z; }
    const Grid3D& get_energy() const { return energy; }

    const Grid3D& get_pressure() const { return pressure; }
    const Grid3D& get_velocity_x() const { return velocity_x; }
    const Grid3D& get_velocity_y() const { return velocity_y; }
    const Grid3D& get_velocity_z() const { return velocity_z; }

    const Grid3D& get_internal_energy() const { return internal_energy; }
    const Grid3D& get_metal_density() const { return metal_density; }

    size_t get_cooling_failed_cells() const { return cooling_failed_cells; }
    size_t get_cooling_total_cycles() const { return cooling_total_cycles; }
    double get_accumulated_radiated_energy() const {
        return accumulated_radiated_energy;
    }

    double get_accumulated_photoheating_energy() const {
        return accumulated_photoheating_energy;
    }

    double get_accumulated_gravitational_work() const {
        return accumulated_gravitational_work;
    }
    double get_accumulated_expansion_work() const {
        return accumulated_expansion_work;
    }
    void add_gravitational_work(double w) {
        accumulated_gravitational_work += w;
    }
    void add_expansion_work(double w) { accumulated_expansion_work += w; }

    Grid3D compute_thermal_timescale(double a) const;
};