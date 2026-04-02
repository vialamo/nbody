#pragma once
#include "types.h"
#include "config.h"
#include <limits>

struct HydroScratchPad {
    Grid density, mom_n, mom_t, energy, v_n, v_t, pressure;
    Grid rho_L, p_L, vn_L, vt_L, E_L, mom_n_L, mom_t_L;
    Grid rho_R, p_R, vn_R, vt_R, E_R, mom_n_R, mom_t_R;
    Grid cs_L, cs_R, S_L, S_R, S_R_minus_S_L;
    Grid F_dens_L, F_dens_R, F_momn_L, F_momn_R, F_momt_L, F_momt_R, F_en_L, F_en_R;
    Grid flux_density, flux_mom_n, flux_mom_t, flux_energy;

    HydroScratchPad( int mesh_size );
};

Grid roll( const Grid& m, int shift, int axis );

struct GasGrid {
    Grid density, momentum_x, momentum_y, energy;
    Grid pressure, velocity_x, velocity_y;
    Grid accel_x, accel_y;

    HydroScratchPad scratch;
    const Config& config;

    GasGrid( const Config& conf );

    void update_primitive_variables();
    void calculate_fluxes( int axis );
    void hydro_step( double dt );

    double get_cfl_timestep() const;
};
