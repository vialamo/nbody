#pragma once
#include "types.h"
#include "config.h"
#include <limits>

struct HydroScratchPad {
    Grid3D density, mom_n, mom_t1, mom_t2, energy, v_n, v_t1, v_t2, pressure;
    Grid3D rho_L, p_L, vn_L, vt1_L, vt2_L, E_L, mom_n_L, mom_t1_L, mom_t2_L;
    Grid3D rho_R, p_R, vn_R, vt1_R, vt2_R, E_R, mom_n_R, mom_t1_R, mom_t2_R;
    Grid3D cs_L, cs_R, S_L, S_R, S_R_minus_S_L;
    Grid3D F_dens_L, F_dens_R, F_momn_L, F_momn_R, F_momt1_L, F_momt1_R, F_momt2_L, F_momt2_R, F_en_L, F_en_R;
    Grid3D flux_density, flux_mom_n, flux_mom_t1, flux_mom_t2, flux_energy;

    HydroScratchPad( int mesh_size );
};

Grid3D roll( const Grid3D& m, int shift, int axis );

struct GasGrid {
    Grid3D density, momentum_x, momentum_y, momentum_z, energy;
    Grid3D pressure, velocity_x, velocity_y, velocity_z;
    Grid3D accel_x, accel_y, accel_z;

    HydroScratchPad scratch;
    const Config& config;

    GasGrid( const Config& conf );

    void update_primitive_variables();
    void calculate_fluxes( int axis );
    void hydro_step( double dt );

    double get_cfl_timestep() const;
};