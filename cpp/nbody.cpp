#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <complex>
#include <chrono>
#include <iomanip>

#include <Eigen/Dense>
#include <SFML/Graphics.hpp>
#include <H5Cpp.h>

// Header-only pocketfft library
#include "pocketfft_hdronly.h"

#include "ini.h"

// Define M_PI if it's not provided by the compiler (it's non-standard)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ------------------------
// Simulation Config Struct
// ------------------------
using Grid = Eigen::MatrixXd;

struct Config {
    // [domain]
    double DOMAIN_SIZE;
    int MESH_SIZE;

    // [matter]
    double OMEGA_BARYON;
    int N_PER_SIDE;
    bool USE_HYDRO;

    // [physics]
    double G;
    double GAMMA;
    bool STANDING_PARTICLES;
    bool EXPANDING_UNIVERSE;
    double EXPANSION_START_T;
    double INITIAL_POWER_SPECTRUM_INDEX;

    // [p3m]
    bool USE_PM;
    bool USE_PP;
    double CUTOFF_RADIUS_CELLS;
    double CUTOFF_TRANSITION_WIDTH_FACTOR;

    // [time]
    double DT_FACTOR;

    // [output]
    int SAVE_HDF5_EVERY_CYCLES;
    int DEBUG_INFO_EVERY_CYCLES;
    int SAVE_RENDER_EVERY_CYCLES;
    int SEED;

    // [visualization]
    unsigned int RENDER_SIZE;

    // Derived Parameters
    double CELL_SIZE;
    double OMEGA_DM;
    double CUTOFF_RADIUS;
    double CUTOFF_RADIUS_SQUARED;
    double CUTOFF_TRANSITION_WIDTH;
    double R_SWITCH_START;
    double R_SWITCH_START_SQ;
    int NUM_DM_PARTICLES;
    double TOTAL_MASS; // Assuming 1.0
    double DM_PARTICLE_MASS;
    double GAS_TOTAL_MASS;
    double MEAN_INTERPARTICLE_SPACING;
    double SOFTENING_SQUARED;
    double DYNAMICAL_TIME;
    double DT;
    double RENDER_SCALE;

    // Global Cosmological Variables
    double scale_factor;
    double hubble_param;

    // Load and calculate all parameters
    void load( const std::string& filename ) {
        ini::Ini config_file;
        config_file.load( filename );

        // Load from INI
        DOMAIN_SIZE = config_file.get_double( "domain", "domain_size", 1.0 );
        MESH_SIZE = config_file.get_int( "domain", "mesh_size", 64 );

        OMEGA_BARYON = config_file.get_double( "matter", "omega_baryon", 0.15 );
        N_PER_SIDE = config_file.get_int( "matter", "n_per_side", 30 );
        USE_HYDRO = config_file.get_bool( "matter", "use_hydro", true );

        G = config_file.get_double( "physics", "g_const", 1.0 / ( 6.0 * M_PI ) );
        GAMMA = config_file.get_double( "physics", "gamma", 5.0 / 3.0 );
        STANDING_PARTICLES = config_file.get_bool( "physics", "standing_particles", false );
        EXPANDING_UNIVERSE = config_file.get_bool( "physics", "expanding_universe", true );
        EXPANSION_START_T = config_file.get_double( "physics", "expansion_start_t", 0.1 );
        INITIAL_POWER_SPECTRUM_INDEX = config_file.get_double( "physics", "initial_power_spectrum_index", -2.0 );

        USE_PM = config_file.get_bool( "p3m", "use_pm", true );
        USE_PP = config_file.get_bool( "p3m", "use_pp", true );
        CUTOFF_RADIUS_CELLS = config_file.get_double( "p3m", "cutoff_radius_cells", 2.5 );
        CUTOFF_TRANSITION_WIDTH_FACTOR = config_file.get_double( "p3m", "cutoff_transition_width_factor", 0.2 );

        DT_FACTOR = config_file.get_double( "time", "dt_factor", 1e-6 );

        SAVE_HDF5_EVERY_CYCLES = config_file.get_int( "output", "save_hdf5_every_cycles", 100 );
        DEBUG_INFO_EVERY_CYCLES = config_file.get_int( "output", "debug_info_every_cycles", 40 );
        SAVE_RENDER_EVERY_CYCLES = config_file.get_int( "output", "save_render_every_cycles", 0 );
        SEED = config_file.get_int( "output", "seed", 42 );

        RENDER_SIZE = config_file.get_int( "visualization", "render_size", 800 );

        // Calculate Derived Parameters
        CELL_SIZE = DOMAIN_SIZE / MESH_SIZE;
        OMEGA_DM = 1.0 - OMEGA_BARYON;

        CUTOFF_RADIUS = CUTOFF_RADIUS_CELLS * CELL_SIZE;
        CUTOFF_RADIUS_SQUARED = CUTOFF_RADIUS * CUTOFF_RADIUS;
        CUTOFF_TRANSITION_WIDTH = CUTOFF_TRANSITION_WIDTH_FACTOR * CUTOFF_RADIUS;
        R_SWITCH_START = CUTOFF_RADIUS - CUTOFF_TRANSITION_WIDTH;
        R_SWITCH_START_SQ = R_SWITCH_START * R_SWITCH_START;

        NUM_DM_PARTICLES = N_PER_SIDE * N_PER_SIDE;
        TOTAL_MASS = 1.0;
        DM_PARTICLE_MASS = ( TOTAL_MASS * OMEGA_DM ) / NUM_DM_PARTICLES;
        GAS_TOTAL_MASS = TOTAL_MASS * OMEGA_BARYON;

        MEAN_INTERPARTICLE_SPACING = DOMAIN_SIZE / sqrt( (double)NUM_DM_PARTICLES );
        SOFTENING_SQUARED = pow( MEAN_INTERPARTICLE_SPACING / 50.0, 2 );

        DYNAMICAL_TIME = 1.0 / sqrt( G );
        DT = DT_FACTOR * DYNAMICAL_TIME;

        RENDER_SCALE = RENDER_SIZE / DOMAIN_SIZE;

        scale_factor = 1.0;
        hubble_param = 0.0;
    }
};

// ------------------------
// Vector and Particle Structs
// ------------------------
struct Vec2 {
    double x = 0.0, y = 0.0;
};

struct Particle {
    Vec2 pos;
    Vec2 vel;
    Vec2 acc;
    double mass;
};


// Helper: Rolls an Eigen matrix
Grid roll( const Grid& m, int shift, int axis ) {
    Grid result( m.rows(), m.cols() );
    if( axis == 0 ) {
        int rows = static_cast< int >( m.rows() );
        shift = ( shift % rows + rows ) % rows; // Handle negative shifts
        result.bottomRows( rows - shift ) = m.topRows( rows - shift );
        result.topRows( shift ) = m.bottomRows( shift );
    }
    else {
        int cols = static_cast< int >( m.cols() );
        shift = ( shift % cols + cols ) % cols; // Handle negative shifts
        result.rightCols( cols - shift ) = m.leftCols( cols - shift );
        result.leftCols( shift ) = m.rightCols( shift );
    }
    return result;
}

struct HydroScratchPad {
    Grid density, mom_n, mom_t, energy, v_n, v_t, pressure;
    Grid rho_L, p_L, vn_L, vt_L, E_L, mom_n_L, mom_t_L;
    Grid rho_R, p_R, vn_R, vt_R, E_R, mom_n_R, mom_t_R;
    Grid cs_L, cs_R, S_L, S_R, S_R_minus_S_L;
    Grid F_dens_L, F_dens_R, F_momn_L, F_momn_R, F_momt_L, F_momt_R, F_en_L, F_en_R;
    Grid flux_density, flux_mom_n, flux_mom_t, flux_energy;

    // Initialize all grids to their size
    HydroScratchPad( int mesh_size ) :
        density( mesh_size, mesh_size ), mom_n( mesh_size, mesh_size ), mom_t( mesh_size, mesh_size ),
        energy( mesh_size, mesh_size ), v_n( mesh_size, mesh_size ), v_t( mesh_size, mesh_size ),
        pressure( mesh_size, mesh_size ), rho_L( mesh_size, mesh_size ), p_L( mesh_size, mesh_size ),
        vn_L( mesh_size, mesh_size ), vt_L( mesh_size, mesh_size ), E_L( mesh_size, mesh_size ),
        mom_n_L( mesh_size, mesh_size ), mom_t_L( mesh_size, mesh_size ), rho_R( mesh_size, mesh_size ),
        p_R( mesh_size, mesh_size ), vn_R( mesh_size, mesh_size ), vt_R( mesh_size, mesh_size ),
        E_R( mesh_size, mesh_size ), mom_n_R( mesh_size, mesh_size ), cs_L( mesh_size, mesh_size ),
        cs_R( mesh_size, mesh_size ), S_L( mesh_size, mesh_size ), S_R( mesh_size, mesh_size ),
        S_R_minus_S_L( mesh_size, mesh_size ), F_dens_L( mesh_size, mesh_size ), F_dens_R( mesh_size, mesh_size ),
        F_momn_L( mesh_size, mesh_size ), F_momn_R( mesh_size, mesh_size ), F_momt_L( mesh_size, mesh_size ),
        F_momt_R( mesh_size, mesh_size ), F_en_L( mesh_size, mesh_size ), F_en_R( mesh_size, mesh_size ),
        flux_density( mesh_size, mesh_size ), flux_mom_n( mesh_size, mesh_size ),
        flux_mom_t( mesh_size, mesh_size ), flux_energy( mesh_size, mesh_size ) {
    }
};

// Global grids, initialized in main() after config is read
Grid g_grid_x;
Grid g_grid_y;

struct GasGrid {
    // Conservative variables
    Grid density;
    Grid momentum_x;
    Grid momentum_y;
    Grid energy;

    // Primitive variables
    Grid pressure;
    Grid velocity_x;
    Grid velocity_y;

    HydroScratchPad scratch;
    const Config& config; // Reference to config

    GasGrid( const Config& conf ) :
        density( conf.MESH_SIZE, conf.MESH_SIZE ),
        momentum_x( conf.MESH_SIZE, conf.MESH_SIZE ),
        momentum_y( conf.MESH_SIZE, conf.MESH_SIZE ),
        energy( conf.MESH_SIZE, conf.MESH_SIZE ),
        pressure( conf.MESH_SIZE, conf.MESH_SIZE ),
        velocity_x( conf.MESH_SIZE, conf.MESH_SIZE ),
        velocity_y( conf.MESH_SIZE, conf.MESH_SIZE ),
        scratch( conf.MESH_SIZE ),
        config( conf )
    {
        density.setZero();
        momentum_x.setZero();
        momentum_y.setZero();
        energy.setZero();
        pressure.setZero();
        velocity_x.setZero();
        velocity_y.setZero();
    }

    void update_primitive_variables() {
        for( int i = 0; i < config.MESH_SIZE; ++i ) {
            for( int j = 0; j < config.MESH_SIZE; ++j ) {
                if( density( i, j ) > 1e-12 ) {
                    velocity_x( i, j ) = momentum_x( i, j ) / density( i, j );
                    velocity_y( i, j ) = momentum_y( i, j ) / density( i, j );
                }
                else {
                    velocity_x( i, j ) = 0.0;
                    velocity_y( i, j ) = 0.0;
                }
            }
        }

        Grid kinetic_energy_density = 0.5 * ( momentum_x.array().square() + momentum_y.array().square() );
        for( int i = 0; i < config.MESH_SIZE; ++i ) {
            for( int j = 0; j < config.MESH_SIZE; ++j ) {
                if( density( i, j ) > 1e-12 ) {
                    kinetic_energy_density( i, j ) /= density( i, j );
                }
                else {
                    kinetic_energy_density( i, j ) = 0.0;
                }
            }
        }

        Grid internal_energy_density = energy - kinetic_energy_density;
        pressure = ( config.GAMMA - 1.0 ) * internal_energy_density;

        // Apply pressure floor
        pressure = ( pressure.array() < 1e-12 ).select( 1e-12, pressure );
    }

    // Helper: Calculates HLL fluxes for one axis
    void calculate_fluxes( int axis ) {
        if( axis == 1 ) { // Permute axes for y-direction
            scratch.density = density.transpose();
            scratch.mom_n = momentum_y.transpose();
            scratch.mom_t = momentum_x.transpose();
            scratch.energy = energy.transpose();
            scratch.v_n = velocity_y.transpose();
            scratch.v_t = velocity_x.transpose();
            scratch.pressure = pressure.transpose();
        }
        else {
            scratch.density = density;
            scratch.mom_n = momentum_x;
            scratch.mom_t = momentum_y;
            scratch.energy = energy;
            scratch.v_n = velocity_x;
            scratch.v_t = velocity_y;
            scratch.pressure = pressure;
        }

        const int rollDir = -1;
        scratch.rho_L = scratch.density;
        scratch.p_L = scratch.pressure;
        scratch.vn_L = scratch.v_n;
        scratch.vt_L = scratch.v_t;
        scratch.E_L = scratch.energy;
        scratch.mom_n_L = scratch.mom_n;
        scratch.mom_t_L = scratch.mom_t;
        scratch.rho_R = roll( scratch.rho_L, rollDir, 0 );
        scratch.p_R = roll( scratch.p_L, rollDir, 0 );
        scratch.vn_R = roll( scratch.vn_L, rollDir, 0 );
        scratch.vt_R = roll( scratch.vt_L, rollDir, 0 );
        scratch.E_R = roll( scratch.E_L, rollDir, 0 );
        scratch.mom_n_R = roll( scratch.mom_n_L, rollDir, 0 );
        scratch.mom_t_R = roll( scratch.mom_t_L, rollDir, 0 );

        scratch.cs_L = ( config.GAMMA * scratch.p_L.array() / scratch.rho_L.array() ).sqrt();
        scratch.cs_R = ( config.GAMMA * scratch.p_R.array() / scratch.rho_R.array() ).sqrt();

        scratch.cs_L = ( scratch.rho_L.array() > 1e-12 ).select( scratch.cs_L, 0.0 );
        scratch.cs_R = ( scratch.rho_R.array() > 1e-12 ).select( scratch.cs_R, 0.0 );
        // Handle NaNs from sqrt(negative pressure)
        scratch.cs_L = ( scratch.cs_L.array() == scratch.cs_L.array() ).select( scratch.cs_L, 0.0 );
        scratch.cs_R = ( scratch.cs_R.array() == scratch.cs_R.array() ).select( scratch.cs_R, 0.0 );


        scratch.S_L = ( scratch.vn_L.array() - scratch.cs_L.array() ).cwiseMin( scratch.vn_R.array() - scratch.cs_R.array() );
        scratch.S_R = ( scratch.vn_L.array() + scratch.cs_L.array() ).cwiseMax( scratch.vn_R.array() + scratch.cs_R.array() );

        scratch.F_dens_L = scratch.rho_L.array() * scratch.vn_L.array();
        scratch.F_dens_R = scratch.rho_R.array() * scratch.vn_R.array();
        scratch.F_momn_L = scratch.rho_L.array() * scratch.vn_L.array().square() + scratch.p_L.array();
        scratch.F_momn_R = scratch.rho_R.array() * scratch.vn_R.array().square() + scratch.p_R.array();
        scratch.F_momt_L = scratch.rho_L.array() * scratch.vn_L.array() * scratch.vt_L.array();
        scratch.F_momt_R = scratch.rho_R.array() * scratch.vn_R.array() * scratch.vt_R.array();
        scratch.F_en_L = ( scratch.E_L.array() + scratch.p_L.array() ) * scratch.vn_L.array();
        scratch.F_en_R = ( scratch.E_R.array() + scratch.p_R.array() ) * scratch.vn_R.array();

        scratch.S_R_minus_S_L = scratch.S_R - scratch.S_L;
        scratch.S_R_minus_S_L = ( scratch.S_R_minus_S_L.array().abs() < 1e-9 ).select( 1e-9, scratch.S_R_minus_S_L );

        scratch.flux_density = ( scratch.S_L.array() >= 0 ).select( scratch.F_dens_L,
            ( scratch.S_R.array() <= 0 ).select( scratch.F_dens_R,
                ( scratch.S_R.array() * scratch.F_dens_L.array() - scratch.S_L.array() * scratch.F_dens_R.array() + scratch.S_L.array() * scratch.S_R.array() * ( scratch.rho_R.array() - scratch.rho_L.array() ) ) / scratch.S_R_minus_S_L.array() ) );
        scratch.flux_mom_n = ( scratch.S_L.array() >= 0 ).select( scratch.F_momn_L,
            ( scratch.S_R.array() <= 0 ).select( scratch.F_momn_R,
                ( scratch.S_R.array() * scratch.F_momn_L.array() - scratch.S_L.array() * scratch.F_momn_R.array() + scratch.S_L.array() * scratch.S_R.array() * ( scratch.mom_n_R.array() - scratch.mom_n_L.array() ) ) / scratch.S_R_minus_S_L.array() ) );
        scratch.flux_mom_t = ( scratch.S_L.array() >= 0 ).select( scratch.F_momt_L,
            ( scratch.S_R.array() <= 0 ).select( scratch.F_momt_R,
                ( scratch.S_R.array() * scratch.F_momt_L.array() - scratch.S_L.array() * scratch.F_momt_R.array() + scratch.S_L.array() * scratch.S_R.array() * ( scratch.mom_t_R.array() - scratch.mom_t_L.array() ) ) / scratch.S_R_minus_S_L.array() ) );
        scratch.flux_energy = ( scratch.S_L.array() >= 0 ).select( scratch.F_en_L,
            ( scratch.S_R.array() <= 0 ).select( scratch.F_en_R,
                ( scratch.S_R.array() * scratch.F_en_L.array() - scratch.S_L.array() * scratch.F_en_R.array() + scratch.S_L.array() * scratch.S_R.array() * ( scratch.E_R.array() - scratch.E_L.array() ) ) / scratch.S_R_minus_S_L.array() ) );

        if( axis == 1 ) {
            scratch.flux_density.transposeInPlace();
            scratch.flux_mom_t.transposeInPlace();
            scratch.flux_mom_n.transposeInPlace();
            scratch.flux_energy.transposeInPlace();
        }
    }

    void hydro_step( double dt ) {
        update_primitive_variables();

        calculate_fluxes( 0 );
        double factor = dt / config.CELL_SIZE;
        density -= factor * ( scratch.flux_density - roll( scratch.flux_density, 1, 0 ) );
        momentum_x -= factor * ( scratch.flux_mom_n - roll( scratch.flux_mom_n, 1, 0 ) );
        momentum_y -= factor * ( scratch.flux_mom_t - roll( scratch.flux_mom_t, 1, 0 ) );
        energy -= factor * ( scratch.flux_energy - roll( scratch.flux_energy, 1, 0 ) );
        update_primitive_variables();

        calculate_fluxes( 1 );
        density -= factor * ( scratch.flux_density - roll( scratch.flux_density, 1, 1 ) );
        momentum_x -= factor * ( scratch.flux_mom_t - roll( scratch.flux_mom_t, 1, 1 ) );
        momentum_y -= factor * ( scratch.flux_mom_n - roll( scratch.flux_mom_n, 1, 1 ) );
        energy -= factor * ( scratch.flux_energy - roll( scratch.flux_energy, 1, 1 ) );
        update_primitive_variables();
    }
};

struct CIC_Data {
    int ix, iy;
    double w1, w2, w3, w4;
};


// ------------------------
// Helper Functions
// ------------------------
double displacement( double dx, const Config& config ) {
    return fmod( dx + 0.5 * config.DOMAIN_SIZE + config.DOMAIN_SIZE, config.DOMAIN_SIZE ) - 0.5 * config.DOMAIN_SIZE;
}

void update_cosmology( Config& config, double total_time ) {
    if( config.EXPANDING_UNIVERSE ) {
        double expansion_time = config.EXPANSION_START_T + total_time;
        config.scale_factor = pow( expansion_time, 2.0 / 3.0 );
        config.hubble_param = ( 2.0 / 3.0 ) / expansion_time;
    }
    else {
        config.scale_factor = 1.0;
        config.hubble_param = 0.0;
    }
}

void cic_mass_assignment( const std::vector<Particle>& particles, Grid& dm_rho, std::vector<CIC_Data>& cic_data, const Config& config ) {
    dm_rho.setZero();
    cic_data.assign( particles.size(), {} );

    for( size_t i = 0; i < particles.size(); ++i ) {
        const auto& p = particles[i];
        int ix = static_cast< int >( p.pos.x / config.CELL_SIZE );
        int iy = static_cast< int >( p.pos.y / config.CELL_SIZE );
        double frac_x = ( p.pos.x / config.CELL_SIZE ) - ix;
        double frac_y = ( p.pos.y / config.CELL_SIZE ) - iy;
        double w1 = ( 1 - frac_x ) * ( 1 - frac_y ), w2 = frac_x * ( 1 - frac_y );
        double w3 = ( 1 - frac_x ) * frac_y, w4 = frac_x * frac_y;
        cic_data[i] = { ix, iy, w1, w2, w3, w4 };

        dm_rho( (ix + config.MESH_SIZE) % config.MESH_SIZE, (iy + config.MESH_SIZE) % config.MESH_SIZE ) += p.mass * w1;
        dm_rho( ( ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, (iy + config.MESH_SIZE) % config.MESH_SIZE ) += p.mass * w2;
        dm_rho( (ix + config.MESH_SIZE) % config.MESH_SIZE, ( iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) += p.mass * w3;
        dm_rho( ( ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) += p.mass * w4;
    }
    dm_rho /= ( config.CELL_SIZE * config.CELL_SIZE );
}

void cic_force_interpolation( const std::vector<Particle>& particles, const Grid& ax_grid, const Grid& ay_grid, const std::vector<CIC_Data>& cic_data, std::vector<Vec2>& forces, const Config& config ) {
    forces.assign( particles.size(), { 0.0, 0.0 } );
    for( size_t i = 0; i < particles.size(); ++i ) {
        const auto& p = particles[i];
        const auto& cd = cic_data[i];
        double ax = ( ax_grid( (cd.ix + config.MESH_SIZE) % config.MESH_SIZE, (cd.iy + config.MESH_SIZE) % config.MESH_SIZE ) * cd.w1 +
            ax_grid( ( cd.ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, (cd.iy + config.MESH_SIZE) % config.MESH_SIZE ) * cd.w2 +
            ax_grid( (cd.ix + config.MESH_SIZE) % config.MESH_SIZE, ( cd.iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w3 +
            ax_grid( ( cd.ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w4 );
        double ay = ( ay_grid( (cd.ix + config.MESH_SIZE) % config.MESH_SIZE, (cd.iy + config.MESH_SIZE) % config.MESH_SIZE ) * cd.w1 +
            ay_grid( ( cd.ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, (cd.iy + config.MESH_SIZE) % config.MESH_SIZE ) * cd.w2 +
            ay_grid( (cd.ix + config.MESH_SIZE) % config.MESH_SIZE, ( cd.iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w3 +
            ay_grid( ( cd.ix + 1 + config.MESH_SIZE ) % config.MESH_SIZE, ( cd.iy + 1 + config.MESH_SIZE ) % config.MESH_SIZE ) * cd.w4 );
        forces[i].x = ax * p.mass;
        forces[i].y = ay * p.mass;
    }
}

void compute_gravitational_acceleration( std::vector<Particle>& particles, GasGrid& gas, Grid& g_grid_x, Grid& g_grid_y, std::vector<CIC_Data>& cic_data, const Config& config ) {
    Grid dm_rho( config.MESH_SIZE, config.MESH_SIZE );
    cic_mass_assignment( particles, dm_rho, cic_data, config );
    Grid total_rho = dm_rho + ( config.USE_HYDRO ? gas.density : Grid::Zero( config.MESH_SIZE, config.MESH_SIZE ) );

    pocketfft::shape_t shape = { ( size_t )config.MESH_SIZE, ( size_t )config.MESH_SIZE };
    pocketfft::stride_t stride_r = { static_cast< ptrdiff_t >( sizeof( double ) * config.MESH_SIZE ), static_cast< ptrdiff_t >( sizeof( double ) ) };
    pocketfft::stride_t stride_c = { static_cast< ptrdiff_t >( sizeof( std::complex<double> ) * ( config.MESH_SIZE / 2 + 1 ) ), static_cast< ptrdiff_t >( sizeof( std::complex<double> ) ) };

    std::vector<std::complex<double>> rho_k( config.MESH_SIZE * ( config.MESH_SIZE / 2 + 1 ) );
    pocketfft::r2c( shape, stride_r, stride_c, { 0, 1 }, true, total_rho.data(), rho_k.data(), 1.0 );

    std::vector<std::complex<double>> phi_k( config.MESH_SIZE * ( config.MESH_SIZE / 2 + 1 ) );
    for( int i = 0; i < config.MESH_SIZE; ++i ) {
        for( int j = 0; j < config.MESH_SIZE / 2 + 1; ++j ) {
            if( i == 0 && j == 0 ) {
                phi_k[0] = { 0.0, 0.0 };
                continue;
            }
            double kx_freq = ( i < config.MESH_SIZE / 2 ) ? i : ( i - config.MESH_SIZE );
            double kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
            double ky = ( double )j * 2.0 * M_PI / config.DOMAIN_SIZE;
            double k2 = kx * kx + ky * ky;

            phi_k[i * ( config.MESH_SIZE / 2 + 1 ) + j] = rho_k[i * ( config.MESH_SIZE / 2 + 1 ) + j] * ( -2.0 * M_PI * config.G / sqrt( k2 ) );
        }
    }

    std::vector<double> phi_flat( config.MESH_SIZE * config.MESH_SIZE, 0.0 );
    pocketfft::c2r( shape, stride_c, stride_r, { 0, 1 }, false, phi_k.data(), phi_flat.data(), 1.0 );

    Grid phi = Eigen::Map<const Grid>( phi_flat.data(), config.MESH_SIZE, config.MESH_SIZE );
    double norm = 1.0 / ( config.MESH_SIZE * config.MESH_SIZE );

    g_grid_x = ( roll( phi, 1, 0 ) - roll( phi, -1, 0 ) ) * norm / ( 2.0 * config.CELL_SIZE );
    g_grid_y = ( roll( phi, 1, 1 ) - roll( phi, -1, 1 ) ) * norm / ( 2.0 * config.CELL_SIZE );
}

void compute_PP_forces( std::vector<Particle>& particles, std::vector<Vec2>& pp_forces, const Config& config ) {
    pp_forces.assign( particles.size(), { 0.0, 0.0 } );

    std::vector<std::vector<int>> cell_grid( config.MESH_SIZE * config.MESH_SIZE );
    for( size_t i = 0; i < particles.size(); ++i ) {
        int ix = static_cast< int >( particles[i].pos.x / config.CELL_SIZE );
        int iy = static_cast< int >( particles[i].pos.y / config.CELL_SIZE );
        ix = ( ix % config.MESH_SIZE + config.MESH_SIZE ) % config.MESH_SIZE;
        iy = ( iy % config.MESH_SIZE + config.MESH_SIZE ) % config.MESH_SIZE;
        int cell_index = iy * config.MESH_SIZE + ix;
        cell_grid[cell_index].push_back( static_cast< int >( i ) );
    }

    const int search_radius_cells = static_cast< int >( ceil( config.CUTOFF_RADIUS / config.CELL_SIZE ) );

    for( size_t i = 0; i < particles.size(); ++i ) {
        auto& p1 = particles[i];
        int ix = static_cast< int >( p1.pos.x / config.CELL_SIZE );
        int iy = static_cast< int >( p1.pos.y / config.CELL_SIZE );
        ix = ( ix % config.MESH_SIZE + config.MESH_SIZE ) % config.MESH_SIZE;
        iy = ( iy % config.MESH_SIZE + config.MESH_SIZE ) % config.MESH_SIZE;

        for( int dx_cell = -search_radius_cells; dx_cell <= search_radius_cells; ++dx_cell ) {
            for( int dy_cell = -search_radius_cells; dy_cell <= search_radius_cells; ++dy_cell ) {
                int neighbor_ix = ( ix + dx_cell + config.MESH_SIZE ) % config.MESH_SIZE;
                int neighbor_iy = ( iy + dy_cell + config.MESH_SIZE ) % config.MESH_SIZE;
                int neighbor_cell_index = neighbor_iy * config.MESH_SIZE + neighbor_ix;

                for( int j : cell_grid[neighbor_cell_index] ) {
                    if( i >= j ) continue;
                    auto& p2 = particles[j];

                    double dx = displacement( p2.pos.x - p1.pos.x, config );
                    double dy = displacement( p2.pos.y - p1.pos.y, config );
                    double dist_sq = dx * dx + dy * dy;

                    if( config.USE_PM && dist_sq > config.CUTOFF_RADIUS_SQUARED ) continue;

                    double S = 1.0;
                    if( config.USE_PM && dist_sq > config.R_SWITCH_START_SQ ) {
                        double dist = sqrt( dist_sq );
                        double x = ( dist - config.R_SWITCH_START ) / config.CUTOFF_TRANSITION_WIDTH;
                        S = 2 * pow( x, 3 ) - 3 * pow( x, 2 ) + 1;
                    }

                    double soft_dist_sq = dist_sq + pow( 0.5 * config.CELL_SIZE, 2 );
                    double f_pm_short = config.G * p1.mass * p2.mass / soft_dist_sq;
                    double soft_dist = sqrt( soft_dist_sq );
                    Vec2 f_pm_short_vec = { f_pm_short * dx / soft_dist, f_pm_short * dy / soft_dist };

                    double pp_dist_sq = dist_sq + config.SOFTENING_SQUARED;
                    double f_pp = config.G * p1.mass * p2.mass / pp_dist_sq;
                    double pp_dist = sqrt( pp_dist_sq );
                    Vec2 f_pp_vec = { f_pp * dx / pp_dist, f_pp * dy / pp_dist };

                    Vec2 correction_f;
                    if( !config.USE_PM ) {
                        f_pm_short_vec.x = 0;
                        f_pm_short_vec.y = 0;
                    }
                    correction_f.x = S * ( f_pp_vec.x - f_pm_short_vec.x );
                    correction_f.y = S * ( f_pp_vec.y - f_pm_short_vec.y );

                    pp_forces[i].x += correction_f.x;
                    pp_forces[i].y += correction_f.y;
                    pp_forces[j].x -= correction_f.x;
                    pp_forces[j].y -= correction_f.y;
                }
            }
        }
    }
}

void create_zeldovich_ics( std::vector<Particle>& particles, const Config& config ) {
    int n_per_side = config.N_PER_SIDE;
    int ic_mesh_size = n_per_side;
    double cell_size = config.DOMAIN_SIZE / ic_mesh_size;

    std::vector<double> real_space_random_field( ic_mesh_size * ic_mesh_size );
    std::default_random_engine generator( config.SEED );
    std::normal_distribution<double> distribution( 0.0, 1.0 );
    for( auto& val : real_space_random_field ) {
        val = distribution( generator );
    }

    pocketfft::shape_t shape_ic = { ( size_t )ic_mesh_size, ( size_t )ic_mesh_size };
    pocketfft::stride_t stride_r_ic = { static_cast< ptrdiff_t >( sizeof( double ) * ic_mesh_size ), sizeof( double ) };
    pocketfft::stride_t stride_c_ic = { static_cast< ptrdiff_t >( sizeof( std::complex<double> ) ) * ( ic_mesh_size / 2 + 1 ), sizeof( std::complex<double> ) };

    std::vector<std::complex<double>> random_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );
    pocketfft::r2c( shape_ic, stride_r_ic, stride_c_ic, { 0, 1 }, true, real_space_random_field.data(), random_k.data(), 1.0 );

    std::vector<std::complex<double>> disp_x_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );
    std::vector<std::complex<double>> disp_y_k( ic_mesh_size * ( ic_mesh_size / 2 + 1 ) );

    for( int i = 0; i < ic_mesh_size; ++i ) {
        for( int j = 0; j < ic_mesh_size / 2 + 1; ++j ) {
            double kx_freq = ( i < ic_mesh_size / 2 ) ? i : ( i - ic_mesh_size );
            double ky_freq = j;
            double kx = kx_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
            double ky = ky_freq * 2.0 * M_PI / config.DOMAIN_SIZE;
            double k2 = kx * kx + ky * ky;
            int idx = i * ( ic_mesh_size / 2 + 1 ) + j;
            if( i == 0 && j == 0 ) {
                disp_x_k[idx] = { 0, 0 };
                disp_y_k[idx] = { 0, 0 };
                continue;
            }
            double power_spectrum_sqrt = sqrt( pow( k2, config.INITIAL_POWER_SPECTRUM_INDEX / 2.0 ) );
            std::complex<double> delta_k = random_k[idx] * power_spectrum_sqrt;
            std::complex<double> phi_k = -delta_k / k2;
            disp_x_k[idx] = std::complex<double>( 0, -1 ) * kx * phi_k;
            disp_y_k[idx] = std::complex<double>( 0, -1 ) * ky * phi_k;
        }
    }

    std::vector<double> disp_x_real( ic_mesh_size * ic_mesh_size );
    std::vector<double> disp_y_real( ic_mesh_size * ic_mesh_size );
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1 }, false, disp_x_k.data(), disp_x_real.data(), 1.0 );
    pocketfft::c2r( shape_ic, stride_c_ic, stride_r_ic, { 0, 1 }, false, disp_y_k.data(), disp_y_real.data(), 1.0 );

    double norm_ic = 1.0 / ( ic_mesh_size * ic_mesh_size );
    double std_x = 0, std_y = 0;
    for( size_t i = 0; i < disp_x_real.size(); ++i ) {
        disp_x_real[i] *= norm_ic;
        disp_y_real[i] *= norm_ic;
        std_x += disp_x_real[i] * disp_x_real[i];
        std_y += disp_y_real[i] * disp_y_real[i];
    }
    std_x = sqrt( std_x / disp_x_real.size() );
    std_y = sqrt( std_y / disp_y_real.size() );

    double initial_scale_factor = config.EXPANDING_UNIVERSE ? pow( config.EXPANSION_START_T, 2.0 / 3.0 ) : 0.5;
    double initial_hubble_param = config.EXPANDING_UNIVERSE ? ( 2.0 / 3.0 ) / config.EXPANSION_START_T : 0.0;

    particles.clear();
    double spacing = config.DOMAIN_SIZE / n_per_side;
    for( int i = 0; i < n_per_side; ++i ) {
        for( int j = 0; j < n_per_side; ++j ) {
            double qx = ( i + 0.5 ) * spacing;
            double qy = ( j + 0.5 ) * spacing;
            double dx = ( disp_x_real[i * n_per_side + j] / std_x ) * initial_scale_factor * cell_size;
            double dy = ( disp_y_real[i * n_per_side + j] / std_y ) * initial_scale_factor * cell_size;
            Particle p;
            p.pos.x = fmod( qx + dx + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
            p.pos.y = fmod( qy + dy + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
            p.vel.x = config.STANDING_PARTICLES ? 0 : initial_hubble_param * dx;
            p.vel.y = config.STANDING_PARTICLES ? 0 : initial_hubble_param * dy;
            p.mass = config.DM_PARTICLE_MASS;
            particles.push_back( p );
        }
    }
}


void KDK_step( double total_time, std::vector<Particle>& particles, GasGrid& gas, Config& config ) {
    update_cosmology( config, total_time );
    double a = config.scale_factor;
    double H = config.hubble_param;
    double a3 = a * a * a;
    double dt = config.DT;

    // KICK 1 (t -> t + dt/2)
    gas.update_primitive_variables();
    Grid total_ax_gas = ( g_grid_x.array() / a3 ) - ( 2 * H * gas.velocity_x.array() );
    Grid total_ay_gas = ( g_grid_y.array() / a3 ) - ( 2 * H * gas.velocity_y.array() );
    Grid g_mom_x_source = gas.density.array() * total_ax_gas.array();
    Grid g_mom_y_source = gas.density.array() * total_ay_gas.array();
    Grid power_density = gas.velocity_x.array() * g_mom_x_source.array() + gas.velocity_y.array() * g_mom_y_source.array();

    gas.momentum_x += g_mom_x_source * ( dt / 2.0 );
    gas.momentum_y += g_mom_y_source * ( dt / 2.0 );
    gas.energy += power_density * ( dt / 2.0 );

    for( auto& p : particles ) {
        double total_ax_p = ( p.acc.x / a3 ) - ( 2 * H * p.vel.x );
        double total_ay_p = ( p.acc.y / a3 ) - ( 2 * H * p.vel.y );
        p.vel.x += total_ax_p * dt / 2.0;
        p.vel.y += total_ay_p * dt / 2.0;
    }

    // DRIFT (t -> t + dt)
    for( auto& p : particles ) {
        p.pos.x = fmod( p.pos.x + p.vel.x * dt + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
        p.pos.y = fmod( p.pos.y + p.vel.y * dt + config.DOMAIN_SIZE, config.DOMAIN_SIZE );
    }

    if( config.USE_HYDRO ) {
        gas.hydro_step( dt );
    }

    // UPDATE COSMOLOGY to t + dt
    update_cosmology( config, total_time + dt );
    a = config.scale_factor;
    H = config.hubble_param;
    a3 = a * a * a;

    // COMPUTE FORCES at t + dt
    static Grid g_grid_x_new( config.MESH_SIZE, config.MESH_SIZE ), g_grid_y_new( config.MESH_SIZE, config.MESH_SIZE );
    std::vector<CIC_Data> cic_data;
    compute_gravitational_acceleration( particles, gas, g_grid_x_new, g_grid_y_new, cic_data, config );

    std::vector<Vec2> pp_forces, pm_forces;
    if( config.USE_PP ) { compute_PP_forces( particles, pp_forces, config ); }
    else { pp_forces.assign( particles.size(), { 0.0, 0.0 } ); }
    if( config.USE_PM ) { cic_force_interpolation( particles, g_grid_x_new, g_grid_y_new, cic_data, pm_forces, config ); }
    else { pm_forces.assign( particles.size(), { 0.0, 0.0 } ); }

    // Store new gravity field for next step's kick
    g_grid_x = g_grid_x_new;
    g_grid_y = g_grid_y_new;

    // KICK 2 (t + dt/2 -> t + dt)
    gas.update_primitive_variables();
    total_ax_gas = ( g_grid_x.array() / a3 ) - ( 2 * H * gas.velocity_x.array() );
    total_ay_gas = ( g_grid_y.array() / a3 ) - ( 2 * H * gas.velocity_y.array() );
    g_mom_x_source = gas.density.array() * total_ax_gas.array();
    g_mom_y_source = gas.density.array() * total_ay_gas.array();
    power_density = gas.velocity_x.array() * g_mom_x_source.array() + gas.velocity_y.array() * g_mom_y_source.array();

    gas.momentum_x += g_mom_x_source * ( dt / 2.0 );
    gas.momentum_y += g_mom_y_source * ( dt / 2.0 );
    gas.energy += power_density * ( dt / 2.0 );

    for( size_t i = 0; i < particles.size(); ++i ) {
        auto& p = particles[i];
        Vec2 f = { pp_forces[i].x + pm_forces[i].x, pp_forces[i].y + pm_forces[i].y };
        p.acc.x = config.STANDING_PARTICLES ? 0 : f.x / p.mass;
        p.acc.y = config.STANDING_PARTICLES ? 0 : f.y / p.mass;
        double total_ax_p = ( p.acc.x / a3 ) - ( 2 * H * p.vel.x );
        double total_ay_p = ( p.acc.y / a3 ) - ( 2 * H * p.vel.y );
        p.vel.x += total_ax_p * dt / 2.0;
        p.vel.y += total_ay_p * dt / 2.0;
    }
}

double calculate_total_energy( const std::vector<Particle>& particles, double a, const Config& config ) {
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;
    for( const auto& p : particles ) {
        double proper_vel_sq = ( a * p.vel.x ) * ( a * p.vel.x ) + ( a * p.vel.y ) * ( a * p.vel.y );
        kinetic_energy += 0.5 * p.mass * proper_vel_sq;
    }
    for( size_t i = 0; i < particles.size(); ++i ) {
        for( size_t j = i + 1; j < particles.size(); ++j ) {
            const auto& p1 = particles[i];
            const auto& p2 = particles[j];
            double dx = displacement( p2.pos.x - p1.pos.x, config );
            double dy = displacement( p2.pos.y - p1.pos.y, config );
            double proper_dist_sq = ( a * a ) * ( dx * dx + dy * dy + config.SOFTENING_SQUARED );
            if( proper_dist_sq > 0 ) {
                potential_energy -= config.G * p1.mass * p2.mass / sqrt( proper_dist_sq );
            }
        }
    }
    return kinetic_energy + potential_energy;
}

void save_frame( sf::RenderWindow& window, int frame_num ) {
    char filename[100];
    snprintf( filename, sizeof( filename ), "frame_%03d.png", frame_num );
    
    // --- SFML 2.x Change ---
    // sf::Texture texture( window.getSize() ); // sf::Texture constructor from sf::Vector2u is SFML 3
    sf::Texture texture;
    texture.create(window.getSize().x, window.getSize().y); // SFML 2.x way
    // --- End of Change ---

    texture.update( window );
    if( texture.copyToImage().saveToFile( filename ) ) {
        std::cout << "Saved " << filename << std::endl;
    }
}

sf::Color lerpColor( sf::Color c1, sf::Color c2, double t ) {
    t = std::max( 0.0, std::min( 1.0, t ) );
    return sf::Color(
        static_cast< uint8_t >( c1.r + ( c2.r - c1.r ) * t ),
        static_cast< uint8_t >( c1.g + ( c2.g - c1.g ) * t ),
        static_cast< uint8_t >( c1.b + ( c2.b - c1.b ) * t )
    );
}

sf::Color getPlasmaColor( double value ) {
    if( value < 0.25 ) return lerpColor( sf::Color( 0, 0, 0 ), sf::Color( 84, 2, 163 ), value / 0.25 );
    if( value < 0.5 ) return lerpColor( sf::Color( 84, 2, 163 ), sf::Color( 158, 28, 133 ), ( value - 0.25 ) / 0.25 );
    if( value < 0.75 ) return lerpColor( sf::Color( 158, 28, 133 ), sf::Color( 218, 107, 49 ), ( value - 0.5 ) / 0.25 );
    return lerpColor( sf::Color( 218, 107, 49 ), sf::Color( 255, 237, 200 ), ( value - 0.75 ) / 0.25 );
}

void render_field( sf::RenderWindow& window, const Grid& field, const Config& config ) {
    double min_f = field.minCoeff();
    double max_f = field.maxCoeff();
    double range_f = max_f - min_f;
    if( range_f < 1e-9 ) range_f = 1.0;

    sf::VertexArray field_mesh( sf::PrimitiveType::TriangleStrip, config.MESH_SIZE * config.MESH_SIZE * 4 );

    for( int i = 0; i < config.MESH_SIZE; ++i ) {
        for( int j = 0; j < config.MESH_SIZE; ++j ) {
            int i_next = ( i + 1 ); // No wrap for drawing
            int j_next = ( j + 1 ); // No wrap for drawing

            double val1 = ( field( i, j ) - min_f ) / range_f;
            sf::Color c1 = getPlasmaColor( val1 );
            sf::Vector2f pos1( static_cast< float >( i * config.CELL_SIZE * config.RENDER_SCALE ), static_cast< float >( j * config.CELL_SIZE * config.RENDER_SCALE ) );

            double val2 = ( j_next < config.MESH_SIZE ) ? ( field( i, j_next ) - min_f ) / range_f : val1;
            sf::Color c2 = ( j_next < config.MESH_SIZE ) ? getPlasmaColor( val2 ) : c1;
            sf::Vector2f pos2( static_cast< float >( i * config.CELL_SIZE * config.RENDER_SCALE ), static_cast< float >( j_next * config.CELL_SIZE * config.RENDER_SCALE ) );

            double val3 = ( i_next < config.MESH_SIZE && j_next < config.MESH_SIZE ) ? ( field( i_next, j_next ) - min_f ) / range_f : val1;
            sf::Color c3 = ( i_next < config.MESH_SIZE && j_next < config.MESH_SIZE ) ? getPlasmaColor( val3 ) : c1;
            sf::Vector2f pos3( static_cast< float >( i_next * config.CELL_SIZE * config.RENDER_SCALE ), static_cast< float >( j_next * config.CELL_SIZE * config.RENDER_SCALE ) );

            double val4 = ( i_next < config.MESH_SIZE ) ? ( field( i_next, j ) - min_f ) / range_f : val1;
            sf::Color c4 = ( i_next < config.MESH_SIZE ) ? getPlasmaColor( val4 ) : c1;
            sf::Vector2f pos4( static_cast< float >( i_next * config.CELL_SIZE * config.RENDER_SCALE ), static_cast< float >( j * config.CELL_SIZE * config.RENDER_SCALE ) );

            size_t idx = ( i * config.MESH_SIZE + j ) * 4;
            field_mesh[idx + 0].position = pos1; field_mesh[idx + 0].color = c1;
            field_mesh[idx + 1].position = pos2; field_mesh[idx + 1].color = c2;
            field_mesh[idx + 2].position = pos4; field_mesh[idx + 2].color = c4;
            field_mesh[idx + 3].position = pos3; field_mesh[idx + 3].color = c3;
        }
    }
    window.draw( field_mesh );
}

void render( sf::RenderWindow& window, const std::vector<Particle>& particles, const GasGrid& gas, const Config& config ) {
    window.clear( sf::Color::Black );
    if( config.USE_HYDRO ) {
        render_field( window, gas.pressure, config );
    }
    sf::CircleShape particle_shape( 1.0f );
    particle_shape.setFillColor( sf::Color::White );
    for( const auto& p : particles ) {
        particle_shape.setPosition( { static_cast< float >( p.pos.x * config.RENDER_SCALE ),
            static_cast< float >( p.pos.y * config.RENDER_SCALE ) } );
        window.draw( particle_shape );
    }
    window.display();
}

// HDF5 Helper Functions
void set_hdf5_attr_double( H5::H5Object& obj, const char* attr_name, double value ) {
    H5::DataSpace scalar_space( H5S_SCALAR );
    H5::Attribute attr = obj.createAttribute( attr_name, H5::PredType::NATIVE_DOUBLE, scalar_space );
    attr.write( H5::PredType::NATIVE_DOUBLE, &value );
    attr.close();
}
void set_hdf5_attr_int( H5::H5Object& obj, const char* attr_name, int value ) {
    H5::DataSpace scalar_space( H5S_SCALAR );
    H5::Attribute attr = obj.createAttribute( attr_name, H5::PredType::NATIVE_INT, scalar_space );
    attr.write( H5::PredType::NATIVE_INT, &value );
    attr.close();
}
void set_hdf5_attr_bool( H5::H5Object& obj, const char* attr_name, bool value ) {
    int int_val = value ? 1 : 0;
    set_hdf5_attr_int( obj, attr_name, int_val );
}

void write_grid_to_hdf5( H5::Group& group, const char* dataset_name, const Grid& grid ) {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row_major_grid = grid;
    hsize_t dims[2] = { ( hsize_t )grid.rows(), ( hsize_t )grid.cols() };
    H5::DataSpace dataspace( 2, dims );
    H5::DataSet dataset = group.createDataSet( dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace );
    dataset.write( row_major_grid.data(), H5::PredType::NATIVE_DOUBLE );
    dataset.close();
}

void write_particle_vec_to_hdf5( H5::Group& group, const char* dataset_name, const std::vector<double>& vec ) {
    hsize_t dims[1] = { vec.size() };
    H5::DataSpace dataspace( 1, dims );
    H5::DataSet dataset = group.createDataSet( dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace );
    dataset.write( vec.data(), H5::PredType::NATIVE_DOUBLE );
    dataset.close();
}

void save_hdf5_snapshot( H5::H5File& file, int snapshot_index, double sim_time, double a,
    const std::vector<Particle>& particles, const GasGrid& gas, const Config& config ) {
    char group_name[100];
    snprintf( group_name, sizeof( group_name ), "snapshot_%04d", snapshot_index );
    H5::Group snapshot_group = file.createGroup( group_name );

    set_hdf5_attr_double( snapshot_group, "simulation_time", sim_time );
    set_hdf5_attr_double( snapshot_group, "scale_factor", a );

    H5::Group particle_group = snapshot_group.createGroup( "particles" );
    size_t n_particles = particles.size();
    std::vector<double> pos_x( n_particles ), pos_y( n_particles ), vel_x( n_particles ), vel_y( n_particles ), acc_x( n_particles ), acc_y( n_particles ), mass( n_particles );
    for( size_t i = 0; i < n_particles; ++i ) {
        pos_x[i] = particles[i].pos.x; pos_y[i] = particles[i].pos.y;
        vel_x[i] = particles[i].vel.x; vel_y[i] = particles[i].vel.y;
        acc_x[i] = particles[i].acc.x; acc_y[i] = particles[i].acc.y;
        mass[i] = particles[i].mass;
    }
    write_particle_vec_to_hdf5( particle_group, "position_x", pos_x );
    write_particle_vec_to_hdf5( particle_group, "position_y", pos_y );
    write_particle_vec_to_hdf5( particle_group, "velocity_x", vel_x );
    write_particle_vec_to_hdf5( particle_group, "velocity_y", vel_y );
    write_particle_vec_to_hdf5( particle_group, "acceleration_x", acc_x );
    write_particle_vec_to_hdf5( particle_group, "acceleration_y", acc_y );
    write_particle_vec_to_hdf5( particle_group, "mass", mass );
    particle_group.close();

    if( config.USE_HYDRO ) {
        H5::Group gas_group = snapshot_group.createGroup( "gas" );
        write_grid_to_hdf5( gas_group, "density", gas.density );
        write_grid_to_hdf5( gas_group, "momentum_x", gas.momentum_x );
        write_grid_to_hdf5( gas_group, "momentum_y", gas.momentum_y );
        write_grid_to_hdf5( gas_group, "energy", gas.energy );
        write_grid_to_hdf5( gas_group, "pressure", gas.pressure );
        gas_group.close();
    }
    snapshot_group.close();
}

void initialize_hdf5_file( H5::H5File& file, const Config& config ) {
    std::time_t now = std::time( nullptr );
    std::tm ltm;
#ifdef _MSC_VER
    localtime_s( &ltm, &now );
#else
    ltm = *std::localtime( &now );
#endif
    std::stringstream ss;
    ss << "sim_" << std::put_time( &ltm, "%Y-%m-%d_%H-%M-%S" ) << ".hdf5";
    std::string filename_str = ss.str();
    const char* filename = filename_str.c_str();
    std::cout << "Creating HDF5 output file: " << filename << std::endl;

    try {
        file = H5::H5File( filename, H5F_ACC_TRUNC );
        H5::Group root_group = file.openGroup( "/" );

        // Save all config parameters
        set_hdf5_attr_double( root_group, "domain_size", config.DOMAIN_SIZE );
        set_hdf5_attr_int( root_group, "mesh_size", config.MESH_SIZE );
        set_hdf5_attr_double( root_group, "omega_baryon", config.OMEGA_BARYON );
        set_hdf5_attr_int( root_group, "n_per_side", config.N_PER_SIDE );
        set_hdf5_attr_bool( root_group, "use_hydro", config.USE_HYDRO );
        set_hdf5_attr_double( root_group, "g_const", config.G );
        set_hdf5_attr_double( root_group, "gamma", config.GAMMA );
        set_hdf5_attr_bool( root_group, "standing_particles", config.STANDING_PARTICLES );
        set_hdf5_attr_bool( root_group, "expanding_universe", config.EXPANDING_UNIVERSE );
        set_hdf5_attr_double( root_group, "expansion_start_t", config.EXPANSION_START_T );
        set_hdf5_attr_double( root_group, "initial_power_spectrum_index", config.INITIAL_POWER_SPECTRUM_INDEX );
        set_hdf5_attr_bool( root_group, "use_pm", config.USE_PM );
        set_hdf5_attr_bool( root_group, "use_pp", config.USE_PP );
        set_hdf5_attr_double( root_group, "cutoff_radius_cells", config.CUTOFF_RADIUS_CELLS );
        set_hdf5_attr_double( root_group, "dt_factor", config.DT_FACTOR );
        set_hdf5_attr_int( root_group, "seed", config.SEED );

        root_group.close();
    }
    catch( H5::Exception& e ) {
        std::cerr << "Error: Could not create HDF5 file." << std::endl;
        e.printErrorStack();
    }
}


int main() {
    // Load Config
    Config config;
    try {
        config.load( "simulation.ini" );
    }
    catch( const std::exception& e ) {
        std::cerr << "Error loading simulation.ini: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Successfully loaded simulation.ini" << std::endl;
    std::cout << "Mesh size: " << config.MESH_SIZE << "x" << config.MESH_SIZE << std::endl;
    std::cout << "Particle count: " << config.NUM_DM_PARTICLES << std::endl;

    // Initialize Global Grids
    g_grid_x = Grid( config.MESH_SIZE, config.MESH_SIZE );
    g_grid_y = Grid( config.MESH_SIZE, config.MESH_SIZE );
    g_grid_x.setZero();
    g_grid_y.setZero();

    // Initialize SFML Window
    sf::RenderWindow window( sf::VideoMode( config.RENDER_SIZE, config.RENDER_SIZE ), "N-Body + Hydro Simulation" );

    // Initialize Simulation State
    std::vector<Particle> particles;
    create_zeldovich_ics( particles, config );

    GasGrid gas( config );
    if( config.USE_HYDRO ) {
        gas.density.fill( config.GAS_TOTAL_MASS / ( config.DOMAIN_SIZE * config.DOMAIN_SIZE ) );
        double initial_internal_energy = 1e-6;
        gas.energy = gas.density.array() * initial_internal_energy;
        gas.update_primitive_variables();
    }

    std::vector<CIC_Data> cic_data;
    std::vector<Vec2> pp_forces, pm_forces;

    // Compute Initial Forces
    compute_gravitational_acceleration( particles, gas, g_grid_x, g_grid_y, cic_data, config );
    if( config.USE_PM ) { cic_force_interpolation( particles, g_grid_x, g_grid_y, cic_data, pm_forces, config ); }
    else { pm_forces.assign( particles.size(), { 0.0, 0.0 } ); }
    if( config.USE_PP ) { compute_PP_forces( particles, pp_forces, config ); }
    else { pp_forces.assign( particles.size(), { 0.0, 0.0 } ); }

    for( size_t i = 0; i < particles.size(); ++i ) {
        particles[i].acc.x = ( pp_forces[i].x + pm_forces[i].x ) / particles[i].mass;
        particles[i].acc.y = ( pp_forces[i].y + pm_forces[i].y ) / particles[i].mass;
    }

    // Initialize HDF5 File
    H5::H5File hdf5_file;
    initialize_hdf5_file( hdf5_file, config );
    int hdf5_snapshot_count = 0;

    // Main Loop
    int cycle_count = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    update_cosmology( config, 0 );

    double initial_energy = calculate_total_energy( particles, config.scale_factor, config );

    std::cout << "Starting simulation..." << std::endl;
    while( window.isOpen() ) {
        
        sf::Event event;
        while( window.pollEvent(event) ) {
            if( event.type == sf::Event::Closed ) {
                window.close();
            }
        }

        double total_simulation_time = config.DT * cycle_count;
        KDK_step( total_simulation_time, particles, gas, config );

        render( window, particles, gas, config );

        if( config.DEBUG_INFO_EVERY_CYCLES > 0 && cycle_count % config.DEBUG_INFO_EVERY_CYCLES == 0 ) {
            update_cosmology( config, total_simulation_time );
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed_seconds = std::chrono::duration_cast< std::chrono::duration<double> >( now - start_time ).count();
            double energy = calculate_total_energy( particles, config.scale_factor, config );
            double energy_error = 100.0 * fabs( energy - initial_energy ) / fabs( initial_energy );

            std::cout << "Cycles: " << cycle_count
                << " | SimTime: " << std::fixed << std::setprecision( 4 ) << total_simulation_time
                << " | ScaleFactor: " << std::fixed << std::setprecision( 3 ) << config.scale_factor
                << " | Energy Err: " << std::fixed << std::setprecision( 2 ) << energy_error << "%"
                << " | Wall Time: " << static_cast< int >( elapsed_seconds ) << "s"
                << std::endl;
        }

        if( config.SAVE_RENDER_EVERY_CYCLES > 0 && cycle_count > 0 && cycle_count % config.SAVE_RENDER_EVERY_CYCLES == 0 ) {
            int frame_num = cycle_count / config.SAVE_RENDER_EVERY_CYCLES;
            save_frame( window, frame_num );
        }

        if( config.SAVE_HDF5_EVERY_CYCLES > 0 && cycle_count > 0 && cycle_count % config.SAVE_HDF5_EVERY_CYCLES == 0 ) {
            try {
                save_hdf5_snapshot( hdf5_file, hdf5_snapshot_count, total_simulation_time, config.scale_factor, particles, gas, config );
                hdf5_snapshot_count++;
                std::cout << "Saved HDF5 snapshot " << hdf5_snapshot_count << std::endl;
            }
            catch( H5::Exception& e ) {
                std::cerr << "Error: Could not save HDF5 snapshot." << std::endl;
                e.printErrorStack();
            }
        }

        cycle_count++;
    }

    hdf5_file.close();

    return 0;
}


