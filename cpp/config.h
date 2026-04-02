#pragma once
#include <string>
#include <cmath>
#include "ini.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Config {
    double DOMAIN_SIZE;
    int MESH_SIZE;
    double OMEGA_BARYON;
    int N_PER_SIDE;
    bool USE_HYDRO;
    double G;
    double GAMMA;
    bool STANDING_PARTICLES;
    bool EXPANDING_UNIVERSE;
    double EXPANSION_START_T;
    double INITIAL_POWER_SPECTRUM_INDEX;
    bool USE_PM;
    bool USE_PP;
    double CUTOFF_RADIUS_CELLS;
    double CUTOFF_TRANSITION_WIDTH_FACTOR;
    double DT_FACTOR;
    double CFL_SAFETY_FACTOR;
    double GRAVITY_DT_FACTOR;
    bool USE_ADAPTIVE_DT;
    int SAVE_HDF5_EVERY_CYCLES;
    int DEBUG_INFO_EVERY_CYCLES;
    int SAVE_RENDER_EVERY_CYCLES;
    int SEED;
    unsigned int RENDER_SIZE;
    bool ENABLE_VISUALIZATION;
    int MAX_CYCLES;

    // Derived Parameters
    double CELL_SIZE, CELL_VOLUME, OMEGA_DM;
    double CUTOFF_RADIUS, CUTOFF_RADIUS_SQUARED;
    double CUTOFF_TRANSITION_WIDTH, R_SWITCH_START, R_SWITCH_START_SQ;
    int NUM_DM_PARTICLES;
    double TOTAL_MASS, DM_PARTICLE_MASS, GAS_TOTAL_MASS;
    double MEAN_INTERPARTICLE_SPACING, SOFTENING_SQUARED;
    double DYNAMICAL_TIME, FIXED_DT, RENDER_SCALE;

    void load( const std::string& filename ) {
        ini::Ini config_file;
        config_file.load( filename );

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
        CFL_SAFETY_FACTOR = config_file.get_double( "time", "cfl_safety_factor", 0.5 );
        GRAVITY_DT_FACTOR = config_file.get_double( "time", "gravity_dt_factor", 0.2 );
        USE_ADAPTIVE_DT = config_file.get_bool( "time", "use_adaptive_dt", true );
        MAX_CYCLES = config_file.get_int( "time", "max_cycles", 1000000000 );
        SAVE_HDF5_EVERY_CYCLES = config_file.get_int( "output", "save_hdf5_every_cycles", 100 );
        DEBUG_INFO_EVERY_CYCLES = config_file.get_int( "output", "debug_info_every_cycles", 40 );
        SAVE_RENDER_EVERY_CYCLES = config_file.get_int( "output", "save_render_every_cycles", 0 );
        SEED = config_file.get_int( "output", "seed", 42 );
        RENDER_SIZE = config_file.get_int( "visualization", "render_size", 800 );
        ENABLE_VISUALIZATION = config_file.get_bool( "visualization", "enable_visualization", true );

        CELL_SIZE = DOMAIN_SIZE / MESH_SIZE;
        CELL_VOLUME = CELL_SIZE * CELL_SIZE;
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
        MEAN_INTERPARTICLE_SPACING = DOMAIN_SIZE / sqrt( ( double )NUM_DM_PARTICLES );
        SOFTENING_SQUARED = pow( MEAN_INTERPARTICLE_SPACING / 50.0, 2 );
        DYNAMICAL_TIME = 1.0 / sqrt( G );
        FIXED_DT = DT_FACTOR * DYNAMICAL_TIME;
        RENDER_SCALE = ( double )RENDER_SIZE / DOMAIN_SIZE;
    }
};
