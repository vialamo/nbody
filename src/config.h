#pragma once
#include <string>

struct Config {
    double DOMAIN_SIZE = 1.0;
    int MESH_SIZE = 32;
    double BOX_SIZE_MPC = 16.0;

    double OMEGA_BARYON = 0.045;
    double OMEGA_M = 0.3;
    double OMEGA_LAMBDA = 0.7;
    double HUBBLE_PARAM = 0.7;
    bool EXPANDING_UNIVERSE = true;

    int N_PER_SIDE = 32;
    bool STANDING_PARTICLES = false;
    double START_A = 0.02;
    double SIGMA_8 = 0.81;
    double SPECTRAL_INDEX = 0.96;
    double INITIAL_GAS_TEMPERATURE_K = 50.0;
    int SEED = 42;

    bool USE_HYDRO = true;
    double GAMMA = 5.0 / 3.0;

    bool USE_PM = true;
    bool USE_PP = true;
    double CUTOFF_RADIUS_CELLS = 2.5;
    double CUTOFF_TRANSITION_WIDTH_FACTOR = 0.2;

    double DT_FACTOR = 1e-3;
    double MAX_SCALE_FACTOR = 1.0;
    double CFL_SAFETY_FACTOR = 0.5;
    double GRAVITY_DT_FACTOR = 0.2;
    bool USE_ADAPTIVE_DT = true;
    int MAX_CYCLES = 1000000000;

    double SAVE_HDF5_EVERY_DELTA_A = 0.005;
    int DEBUG_INFO_EVERY_CYCLES = 40;
    bool ENABLE_ENERGY_DIAGNOSTICS = true;

    // Derived Parameters
    double G = 0.0;
    double CELL_SIZE = 0.0, CELL_VOLUME = 0.0, OMEGA_DM = 0.0;
    double CUTOFF_RADIUS = 0.0, CUTOFF_RADIUS_SQUARED = 0.0;
    double CUTOFF_TRANSITION_WIDTH = 0.0, R_SWITCH_START = 0.0, R_SWITCH_START_SQ = 0.0;
    int NUM_DM_PARTICLES = 0;
    double DM_PARTICLE_MASS = 0.0, GAS_TOTAL_MASS = 0.0;
    double SOFTENING_SQUARED = 0.0;
    double FIXED_DT = 0.0;

    Config();
    void load(const std::string& filename);
    void compute_derived_data();
};