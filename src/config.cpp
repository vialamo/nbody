#include "config.h"

#include <cmath>

#include "ini.h"
#include "math_utils.h"

Config::Config() { compute_derived_data(); }

void Config::load(const std::string& filename) {
    ini::Ini config_file;
    config_file.load(filename);

    DOMAIN_SIZE = config_file.get_double("domain", "domain_size", DOMAIN_SIZE);
    MESH_SIZE = config_file.get_int("domain", "mesh_size", MESH_SIZE);
    BOX_SIZE_MPC =
        config_file.get_double("domain", "box_size_mpc", BOX_SIZE_MPC);

    OMEGA_BARYON =
        config_file.get_double("cosmology", "omega_baryon", OMEGA_BARYON);
    OMEGA_M = config_file.get_double("cosmology", "omega_M", OMEGA_M);
    OMEGA_LAMBDA =
        config_file.get_double("cosmology", "omega_lambda", OMEGA_LAMBDA);
    HUBBLE_PARAM =
        config_file.get_double("cosmology", "hubble_param", HUBBLE_PARAM);
    EXPANDING_UNIVERSE = config_file.get_bool("cosmology", "expanding_universe",
                                              EXPANDING_UNIVERSE);

    N_PER_SIDE =
        config_file.get_int("initial_conditions", "n_per_side", N_PER_SIDE);
    STANDING_PARTICLES = config_file.get_bool(
        "initial_conditions", "standing_particles", STANDING_PARTICLES);
    START_A = config_file.get_double("initial_conditions", "start_a", START_A);
    SIGMA_8 = config_file.get_double("initial_conditions", "sigma_8", SIGMA_8);
    SPECTRAL_INDEX = config_file.get_double("initial_conditions",
                                            "spectral_index", SPECTRAL_INDEX);
    INITIAL_GAS_TEMPERATURE_K = config_file.get_double(
        "initial_conditions", "initial_gas_temp_k", INITIAL_GAS_TEMPERATURE_K);
    SEED = config_file.get_int("initial_conditions", "seed", SEED);

    USE_HYDRO = config_file.get_bool("hydro", "use_hydro", USE_HYDRO);
    GAMMA = config_file.get_double("hydro", "gamma", GAMMA);

    USE_PM = config_file.get_bool("p3m", "use_pm", USE_PM);
    USE_PP = config_file.get_bool("p3m", "use_pp", USE_PP);
    CUTOFF_RADIUS_CELLS = config_file.get_double("p3m", "cutoff_radius_cells",
                                                 CUTOFF_RADIUS_CELLS);
    CUTOFF_TRANSITION_WIDTH_FACTOR =
        config_file.get_double("p3m", "cutoff_transition_width_factor",
                               CUTOFF_TRANSITION_WIDTH_FACTOR);

    DT_FACTOR = config_file.get_double("time", "dt_factor", DT_FACTOR);
    MAX_SCALE_FACTOR =
        config_file.get_double("time", "max_scale_factor", MAX_SCALE_FACTOR);
    CFL_SAFETY_FACTOR =
        config_file.get_double("time", "cfl_safety_factor", CFL_SAFETY_FACTOR);
    GRAVITY_DT_FACTOR =
        config_file.get_double("time", "gravity_dt_factor", GRAVITY_DT_FACTOR);
    USE_ADAPTIVE_DT =
        config_file.get_bool("time", "use_adaptive_dt", USE_ADAPTIVE_DT);
    MAX_CYCLES = config_file.get_int("time", "max_cycles", MAX_CYCLES);

    SAVE_HDF5_EVERY_DELTA_A = config_file.get_double(
        "output", "save_hdf5_every_delta_a", SAVE_HDF5_EVERY_DELTA_A);
    DEBUG_INFO_EVERY_CYCLES = config_file.get_int(
        "output", "debug_info_every_cycles", DEBUG_INFO_EVERY_CYCLES);
    ENABLE_ENERGY_DIAGNOSTICS = config_file.get_bool(
        "output", "enable_energy_diagnostics", ENABLE_ENERGY_DIAGNOSTICS);

    compute_derived_data();
}

void Config::compute_derived_data() {
    CELL_SIZE = DOMAIN_SIZE / MESH_SIZE;
    CELL_VOLUME = CELL_SIZE * CELL_SIZE * CELL_SIZE;
    OMEGA_DM = OMEGA_M - OMEGA_BARYON;
    CUTOFF_RADIUS = CUTOFF_RADIUS_CELLS * CELL_SIZE;
    CUTOFF_RADIUS_SQUARED = CUTOFF_RADIUS * CUTOFF_RADIUS;
    CUTOFF_TRANSITION_WIDTH = CUTOFF_TRANSITION_WIDTH_FACTOR * CUTOFF_RADIUS;
    R_SWITCH_START = CUTOFF_RADIUS - CUTOFF_TRANSITION_WIDTH;
    R_SWITCH_START_SQ = R_SWITCH_START * R_SWITCH_START;

    NUM_DM_PARTICLES = N_PER_SIDE * N_PER_SIDE * N_PER_SIDE;

    const double total_mass = 1.0;
    G = std::pow(DOMAIN_SIZE, 3) / (6.0 * M_PI * total_mass);

    const double baryon_fraction = OMEGA_BARYON / OMEGA_M;
    const double dm_total_mass = total_mass * (1.0 - baryon_fraction);
    DM_PARTICLE_MASS = dm_total_mass / NUM_DM_PARTICLES;
    GAS_TOTAL_MASS = total_mass * baryon_fraction;

    double mean_interparticle_spacing =
        DOMAIN_SIZE / std::cbrt(static_cast<double>(NUM_DM_PARTICLES));

    SOFTENING_SQUARED = std::pow(mean_interparticle_spacing / 50.0, 2);
    const double dynamical_time = 1.0 / std::sqrt(G);
    FIXED_DT = DT_FACTOR * dynamical_time;
}