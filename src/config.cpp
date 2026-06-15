#include "config.h"

#include <cmath>

#include "constants.h"
#include "ini.h"
#include "math_utils.h"

static void init_derived_units(Config& config) {
    // Length
    config.UNIT_LENGTH_MPC = config.BOX_SIZE_MPC / config.DOMAIN_SIZE;

    // Time
    // Conversion factor to turn the base Hubble constant (100 km/s/Mpc) into
    // inverse Gigayears (1/Gyr).
    // 100 * (3.15576e16 s / Gyr) / (3.085677e19 km / Mpc) = 0.10227
    constexpr double km100_sMpc_to_gyr_inv = 0.10227;
    // Calculate physical H0 in 1/Gyr and dimensionless code H0
    // (matter-dominated)
    double H0_phys_gyr_inv = km100_sMpc_to_gyr_inv * config.HUBBLE_PARAM;
    double H0_code = 2.0 / (3.0 * std::sqrt(config.OMEGA_M));
    // Set how many physical Gigayears correspond to 1.0 unit of code time
    config.UNIT_TIME_GYR = H0_code / H0_phys_gyr_inv;

    // Velocity
    // Velocity is Length / Time. When dividing UNIT_LENGTH_MPC by
    // UNIT_TIME_GYR, the 0.10227 conversion mathematically cancels out.
    // We are left with the base 100 km/s/Mpc multiplied by the inverse of the
    // 2/3 scaling factor defined in H0_code. (100 * 3/2 = 150.0).
    constexpr double hubble_velocity_scaling = 150;
    config.UNIT_VELOCITY_KMS = config.UNIT_LENGTH_MPC *
                               hubble_velocity_scaling * config.HUBBLE_PARAM *
                               std::sqrt(config.OMEGA_M);
    config.UNIT_VELOCITY_CGS = config.UNIT_VELOCITY_KMS * 1e5;

    // Mass
    // Base critical density of the universe: (3 * H0^2) / (8 * pi * G)
    // Evaluated for a baseline H0 = 100 km/s/Mpc, and converted into
    // units of Solar Masses per cubic Megaparsec (M_sun / Mpc^3).
    constexpr double rho_crit_base_msun_mpc3 = 2.775e11;
    double rho_crit =
        rho_crit_base_msun_mpc3 * (config.HUBBLE_PARAM * config.HUBBLE_PARAM);
    double box_volume_mpc3 = std::pow(config.BOX_SIZE_MPC, 3);
    double total_mass_msun = config.OMEGA_M * rho_crit * box_volume_mpc3;
    config.UNIT_MASS_MSUN = total_mass_msun / Config::TOTAL_MASS;

    // Thermodynamics derived factors
    // T = u * (gamma - 1) * mu * m_p / k_B
    double V_unit_m_s = config.UNIT_VELOCITY_KMS * 1000.0;
    config.FACTOR_U_TO_T = (V_unit_m_s * V_unit_m_s) * (config.GAMMA - 1.0) *
                           config.PRIMORDIAL_MU * constants::M_P_SI /
                           constants::K_B_SI;
    config.FACTOR_T_TO_U = 1.0 / config.FACTOR_U_TO_T;

    // Density conversion: Code -> CGS (g/cm^3)
    double mass_cgs = config.UNIT_MASS_MSUN * constants::MSUN_CGS;
    double length_cgs = config.UNIT_LENGTH_MPC * constants::MPC_CGS;
    double vol_cgs = length_cgs * length_cgs * length_cgs;
    config.UNIT_DENSITY_CGS = mass_cgs / vol_cgs;

    // Cooling rate conversion: (erg/g/s) -> Code (u_code / t_code)
    double u_unit_cgs = config.UNIT_VELOCITY_CGS * config.UNIT_VELOCITY_CGS;
    double t_unit_cgs = config.UNIT_TIME_GYR * constants::GYR_CGS;
    config.COOLING_CONVERSION_FACTOR = t_unit_cgs / u_unit_cgs;
}

void Config::compute_derived_data() {
    CELL_SIZE = DOMAIN_SIZE / MESH_SIZE;
    CELL_VOLUME = CELL_SIZE * CELL_SIZE * CELL_SIZE;
    OMEGA_DM = OMEGA_M - OMEGA_BARYON;
    CUTOFF_RADIUS_CELLS = CUTOFF_RADIUS_FACTOR * PM_SMOOTHING_CELLS;
    CUTOFF_RADIUS = CUTOFF_RADIUS_CELLS * CELL_SIZE;
    CUTOFF_RADIUS_SQUARED = CUTOFF_RADIUS * CUTOFF_RADIUS;

    NUM_DM_PARTICLES = N_PER_SIDE * N_PER_SIDE * N_PER_SIDE;

    G = std::pow(DOMAIN_SIZE, 3) / (6.0 * M_PI * TOTAL_MASS);

    const double baryon_fraction = OMEGA_BARYON / OMEGA_M;
    const double dm_total_mass = TOTAL_MASS * (1.0 - baryon_fraction);
    DM_PARTICLE_MASS = dm_total_mass / NUM_DM_PARTICLES;
    GAS_TOTAL_MASS = TOTAL_MASS * baryon_fraction;

    double mean_interparticle_spacing =
        DOMAIN_SIZE / std::cbrt(static_cast<double>(NUM_DM_PARTICLES));

    SOFTENING_SQUARED = std::pow(mean_interparticle_spacing / 50.0, 2);
    const double dynamical_time = 1.0 / std::sqrt(G);
    FIXED_DT = DT_FACTOR * dynamical_time;

    init_derived_units(*this);
}

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
    ENABLE_COOLING =
        config_file.get_bool("hydro", "enable_cooling", ENABLE_COOLING);
    PRIMORDIAL_MU =
        config_file.get_double("hydro", "primordial_mu", PRIMORDIAL_MU);
    TEMP_FLOOR_KELVIN =
        config_file.get_double("hydro", "temp_floor_k", TEMP_FLOOR_KELVIN);

    USE_PM = config_file.get_bool("p3m", "use_pm", USE_PM);
    USE_PP = config_file.get_bool("p3m", "use_pp", USE_PP);
    CUTOFF_RADIUS_FACTOR = config_file.get_double("p3m", "cutoff_radius_factor",
                                                 CUTOFF_RADIUS_FACTOR);
    PM_SMOOTHING_CELLS =
        config_file.get_double("p3m", "pm_smoothing_cells", PM_SMOOTHING_CELLS);

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
    DEBUG_INFO_EVERY_SECONDS = config_file.get_int(
        "output", "debug_info_every_seconds", DEBUG_INFO_EVERY_SECONDS);
    ENABLE_ENERGY_DIAGNOSTICS = config_file.get_bool(
        "output", "enable_energy_diagnostics", ENABLE_ENERGY_DIAGNOSTICS);

    compute_derived_data();

    if (OMEGA_BARYON > OMEGA_M) {
        throw std::runtime_error(
            "Error: OMEGA_BARYON cannot be larger than total OMEGA_M.");
    }
}