#include "config.h"

#include <cmath>

#include "constants.h"
#include "ini.h"
#include "math_utils.h"

static void init_derived_units(Config& config) {
    // Length
    config.unit_length_mpc = config.box_size_mpc / config.domain_size;

    // Time
    // Conversion factor to turn the base Hubble constant (100 km/s/Mpc) into
    // inverse Gigayears (1/Gyr).
    // 100 * (3.15576e16 s / Gyr) / (3.085677e19 km / Mpc) = 0.10227
    constexpr double km100_sMpc_to_gyr_inv = 0.10227;
    // Calculate physical H0 in 1/Gyr and dimensionless code H0
    // (matter-dominated)
    double H0_phys_gyr_inv = km100_sMpc_to_gyr_inv * config.hubble_h;
    double H0_code = 2.0 / (3.0 * std::sqrt(config.omega_m));
    // Set how many physical Gigayears correspond to 1.0 unit of code time
    config.unit_time_gyr = H0_code / H0_phys_gyr_inv;

    // Velocity
    // Velocity is Length / Time. When dividing unit_length_mpc by
    // unit_time_gyr, the 0.10227 conversion mathematically cancels out.
    // We are left with the base 100 km/s/Mpc multiplied by the inverse of the
    // 2/3 scaling factor defined in H0_code. (100 * 3/2 = 150.0).
    constexpr double hubble_velocity_scaling = 150;
    config.unit_velocity_kms = config.unit_length_mpc *
                               hubble_velocity_scaling * config.hubble_h *
                               std::sqrt(config.omega_m);
    config.unit_velocity_cgs = config.unit_velocity_kms * 1e5;

    // Mass
    // Base critical density of the universe: (3 * H0^2) / (8 * pi * G)
    // Evaluated for a baseline H0 = 100 km/s/Mpc, and converted into
    // units of Solar Masses per cubic Megaparsec (M_sun / Mpc^3).
    constexpr double rho_crit_base_msun_mpc3 = 2.775e11;
    double rho_crit =
        rho_crit_base_msun_mpc3 * (config.hubble_h * config.hubble_h);
    double box_volume_mpc3 = std::pow(config.box_size_mpc, 3);
    double total_mass_msun = config.omega_m * rho_crit * box_volume_mpc3;
    config.unit_mass_msun = total_mass_msun / Config::total_mass;

    // Thermodynamics derived factors
    // T = u * (gamma - 1) * mu * m_p / k_B
    double V_unit_m_s = config.unit_velocity_kms * 1000.0;
    config.factor_u_to_t = (V_unit_m_s * V_unit_m_s) * (config.gamma - 1.0) *
                           config.primordial_mu * constants::M_P_SI /
                           constants::K_B_SI;
    config.factor_t_to_u = 1.0 / config.factor_u_to_t;

    // Density conversion: Code -> CGS (g/cm^3)
    double mass_cgs = config.unit_mass_msun * constants::MSUN_CGS;
    double length_cgs = config.unit_length_mpc * constants::MPC_CGS;
    double vol_cgs = length_cgs * length_cgs * length_cgs;
    config.unit_density_cgs = mass_cgs / vol_cgs;

    // Cooling rate conversion: (erg/g/s) -> Code (u_code / t_code)
    double u_unit_cgs = config.unit_velocity_cgs * config.unit_velocity_cgs;
    double t_unit_cgs = config.unit_time_gyr * constants::GYR_CGS;
    config.cooling_conversion_factor = t_unit_cgs / u_unit_cgs;
}

void Config::compute_derived_data() {
    if (hubble_h > 10.0) {
        hubble_h /= 100.0;
    }

    cell_size = domain_size / mesh_size;
    cell_volume = cell_size * cell_size * cell_size;
    omega_dm = omega_m - omega_baryon;
    cutoff_radius_cells = cutoff_radius_factor * PM_smoothing_cells;
    cutoff_radius = cutoff_radius_cells * cell_size;
    cutoff_radius_squared = cutoff_radius * cutoff_radius;

    num_dm_particles = num_particles_1d * num_particles_1d * num_particles_1d;

    G = std::pow(domain_size, 3) / (6.0 * M_PI * total_mass);

    const double baryon_fraction = omega_baryon / omega_m;
    const double dm_total_mass = total_mass * (1.0 - baryon_fraction);
    dm_particle_mass = dm_total_mass / num_dm_particles;
    gas_total_mass = total_mass * baryon_fraction;

    double mean_interparticle_spacing =
        domain_size / std::cbrt(static_cast<double>(num_dm_particles));

    double actual_softening =
        mean_interparticle_spacing * comoving_softening_factor;

    // Automatically detect an oversampled gas grid
    if (use_hydro && mesh_size > num_particles_1d) {
        // Protect gas cells from DM "bullets"
        double min_softening_for_gas = 1.5 * cell_size;
        if (min_softening_for_gas > actual_softening) {
            actual_softening = min_softening_for_gas;
        }

        // Automatically blur the sparse DM mass on the fine grid
        double oversample_ratio =
            static_cast<double>(mesh_size) / num_particles_1d;

        // Dynamically scale it up
        if (PM_smoothing_cells < 1.25 * oversample_ratio) {
            PM_smoothing_cells = 1.25 * oversample_ratio;
        }
    }

    base_comoving_softening = actual_softening;
    softening_squared = actual_softening * actual_softening;

    const double dynamical_time = 1.0 / std::sqrt(G);
    fixed_dt = max_dt_dynamical_factor * dynamical_time;

    init_derived_units(*this);
}

Config::Config() { compute_derived_data(); }

void Config::load(const std::string& filename) {
    ini::Ini config_file;
    config_file.load(filename);

    mesh_size = config_file.get_int("domain", "mesh_size_1d", mesh_size);
    num_particles_1d =
        config_file.get_int("domain", "num_particles_1d", num_particles_1d);
    box_size_mpc =
        config_file.get_double("domain", "box_size_mpc", box_size_mpc);

    omega_baryon =
        config_file.get_double("cosmology", "omega_baryon", omega_baryon);
    omega_m = config_file.get_double("cosmology", "omega_M", omega_m);
    omega_lambda =
        config_file.get_double("cosmology", "omega_lambda", omega_lambda);
    hubble_h = config_file.get_double("cosmology", "Hubble_h", hubble_h);
    expanding_universe = config_file.get_bool("cosmology", "expanding_universe",
                                              expanding_universe);
    sigma_8 = config_file.get_double("cosmology", "sigma_8", sigma_8);
    spectral_index =
        config_file.get_double("cosmology", "spectral_index", spectral_index);

    comoving_softening_factor = config_file.get_double(
        "gravity", "comoving_softening_factor", comoving_softening_factor);
    physical_softening_cap_a = config_file.get_double(
        "gravity", "softening_cap_scale_factor", physical_softening_cap_a);

    standing_particles = config_file.get_bool(
        "initial_conditions", "standing_particles", standing_particles);
    initial_gas_temperature_k = config_file.get_double(
        "initial_conditions", "initial_gas_temp_k", initial_gas_temperature_k);
    seed = config_file.get_int("initial_conditions", "seed", seed);

    use_hydro = config_file.get_bool("hydro", "use_hydro", use_hydro);
    gamma = config_file.get_double("hydro", "gamma", gamma);
    enable_cooling =
        config_file.get_bool("hydro", "enable_cooling", enable_cooling);
    primordial_mu =
        config_file.get_double("hydro", "primordial_mu", primordial_mu);
    temp_floor_k =
        config_file.get_double("hydro", "temp_floor_k", temp_floor_k);
    cooling_cutoff_k =
        config_file.get_double("hydro", "cooling_cutoff_k", cooling_cutoff_k);

    enable_subgrid_gas_gravity = config_file.get_bool(
        "subgrid", "enable_subgrid_gravity", enable_subgrid_gas_gravity);
    enable_subgrid_clumping = config_file.get_bool(
        "subgrid", "enable_subgrid_clumping", enable_subgrid_clumping);
    subgrid_clumping_amplitude = config_file.get_double(
        "subgrid", "subgrid_clumping_amplitude", subgrid_clumping_amplitude);

    use_PM = config_file.get_bool("p3m", "use_pm", use_PM);
    use_PP = config_file.get_bool("p3m", "use_pp", use_PP);
    cutoff_radius_factor = config_file.get_double("p3m", "cutoff_radius_factor",
                                                  cutoff_radius_factor);
    PM_smoothing_cells =
        config_file.get_double("p3m", "pm_smoothing_cells", PM_smoothing_cells);

    max_dt_dynamical_factor = config_file.get_double(
        "time", "max_dt_dynamical_factor", max_dt_dynamical_factor);
    a_start = config_file.get_double("time", "a_start", a_start);
    a_end = config_file.get_double("time", "a_end", a_end);
    hydro_courant_factor = config_file.get_double(
        "time", "hydro_courant_factor", hydro_courant_factor);
    gravity_accuracy_eta = config_file.get_double(
        "time", "gravity_accuracy_eta", gravity_accuracy_eta);
    use_adaptive_dt =
        config_file.get_bool("time", "use_adaptive_dt", use_adaptive_dt);
    max_cycles = config_file.get_int("time", "max_cycles", max_cycles);

    save_HDF5_every_delta_a = config_file.get_double(
        "output", "save_hdf5_every_delta_a", save_HDF5_every_delta_a);
    debug_info_every_seconds = config_file.get_int(
        "output", "debug_info_every_seconds", debug_info_every_seconds);
    enable_energy_diagnostics = config_file.get_bool(
        "output", "enable_energy_diagnostics", enable_energy_diagnostics);

    num_threads = config_file.get_int("HPC", "num_threads", num_threads);
    enable_GPU = config_file.get_bool("HPC", "use_gpu", enable_GPU);

    compute_derived_data();

    if (omega_baryon > omega_m) {
        throw std::runtime_error(
            "Error: omega_baryon cannot be larger than total omega_m.");
    }
}
