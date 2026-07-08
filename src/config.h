#pragma once
#include <string>

struct Config {
    // Domain
    int mesh_size = 32;
    int num_particles_1d = 32;
    double box_size_mpc = 16.0;

    // Cosmology
    double omega_baryon = 0.045;
    double omega_m = 0.3;
    double omega_lambda = 0.7;
    double hubble_h = 0.7;
    double sigma_8 = 0.81;
    double spectral_index = 0.96;
    bool expanding_universe = true;

    // Gravity
    double comoving_softening_factor = 0.0334;
    double physical_softening_cap_a = 0.3;

    // Initial conditions
    bool standing_particles = false;
    double initial_gas_temperature_k = 50.0;
    int seed = 42;

    // Hydro
    bool use_hydro = true;
    double gamma = 5.0 / 3.0;
    bool enable_cooling = true;
    double primordial_mu = 1.22;  // Mean molecular weight for neutral
                                  // primordial gas (76% H, 24% He)
    double temp_floor_k = 10.0;
    double cooling_cutoff_k = 10000.0;

    // Subgrid
    bool enable_subgrid_gas_gravity = false;
    bool enable_subgrid_clumping = true;
    double subgrid_clumping_amplitude = 10.0;

    // P3M
    bool use_PM = true;
    bool use_PP = true;
    double cutoff_radius_factor = 4.5;
    double PM_smoothing_cells = 1.25;

    // Time
    bool use_adaptive_dt = true;
    double max_dt_dynamical_factor = 1e-3;
    double hydro_courant_factor = 0.4;
    double gravity_accuracy_eta = 0.2;
    double a_start = 0.02;
    double a_end = 1.0;
    int max_cycles = 1000000000;

    // Output
    double save_HDF5_every_delta_a = 0.005;
    int debug_info_every_seconds = 30;
    bool enable_energy_diagnostics = true;

    // HPC
    int num_threads = 0;
    bool enable_GPU = true;

    // Derived Parameters
    double G = 0.0;
    double cell_size = 0.0, cell_volume = 0.0, omega_dm = 0.0;
    double cutoff_radius_cells = 0;
    double cutoff_radius = 0.0, cutoff_radius_squared = 0.0;
    int num_dm_particles = 0;
    double dm_particle_mass = 0.0, gas_total_mass = 0.0;
    double softening_squared = 0.0;
    double base_comoving_softening = 0.0;

    double fixed_dt = 0.0;

    // Derived unit conversions (calculated once)
    double unit_length_mpc = 0.0;    // 1 code length = X Mpc
    double unit_time_gyr = 0.0;      // 1 code time = X Gyr
    double unit_velocity_kms = 0.0;  // 1 code velocity = X km/s
    double unit_velocity_cgs = 0.0;  // 1 code velocity = X cm/s
    double unit_mass_msun = 0.0;     // 1 code mass = X Solar Masses

    // Thermodynamics
    double factor_u_to_t = 0.0;  // 1 code specific internal energy = X Kelvin
    double factor_t_to_u = 0.0;  // 1 Kelvin = X code specific internal energy
    double unit_density_cgs = 0.0;  // 1 code density = X g/cm^3
    double cooling_conversion_factor =
        0.0;  // Converts physical cooling rate (erg/g/s) to code units (du/dt)

    // Constant Parameters
    static constexpr double total_mass = 1;
    static constexpr double domain_size = 1.0;

    Config();
    void load(const std::string& filename);
    void compute_derived_data();
};
