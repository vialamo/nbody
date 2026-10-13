#pragma once
#include <algorithm>
#include <cctype>
#include <string>

#include "types.h"

// Simple struct to hold the 3D displacement vectors

namespace EnumUtils {
template <typename E>
struct Traits;

template <typename E>
std::string to_string(E value) {
    size_t idx = static_cast<size_t>(value);
    if (idx < static_cast<size_t>(E::Count)) {
        return Traits<E>::names[idx];
    }
    return Traits<E>::names[static_cast<size_t>(Traits<E>::default_value)];
}

template <typename E>
E from_string(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    for (size_t i = 0; i < static_cast<size_t>(E::Count); ++i) {
        if (str == Traits<E>::names[i]) {
            return static_cast<E>(i);
        }
    }
    return Traits<E>::default_value;
}
}  // namespace EnumUtils

enum class HydroMethod { None, Eulerian, MFM, Count };

enum class InitialSetup {
    Cosmological,
    SodShockTube,
    AdiabaticExpansion,
    SedovBlastwave,
    SoundWave,
    Glass,
    Count
};

template <>
struct EnumUtils::Traits<HydroMethod> {
    static constexpr const char* names[] = {"none", "eulerian", "mfm"};
    static constexpr HydroMethod default_value = HydroMethod::None;
};

template <>
struct EnumUtils::Traits<InitialSetup> {
    static constexpr const char* names[] = {
        "cosmological",    "sod_shock_tube", "adiabatic_expansion",
        "sedov_blastwave", "sound_wave",     "glass"};
    static constexpr InitialSetup default_value = InitialSetup::Cosmological;
};

namespace HydroConfig {
inline std::string to_string(HydroMethod method) {
    return EnumUtils::to_string(method);
}

// Wrapper allows calling without <HydroMethod> template brackets
inline HydroMethod from_string(std::string str) {
    return EnumUtils::from_string<HydroMethod>(std::move(str));
}
}  // namespace HydroConfig

namespace InitialConfig {
inline std::string to_string(InitialSetup setup) {
    return EnumUtils::to_string(setup);
}

inline InitialSetup from_string(std::string str) {
    return EnumUtils::from_string<InitialSetup>(std::move(str));
}
}  // namespace InitialConfig

struct Config {
    // Domain
    double mesh_size_1d = 64;
    int num_particles_1d = 64;
    int num_gas_particles_1d = 64;
    double box_size_mpc_1d = 96.0;

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
    bool use_PM = true;
    bool use_PP = true;
    double cutoff_radius_factor = 4.5;
    double PM_smoothing_cells = 1.25;

    // Initial conditions
    InitialSetup initial_setup = InitialSetup::Cosmological;
    bool use_glass = true;
    bool fixed_ics = false;
    bool invert_phases = false;
    double initial_gas_temperature_k = 50.0;
    double seed_metallicity_solar = 0.0;
    int seed = 42;

    // Hydro
    HydroMethod hydro_method = HydroMethod::Eulerian;
    double gamma = 5.0 / 3.0;
    bool enable_cooling = true;
    double primordial_mu = 1.22;  // Mean molecular weight for neutral
                                  // primordial gas (76% H, 24% He)
    double temp_floor_k = 10.0;
    double cooling_cutoff_k = 10000.0;
    std::string cooling_table_path = "";

    // MFM
    double mfm_target_neighbors = 32.0;
    double mfm_neighbor_tolerance = 0.01;
    int mfm_max_iterations = 50;
    double mfm_min_hsml_fraction = 0.25;
    double dual_energy_threshold = 0.001;
    bool use_exact_riemann = false;

    // Subgrid
    bool enable_subgrid_clumping = true;
    double subgrid_clumping_amplitude = 10.0;

    // Time
    double max_dt_dynamical_factor = 1e-3;
    double hydro_courant_factor = 0.4;
    double gravity_accuracy_eta = 0.2;
    double a_start = 0.02;
    double a_end = 1.0;
    int max_cycles = 1000000000;
    bool individual_particle_timesteps = true;

    // Output
    double save_HDF5_every_delta_a = 0.005;
    double debug_info_every_seconds = 30;

    // HPC
    int num_threads = 0;

    // Derived Parameters
    double G = 0.0;
    double cell_size = 0.0;
    double cell_volume = 0.0, omega_dm = 0.0;
    double cutoff_radius_cells = 0;
    double cutoff_radius = 0.0, cutoff_radius_squared = 0.0;
    int num_dm_particles = 0;
    int num_gas_particles = 0;
    double dm_particle_mass = 0.0, gas_total_mass = 0.0;
    double softening = 0.0;
    double base_comoving_softening = 0.0;

    Vec3 mesh_size;
    Vec3 box_size_mpc;

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
    double total_mass = 1;
    static constexpr double domain_size_1d = 1.0;
    Vec3 domain_size;
#ifdef USE_GPU
    static constexpr bool using_GPU = true;
#else
    static constexpr bool using_GPU = false;
#endif

    Config();
    void load(const std::string& filename);
    void compute_derived_data();

   private:
    void init_derived_units();
};
