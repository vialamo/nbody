#pragma once
#include <Eigen/Dense>
#include <array>
#include <chrono>

struct Config;
struct SimState;
struct TimestepInfo;
struct Grid3D;
class GasGrid;
class GasParticleSystem;
class ParticleSystem;

enum class TimerRegion {
    Step,
    PM,
    PP,
    Hydro,
    Cool,
    IO,
    Tree,
    Density,
    PP_GPU,
    NUM_REGIONS
};
enum class ProfRegion { Transf, Ret, Compute, Tree, NUM_PROF_REGIONS };

class Diagnostics {
   private:
    // Physics State
    int cycle = 0;
    double sim_time = 0.0;
    double scale_factor = 1.0;
    double total_mass_gas = 0.0;
    double total_mass_dm = 0.0;
    double initial_mass = 0.0;
    double mass_err = 0.0;
    Eigen::Vector3d total_momentum_gas = {0.0, 0.0, 0.0};
    Eigen::Vector3d total_momentum_dm = {0.0, 0.0, 0.0};
    double total_momentum = 0.0;
    double ke_gas = 0.0, ke_dm = 0.0, pe_total = 0.0, ie_gas = 0.0;
    double dt_cfl = 0.0, dt_gravity = 0.0, dt_cool = 0.0, dt_final = 0.0;
    double max_gas_density = 0.0, max_gas_pressure = 0.0,
           max_gas_velocity = 0.0;
    double total_radiated_energy = 0.0;
    double total_heated_energy = 0.0;
    double initial_energy = 0.0;
    size_t non_converged_cooling_cells = 0;
    double energy_err = 0.0;
    double initial_gas_energy = 0.0;
    double dm_energy_err = 0.0;
    double initial_dm_energy = 0.0;
    bool energy_initialized = false;
    bool dm_energy_initialized = false;
    size_t MFM_h_non_converged_cases = 0;
    size_t MFM_not_enough_nbg_cases = 0;
    size_t MFM_sph_fallback_cases = 0;
    size_t MFM_total_particles = 0;
    size_t MFM_riemann_failed_cases = 0;
    double MFM_ill_percent = 0.0;
    double MFM_clamp_percent = 0.0;
    double MFM_inc_percent = 0.0;
    double MFM_flux_energy_error = 0.0;
    double MFM_flux_momentum_error = 0.0;
    double percent_particles_updated = 0;
    double volume_error = 0.0;

    // Performance State
    int accumulated_cycles = 0;
    std::array<double, static_cast<size_t>(TimerRegion::NUM_REGIONS)>
        accumulated_times{};
    std::array<double, static_cast<size_t>(ProfRegion::NUM_PROF_REGIONS)>
        accumulated_prof{};

    void compute_dm_metrics(const ParticleSystem& dm);
    double calculate_hydro_energy_error(double current_energy, double w_grav,
                                        double w_exp, double e_rad,
                                        double e_heat);
    void compute_eulerian_metrics(const GasGrid& gas, const Config& config);
    void compute_mfm_metrics(GasParticleSystem& gas, const Config& config);
    double compute_potential_energy(const Grid3D& total_rho, const Grid3D& phi,
                                    const Config& config);

    friend class Logger;

   public:
    Diagnostics() = default;
    ~Diagnostics() = default;

    void add_prof_time(ProfRegion region, double time_sec);
    void add_time(TimerRegion region, double time_sec);
    void increment_cycle();
    double get_average(TimerRegion region) const;
    double get_prof_average(ProfRegion region) const;
    void reset_accumulators();
    double get_average_overhead() const;
    double get_io_time() const;

    void update_physics(const SimState& state, const TimestepInfo& ts,
                        const Config& config);
};

class ScopedTimer {
   private:
    Diagnostics& diag;
    TimerRegion region;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

   public:
    ScopedTimer(Diagnostics& d, TimerRegion r);
    ~ScopedTimer();
};
