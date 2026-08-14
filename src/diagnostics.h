#pragma once
#include <Eigen/Dense>
#include <array>
#include <chrono>

#include "config.h"
#include "state.h"

enum class TimerRegion { Step, PM, PP, Hydro, Cool, IO, NUM_REGIONS };
enum class SubstepCounter { Hydro, Gravity, Cool, NUM_COUNTERS };
enum class ProfRegion { Transf, Ret, Compute, NUM_PROF_REGIONS };

class Diagnostics {
   private:
    // Physics State
    int cycle = 0;
    double sim_time = 0.0;
    double scale_factor = 1.0;
    double total_mass_gas = 0.0;
    double total_mass_dm = 0.0;
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

    // Performance State
    int accumulated_cycles = 0;
    std::array<double, static_cast<size_t>(TimerRegion::NUM_REGIONS)>
        accumulated_times{};
    std::array<int, static_cast<size_t>(SubstepCounter::NUM_COUNTERS)>
        accumulated_substeps{};
    std::array<double, static_cast<size_t>(ProfRegion::NUM_PROF_REGIONS)>
        accumulated_prof{};

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

    void add_substeps(SubstepCounter counter, int count = 1) {
        accumulated_substeps[static_cast<size_t>(counter)] += count;
    }

    double get_average_substeps(SubstepCounter counter) const {
        if (accumulated_cycles == 0) return 0.0;
        return static_cast<double>(
                   accumulated_substeps[static_cast<size_t>(counter)]) /
               accumulated_cycles;
    }
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