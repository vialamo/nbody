#pragma once
#include <Eigen/Dense>
#include <array>
#include <chrono>

#include "config.h"
#include "state.h"

enum class TimerRegion { Step, PM, PP, Hydro, IO, NUM_REGIONS };

class Diagnostics {
   private:
    // Physics State
    int cycle = 0;
    double sim_time = 0.0;
    double scale_factor = 1.0;
    double total_mass_gas = 0.0;
    double total_mass_dm = 0.0;
    Eigen::Vector3d total_momentum_gas = {0.0, 0.0, 0.0};
    Eigen::Vector3d total_momentum_dm = {0.0, 0.0, 0.0};
    double ke_gas = 0.0, ke_dm = 0.0, pe_total = 0.0, ie_gas = 0.0;
    double dt_cfl = 0.0, dt_gravity = 0.0, dt_final = 0.0;
    double max_gas_density = 0.0, max_gas_pressure = 0.0,
           max_gas_velocity = 0.0;

    // Performance State
    std::array<double, static_cast<size_t>(TimerRegion::NUM_REGIONS)>
        accumulated_times{};
    int accumulated_cycles = 0;

    friend class Logger;

   public:
    Diagnostics() = default;
    ~Diagnostics() = default;

    void add_time(TimerRegion region, double time_sec);
    void increment_cycle();
    double get_average(TimerRegion region) const;
    void reset_accumulators();
    double get_average_overhead() const;
    double get_io_time() const;

    double total_mass() const;
    double total_energy() const;

    void update_physics(const SimState& state, double dt, const Config& config);
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