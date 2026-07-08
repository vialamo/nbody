#include "engine.h"

#include <algorithm>

#include "gas.h"
#include "ics.h"
#include "integrator.h"
#include "particles.h"

SimulationEngine::SimulationEngine(Config& conf, Logger& log, HDF5Writer& h5,
                                   Diagnostics& diag)
    : config(conf),
      logger(log),
      h5_writer(h5),
      diagnostics(diag),
      state(initialize_state(conf)) {
    // Determine initial timestep
    current_ts = get_timestep();

    // Log initial state (Cycle 0)
    diagnostics.update_physics(state, current_ts, config);
    last_debug_time = std::chrono::high_resolution_clock::now();

    logger.log(diag);

    if (config.save_HDF5_every_delta_a > 0.0) {
        h5_writer.save_snapshot(snapshot_count, cycle_count, state, config);
        snapshot_count++;
        next_output_a = config.a_start + config.save_HDF5_every_delta_a;
    }
}

TimestepInfo SimulationEngine::get_timestep() const {
    TimestepInfo ts = {0};
    ts.subcycle_hydro = false;
    ts.subcycle_grav = false;

    if (!config.use_adaptive_dt) {
        ts.dt_macro = config.fixed_dt;
        return ts;
    }

    // Get raw physics constraints
    ts.dt_hydro = state.gas.get_cfl_timestep();
    ts.dt_grav = state.dm.get_gravity_timestep(config);

    // We calculate cooling just for logging/diagnostics.
    // It does NOT influence the macro-step
    ts.dt_cool = state.gas.get_cooling_timestep(state.scale_factor);

    // The Macro Step is bounded by the SLOWER of the two primary physics
    double base_macro = std::max(ts.dt_hydro, ts.dt_grav);

    // Cosmological Expansion Limiter
    double dt_expansion = std::numeric_limits<double>::infinity();
    if (config.expanding_universe && state.hubble_param > 0.0) {
        dt_expansion = 0.01 / state.hubble_param;
    }

    // Apply safety caps
    ts.dt_macro = std::min({base_macro, dt_expansion, config.fixed_dt});

    // Force the simulation to land on the next output target
    if (config.save_HDF5_every_delta_a > 0.0 && config.expanding_universe) {
        double t_start = get_time_from_scale_factor(config.a_start, config);
        double current_t = t_start + state.total_time;
        double target_t = get_time_from_scale_factor(next_output_a, config);

        double dt_snapshot = target_t - current_t;

        const double MIN_VALID_DT = 1e-10;

        if (dt_snapshot > MIN_VALID_DT && dt_snapshot < ts.dt_macro) {
            ts.dt_macro = dt_snapshot;
        }
    }

    // Determine who subcycles
    if (ts.dt_hydro < ts.dt_grav) {
        ts.subcycle_hydro = true;
        ts.subcycle_grav = false;
    } else if (ts.dt_hydro > ts.dt_grav) {
        ts.subcycle_hydro = false;
        ts.subcycle_grav = true;
    }

    return ts;
}

void SimulationEngine::step() {
    {
        ScopedTimer step_timer(diagnostics, TimerRegion::Step);
        KDK_step(state, current_ts, config, diagnostics);
    }

    diagnostics.increment_cycle();

    cycle_count++;

    // Update Timestep for next cycle
    current_ts = get_timestep();

    needs_more_cycles =
        state.scale_factor < config.a_end && cycle_count < config.max_cycles;

    // I/O and Logging
    const double TOLERANCE = 1e-7;
    bool must_save_snapshot = state.scale_factor >= (next_output_a - TOLERANCE);
    if (config.save_HDF5_every_delta_a > 0.0 &&
        (!needs_more_cycles || must_save_snapshot)) {
        ScopedTimer io_timer(diagnostics, TimerRegion::IO);
        h5_writer.save_snapshot(snapshot_count, cycle_count, state, config);

        snapshot_count++;
        while (next_output_a <= state.scale_factor) {
            next_output_a += config.save_HDF5_every_delta_a;
        }
    }

    // Replace the old DEBUG_INFO_EVERY_CYCLES check with this:
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = now - last_debug_time;

    if (config.debug_info_every_seconds > 0.0 &&
        (!needs_more_cycles ||
         elapsed.count() >= config.debug_info_every_seconds)) {
        diagnostics.update_physics(state, current_ts, config);
        logger.log(diagnostics);
        diagnostics.reset_accumulators();

        // Reset the timer for the next interval
        last_debug_time = std::chrono::high_resolution_clock::now();
    }
}

ExitStatus SimulationEngine::run() {
    while (!stop_requested && needs_more_cycles) {
        step();
    }

    if (stop_requested) {
        return ExitStatus::UserAborted;
    }
    if (cycle_count >= config.max_cycles) {
        return ExitStatus::ReachedMaxCycles;
    }

    return ExitStatus::ReachedMaxScaleFactor;
}