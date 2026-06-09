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
    current_dt = get_timestep();

    // Log initial state (Cycle 0)
    diagnostics.update_physics(state, current_dt, config);

    logger.log(diag);

    if (config.SAVE_HDF5_EVERY_DELTA_A > 0.0) {
        h5_writer.save_snapshot(snapshot_count, cycle_count, state, config);
        snapshot_count++;
        next_output_a = config.START_A + config.SAVE_HDF5_EVERY_DELTA_A;
    }
}

double SimulationEngine::get_timestep() const {
    double dt = config.FIXED_DT;

    if (config.USE_ADAPTIVE_DT) {
        double dt_cfl = state.gas.get_cfl_timestep();
        double dt_grav = state.dm.get_gravity_timestep(config);
        double dt_cool = state.gas.get_cooling_timestep(state.scale_factor);

        // Cosmological Expansion Limiter
        double dt_expansion = std::numeric_limits<double>::infinity();
        if (config.EXPANDING_UNIVERSE && state.hubble_param > 0.0) {
            // Restrict timestep so the universe expands by at most 1% per step
            dt_expansion = 0.01 / state.hubble_param;
        }

        dt =
            std::min({dt_cfl, dt_grav, dt_cool, dt_expansion, config.FIXED_DT});
    }

    // Force the simulation to land exactly on the next output target
    if (config.SAVE_HDF5_EVERY_DELTA_A > 0.0 && config.EXPANDING_UNIVERSE) {
        double t_start = get_time_from_scale_factor(config.START_A, config);
        double current_t = t_start + state.total_time;

        // Target absolute code time for the snapshot
        double target_t = get_time_from_scale_factor(next_output_a, config);

        double dt_snapshot = target_t - current_t;

        // If the snapshot is the most imminent event, shrink dt to match it
        if (dt_snapshot > 0.0 && dt_snapshot < dt) {
            dt = dt_snapshot;
        }
    }

    return dt;
}

void SimulationEngine::step() {
    {
        ScopedTimer step_timer(diagnostics, TimerRegion::Step);
        KDK_step(state, current_dt, config, diagnostics);
    }

    diagnostics.increment_cycle();

    cycle_count++;

    // Update Timestep for next cycle
    current_dt = get_timestep();

    needs_more_cycles = state.scale_factor < config.MAX_SCALE_FACTOR &&
                        cycle_count < config.MAX_CYCLES;

    // I/O and Logging
    bool must_save_snapshot = state.scale_factor >= next_output_a;
    if (config.SAVE_HDF5_EVERY_DELTA_A > 0.0 &&
        (!needs_more_cycles || must_save_snapshot)) {
        ScopedTimer io_timer(diagnostics, TimerRegion::IO);
        h5_writer.save_snapshot(snapshot_count, cycle_count, state, config);

        snapshot_count++;
        while (next_output_a <= state.scale_factor) {
            next_output_a += config.SAVE_HDF5_EVERY_DELTA_A;
        }
    }

    if (config.DEBUG_INFO_EVERY_CYCLES > 0 &&
        (!needs_more_cycles ||
         cycle_count % config.DEBUG_INFO_EVERY_CYCLES == 0)) {
        diagnostics.update_physics(state, current_dt, config);
        logger.log(diagnostics);
        diagnostics.reset_accumulators();
    }
}

ExitStatus SimulationEngine::run() {
    while (!stop_requested && needs_more_cycles) {
        step();
    }

    if (stop_requested) {
        return ExitStatus::UserAborted;
    }
    if (cycle_count >= config.MAX_CYCLES) {
        return ExitStatus::ReachedMaxCycles;
    }

    return ExitStatus::ReachedMaxScaleFactor;
}