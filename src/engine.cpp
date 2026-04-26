#include "engine.h"

#include <algorithm>

#include "gas.h"
#include "ics.h"
#include "integrator.h"
#include "particles.h"

SimulationEngine::SimulationEngine(Config& conf, Logger& log, HDF5Writer& h5)
    : config(conf), logger(log), h5_writer(h5), state(initialize_state(conf)) {
    // Determine initial timestep
    double dt_cfl = state.gas.get_cfl_timestep();
    double dt_grav = state.dm.get_gravity_timestep(config);
    current_dt = config.USE_ADAPTIVE_DT
                     ? std::min({dt_cfl, dt_grav, config.FIXED_DT})
                     : config.FIXED_DT;

    // Log initial state (Cycle 0)
    std::map<std::string, double> initial_timings;
    Diagnostics diag = calculate_diagnostics(state, initial_timings, current_dt,
                                             cycle_count, config);
    logger.log(diag);

    if (config.SAVE_HDF5_EVERY_DELTA_A > 0.0) {
        h5_writer.save_snapshot(snapshot_count, cycle_count, state, config);
        snapshot_count++;
        next_output_a = config.START_A + config.SAVE_HDF5_EVERY_DELTA_A;
    }
}

void SimulationEngine::step() {
    // Advance Physics
    std::map<std::string, double> timings = KDK_step(state, current_dt, config);
    cycle_count++;

    // Update Timestep for next cycle
    if (config.USE_ADAPTIVE_DT) {
        double dt_cfl = state.gas.get_cfl_timestep();
        double dt_grav = state.dm.get_gravity_timestep(config);
        current_dt = std::min({dt_cfl, dt_grav, config.FIXED_DT});
    }

    needs_more_cycles = state.scale_factor < config.MAX_SCALE_FACTOR &&
                        cycle_count < config.MAX_CYCLES;

    // I/O and Logging
    double io_time = 0.0;
    bool must_save_snapshot = state.scale_factor >= next_output_a;
    if (config.SAVE_HDF5_EVERY_DELTA_A > 0.0 &&
        (!needs_more_cycles || must_save_snapshot)) {
        io_time =
            h5_writer.save_snapshot(snapshot_count, cycle_count, state, config);
        snapshot_count++;
        while (next_output_a <= state.scale_factor) {
            next_output_a += config.SAVE_HDF5_EVERY_DELTA_A;
        }
    }
    timings["io"] = io_time;

    if (config.DEBUG_INFO_EVERY_CYCLES > 0 &&
        (!needs_more_cycles ||
         cycle_count % config.DEBUG_INFO_EVERY_CYCLES == 0)) {
        Diagnostics diag = calculate_diagnostics(state, timings, current_dt,
                                                 cycle_count, config);
        logger.log(diag);
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