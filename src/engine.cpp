#include "engine.h"

#include "cosmology.h"
#include "diagnostics.h"
#include "gravity.h"
#include "hdf5_writer.h"
#include "hydro.h"
#include "ics.h"
#include "integrator.h"
#include "logger.h"

namespace {
uint64_t get_next_sync_ti(const std::vector<uint64_t>& ti_end,
                          size_t num_particles, uint64_t ti_global_time) {
    uint64_t next_ti = std::numeric_limits<uint64_t>::max();
    for (size_t i = 0; i < num_particles; ++i) {
        if (ti_end[i] > ti_global_time) {
            next_ti = std::min(next_ti, ti_end[i]);
        }
    }
    return next_ti - ti_global_time;
}
}  // namespace

SimulationEngine::SimulationEngine(Config& conf, Logger& log, HDF5Writer& h5,
                                   Diagnostics& diag)
    : config(conf), logger(log), h5_writer(h5), diagnostics(diag), state(conf) {
    initialize_state(state, conf);
    integrator.initialize(state, config);

    // Determine initial timestep
    current_ts = get_timestep();

    // Log initial state (Cycle 0)
    diagnostics.update_physics(state, current_ts, config);
    last_debug_time = std::chrono::high_resolution_clock::now();

    logger.log(diag, conf);

    if (config.save_HDF5_every_delta_a > 0.0) {
        h5_writer.save_snapshot(snapshot_count, cycle_count, state, config,
                                current_ts);
        snapshot_count++;
        if (config.expanding_universe) {
            next_output_a = config.a_start + config.save_HDF5_every_delta_a;
        } else {
            next_output_a = config.save_HDF5_every_delta_a;
        }
    }
}

TimestepInfo SimulationEngine::get_timestep() const {
    TimestepInfo ts = {0};

    if (config.individual_particle_timesteps) {
        uint64_t ti_gas =
            get_next_sync_ti(state.gas.ti_end, state.gas.num_particles,
                             state.gas.ti_global_time);
        uint64_t ti_dm = get_next_sync_ti(
            state.dm.ti_end, state.dm.num_particles, state.dm.ti_global_time);

        ts.ti_macro = std::min(ti_gas, ti_dm);
        ts.dt_macro = ts.ti_macro * config.timebase_to_float;
        return ts;
    }

    // Get raw physics constraints
    ts.dt_hydro = HydroSolver::get_timestep(state.gas, config);

    if (config.use_PM || config.use_PP) {
        ts.dt_grav = state.gravity.get_timestep(state.dm, state.gas, config);
    } else {
        ts.dt_grav = std::numeric_limits<double>::infinity();
    }

    // Cosmological Expansion Limiter
    double dt_expansion = std::numeric_limits<double>::infinity();
    if (config.expanding_universe && state.hubble_param > 0.0) {
        dt_expansion = 0.01 / state.hubble_param;
    }

    ts.dt_macro =
        std::min({ts.dt_hydro, ts.dt_grav, dt_expansion, config.fixed_dt});

    // Force the simulation to land on the next output target
    if (config.save_HDF5_every_delta_a > 0.0 && config.expanding_universe) {
        double t_start =
            Cosmology::get_time_from_scale_factor(config.a_start, config);
        double current_t = t_start + state.total_time;
        double target_t =
            Cosmology::get_time_from_scale_factor(next_output_a, config);

        double dt_snapshot = target_t - current_t;
        if (dt_snapshot >= config.timebase_to_float &&
            dt_snapshot < ts.dt_macro) {
            ts.dt_macro = dt_snapshot;
        }
    }

    // Quantize the continuous physics limit to the nearest integer tick BELOW
    // the limit. We truncate (static_cast) rather than round
    ts.ti_macro = static_cast<uint64_t>(ts.dt_macro / config.timebase_to_float);
    if (ts.ti_macro == 0) ts.ti_macro = 1;
    ts.dt_macro = ts.ti_macro * config.timebase_to_float;

    return ts;
}

void SimulationEngine::step() {
    {
        ScopedTimer step_timer(diagnostics, TimerRegion::Step);
        integrator.KDK_step(state, current_ts.ti_macro, config, diagnostics);
    }

    diagnostics.increment_cycle();

    cycle_count++;

    // Update Timestep for next cycle
    current_ts = get_timestep();

    // I/O and Logging
    bool is_synchronized = true;
    if (config.individual_particle_timesteps) {
        // Synchronized when the integer total time is a multiple of the
        // TIMEBASE
        is_synchronized = (state.ti_total_time % TimeConstants::TIMEBASE == 0);
    }

    needs_more_cycles =
        !is_synchronized ||
        (cycle_count < config.max_cycles &&
         ((config.expanding_universe && state.scale_factor < config.a_end) ||
          (!config.expanding_universe && state.total_time < config.a_end)));

    double next_out = next_output_a - (config.timebase_to_float / 2.0);
    bool must_save_snapshot =
        is_synchronized &&
        ((config.expanding_universe && state.scale_factor >= next_out) ||
         (!config.expanding_universe && state.total_time >= next_out));

    if (config.save_HDF5_every_delta_a > 0.0 &&
        (!needs_more_cycles || must_save_snapshot)) {
        ScopedTimer io_timer(diagnostics, TimerRegion::IO);
        h5_writer.save_snapshot(snapshot_count, cycle_count, state, config,
                                current_ts);

        snapshot_count++;

        if (config.expanding_universe) {
            while (next_output_a <= state.scale_factor) {
                next_output_a += config.save_HDF5_every_delta_a;
            }
        } else {
            while (next_output_a <= state.total_time) {
                next_output_a += config.save_HDF5_every_delta_a;
            }
        }
    }

    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = now - last_debug_time;

    if (config.debug_info_every_seconds > 0.0 &&
        (!needs_more_cycles ||
         elapsed.count() >= config.debug_info_every_seconds)) {
        // Update dt_cool for diagnostics
        current_ts.dt_cool = std::numeric_limits<double>::infinity();

        diagnostics.update_physics(state, current_ts, config);
        logger.log(diagnostics, config);
        diagnostics.reset_accumulators();
        state.gas.reset_diagnostic_counters();

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
