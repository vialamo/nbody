#include "engine.h"

#include "cosmology.h"
#include "diagnostics.h"
#include "gas.h"
#include "gravity.h"
#include "hdf5_writer.h"
#include "ics.h"
#include "integrator.h"
#include "logger.h"

static double get_next_sync_dt(const std::vector<double>& t_end,
                               size_t num_particles, double global_time) {
    double next_time = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < num_particles; ++i) {
        if (t_end[i] > global_time + 1e-10) {
            next_time = std::min(next_time, t_end[i]);
        }
    }
    return (next_time - global_time > 1e-10)
               ? (next_time - global_time)
               : std::numeric_limits<double>::infinity();
}

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
        double dt_gas = std::numeric_limits<double>::infinity();
        if (config.hydro_method == HydroMethod::MFM) {
            dt_gas = get_next_sync_dt(state.mfm_gas->t_end,
                                      state.mfm_gas->num_particles,
                                      state.total_time);
        }
        double dt_dm = get_next_sync_dt(state.dm.t_end, state.dm.num_particles,
                                        state.total_time);
        ts.dt_macro = std::min(dt_gas, dt_dm);
        return ts;
    }

    // Get raw physics constraints
    if (config.hydro_method == HydroMethod::Eulerian) {
        ts.dt_hydro = state.gas->get_cfl_timestep();
    } else if (config.hydro_method == HydroMethod::MFM) {
        ts.dt_hydro =
            state.mfm_gas->get_cfl_timestep(state.scale_factor, config);
    } else {
        ts.dt_hydro = std::numeric_limits<double>::infinity();
    }

    if (config.use_PM || config.use_PP) {
        double dt_grav_dm = GravitySolver::get_timestep(state.dm, config);
        double dt_grav_mfm = std::numeric_limits<double>::infinity();
        if (config.hydro_method == HydroMethod::MFM) {
            dt_grav_mfm = GravitySolver::get_timestep(*state.mfm_gas, config);
        }
        ts.dt_grav = std::min(dt_grav_dm, dt_grav_mfm);
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
        const double MIN_VALID_DT = 1e-10;
        if (dt_snapshot > MIN_VALID_DT && dt_snapshot < ts.dt_macro) {
            ts.dt_macro = dt_snapshot;
        }
    }

    return ts;
}

void SimulationEngine::step() {
    {
        ScopedTimer step_timer(diagnostics, TimerRegion::Step);
        Integrator::KDK_step(state, current_ts, config, diagnostics);
    }

    diagnostics.increment_cycle();

    cycle_count++;

    // Update Timestep for next cycle
    current_ts = get_timestep();

    needs_more_cycles =
        ((config.expanding_universe && state.scale_factor < config.a_end) ||
         (!config.expanding_universe && state.total_time < config.a_end)) &&
        cycle_count < config.max_cycles;

    // I/O and Logging
    bool is_synchronized = true;
    if (config.individual_particle_timesteps) {
        // At a true macrostep boundary, all particles have reached their t_end
        // and are flagged as active for the upcoming cycle
        double ratio = state.total_time / config.fixed_dt;
        is_synchronized = (std::abs(ratio - std::round(ratio)) < 1e-7);
    }

    constexpr double TOLERANCE = 1e-7;
    double next_out = next_output_a - TOLERANCE;
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
        if (config.hydro_method == HydroMethod::Eulerian) {
            current_ts.dt_cool = state.gas->get_cooling_timestep(
                state.scale_factor, state.cooling);
        } else {
            current_ts.dt_cool = std::numeric_limits<double>::infinity();
        }
        diagnostics.update_physics(state, current_ts, config);
        logger.log(diagnostics, config);
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
