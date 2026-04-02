#include "engine.h"
#include "ics.h"
#include "integrator.h"
#include "gas.h"
#include "particles.h"
#include <algorithm>

SimulationEngine::SimulationEngine( Config& conf, Logger& log, HDF5Writer& h5 )
    : config( conf ), logger( log ), h5_writer( h5 ), state( initialize_state( conf ) )
{
    // Determine initial timestep
    dt_cfl = state.gas.get_cfl_timestep();
    dt_grav = state.dm.get_gravity_timestep( config );
    current_dt = config.USE_ADAPTIVE_DT ? std::min( { dt_cfl, dt_grav, config.FIXED_DT } ) : config.FIXED_DT;

    // Log initial state (Cycle 0)
    std::map<std::string, double> initial_timings;
    Diagnostics diag = calculate_diagnostics( state, initial_timings, current_dt, cycle_count, config );
    logger.log( diag );
}

void SimulationEngine::step() {
    // Advance Physics
    std::map<std::string, double> timings = KDK_step( state, current_dt, config );
    cycle_count++;

    // Update Timestep for next cycle
    if( config.USE_ADAPTIVE_DT ) {
        dt_cfl = state.gas.get_cfl_timestep();
        dt_grav = state.dm.get_gravity_timestep( config );
        current_dt = std::min( { dt_cfl, dt_grav, config.FIXED_DT } );
    }
    else {
        current_dt = config.FIXED_DT;
    }

    // I/O and Logging
    double io_time = 0.0;
    if( config.SAVE_HDF5_EVERY_CYCLES > 0 && cycle_count % config.SAVE_HDF5_EVERY_CYCLES == 0 ) {
        io_time = h5_writer.save_snapshot( snapshot_count, state, config );
        snapshot_count++;
    }
    timings["io"] = io_time;

    if( config.DEBUG_INFO_EVERY_CYCLES > 0 && cycle_count % config.DEBUG_INFO_EVERY_CYCLES == 0 ) {
        Diagnostics diag = calculate_diagnostics( state, timings, current_dt, cycle_count, config );
        logger.log( diag );
    }
}