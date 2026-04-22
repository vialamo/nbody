#pragma once
#include "state.h"
#include "config.h"
#include "utils.h"

class SimulationEngine {
public:
    Config& config;
    Logger& logger;
    HDF5Writer& h5_writer;

    SimState state;

    int cycle_count = 0;
    int snapshot_count = 0;
    double current_dt = 0.0;
    double dt_cfl = 0.0;
    double dt_grav = 0.0;

    SimulationEngine( Config& conf, Logger& log, HDF5Writer& h5 );
    void step();
};
