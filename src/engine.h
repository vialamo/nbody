#pragma once
#include "config.h"
#include "state.h"
#include "logger.h"
#include "diagnostics.h"
#include "hdf5_writter.h"

enum class ExitStatus { ReachedMaxScaleFactor, ReachedMaxCycles, UserAborted };

class SimulationEngine {
   private:
    std::atomic<bool> stop_requested{false};
    bool needs_more_cycles = true;
    int cycle_count = 0;
    int snapshot_count = 0;
    double current_dt = 0.0;
    double next_output_a = 0.0;

    Config& config;
    Logger& logger;
    Diagnostics& diagnostics;
    HDF5Writer& h5_writer;

    SimState state;

    void step();

   public:

    SimulationEngine(Config& conf, Logger& log, HDF5Writer& h5, Diagnostics& diag);

    ExitStatus run();
    void request_stop() { stop_requested = true; }
};
