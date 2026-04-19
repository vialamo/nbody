#include <iostream>
#include <csignal>
#include <atomic>
#include "config.h"
#include "utils.h"
#include "engine.h"

std::atomic<bool> keep_running{ true };

void signal_handler( int signal ) {
    if( signal == SIGINT ) {
        std::cout << "\n[Ctrl+C Detected] Finishing the current cycle and shutting down safely..." << std::endl;
        keep_running = false;
    }
}

int main( int argc, char* argv[] ) {
    std::signal( SIGINT, signal_handler );

    std::string config_filename = "simulation.ini";
    if( argc >= 2 ) {
        config_filename = argv[1];
    }

    Config config;
    try {
        config.load( config_filename );
    }
    catch( const std::exception& e ) {
        std::cerr << "Error loading " << config_filename << ": " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Successfully loaded " << config_filename << std::endl;

    HDF5Writer h5_writer( config );
    Logger logger( h5_writer.get_base_name() + "_diagnostics.csv" );

    SimulationEngine engine( config, logger, h5_writer );

    try {
        while( keep_running && engine.cycle_count < config.MAX_CYCLES && 
               engine.state.scale_factor < config.MAX_SCALE_FACTOR ) {
            engine.step();
        }
    }
    catch( const std::exception& e ) {
        std::cerr << "\nSimulation crashed: " << e.what() << std::endl;
    }

    if( !keep_running ) {
        std::cout << "Simulation aborted by user." << std::endl;
    }
    else if (engine.state.scale_factor >= config.MAX_SCALE_FACTOR) {
        std::cout << "Simulation successfully reached a = " 
                  << config.MAX_SCALE_FACTOR << "." << std::endl;
    }
    else {
        std::cout << "Simulation reached MAX_CYCLES." << std::endl;
    }

    return 0;
}