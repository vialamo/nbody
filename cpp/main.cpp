#include <iostream>
#include "config.h"
#include "utils.h"
#include "engine.h"
#include "visualization.h"

int main( int argc, char* argv[] ) {
    std::string config_filename = "simulation.ini";
    if( argc >= 2 ) {
        config_filename = argv[1];
    }

    // Load Configuration
    Config config;
    try {
        config.load( config_filename );
    }
    catch( const std::exception& e ) {
        std::cerr << "Error loading " << config_filename << ": " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Successfully loaded " << config_filename << std::endl;

    // Initialize Output Files
    HDF5Writer h5_writer( config );
    Logger logger( h5_writer.get_base_name() + "_diagnostics.csv");

    // Create the Physics Engine
    SimulationEngine engine( config, logger, h5_writer );

    // Run the Simulation
    try {
        if( config.ENABLE_VISUALIZATION ) {
            SimulationApp sim_app( engine );
            sim_app.run(); // Blocks until window is closed
        }
        else {
            std::cout << "Starting headless simulation loop..." << std::endl;
            while( engine.cycle_count < config.MAX_CYCLES ) {
                engine.step();
            }
        }
    }
    catch( const std::exception& e ) {
        std::cerr << "\nSimulation interrupted: " << e.what() << std::endl;
    }

    // Cleanup
    std::cout << "Simulation finished" << std::endl;

    return 0;
}