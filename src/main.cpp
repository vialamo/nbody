#include <atomic>
#include <csignal>
#include <filesystem>
#include <iostream>
#include <string>

#include "config.h"
#include "engine.h"
#include "utils.h"

SimulationEngine* g_engine = nullptr;

void signal_handler(int signal) {
    if (signal == SIGINT) {
        std::cout << "\n[Ctrl+C Detected] Finishing the current cycle and "
                     "shutting down safely..."
                  << std::endl;
        if (g_engine != nullptr) {
            g_engine->request_stop();
        }
    }
}

int main(int argc, char* argv[]) {
    std::signal(SIGINT, signal_handler);

    std::string config_filename = "simulation.ini";
    if (argc >= 2) {
        config_filename = argv[1];
    }

    Config config;
    try {
        config.load(config_filename);

        if (config.OMEGA_BARYON > config.OMEGA_M) {
            throw std::runtime_error(
                "Error: OMEGA_BARYON cannot be larger than total OMEGA_M.");
        }

    } catch (const std::exception& e) {
        std::cerr << "Error loading " << config_filename << ": " << e.what()
                  << std::endl;
        return 1;
    }
    std::cout << "Successfully loaded " << config_filename << std::endl;

    std::cout << "N(P,G): (" << config.N_PER_SIDE << "³," << config.MESH_SIZE
              << "³) G: " << config.G << " fixed_dt: " << config.FIXED_DT
              << std::endl;

    // Create the output directories
    std::string timestamp = utils::get_timestamp();
    std::string run_dir = "outputs/run_" + timestamp;
    std::filesystem::create_directories(run_dir);

    // Copy the config file into the run directory for reproducibility
    std::filesystem::copy_file(
        "simulation.ini", run_dir + "/simulation.ini",
        std::filesystem::copy_options::overwrite_existing);

    HDF5Writer h5_writer(run_dir, config);
    Logger logger(run_dir);
    Diagnostics diagnostics;

    SimulationEngine engine(config, logger, h5_writer, diagnostics);
    g_engine = &engine;

    try {
        ExitStatus status = engine.run();
        switch (status) {
            case ExitStatus::UserAborted:
                std::cout << "\nSimulation aborted by user (Ctrl+C)."
                          << std::endl;
                break;
            case ExitStatus::ReachedMaxScaleFactor:
                std::cout << "\nSimulation successfully reached target a = "
                          << config.MAX_SCALE_FACTOR << "." << std::endl;
                break;
            case ExitStatus::ReachedMaxCycles:
                std::cout << "\nSimulation reached maximum allowed cycles ("
                          << config.MAX_CYCLES << ")" << std::endl;
                break;
        }
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL] Simulation crashed: " << e.what() << std::endl;
    }

    return 0;
}