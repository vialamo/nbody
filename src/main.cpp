#include <csignal>
#include <filesystem>
#include <iostream>
#include <string>

#include "config.h"
#include "diagnostics.h"
#include "engine.h"
#include "hdf5_writer.h"
#include "logger.h"
#include "utils.h"

SimulationEngine* g_engine = nullptr;

void signal_handler(int signal) {
    if (signal == SIGINT) {
        std::cout << "\n[Ctrl+C Detected] Finishing the current cycle and "
                     "shutting down..."
                  << std::endl;
        if (g_engine != nullptr) {
            g_engine->request_stop();
        }
    }
}

std::string generate_directory_name(const Config& config,
                                    const std::string& timestamp) {
    std::string base = InitialConfig::to_string(config.initial_setup);
    std::string run_dir =
        config.outputs_folder_path + "/" + base + "_" + timestamp;
    return run_dir;
}

int main(int argc, char* argv[]) {
    std::signal(SIGINT, signal_handler);

    std::string config_filename = "config.ini";
    if (argc >= 2) {
        config_filename = argv[1];
    }

    Config config;
    try {
        config.load(config_filename);
    } catch (const std::exception& e) {
        std::cerr << "Error loading " << config_filename << ": " << e.what()
                  << std::endl;
        return 1;
    }
    std::cout << "Successfully loaded " << config_filename << std::endl;

    if (config.num_threads > 0) {
        omp_set_num_threads(config.num_threads);
    }

    Logger::print_startup_info(config);

    // Create the output directories
    std::string timestamp = utils::get_timestamp();
    std::string run_dir = generate_directory_name(config, timestamp);
    std::filesystem::create_directories(run_dir);

    // Copy the config file into the run directory for reproducibility
    std::filesystem::copy_file(
        config_filename, run_dir + "/config.ini",
        std::filesystem::copy_options::overwrite_existing);

    HDF5Writer h5_writer(run_dir, config);
    Logger logger(run_dir);
    Diagnostics diagnostics;

    SimulationEngine engine(config, logger, h5_writer, diagnostics);
    g_engine = &engine;

    std::cout << "\nSimulation loop started: " << timestamp << "\n"
              << std::endl;
    try {
        ExitStatus status = engine.run();
        switch (status) {
            case ExitStatus::UserAborted:
                std::cout << "\nSimulation aborted by user (Ctrl+C)."
                          << std::endl;
                break;
            case ExitStatus::ReachedMaxScaleFactor:
                std::cout << "\nSimulation successfully reached target a = "
                          << config.a_end << "." << std::endl;
                break;
            case ExitStatus::ReachedMaxCycles:
                std::cout << "\nSimulation reached maximum allowed cycles ("
                          << config.max_cycles << ")" << std::endl;
                break;
        }
    } catch (const std::exception& e) {
        std::cerr << "\n[FATAL] Simulation crashed: " << e.what() << std::endl;
    }

    return 0;
}
