#include <atomic>
#include <csignal>
#include <filesystem>
#include <iostream>
#include <string>

#include "config.h"
#include "engine.h"
#include "utils.h"

SimulationEngine* g_engine = nullptr;

static void print_info(const Config& config) {
    std::cout << "\nInternal Unit System\n";
    std::cout << "1.0 Code Length   = " << config.unit_length_mpc
              << " comoving Mpc\n";
    std::cout << "1.0 Code Time     = " << config.unit_time_gyr << " Gyr\n";
    std::cout << "1.0 Code Velocity = " << config.unit_velocity_kms
              << " km/s\n";
    std::cout << "1.0 Code Mass     = " << config.unit_mass_msun
              << " Solar Masses\n";
    std::cout << "G = " << config.G << "\n" << std::endl;

    std::cout << "Mesh Size = " << config.mesh_size << "\n";
    std::cout << "Hydro Mode = " << HydroConfig::to_string(config.hydro_method)
              << "\n";
    std::cout << "Initial Setup = "
              << InitialConfig::to_string(config.initial_setup) << "\n";

    std::cout << "Threads: " << omp_get_max_threads() << " | ";
    int num_devices = omp_get_num_devices();
    if (num_devices == 0) {
        std::cout << "GPU: not available\n" << std::endl;
    } else {
        if (config.enable_GPU) {
            std::cout << "GPU: enabled\n" << std::endl;
        } else {
            std::cout << "GPU: disabled\n" << std::endl;
        }
    }
}

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
    } catch (const std::exception& e) {
        std::cerr << "Error loading " << config_filename << ": " << e.what()
                  << std::endl;
        return 1;
    }
    std::cout << "Successfully loaded " << config_filename << std::endl;

    if (config.num_threads > 0) {
        omp_set_num_threads(config.num_threads);
    }
    if (omp_get_num_devices() == 0) {
        config.enable_GPU = false;
    }

    print_info(config);

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
