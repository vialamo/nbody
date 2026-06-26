#include <atomic>
#include <csignal>
#include <filesystem>
#include <iostream>
#include <string>

#include "config.h"
#include "engine.h"
#include "utils.h"

SimulationEngine* g_engine = nullptr;

void print_info(const Config& config) {
    std::cout << "Mesh Size = " << config.MESH_SIZE << "\n\n";
    std::cout << "Internal Unit System\n";
    std::cout << "1.0 Code Length   = " << config.UNIT_LENGTH_MPC
              << " comoving Mpc\n";
    std::cout << "1.0 Code Time     = " << config.UNIT_TIME_GYR << " Gyr\n";
    std::cout << "1.0 Code Velocity = " << config.UNIT_VELOCITY_KMS
              << " km/s\n";
    std::cout << "1.0 Code Mass     = " << config.UNIT_MASS_MSUN
              << " Solar Masses\n";
    std::cout << "G = " << config.G << "\n" << std::endl;

    constexpr double sigma_safe = 0.3;
    double a_safe_max =
        (sigma_safe / config.SIGMA_8) * pow((config.BOX_SIZE_MPC / 8.0), 0.9);
    if (a_safe_max < config.MAX_SCALE_FACTOR) {
        std::cout << "Warning: Fundamental mode collapse risk when a > "
                  << a_safe_max << std::endl;
    }

    int num_devices = omp_get_num_devices();
    if (num_devices == 0) {
        std::cout << "Warning: OpenMP sees NO GPUs. Falling back to CPU...\n"
                  << std::endl;
    } else {
        std::cout << "OpenMP GPU available\n" << std::endl;
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