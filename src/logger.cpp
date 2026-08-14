#include "logger.h"

#include <sys/resource.h>
#include <unistd.h>  // For sysconf and _SC_PAGESIZE

#include <iomanip>
#include <iostream>
#include <sstream>

#include "config.h"
#include "diagnostics.h"

// Returns the Peak Resident Set Size (RSS) in Megabytes
static double get_peak_memory_usage_MB() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    // On Linux, ru_maxrss is returned in Kilobytes.
    // Divide by 1024 to get Megabytes.
    return static_cast<double>(usage.ru_maxrss) / 1024.0;
}

// Returns the Current Resident Set Size (RSS) in Megabytes
static double get_current_memory_usage_MB() {
    // /proc/self/statm contains memory stats for the current process
    std::ifstream statm("/proc/self/statm");
    if (!statm.is_open()) {
        return 0.0;  // Fallback in case of an error
    }

    long virtual_size, rss_pages;
    // The first value is total virtual memory, the second is current RSS in
    // pages
    statm >> virtual_size >> rss_pages;
    statm.close();

    // Linux allocates memory in "pages" (usually 4096 bytes).
    // We multiply the page count by the page size to get total bytes.
    long page_size_bytes = sysconf(_SC_PAGESIZE);

    // Convert bytes to Megabytes
    return static_cast<double>(rss_pages * page_size_bytes) / (1024.0 * 1024.0);
}

std::string Logger::format_double(double val, int precision, bool scientific) {
    std::stringstream ss;
    if (scientific) {
        ss << std::scientific;
    } else {
        ss << std::fixed;
    }
    ss << std::setprecision(precision) << val;
    return ss.str();
}

Logger::Logger(const std::string& run_dir)
    : start_time(std::chrono::high_resolution_clock::now()) {
    std::string log_filename = run_dir + "/diagnostics.csv";
    log_file.open(log_filename);
    if (!log_file.is_open()) {
        std::cerr << "Warning: Could not open log file: " << log_filename
                  << "\n";
    }
}

Logger::~Logger() {
    if (log_file.is_open()) {
        log_file.close();
    }
}

void Logger::write_header() {
    if (!log_file.is_open() || header_written) return;

    log_file
        << "cycle,sim_time,scale_factor,"
        << "total_mass_gas,total_mass_dm,total_momentum_gas_x,total_"
           "momentum_gas_y,total_momentum_gas_z,total_momentum_dm_x,total_"
           "momentum_dm_y,total_momentum_dm_z,"
        << "ke_gas,ke_dm,pe_dm,ie_gas,radiated_energy,heated_energy,"
        << "dt_cfl,dt_gravity,dt_cool,dt_final,"
        << "max_gas_density,max_gas_pressure,max_gas_velocity,"
        << "wall_time_total,wall_time_pm,wall_time_pp,wall_time_hydro,"
           "wall_time_cooling,wall_time_io,cumulative_wall_time,"
           "wall_time_GPU_transf,wall_time_GPU_back,wall_time_GPU_compute,"
           "memory_peak,memory_current,"
        << "avg_substeps_hydro,avg_substeps_gravity,avg_substeps_cooling\n";
    header_written = true;
}

void Logger::log(const Diagnostics& diag, const Config& conf) {
    if (!header_written) {
        write_header();
    }

    auto now = std::chrono::high_resolution_clock::now();
    double wall_time_s =
        std::chrono::duration_cast<std::chrono::duration<double> >(now -
                                                                   start_time)
            .count();

    int peak_mem = get_peak_memory_usage_MB();
    int curr_mem = get_current_memory_usage_MB();

    // Write to CSV file
    if (log_file.is_open()) {
        log_file << diag.cycle << "," << diag.sim_time << ","
                 << diag.scale_factor << "," << diag.total_mass_gas << ","
                 << diag.total_mass_dm << "," << diag.total_momentum_gas.x()
                 << "," << diag.total_momentum_gas.y() << ","
                 << diag.total_momentum_gas.z() << ","
                 << diag.total_momentum_dm.x() << ","
                 << diag.total_momentum_dm.y() << ","
                 << diag.total_momentum_dm.z() << "," << diag.ke_gas << ","
                 << diag.ke_dm << "," << diag.pe_total << "," << diag.ie_gas
                 << "," << diag.total_radiated_energy << ","
                 << diag.total_heated_energy << "," << diag.dt_cfl << ","
                 << diag.dt_gravity << "," << diag.dt_cool << ","
                 << diag.dt_final << "," << diag.max_gas_density << ","
                 << diag.max_gas_pressure << "," << diag.max_gas_velocity << ","
                 << diag.get_average(TimerRegion::Step) << ","
                 << diag.get_average(TimerRegion::PM) << ","
                 << diag.get_average(TimerRegion::PP) << ","
                 << diag.get_average(TimerRegion::Hydro) << ","
                 << diag.get_average(TimerRegion::Cool) << ","
                 << diag.get_io_time() << "," << wall_time_s << ","
                 << diag.get_prof_average(ProfRegion::Transf) / 1000.0 << ","
                 << diag.get_prof_average(ProfRegion::Ret) / 1000.0 << ","
                 << diag.get_prof_average(ProfRegion::Compute) / 1000.0 << ","
                 << peak_mem << "," << curr_mem << ","
                 << diag.get_average_substeps(SubstepCounter::Hydro) << ","
                 << diag.get_average_substeps(SubstepCounter::Gravity) << ","
                 << diag.get_average_substeps(SubstepCounter::Cool) << "\n";
    }

    // Write to stdout (console)
    std::cout << "[Cycle " << diag.cycle << "] "
              << "Mem (peak/curr MB): " << peak_mem << ", " << curr_mem << " | "
              << "Wall: " << format_double(wall_time_s, 1) << "s | "
              << "SimTime: " << format_double(diag.sim_time, 3) << " | "
              << "a: " << format_double(diag.scale_factor, 3) << "\n";

    std::cout << "  [Physics]" << "\n";
    std::cout << "    - Mass Err: " << format_double(diag.mass_err, 1, true)
              << " | Momentum: " << format_double(diag.total_momentum, 1, true)
              << " | DM Err: " << format_double(diag.dm_energy_err, 1, true)
              << " | Gas Err: " << format_double(diag.energy_err, 1, true)
              << "\n";

    std::cout << "  [Stability]\n";
    std::cout << "    - Timestep (Macro): "
              << format_double(diag.dt_final, 1, true)
              << " | Grav: " << format_double(diag.dt_gravity, 1, true)
              << " | CFL: " << format_double(diag.dt_cfl, 1, true)
              << " | Cool: " << format_double(diag.dt_cool, 1, true) << "\n";
    std::cout
        << "    - Subcycles/step:   Hydro: "
        << format_double(diag.get_average_substeps(SubstepCounter::Hydro), 1)
        << " | Gravity: "
        << format_double(diag.get_average_substeps(SubstepCounter::Gravity), 1)
        << " | Cooling: "
        << format_double(diag.get_average_substeps(SubstepCounter::Cool), 1)
        << "\n";

    std::cout << "  [Performance (ms / cycle)]" << "\n";
    std::cout << "    - Step: "
              << format_double(diag.get_average(TimerRegion::Step) * 1000.0, 0)
              << " | Overhead: "
              << format_double(diag.get_average_overhead() * 1000.0, 0)
              << " | PM: "
              << format_double(diag.get_average(TimerRegion::PM) * 1000.0, 0)
              << " | PP: "
              << format_double(diag.get_average(TimerRegion::PP) * 1000.0, 0)
              << " | Hydro: "
              << format_double(diag.get_average(TimerRegion::Hydro) * 1000.0, 0)
              << " | Cooling: "
              << format_double(diag.get_average(TimerRegion::Cool) * 1000.0, 0)
              << "\n"
              << "    - I/O Spike: "
              << format_double(diag.get_io_time() * 1000.0, 0) << "\n";
    if (diag.get_prof_average(ProfRegion::Compute) > 0.0) {
        std::cout << "    - GPU Transfer: "
                  << format_double(
                         diag.get_prof_average(ProfRegion::Transf) / 1000.0, 1)
                  << " | Return: "
                  << format_double(
                         diag.get_prof_average(ProfRegion::Ret) / 1000.0, 1)
                  << " | Compute: "
                  << format_double(
                         diag.get_prof_average(ProfRegion::Compute) / 1000.0, 1)
                  << "\n";
    }
    std::cout << "-------------------------------------------------------------"
                 "---------"
              << "\n";

    if (diag.non_converged_cooling_cells > 0) {
        std::cout << "  [WARNING] Cooling solver failed to converge in "
                  << diag.non_converged_cooling_cells << " cells\n";
    }
}
