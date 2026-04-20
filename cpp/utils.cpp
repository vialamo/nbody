#include "utils.h"

#include <sys/resource.h>
#include <unistd.h>  // For sysconf and _SC_PAGESIZE

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "gas.h"        // For get_cfl_timestep
#include "particles.h"  // For calculate_particle_energies and get_gravity_timestep

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

// ---------------------------------------------------------
// Helper: Timestamp Generator
// ---------------------------------------------------------
std::string get_timestamp() {
    std::time_t now = std::time(nullptr);
    std::tm ltm;
#ifdef _MSC_VER
    localtime_s(&ltm, &now);
#else
    ltm = *std::localtime(&now);
#endif
    std::stringstream ss;
    ss << std::put_time(&ltm, "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}

// ---------------------------------------------------------
// Logger Implementation
// ---------------------------------------------------------
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
        << "ke_gas,ke_dm,pe_dm,ie_gas,"
        << "dt_cfl,dt_gravity,dt_final,"
        << "max_gas_density,max_gas_pressure,max_gas_velocity,"
        << "wall_time_total,wall_time_pm,wall_time_pp,wall_time_hydro,"
           "wall_time_io,cumulative_wall_time,memory_peak,memory_current\n";
    header_written = true;
}

void Logger::log(const Diagnostics& diag) {
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
                 << diag.total_mass_dm << "," << diag.total_momentum_gas.x
                 << "," << diag.total_momentum_gas.y << ","
                 << diag.total_momentum_gas.z << "," << diag.total_momentum_dm.x
                 << "," << diag.total_momentum_dm.y << ","
                 << diag.total_momentum_dm.z << "," << diag.ke_gas << ","
                 << diag.ke_dm << "," << diag.pe_dm << "," << diag.ie_gas << ","
                 << diag.dt_cfl << "," << diag.dt_gravity << ","
                 << diag.dt_final << "," << diag.max_gas_density << ","
                 << diag.max_gas_pressure << "," << diag.max_gas_velocity << ","
                 << diag.wall_time_total << "," << diag.wall_time_pm << ","
                 << diag.wall_time_pp << "," << diag.wall_time_hydro << ","
                 << diag.wall_time_io << "," << wall_time_s << "," << peak_mem
                 << "," << curr_mem << "\n";
    }

    // Write to stdout (console)
    double mass_err = diag.total_mass() - 1.0;

    std::cout << "[Cycle " << diag.cycle << "] "
              << "Mem (peak/curr MB): " << peak_mem << ", " << curr_mem << " | "
              << "Wall: " << format_double(wall_time_s, 1) << "s | "
              << "SimTime: " << format_double(diag.sim_time, 3) << " | "
              << "a: " << format_double(diag.scale_factor, 3) << "\n";

    std::cout << "  [Physics]" << "\n";
    std::cout << "    - Mass (P/G/T):   "
              << format_double(diag.total_mass_dm, 4) << " | "
              << format_double(diag.total_mass_gas, 4) << " | "
              << format_double(diag.total_mass(), 4)
              << " (Err: " << std::scientific << std::setprecision(1)
              << mass_err << std::fixed << ")" << "\n";

    std::cout << "    - Momentum (P/G): ("
              << format_double(diag.total_momentum_dm.x, 1, true) << ", "
              << format_double(diag.total_momentum_dm.y, 1, true) << ", "
              << format_double(diag.total_momentum_dm.z, 1, true) << ") | ("
              << format_double(diag.total_momentum_gas.x, 1, true) << ", "
              << format_double(diag.total_momentum_gas.y, 1, true) << ", "
              << format_double(diag.total_momentum_gas.z, 1, true) << ")"
              << "\n";

    std::cout << "    - Energy (KE/PE/IE): "
              << format_double(diag.ke_dm + diag.ke_gas, 3, true) << " | "
              << format_double(diag.pe_dm, 3, true) << " | "
              << format_double(diag.ie_gas, 3, true)
              << " (Total: " << format_double(diag.total_energy(), 3, true)
              << ")" << "\n";

    std::cout << "  [Stability]" << "\n";
    std::cout << "    - Timestep (CFL): " << format_double(diag.dt_cfl, 2, true)
              << " | (Grav): " << format_double(diag.dt_gravity, 2, true)
              << " | (Final): " << format_double(diag.dt_final, 2, true)
              << "\n";
    std::cout << "    - Max(rho): "
              << format_double(diag.max_gas_density, 2, true)
              << " | Max(press): "
              << format_double(diag.max_gas_pressure, 2, true)
              << " | Max(vel): "
              << format_double(diag.max_gas_velocity, 2, true) << "\n";

    std::cout << "  [Performance (ms)]" << "\n";
    std::cout << "    - PM: " << format_double(diag.wall_time_pm * 1000.0, 1)
              << " | PP: " << format_double(diag.wall_time_pp * 1000.0, 1)
              << " | Hydro: " << format_double(diag.wall_time_hydro * 1000.0, 1)
              << " | I/O: " << format_double(diag.wall_time_io * 1000.0, 1)
              << " | Total: " << format_double(diag.wall_time_total * 1000.0, 1)
              << "\n";
    std::cout << "-------------------------------------------------------------"
                 "---------"
              << "\n";
}

// ---------------------------------------------------------
// HDF5Writer Implementation
// ---------------------------------------------------------
void HDF5Writer::set_attr_double(H5::H5Object& obj, const char* attr_name,
                                 double value) {
    H5::DataSpace scalar_space(H5S_SCALAR);
    H5::Attribute attr = obj.createAttribute(
        attr_name, H5::PredType::NATIVE_DOUBLE, scalar_space);
    attr.write(H5::PredType::NATIVE_DOUBLE, &value);
    attr.close();
}

void HDF5Writer::set_attr_int(H5::H5Object& obj, const char* attr_name,
                              int value) {
    H5::DataSpace scalar_space(H5S_SCALAR);
    H5::Attribute attr =
        obj.createAttribute(attr_name, H5::PredType::NATIVE_INT, scalar_space);
    attr.write(H5::PredType::NATIVE_INT, &value);
    attr.close();
}

void HDF5Writer::set_attr_bool(H5::H5Object& obj, const char* attr_name,
                               bool value) {
    int int_val = value ? 1 : 0;
    set_attr_int(obj, attr_name, int_val);
}

void HDF5Writer::write_grid(H5::Group& group, const char* dataset_name,
                            const Grid3D& grid) {
    hsize_t N = grid.n;
    hsize_t dims[3] = {N, N, N};

    H5::DataSpace dataspace(3, dims);

    H5::DataSet dataset = group.createDataSet(
        dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(grid.raw_data(), H5::PredType::NATIVE_DOUBLE);
    dataset.close();
}

void HDF5Writer::write_particle_vec(H5::Group& group, const char* dataset_name,
                                    const std::vector<double>& vec) {
    hsize_t dims[1] = {vec.size()};
    H5::DataSpace dataspace(1, dims);
    H5::DataSet dataset = group.createDataSet(
        dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(vec.data(), H5::PredType::NATIVE_DOUBLE);
    dataset.close();
}

HDF5Writer::HDF5Writer(const std::string& run_dir, const Config& config)
    : output_directory(run_dir) {}

HDF5Writer::~HDF5Writer() {}

double HDF5Writer::save_snapshot(int snapshot_index, const SimState& state,
                                 const Config& config) {
    auto start_time = std::chrono::high_resolution_clock::now();

    char filename[256];
    snprintf(filename, sizeof(filename), "%s/snapshot_%04d.hdf5",
             output_directory.c_str(), snapshot_index);

    try {
        H5::H5File file(filename, H5F_ACC_TRUNC | H5F_ACC_SWMR_WRITE);

        H5::Group root_group = file.openGroup("/");

        set_attr_double(root_group, "domain_size", config.DOMAIN_SIZE);
        set_attr_int(root_group, "mesh_size", config.MESH_SIZE);
        set_attr_double(root_group, "omega_baryon", config.OMEGA_BARYON);
        set_attr_int(root_group, "n_per_side", config.N_PER_SIDE);
        set_attr_bool(root_group, "use_hydro", config.USE_HYDRO);
        set_attr_double(root_group, "g_const", config.G);
        set_attr_double(root_group, "gamma", config.GAMMA);
        set_attr_bool(root_group, "standing_particles",
                      config.STANDING_PARTICLES);
        set_attr_bool(root_group, "expanding_universe",
                      config.EXPANDING_UNIVERSE);
        set_attr_double(root_group, "expansion_start_t",
                        config.EXPANSION_START_T);
        set_attr_double(root_group, "initial_power_spectrum_index",
                        config.INITIAL_POWER_SPECTRUM_INDEX);
        set_attr_bool(root_group, "use_pm", config.USE_PM);
        set_attr_bool(root_group, "use_pp", config.USE_PP);
        set_attr_double(root_group, "cutoff_radius_cells",
                        config.CUTOFF_RADIUS_CELLS);
        set_attr_double(root_group, "dt_factor", config.DT_FACTOR);
        set_attr_double(root_group, "cfl_safety_factor",
                        config.CFL_SAFETY_FACTOR);
        set_attr_double(root_group, "gravity_dt_factor",
                        config.GRAVITY_DT_FACTOR);
        set_attr_bool(root_group, "use_adaptive_dt", config.USE_ADAPTIVE_DT);
        set_attr_int(root_group, "seed", config.SEED);

        set_attr_double(root_group, "simulation_time", state.total_time);
        set_attr_double(root_group, "scale_factor", state.scale_factor);

        H5::Group particle_group = root_group.createGroup("particles");
        size_t n_particles = state.dm.particles.size();

        std::vector<double> pos_x(n_particles), pos_y(n_particles),
            pos_z(n_particles);
        std::vector<double> vel_x(n_particles), vel_y(n_particles),
            vel_z(n_particles);
        std::vector<double> acc_x(n_particles), acc_y(n_particles),
            acc_z(n_particles);
        std::vector<double> mass(n_particles);

        for (size_t i = 0; i < n_particles; ++i) {
            pos_x[i] = state.dm.particles[i].pos.x;
            pos_y[i] = state.dm.particles[i].pos.y;
            pos_z[i] = state.dm.particles[i].pos.z;

            vel_x[i] = state.dm.particles[i].vel.x;
            vel_y[i] = state.dm.particles[i].vel.y;
            vel_z[i] = state.dm.particles[i].vel.z;

            acc_x[i] = state.dm.particles[i].acc.x;
            acc_y[i] = state.dm.particles[i].acc.y;
            acc_z[i] = state.dm.particles[i].acc.z;

            mass[i] = state.dm.particles[i].mass;
        }

        write_particle_vec(particle_group, "position_x", pos_x);
        write_particle_vec(particle_group, "position_y", pos_y);
        write_particle_vec(particle_group, "position_z", pos_z);

        write_particle_vec(particle_group, "velocity_x", vel_x);
        write_particle_vec(particle_group, "velocity_y", vel_y);
        write_particle_vec(particle_group, "velocity_z", vel_z);

        write_particle_vec(particle_group, "acceleration_x", acc_x);
        write_particle_vec(particle_group, "acceleration_y", acc_y);
        write_particle_vec(particle_group, "acceleration_z", acc_z);

        write_particle_vec(particle_group, "mass", mass);
        particle_group.close();

        if (config.USE_HYDRO) {
            H5::Group gas_group = root_group.createGroup("gas");
            write_grid(gas_group, "density", state.gas.get_density());
            write_grid(gas_group, "momentum_x", state.gas.get_momentum_x());
            write_grid(gas_group, "momentum_y", state.gas.get_momentum_y());
            write_grid(gas_group, "momentum_z", state.gas.get_momentum_z());
            write_grid(gas_group, "energy", state.gas.get_energy());
            write_grid(gas_group, "pressure", state.gas.get_pressure());
            gas_group.close();
        }
        root_group.close();
        file.close();
    } catch (H5::Exception& e) {
        std::cerr << "Error: Could not save HDF5 snapshot." << "\n";
        e.printErrorStack();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double> >(
               end_time - start_time)
        .count();
}

// ---------------------------------------------------------
// Diagnostics Calculation
// ---------------------------------------------------------
Diagnostics calculate_diagnostics(const SimState& state,
                                  const std::map<std::string, double>& timings,
                                  double dt, int cycle, const Config& config) {
    Diagnostics diag;
    diag.cycle = cycle;
    diag.sim_time = state.total_time;
    diag.scale_factor = state.scale_factor;

    // Performance
    diag.wall_time_pm = timings.count("pm") ? timings.at("pm") : 0.0;
    diag.wall_time_pp = timings.count("pp") ? timings.at("pp") : 0.0;
    diag.wall_time_hydro = timings.count("hydro") ? timings.at("hydro") : 0.0;
    diag.wall_time_io = timings.count("io") ? timings.at("io") : 0.0;
    diag.wall_time_total = diag.wall_time_pm + diag.wall_time_pp +
                           diag.wall_time_hydro + diag.wall_time_io;

    // Stability
    diag.dt_cfl = state.gas.get_cfl_timestep();
    diag.dt_gravity = state.dm.get_gravity_timestep(config);
    diag.dt_final = dt;

    if (config.USE_HYDRO) {
        diag.max_gas_density = state.gas.get_density().maxCoeff();
        diag.max_gas_pressure = state.gas.get_pressure().maxCoeff();
        diag.max_gas_velocity = (state.gas.get_velocity_x().array().square() +
                                 state.gas.get_velocity_y().array().square() +
                                 state.gas.get_velocity_z().array().square())
                                    .sqrt()
                                    .maxCoeff();
    }

    // Conservation (Particles)
    for (const auto& p : state.dm.particles) {
        diag.total_mass_dm += p.mass;
        diag.total_momentum_dm.x += p.mass * p.vel.x;
        diag.total_momentum_dm.y += p.mass * p.vel.y;
        diag.total_momentum_dm.z += p.mass * p.vel.z;
    }
    auto energies = state.dm.calculate_energies(state.scale_factor, config);
    diag.ke_dm = energies.first;
    diag.pe_dm = energies.second;

    // Conservation (Gas)
    if (config.USE_HYDRO) {
        diag.total_mass_gas =
            state.gas.get_density().sum() * config.CELL_VOLUME;
        diag.total_momentum_gas.x = state.gas.get_momentum_x().sum();
        diag.total_momentum_gas.y = state.gas.get_momentum_y().sum();
        diag.total_momentum_gas.z = state.gas.get_momentum_z().sum();

        Grid3D ke_gas_density(config.MESH_SIZE);
        ke_gas_density.data =
            0.5 * (state.gas.get_momentum_x().array().square() +
                   state.gas.get_momentum_y().array().square() +
                   state.gas.get_momentum_z().array().square());

        int total_cells =
            config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;
        for (int i = 0; i < total_cells; ++i) {
            if (state.gas.get_density().data[i] > 1e-12) {
                ke_gas_density.data[i] /= state.gas.get_density().data[i];
            } else {
                ke_gas_density.data[i] = 0.0;
            }
        }
        diag.ke_gas = ke_gas_density.sum() * config.CELL_VOLUME;

        Grid3D internal_energy_density(config.MESH_SIZE);
        internal_energy_density.data =
            state.gas.get_pressure().array() / (config.GAMMA - 1.0);
        diag.ie_gas = internal_energy_density.sum() * config.CELL_VOLUME;
    }

    return diag;
}