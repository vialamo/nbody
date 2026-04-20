#pragma once
#include <chrono>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4251)
#endif
#include <H5Cpp.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "config.h"
#include "state.h"

std::string get_timestamp();

struct Diagnostics {
    int cycle = 0;
    double sim_time = 0.0;
    double scale_factor = 1.0;
    double total_mass_gas = 0.0;
    double total_mass_dm = 0.0;
    Vec3 total_momentum_gas = {0.0, 0.0, 0.0};
    Vec3 total_momentum_dm = {0.0, 0.0, 0.0};
    double ke_gas = 0.0, ke_dm = 0.0, pe_dm = 0.0, ie_gas = 0.0;
    double dt_cfl = 0.0, dt_gravity = 0.0, dt_final = 0.0;
    double max_gas_density = 0.0, max_gas_pressure = 0.0,
           max_gas_velocity = 0.0;
    double wall_time_total = 0.0, wall_time_pm = 0.0, wall_time_pp = 0.0,
           wall_time_hydro = 0.0, wall_time_io = 0.0;

    double total_mass() const { return total_mass_gas + total_mass_dm; }
    Vec3 total_momentum() const {
        return {total_momentum_gas.x + total_momentum_dm.x,
                total_momentum_gas.y + total_momentum_dm.y,
                total_momentum_gas.z + total_momentum_dm.z};
    }
    double total_energy() const { return ke_gas + ke_dm + pe_dm + ie_gas; }
};

class Logger {
    std::ofstream log_file;
    std::chrono::high_resolution_clock::time_point start_time;
    bool header_written = false;
    static std::string format_double(double val, int precision,
                                     bool scientific = false);

   public:
    Logger(const std::string& run_dir);
    ~Logger();
    void write_header();
    void log(const Diagnostics& diag);
};

class HDF5Writer {
   private:
    std::string output_directory;

    void set_attr_double(H5::H5Object& obj, const char* attr_name,
                         double value);
    void set_attr_int(H5::H5Object& obj, const char* attr_name, int value);
    void set_attr_bool(H5::H5Object& obj, const char* attr_name, bool value);
    void write_grid(H5::Group& group, const char* dataset_name,
                    const Grid3D& grid);
    void write_particle_vec(H5::Group& group, const char* dataset_name,
                            const std::vector<double>& vec);

   public:
    HDF5Writer(const std::string& run_dir, const Config& config);
    ~HDF5Writer();

    double save_snapshot(int snapshot_index, const SimState& state,
                         const Config& config);
};

Diagnostics calculate_diagnostics(const SimState& state,
                                  const std::map<std::string, double>& timings,
                                  double dt, int cycle, const Config& config);
