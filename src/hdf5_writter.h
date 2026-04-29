#pragma once
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

    void save_snapshot(int snapshot_index, int cycle_count,
                       const SimState& state, const Config& config);
};