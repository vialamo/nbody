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

// Container to hold the gas arrays from a snapshot
struct MFMGasSnapshot {
    size_t num_particles;
    std::vector<double> pos_x, pos_y, pos_z;
    std::vector<double> vel_x, vel_y, vel_z;
    std::vector<double> mass;
    std::vector<double> internal_energy;
    std::vector<double> smoothing_length;
    std::vector<double> metal_fraction;
};

struct GlassData {
    size_t num_particles;
    std::vector<double> pos_x;
    std::vector<double> pos_y;
    std::vector<double> pos_z;
    std::vector<double> h;
};

class HDF5Reader {
   private:
    H5::H5File file;

   public:
    HDF5Reader(const std::string& filepath);
    ~HDF5Reader();

    // Group management
    H5::Group open_group(const std::string& group_name) const;
    bool group_exists(const std::string& group_name) const;

    // Low-level attribute readers
    double read_attr_double(const H5::Group& group,
                            const char* attr_name) const;
    int read_attr_int(const H5::Group& group, const char* attr_name) const;
    bool read_attr_bool(const H5::Group& group, const char* attr_name) const;
    std::string read_attr_string(const H5::Group& group,
                                 const char* attr_name) const;

    // Low-level dataset readers
    std::vector<double> read_particle_vec(const H5::Group& group,
                                          const char* dataset_name) const;

    // Extensible high-level helpers
    MFMGasSnapshot read_mfm_gas() const;
    static GlassData read_glass_cube(const std::string& filepath);
};