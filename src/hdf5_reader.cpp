#include "hdf5_reader.h"

#include <iostream>
#include <stdexcept>

HDF5Reader::HDF5Reader(const std::string& filepath) {
    try {
        // Open file in read-only mode
        file = H5::H5File(filepath, H5F_ACC_RDONLY);
    } catch (const H5::Exception& e) {
        std::cerr << "Error: Could not open HDF5 file: " << filepath << "\n";
        e.printErrorStack();
        throw std::runtime_error("HDF5 file open failed.");
    }
}

HDF5Reader::~HDF5Reader() { file.close(); }

H5::Group HDF5Reader::open_group(const std::string& group_name) const {
    return file.openGroup(group_name);
}

bool HDF5Reader::group_exists(const std::string& group_name) const {
    // H5Lexists checks if a link (like a group or dataset) exists
    return H5Lexists(file.getId(), group_name.c_str(), H5P_DEFAULT) > 0;
}

double HDF5Reader::read_attr_double(const H5::Group& group,
                                    const char* attr_name) const {
    double value = 0.0;
    H5::Attribute attr = group.openAttribute(attr_name);
    attr.read(H5::PredType::NATIVE_DOUBLE, &value);
    attr.close();
    return value;
}

int HDF5Reader::read_attr_int(const H5::Group& group,
                              const char* attr_name) const {
    int value = 0;
    H5::Attribute attr = group.openAttribute(attr_name);
    attr.read(H5::PredType::NATIVE_INT, &value);
    attr.close();
    return value;
}

bool HDF5Reader::read_attr_bool(const H5::Group& group,
                                const char* attr_name) const {
    return read_attr_int(group, attr_name) != 0;
}

std::string HDF5Reader::read_attr_string(const H5::Group& group,
                                         const char* attr_name) const {
    H5::Attribute attr = group.openAttribute(attr_name);
    H5::StrType str_type = attr.getStrType();
    std::string value;
    attr.read(str_type, value);
    attr.close();
    return value;
}

std::vector<double> HDF5Reader::read_particle_vec(
    const H5::Group& group, const char* dataset_name) const {
    H5::DataSet dataset = group.openDataSet(dataset_name);
    H5::DataSpace dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1) {
        throw std::runtime_error("Expected 1D dataset for particle vector.");
    }

    hsize_t dims_out[1];
    dataspace.getSimpleExtentDims(dims_out, NULL);

    std::vector<double> vec(dims_out[0]);
    dataset.read(vec.data(), H5::PredType::NATIVE_DOUBLE);

    dataset.close();
    dataspace.close();

    return vec;
}

MFMGasSnapshot HDF5Reader::read_mfm_gas() const {
    MFMGasSnapshot gas_data;

    if (!group_exists("/Gas")) {
        throw std::runtime_error("Snapshot does not contain a /Gas group.");
    }

    H5::Group gas_group = open_group("/Gas");

    // Read necessary fields for glass generation/resumption
    gas_data.pos_x = read_particle_vec(gas_group, "position_x");
    gas_data.pos_y = read_particle_vec(gas_group, "position_y");
    gas_data.pos_z = read_particle_vec(gas_group, "position_z");

    gas_data.vel_x = read_particle_vec(gas_group, "velocity_x");
    gas_data.vel_y = read_particle_vec(gas_group, "velocity_y");
    gas_data.vel_z = read_particle_vec(gas_group, "velocity_z");

    gas_data.mass = read_particle_vec(gas_group, "mass");
    gas_data.internal_energy = read_particle_vec(gas_group, "internal_energy");
    gas_data.smoothing_length =
        read_particle_vec(gas_group, "smoothing_length");
    gas_data.metal_fraction = read_particle_vec(gas_group, "metal_fraction");

    gas_data.num_particles = gas_data.pos_x.size();

    gas_group.close();

    return gas_data;
}

GlassData HDF5Reader::read_glass_cube(const std::string& filepath) {
    GlassData data;

    try {
        H5::H5File file(filepath, H5F_ACC_RDONLY);

        // Read Coordinates (Assuming shape [N, 3])
        H5::DataSet pos_dataset = file.openDataSet("/PartType0/Coordinates");
        H5::DataSpace pos_dataspace = pos_dataset.getSpace();
        hsize_t dims[2];
        pos_dataspace.getSimpleExtentDims(dims, NULL);

        data.num_particles = dims[0];
        std::vector<double> pos_buffer(data.num_particles * 3);
        pos_dataset.read(pos_buffer.data(), H5::PredType::NATIVE_DOUBLE);

        data.pos_x.reserve(data.num_particles);
        data.pos_y.reserve(data.num_particles);
        data.pos_z.reserve(data.num_particles);

        for (size_t i = 0; i < data.num_particles; ++i) {
            data.pos_x.push_back(pos_buffer[i * 3 + 0]);
            data.pos_y.push_back(pos_buffer[i * 3 + 1]);
            data.pos_z.push_back(pos_buffer[i * 3 + 2]);
        }

        // Read Smoothing Lengths
        H5::DataSet h_dataset = file.openDataSet("/PartType0/SmoothingLength");
        data.h.resize(data.num_particles);

        // SWIFT files sometimes use float instead of double for properties,
        // read into a float buffer if necessary, then cast.
        if (h_dataset.getFloatType().getSize() == 4) {
            std::vector<float> h_buffer(data.num_particles);
            h_dataset.read(h_buffer.data(), H5::PredType::NATIVE_FLOAT);
            for (size_t i = 0; i < data.num_particles; ++i) {
                // Apply the SWIFT 0.3 scaling factor directly on read
                data.h[i] = static_cast<double>(h_buffer[i]) * 0.3;
            }
        } else {
            h_dataset.read(data.h.data(), H5::PredType::NATIVE_DOUBLE);
            for (size_t i = 0; i < data.num_particles; ++i) {
                data.h[i] *= 0.3;
            }
        }

    } catch (const H5::Exception& e) {
        throw std::runtime_error("Failed to read glass file: " + filepath);
    }

    return data;
}