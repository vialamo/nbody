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