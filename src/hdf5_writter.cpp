#include "hdf5_writter.h"

#include <iostream>

#include "cooling.h"

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

void HDF5Writer::save_snapshot(int snapshot_index, int cycle_count,
                               const SimState& state, const Config& config) {
    char filename[256];
    snprintf(filename, sizeof(filename), "%s/snapshot_%04d.hdf5",
             output_directory.c_str(), snapshot_index);

    try {
        H5::H5File file(filename, H5F_ACC_TRUNC | H5F_ACC_SWMR_WRITE);

        H5::Group root_group = file.openGroup("/");

        // ====================================================================
        // SNAPSHOT HEADER
        // ====================================================================
        H5::Group header_group = root_group.createGroup("Header");
        set_attr_double(header_group, "scale_factor", state.scale_factor);
        set_attr_double(header_group, "simulation_time", state.total_time);
        header_group.close();

        // ====================================================================
        // CONFIGURATION
        // ====================================================================
        H5::Group config_group = root_group.createGroup("Config");
        
        // [domain]
        set_attr_double(config_group, "domain_size", config.DOMAIN_SIZE);
        set_attr_double(config_group, "box_size_mpc", config.BOX_SIZE_MPC);
        set_attr_int(config_group, "mesh_size", config.MESH_SIZE);
        
        // [cosmology]
        set_attr_double(config_group, "omega_baryon", config.OMEGA_BARYON);
        set_attr_double(config_group, "omega_M", config.OMEGA_M);
        set_attr_double(config_group, "omega_lambda", config.OMEGA_LAMBDA);
        set_attr_double(config_group, "hubble_param", config.HUBBLE_PARAM);
        set_attr_bool(config_group, "expanding_universe", config.EXPANDING_UNIVERSE);
        
        // [initial_conditions]
        set_attr_double(config_group, "spectral_index", config.SPECTRAL_INDEX);
        set_attr_double(config_group, "start_a", config.START_A);
        set_attr_double(config_group, "sigma_8", config.SIGMA_8);
        set_attr_int(config_group, "n_per_side", config.N_PER_SIDE);
        set_attr_bool(config_group, "standing_particles", config.STANDING_PARTICLES);
        set_attr_int(config_group, "seed", config.SEED);
        
        // [hydro]
        set_attr_bool(config_group, "use_hydro", config.USE_HYDRO);
        set_attr_double(config_group, "gamma", config.GAMMA);
        set_attr_bool(config_group, "enable_cooling", config.ENABLE_COOLING);
        set_attr_double(config_group, "primordial_mu", config.PRIMORDIAL_MU);
        set_attr_double(config_group, "temp_floor_k", config.TEMP_FLOOR_KELVIN);
        
        // [p3m]
        set_attr_bool(config_group, "use_pm", config.USE_PM);
        set_attr_bool(config_group, "use_pp", config.USE_PP);
        set_attr_double(config_group, "pm_smoothing_cells", config.PM_SMOOTHING_CELLS);
        set_attr_double(config_group, "cutoff_radius_factor", config.CUTOFF_RADIUS_FACTOR);
        
        // [time]
        set_attr_double(config_group, "dt_factor", config.DT_FACTOR);
        set_attr_double(config_group, "cfl_safety_factor", config.CFL_SAFETY_FACTOR);
        set_attr_bool(config_group, "use_adaptive_dt", config.USE_ADAPTIVE_DT);
        
        // Global system physics
        set_attr_double(config_group, "G", config.G);
        config_group.close();

        // ====================================================================
        // PHYSICAL UNITS
        // ====================================================================
        H5::Group units_group = root_group.createGroup("Units");
        set_attr_double(units_group, "unit_length_in_mpc", config.UNIT_LENGTH_MPC);
        set_attr_double(units_group, "unit_time_in_gyr", config.UNIT_TIME_GYR);
        set_attr_double(units_group, "unit_velocity_in_kms", config.UNIT_VELOCITY_KMS);
        set_attr_double(units_group, "unit_velocity_in_cgs", config.UNIT_VELOCITY_CGS);
        set_attr_double(units_group, "unit_mass_in_msun", config.UNIT_MASS_MSUN);
        set_attr_double(units_group, "unit_density_in_cgs", config.UNIT_DENSITY_CGS);
        set_attr_double(units_group, "factor_u_to_t", config.FACTOR_U_TO_T);
        units_group.close();

        H5::Group particle_group = root_group.createGroup("Particles");

        write_particle_vec(particle_group, "position_x", state.dm.pos_x);
        write_particle_vec(particle_group, "position_y", state.dm.pos_y);
        write_particle_vec(particle_group, "position_z", state.dm.pos_z);

        write_particle_vec(particle_group, "velocity_x", state.dm.vel_x);
        write_particle_vec(particle_group, "velocity_y", state.dm.vel_y);
        write_particle_vec(particle_group, "velocity_z", state.dm.vel_z);

        write_particle_vec(particle_group, "acceleration_x", state.dm.acc_x);
        write_particle_vec(particle_group, "acceleration_y", state.dm.acc_y);
        write_particle_vec(particle_group, "acceleration_z", state.dm.acc_z);

        write_particle_vec(particle_group, "mass", state.dm.mass);
        particle_group.close();

        if (config.USE_HYDRO) {
            H5::Group gas_group = root_group.createGroup("Gas");
            set_attr_double(gas_group, "cumulative_radiated_energy",
                            state.gas.get_accumulated_radiated_energy());
            write_grid(gas_group, "density", state.gas.get_density());
            write_grid(gas_group, "momentum_x", state.gas.get_momentum_x());
            write_grid(gas_group, "momentum_y", state.gas.get_momentum_y());
            write_grid(gas_group, "momentum_z", state.gas.get_momentum_z());
            write_grid(gas_group, "energy", state.gas.get_energy());
            write_grid(gas_group, "pressure", state.gas.get_pressure());

            if (config.ENABLE_COOLING) {
                Grid3D temp_grid(config.MESH_SIZE);
                const Grid3D& rho = state.gas.get_density();
                const Grid3D& ie = state.gas.get_internal_energy();

                int total_cells =
                    config.MESH_SIZE * config.MESH_SIZE * config.MESH_SIZE;
                for (int i = 0; i < total_cells; ++i) {
                    if (rho.data[i] > 1e-12) {
                        double u = ie.data[i] / rho.data[i];
                        temp_grid.data[i] =
                            cooling::get_temp_from_internal_energy(
                                u, state.scale_factor, config);
                    } else {
                        temp_grid.data[i] = 10.0;  // Floor
                    }
                }
                write_grid(gas_group, "temperature", temp_grid);
            }

            gas_group.close();
        }
        root_group.close();
        file.close();
    } catch (H5::Exception& e) {
        std::cerr << "Error: Could not save HDF5 snapshot." << "\n";
        e.printErrorStack();
    }
}
