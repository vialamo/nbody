#include "hdf5_writer.h"

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

void HDF5Writer::set_attr_string(H5::H5Object& obj, const char* attr_name,
                                 const std::string& str) {
    H5::DataSpace scalar_space(H5S_SCALAR);
    H5::StrType str_type(H5::PredType::C_S1, str.length() + 1);
    H5::Attribute attr = obj.createAttribute(attr_name, str_type, scalar_space);
    attr.write(str_type, str.c_str());
    attr.close();
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
                               const SimState& state, const Config& config,
                               const TimestepInfo& ts) {
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
        set_attr_double(header_group, "dt_macro", ts.dt_macro);
        set_attr_double(header_group, "dt_grav", ts.dt_grav);
        set_attr_double(header_group, "dt_hydro", ts.dt_hydro);
        set_attr_double(header_group, "dt_cool", ts.dt_cool);
        header_group.close();

        // ====================================================================
        // CONFIGURATION
        // ====================================================================
        H5::Group config_group = root_group.createGroup("Config");

        // [domain]
        set_attr_double(config_group, "domain_size", config.domain_size);
        set_attr_double(config_group, "box_size_mpc", config.box_size_mpc);
        set_attr_int(config_group, "mesh_size_1d", config.mesh_size);
        set_attr_int(config_group, "num_particles_1d", config.num_particles_1d);

        // [cosmology]
        set_attr_double(config_group, "omega_baryon", config.omega_baryon);
        set_attr_double(config_group, "omega_M", config.omega_m);
        set_attr_double(config_group, "omega_lambda", config.omega_lambda);
        set_attr_double(config_group, "Hubble_h", config.hubble_h);
        set_attr_double(config_group, "spectral_index", config.spectral_index);
        set_attr_double(config_group, "sigma_8", config.sigma_8);
        set_attr_bool(config_group, "expanding_universe",
                      config.expanding_universe);

        // [gravity]
        set_attr_double(config_group, "comoving_softening_factor",
                        config.comoving_softening_factor);
        set_attr_double(config_group, "softening_cap_scale_factor",
                        config.physical_softening_cap_a);

        // [initial_conditions]
        set_attr_string(config_group, "setup",
                        InitialConfig::to_string(config.initial_setup));
        set_attr_bool(config_group, "fixed_ics", config.fixed_ics);
        set_attr_bool(config_group, "invert_phases", config.invert_phases);
        set_attr_bool(config_group, "standing_particles",
                      config.standing_particles);
        set_attr_double(config_group, "initial_gas_temp_k",
                        config.initial_gas_temperature_k);
        set_attr_double(config_group, "seed_metallicity_solar",
                        config.seed_metallicity_solar);
        set_attr_int(config_group, "seed", config.seed);

        // [hydro]
        set_attr_string(config_group, "hydro_method",
                        HydroConfig::to_string(config.hydro_method));
        set_attr_double(config_group, "gamma", config.gamma);
        set_attr_bool(config_group, "enable_cooling", config.enable_cooling);
        set_attr_double(config_group, "primordial_mu", config.primordial_mu);
        set_attr_double(config_group, "temp_floor_k", config.temp_floor_k);
        set_attr_double(config_group, "cooling_cutoff_k",
                        config.cooling_cutoff_k);
        set_attr_string(config_group, "cooling_table_path",
                        config.cooling_table_path);

        // [MFM]
        set_attr_double(config_group, "mfm_target_neighbors",
                        config.mfm_target_neighbors);
        set_attr_double(config_group, "mfm_neighbor_tolerance",
                        config.mfm_neighbor_tolerance);
        set_attr_int(config_group, "mfm_max_iterations",
                     config.mfm_max_iterations);

        // [subgrid]
        set_attr_bool(config_group, "enable_subgrid_gravity",
                      config.enable_subgrid_gas_gravity);
        set_attr_bool(config_group, "enable_subgrid_clumping",
                      config.enable_subgrid_clumping);
        set_attr_double(config_group, "subgrid_clumping_amplitude",
                        config.subgrid_clumping_amplitude);

        // [p3m]
        set_attr_bool(config_group, "use_pm", config.use_PM);
        set_attr_bool(config_group, "use_pp", config.use_PP);
        set_attr_double(config_group, "pm_smoothing_cells",
                        config.PM_smoothing_cells);
        set_attr_double(config_group, "cutoff_radius_factor",
                        config.cutoff_radius_factor);

        // [time]
        set_attr_double(config_group, "max_dt_dynamical_factor",
                        config.max_dt_dynamical_factor);
        set_attr_double(config_group, "gravity_accuracy_eta",
                        config.gravity_accuracy_eta);
        set_attr_double(config_group, "hydro_courant_factor",
                        config.hydro_courant_factor);
        set_attr_bool(config_group, "use_adaptive_dt", config.use_adaptive_dt);
        set_attr_double(config_group, "a_start", config.a_start);
        set_attr_double(config_group, "a_end", config.a_end);

        // [HPC]
        set_attr_int(config_group, "num_threads", config.num_threads);
        set_attr_bool(config_group, "use_gpu", config.enable_GPU);

        // Global system physics
        set_attr_double(config_group, "G", config.G);
        config_group.close();

        // ====================================================================
        // PHYSICAL UNITS
        // ====================================================================
        H5::Group units_group = root_group.createGroup("Units");
        set_attr_double(units_group, "unit_length_in_mpc",
                        config.unit_length_mpc);
        set_attr_double(units_group, "unit_time_in_gyr", config.unit_time_gyr);
        set_attr_double(units_group, "unit_velocity_in_kms",
                        config.unit_velocity_kms);
        set_attr_double(units_group, "unit_velocity_in_cgs",
                        config.unit_velocity_cgs);
        set_attr_double(units_group, "unit_mass_in_msun",
                        config.unit_mass_msun);
        set_attr_double(units_group, "unit_density_in_cgs",
                        config.unit_density_cgs);
        set_attr_double(units_group, "factor_u_to_t", config.factor_u_to_t);
        units_group.close();

        H5::Group particle_group = root_group.createGroup("Particles");

        set_attr_double(particle_group, "cumulative_gravitational_work",
                        state.dm.accumulated_gravitational_work);
        set_attr_double(particle_group, "cumulative_expansion_work",
                        state.dm.accumulated_expansion_work);

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

        if (config.hydro_method != HydroMethod::None) {
            H5::Group gas_group = root_group.createGroup("Gas");
            if (config.hydro_method == HydroMethod::Eulerian) {
                GasGrid& gas = *state.gas;
                set_attr_double(gas_group, "cumulative_radiated_energy",
                                gas.get_accumulated_radiated_energy());
                set_attr_double(gas_group, "cumulative_photoheating_energy",
                                gas.get_accumulated_photoheating_energy());
                set_attr_double(gas_group, "cumulative_gravitational_work",
                                gas.get_accumulated_gravitational_work());
                set_attr_double(gas_group, "cumulative_expansion_work",
                                gas.get_accumulated_expansion_work());
                set_attr_double(
                    gas_group, "cumulative_dual_energy_switch_energy",
                    gas.get_accumulated_dual_energy_switch_energy());
                write_grid(gas_group, "density", gas.get_density());
                write_grid(gas_group, "momentum_x", gas.get_momentum_x());
                write_grid(gas_group, "momentum_y", gas.get_momentum_y());
                write_grid(gas_group, "momentum_z", gas.get_momentum_z());
                write_grid(gas_group, "energy", gas.get_energy());
                write_grid(gas_group, "pressure", gas.get_pressure());
                write_grid(gas_group, "metal_density", gas.get_metal_density());
                write_grid(gas_group, "thermal_timescale",
                           gas.compute_thermal_timescale(state.scale_factor,
                                                         state.cooling));

                Grid3D temp_grid(config.mesh_size);
                const Grid3D& rho = gas.get_density();
                const Grid3D& ie = gas.get_internal_energy();

                int total_cells =
                    config.mesh_size * config.mesh_size * config.mesh_size;
                for (int i = 0; i < total_cells; ++i) {
                    if (rho.data[i] > 1e-12) {
                        double u = ie.data[i] / rho.data[i];
                        temp_grid.data[i] =
                            Cooling::get_temp_from_internal_energy(
                                u, state.scale_factor, config);
                    } else {
                        temp_grid.data[i] = 10.0;  // Floor
                    }
                }
                write_grid(gas_group, "temperature", temp_grid);
            } else if (config.hydro_method == HydroMethod::MFM) {
                GasParticleSystem& gas = *state.mfm_gas;
                set_attr_double(gas_group, "cumulative_radiated_energy",
                                gas.accumulated_radiated_energy);
                set_attr_double(gas_group, "cumulative_photoheating_energy",
                                gas.accumulated_photoheating_energy);
                set_attr_double(gas_group, "cumulative_gravitational_work",
                                gas.accumulated_gravitational_work);
                set_attr_double(gas_group, "cumulative_expansion_work",
                                gas.accumulated_expansion_work);
                set_attr_double(gas_group, "cumulative_entropy_switch_energy",
                                gas.accumulated_entropy_switch_energy);

                write_particle_vec(gas_group, "position_x", gas.pos_x);
                write_particle_vec(gas_group, "position_y", gas.pos_y);
                write_particle_vec(gas_group, "position_z", gas.pos_z);
                write_particle_vec(gas_group, "velocity_x", gas.vel_x);
                write_particle_vec(gas_group, "velocity_y", gas.vel_y);
                write_particle_vec(gas_group, "velocity_z", gas.vel_z);
                write_particle_vec(gas_group, "acceleration_x", gas.acc_x);
                write_particle_vec(gas_group, "acceleration_y", gas.acc_y);
                write_particle_vec(gas_group, "acceleration_z", gas.acc_z);

                write_particle_vec(gas_group, "mass", gas.mass);
                write_particle_vec(gas_group, "smoothing_length", gas.h);
                write_particle_vec(gas_group, "density", gas.rho);
                write_particle_vec(gas_group, "pressure", gas.pressure);
                write_particle_vec(gas_group, "internal_energy", gas.u);
                write_particle_vec(gas_group, "metal_fraction", gas.metal_frac);

                std::vector<double> temp_vec(gas.num_particles);
                for (size_t i = 0; i < gas.num_particles; ++i) {
                    temp_vec[i] = Cooling::get_temp_from_internal_energy(
                        gas.u[i], state.scale_factor, config);
                }
                write_particle_vec(gas_group, "temperature", temp_vec);
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
