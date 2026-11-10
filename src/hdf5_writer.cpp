#include "hdf5_writer.h"

#include <iostream>
#include <thread>

#include "config.h"
#include "cooling.h"
#include "state.h"
#include "thermodynamics.h"
#include "types.h"

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
    hsize_t dims[3] = {static_cast<hsize_t>(grid.nx),
                       static_cast<hsize_t>(grid.ny),
                       static_cast<hsize_t>(grid.nz)};

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

void HDF5Writer::write_particle_vec3d(H5::Group& group,
                                      const char* dataset_name,
                                      const std::vector<Eigen::Vector3d>& vec) {
    if (vec.empty()) return;

    // Create a 2D dataset: N particles x 3 spatial components
    hsize_t dims[2] = {vec.size(), 3};
    H5::DataSpace dataspace(2, dims);

    H5::DataSet dataset = group.createDataSet(
        dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);

    // vec.data() returns a pointer to the contiguous block of Eigen::Vector3d
    // structs
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

    constexpr int max_retries = 3;
    bool success = false;

    for (int attempt = 1; attempt <= max_retries; ++attempt) {
        try {
            H5::H5File file(filename, H5F_ACC_TRUNC | H5F_ACC_SWMR_WRITE);

            H5::Group root_group = file.openGroup("/");

            // ====================================================================
            // SNAPSHOT HEADER
            // ====================================================================
            H5::Group header_group = root_group.createGroup("Header");
            set_attr_double(header_group, "scale_factor", state.scale_factor);
            set_attr_double(header_group, "simulation_time", state.total_time);
            set_attr_double(header_group, "softening_length",
                            state.gravity.softening_length);
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
            set_attr_double(config_group, "domain_size", config.domain_size_1d);
            set_attr_double(config_group, "domain_size_x",
                            config.domain_size.x);
            set_attr_double(config_group, "domain_size_y",
                            config.domain_size.y);
            set_attr_double(config_group, "domain_size_z",
                            config.domain_size.z);
            set_attr_double(config_group, "box_size_mpc",
                            config.box_size_mpc_1d);
            set_attr_int(config_group, "mesh_size_1d", config.mesh_size_1d);
            set_attr_int(config_group, "num_dm_particles_1d",
                         config.num_dm_particles_1d);
            set_attr_int(config_group, "num_gas_particles_1d",
                         config.num_gas_particles_1d);

            // [cosmology]
            set_attr_double(config_group, "omega_baryon", config.omega_baryon);
            set_attr_double(config_group, "omega_m", config.omega_m);
            set_attr_double(config_group, "omega_lambda", config.omega_lambda);
            set_attr_double(config_group, "hubble_h", config.hubble_h);
            set_attr_double(config_group, "spectral_index",
                            config.spectral_index);
            set_attr_double(config_group, "sigma_8", config.sigma_8);
            set_attr_bool(config_group, "expanding_universe",
                          config.expanding_universe);

            // [gravity]
            set_attr_double(config_group, "comoving_softening_factor",
                            config.comoving_softening_factor);
            set_attr_double(config_group, "softening_cap_scale_factor",
                            config.softening_cap_a);
            set_attr_bool(config_group, "use_pm", config.use_PM);
            set_attr_bool(config_group, "use_pp", config.use_PP);
            set_attr_double(config_group, "pm_smoothing_cells",
                            config.PM_smoothing_cells);
            set_attr_double(config_group, "cutoff_radius_factor",
                            config.cutoff_radius_factor);

            // [initial_conditions]
            set_attr_string(config_group, "setup",
                            InitialConfig::to_string(config.initial_setup));
            set_attr_bool(config_group, "use_glass", config.use_glass);
            set_attr_bool(config_group, "fixed_ics", config.fixed_ics);
            set_attr_bool(config_group, "invert_phases", config.invert_phases);
            set_attr_double(config_group, "initial_gas_temp_k",
                            config.initial_gas_temperature_k);
            set_attr_double(config_group, "seed_metallicity_solar",
                            config.seed_metallicity_solar);
            set_attr_int(config_group, "seed", config.seed);
            set_attr_string(config_group, "glass_folder_path",
                            config.glass_folder_path);

            // [hydro]
            set_attr_double(config_group, "target_neighbors",
                            config.target_neighbors);
            set_attr_double(config_group, "neighbor_tolerance",
                            config.neighbor_tolerance);
            set_attr_int(config_group, "max_iterations", config.max_iterations);
            set_attr_double(config_group, "min_hsml_fraction",
                            config.min_hsml_fraction);
            set_attr_double(config_group, "max_condition_number",
                            config.max_condition_number);
            set_attr_bool(config_group, "use_exact_riemann",
                          config.use_exact_riemann);
            set_attr_double(config_group, "gradient_limiter_beta_max",
                            config.gradient_limiter_beta_max);
            set_attr_double(config_group, "gradient_limiter_beta_min",
                            config.gradient_limiter_beta_min);
            set_attr_double(config_group, "pairwise_limiter_psi1",
                            config.pairwise_limiter_psi1);
            set_attr_double(config_group, "pairwise_limiter_psi2",
                            config.pairwise_limiter_psi2);

            // [thermodynamics]
            set_attr_double(config_group, "gamma", config.gamma);
            set_attr_double(config_group, "primordial_mu",
                            config.primordial_mu);
            set_attr_double(config_group, "temp_floor_k", config.temp_floor_k);
            set_attr_double(config_group, "dual_energy_threshold",
                            config.dual_energy_threshold);
            set_attr_double(config_group, "entropy_switch_alpha_kin",
                            config.entropy_switch_akin);
            set_attr_double(config_group, "entropy_switch_alpha_grav",
                            config.entropy_switch_agrav);

            // [cooling]
            set_attr_bool(config_group, "enable_cooling",
                          config.enable_cooling);
            set_attr_double(config_group, "cooling_cutoff_k",
                            config.cooling_cutoff_k);
            set_attr_string(config_group, "cooling_table_path",
                            config.cooling_table_path);

            // [time]
            set_attr_double(config_group, "max_dt_dynamical_factor",
                            config.max_dt_dynamical_factor);
            set_attr_double(config_group, "gravity_accuracy_eta",
                            config.gravity_accuracy_eta);
            set_attr_double(config_group, "hydro_courant_factor",
                            config.hydro_courant_factor);
            set_attr_double(config_group, "a_start", config.a_start);
            set_attr_double(config_group, "a_end", config.a_end);
            set_attr_bool(config_group, "individual_particle_time",
                          config.individual_particle_timesteps);

            // [HPC]
            set_attr_int(config_group, "num_threads", config.num_threads);
            set_attr_bool(config_group, "use_gpu", config.using_GPU);

            // Global system physics
            set_attr_double(config_group, "G", config.G);
            config_group.close();

            // ====================================================================
            // PHYSICAL UNITS
            // ====================================================================
            H5::Group units_group = root_group.createGroup("Units");
            set_attr_double(units_group, "unit_length_in_mpc",
                            config.unit_length_mpc);
            set_attr_double(units_group, "unit_time_in_gyr",
                            config.unit_time_gyr);
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
                            state.dm.budget.gravitational_work);
            set_attr_double(particle_group, "cumulative_expansion_work",
                            state.dm.budget.expansion_work);

            write_particle_vec(particle_group, "position_x", state.dm.pos_x);
            write_particle_vec(particle_group, "position_y", state.dm.pos_y);
            write_particle_vec(particle_group, "position_z", state.dm.pos_z);

            write_particle_vec(particle_group, "velocity_x", state.dm.vel_x);
            write_particle_vec(particle_group, "velocity_y", state.dm.vel_y);
            write_particle_vec(particle_group, "velocity_z", state.dm.vel_z);

            write_particle_vec(particle_group, "acceleration_x",
                               state.dm.acc_x);
            write_particle_vec(particle_group, "acceleration_y",
                               state.dm.acc_y);
            write_particle_vec(particle_group, "acceleration_z",
                               state.dm.acc_z);

            write_particle_vec(particle_group, "mass", state.dm.mass);
            particle_group.close();

            H5::Group gas_group = root_group.createGroup("Gas");
            const GasParticleSystem& gas = state.gas;
            set_attr_double(gas_group, "cumulative_radiated_energy",
                            gas.budget.radiated);
            set_attr_double(gas_group, "cumulative_photoheating_energy",
                            gas.budget.photoheated);
            set_attr_double(gas_group, "cumulative_gravitational_work",
                            gas.budget.gravitational_work);
            set_attr_double(gas_group, "cumulative_expansion_work",
                            gas.budget.expansion_work);
            set_attr_double(gas_group, "cumulative_entropy_switch_energy",
                            gas.budget.entropy_switch);
            set_attr_double(gas_group, "cumulative_floor_injected_energy",
                            gas.budget.floor_injection);

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
            write_particle_vec(gas_group, "entropy", gas.entropy);
            write_particle_vec(gas_group, "total_energy", gas.total_energy);
            write_particle_vec(gas_group, "metal_fraction", gas.metal_frac);

            std::vector<double> temp_vec(gas.num_particles);
            for (size_t i = 0; i < gas.num_particles; ++i) {
                temp_vec[i] = Thermodynamics::get_temp_from_internal_energy(
                    gas.u[i], state.scale_factor, config);
            }
            write_particle_vec(gas_group, "temperature", temp_vec);

            // Debugging
            write_particle_vec3d(gas_group, "grad_p", gas.grad_p);
            write_particle_vec(gas_group, "condition_number", gas.cond_num);
            write_particle_vec(gas_group, "n_enc", gas.n_enc_final);

            gas_group.close();

            root_group.close();
            file.close();

            success = true;
            break;
        } catch (H5::Exception& e) {
            std::cerr << "[Warning] HDF5 write failed on attempt " << attempt
                      << " of " << max_retries << ".\n";
            e.printErrorStack();

            if (attempt < max_retries) {
                std::cerr << "Retrying in 5 seconds...\n";
                std::this_thread::sleep_for(std::chrono::seconds(5));
            }
        }
    }

    if (!success) {
        // FATAL: Throw a standard exception to be caught by main()
        throw std::runtime_error(
            "Fatal Error: Could not save HDF5 snapshot after multiple "
            "attempts. Aborting to save compute allocation.");
    }
}
