#include "cooling.h"

#include <H5Cpp.h>

#include <algorithm>  // For std::max, std::min
#include <atomic>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "constants.h"
#include "math_utils.h"  // Needed for QuadraticInterpolator

// Helper function to read an entire HDF5 dataset into a flattened 1D
// std::vector
static std::vector<double> read_hdf5_dataset(const H5::H5File& file,
                                             const std::string& dataset_name) {
    H5::DataSet dataset = file.openDataSet(dataset_name);
    H5::DataSpace dataspace = dataset.getSpace();

    hsize_t num_elements = dataspace.getSimpleExtentNpoints();
    std::vector<double> data(num_elements);

    dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE, H5::DataSpace::ALL,
                 H5::DataSpace::ALL);
    return data;
}

// Helper function to read a 1D HDF5 Attribute attached to a dataset
static std::vector<double> read_hdf5_attribute(const H5::H5File& file,
                                               const std::string& dataset_name,
                                               const std::string& attr_name) {
    H5::DataSet dataset = file.openDataSet(dataset_name);
    H5::Attribute attr = dataset.openAttribute(attr_name);
    H5::DataSpace dataspace = attr.getSpace();

    hsize_t num_elements = dataspace.getSimpleExtentNpoints();
    std::vector<double> data(num_elements);

    attr.read(H5::PredType::NATIVE_DOUBLE, data.data());
    return data;
}

void Cooling::load_cooling_tables(const Config& config) {
    if (config.cooling_table_path.empty()) {
        use_table = false;
        return;
    }

    try {
        H5::H5File file(config.cooling_table_path, H5F_ACC_RDONLY);

        std::string base_path = "/CoolingRates/Metals/";

        // Read the 1D axes from the Attributes of the Cooling dataset
        cooling_data.density =
            read_hdf5_attribute(file, base_path + "Cooling", "Parameter1");
        cooling_data.redshift =
            read_hdf5_attribute(file, base_path + "Cooling", "Parameter2");
        cooling_data.temperature =
            read_hdf5_attribute(file, base_path + "Cooling", "Temperature");

        for (double& t : cooling_data.temperature) {
            t = std::log10(t);
        }

        // Store the dimensions
        cooling_data.Nz = cooling_data.redshift.size();
        cooling_data.Nd = cooling_data.density.size();
        cooling_data.Nt = cooling_data.temperature.size();

        // Read the 3D data tables
        cooling_data.cooling_metal =
            read_hdf5_dataset(file, base_path + "Cooling");
        cooling_data.heating_metal =
            read_hdf5_dataset(file, base_path + "Heating");

        cooling_data.cooling_primordial =
            read_hdf5_dataset(file, "/CoolingRates/Primordial/Cooling");
        cooling_data.heating_primordial =
            read_hdf5_dataset(file, "/CoolingRates/Primordial/Heating");

        // Ensure the arrays size match the product of the axes
        size_t expected_size =
            cooling_data.Nz * cooling_data.Nd * cooling_data.Nt;
        if (cooling_data.cooling_metal.size() != expected_size) {
            throw std::runtime_error(
                "Table size mismatch! Expected " +
                std::to_string(expected_size) + " but got " +
                std::to_string(cooling_data.cooling_metal.size()));
        }

        // Precalculate cgs to code unit conversions
        double X_over_mp = constants::X_MASS_FRAC / constants::M_P_CGS;

        // Heating scales with (1 / a^2)
        double heating_factor = X_over_mp * config.cooling_conversion_factor;

        // Cooling scales with (rho_code / a^5)
        double cooling_factor = config.unit_density_cgs *
                                (X_over_mp * X_over_mp) *
                                config.cooling_conversion_factor;

        size_t stride_d = cooling_data.Nz * cooling_data.Nt;

        for (size_t id = 0; id < cooling_data.Nd; ++id) {
            double nH_bin = std::pow(10.0, cooling_data.density[id]);

            for (size_t iz = 0; iz < cooling_data.Nz; ++iz) {
                for (size_t it = 0; it < cooling_data.Nt; ++it) {
                    size_t idx = id * stride_d + iz * cooling_data.Nt + it;

                    // Bake in the nH correction and the code unit conversions
                    // Grackle/Cloudy tables store heating as (Heating / nH^2)
                    // instead of (/nH). So multiply by nH at load time
                    cooling_data.heating_primordial[idx] *=
                        (nH_bin * heating_factor);
                    cooling_data.heating_metal[idx] *=
                        (nH_bin * heating_factor);

                    cooling_data.cooling_primordial[idx] *= cooling_factor;
                    cooling_data.cooling_metal[idx] *= cooling_factor;
                }
            }
        }

        file.close();
        use_table = true;
        std::cout << "Successfully loaded HM2012 cooling tables: "
                  << config.cooling_table_path << std::endl;
        std::cout << "Grid Dimensions [z, rho, T]: [" << cooling_data.Nz << ", "
                  << cooling_data.Nd << ", " << cooling_data.Nt << "]\n"
                  << std::endl;

    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 Error loading cooling tables: " << e.getDetailMsg()
                  << std::endl;
        use_table = false;
    } catch (const std::exception& e) {
        std::cerr << "Error loading cooling tables: " << e.what() << std::endl;
        use_table = false;
    }
}

Cooling::Cooling(const Config& config) {
    load_cooling_tables(config);

    // Empirical fit for subgrid clumping
    // x = resolution, y = A factor
    Eigen::Vector2d p1(0.325, 0.5);
    Eigen::Vector2d p2(0.5, 1.4);
    Eigen::Vector2d p3(1.04, 10.0);

    QuadraticInterpolator A_curve(p1, p2, p3);
    mean_rho_code =
        config.gas_total_mass /
        (config.domain_size * config.domain_size * config.domain_size);
    double resolution = config.box_size_mpc / config.mesh_size;
    subgrid_clumping_A = config.subgrid_clumping_amplitude < 0
                             ? std::max(0.0, (use_table ? 0.12 : 1.0) *
                                                 A_curve.evaluate(resolution))
                             : config.subgrid_clumping_amplitude;
}

inline int Cooling::get_index(double value,
                              const std::vector<double>& axis) const {
    if (value <= axis.front()) return 0;
    if (value >= axis.back()) return axis.size() - 2;

    auto it = std::upper_bound(axis.begin(), axis.end(), value);
    return std::max(0, static_cast<int>(std::distance(axis.begin(), it)) - 1);
}

double Cooling::trilinear_interpolate(const std::vector<double>& data, double z,
                                      double log_nH, double log_T) const {
    // Get the lower-bound indices
    int iz = get_index(z, cooling_data.redshift);
    int id = get_index(log_nH, cooling_data.density);
    int it = get_index(log_T, cooling_data.temperature);

    // Ensure indices never reach the last element, so index + 1 is always safe.
    iz = std::max(0, std::min(iz, static_cast<int>(cooling_data.Nz - 2)));
    id = std::max(0, std::min(id, static_cast<int>(cooling_data.Nd - 2)));
    it = std::max(0, std::min(it, static_cast<int>(cooling_data.Nt - 2)));

    // Calculate the fractional distances (weights) between 0.0 and 1.0 along
    // each axis
    double wz = (z - cooling_data.redshift[iz]) /
                (cooling_data.redshift[iz + 1] - cooling_data.redshift[iz]);
    double wd = (log_nH - cooling_data.density[id]) /
                (cooling_data.density[id + 1] - cooling_data.density[id]);
    double wt =
        (log_T - cooling_data.temperature[it]) /
        (cooling_data.temperature[it + 1] - cooling_data.temperature[it]);

    // Clamp the weights to [0.0, 1.0] in case the input was outside
    // the table limits
    wz = std::max(0.0, std::min(1.0, wz));
    wd = std::max(0.0, std::min(1.0, wd));
    wt = std::max(0.0, std::min(1.0, wt));

    // Compute the 1D array indices for the 8 corners of the 3D box
    // HDF5 Row-Major Stride for [Nd, Nz, Nt]:
    // index = d * (Nz * Nt) + z * Nt + t
    size_t stride_d = cooling_data.Nz * cooling_data.Nt;
    size_t stride_z = cooling_data.Nt;

    size_t c000 = (id)*stride_d + (iz)*stride_z + (it);
    size_t c001 = (id)*stride_d + (iz)*stride_z + (it + 1);
    size_t c010 = (id)*stride_d + (iz + 1) * stride_z + (it);
    size_t c011 = (id)*stride_d + (iz + 1) * stride_z + (it + 1);
    size_t c100 = (id + 1) * stride_d + (iz)*stride_z + (it);
    size_t c101 = (id + 1) * stride_d + (iz)*stride_z + (it + 1);
    size_t c110 = (id + 1) * stride_d + (iz + 1) * stride_z + (it);
    size_t c111 = (id + 1) * stride_d + (iz + 1) * stride_z + (it + 1);

    // Interpolate along the Temperature axis (t) first
    double c00 = data[c000] * (1.0 - wt) + data[c001] * wt;  // id, iz
    double c01 = data[c010] * (1.0 - wt) + data[c011] * wt;  // id, iz+1
    double c10 = data[c100] * (1.0 - wt) + data[c101] * wt;  // id+1, iz
    double c11 = data[c110] * (1.0 - wt) + data[c111] * wt;  // id+1, iz+1

    // Interpolate along the Redshift axis (z)
    double c0 = c00 * (1.0 - wz) + c01 * wz;  // id
    double c1 = c10 * (1.0 - wz) + c11 * wz;  // id+1

    // Interpolate along the Density axis (d)
    double interpolated_value = c0 * (1.0 - wd) + c1 * wd;

    return interpolated_value;
}

double Cooling::compute_tabulated_cooling(double u_code, double rho_code,
                                          double Z, double a,
                                          double clumping_factor,
                                          const Config& config) const {
    // Values for the axes
    double T_kelvin = get_temp_from_internal_energy(u_code, a, config);
    double log_T = std::log10(T_kelvin);

    double rho_cgs = (rho_code * config.unit_density_cgs) / (a * a * a);
    double nH_cgs = (rho_cgs * constants::X_MASS_FRAC) / constants::M_P_CGS;
    double log_nH = std::log10(nH_cgs);

    double z = 1.0 / a - 1.0;

    // Interpolate the pre-converted rates
    double lambda_metal =
        trilinear_interpolate(cooling_data.cooling_metal, z, log_nH, log_T);
    double gamma_metal =
        trilinear_interpolate(cooling_data.heating_metal, z, log_nH, log_T);

    double lambda_prim = trilinear_interpolate(cooling_data.cooling_primordial,
                                               z, log_nH, log_T);
    double gamma_prim = trilinear_interpolate(cooling_data.heating_primordial,
                                              z, log_nH, log_T);

    // Combine rates based on metallicity
    double metal_ratio = Z / constants::Z_SOLAR;

    double lambda_code = lambda_prim + (lambda_metal * metal_ratio);
    double gamma_code = gamma_prim + (gamma_metal * metal_ratio);

    // Apply dynamic scale factor and density dependencies
    double a2 = a * a;
    double a5 = a2 * a2 * a;

    double heating = gamma_code / a2;
    double cooling = (lambda_code * rho_code * clumping_factor) / a5;

    // Return the final volumetric code energy change
    return heating - cooling;
}

double Cooling::compute_clumping_factor(double rho_code, double a,
                                        const Config& config) const {
    double clumping_factor = 1.0;
    if (config.enable_subgrid_clumping) {
        double overdensity = rho_code / mean_rho_code;
        // Apply clumping only in collapsing regions
        if (overdensity > 1.0) {
            constexpr double B = 0.7;
            double A_dynamic = subgrid_clumping_A * a;
            clumping_factor = 1.0 + A_dynamic * std::pow(overdensity, B);
        }
    }
    return clumping_factor;
}

double Cooling::compute_du_dt(double u_code, double rho_code, double Z,
                              double a, const Config& config) const {
    double clumping_factor = compute_clumping_factor(rho_code, a, config);
    if (use_table) {
        return compute_tabulated_cooling(u_code, rho_code, Z, a,
                                         clumping_factor, config);
    } else {
        double T = get_temp_from_internal_energy(u_code, a, config);

        // Convert code density to CGS
        double rho_cgs = (rho_code * config.unit_density_cgs) / (a * a * a);

        // Number density (n_e = n_i = n)
        double n = rho_cgs / (config.primordial_mu * constants::M_P_CGS);

        // Bremsstrahlung specific energy loss rate (erg / g / s)
        double lambda_cgs = 1.4e-27 * std::sqrt(T) * n * n * clumping_factor;
        double du_dt_cgs = -lambda_cgs / rho_cgs;

        // Convert back to code units
        return (du_dt_cgs * config.cooling_conversion_factor) / (a * a);
    }
}

double Cooling::get_u_rad_floor(double a, const Config& config) const {
    double target_floor_k = 0.0;
    if (use_table) {
        // The table handles real physics. Only enforce the CMB and absolute
        // hydro floor.
        double t_cmb = 2.7315 / a;
        target_floor_k = std::max(config.temp_floor_k, t_cmb);
    } else {
        // The analytical model uses the cutoff to fake the UVB
        target_floor_k = std::max(config.cooling_cutoff_k, config.temp_floor_k);
    }

    double u_rad_floor =
        get_internal_energy_from_temp(target_floor_k, a, config);
    return u_rad_floor;
}

double Cooling::solve_cooling_implicit(double u_old, double rho_code, double Z,
                                       double a, double dt, double u_rad_floor,
                                       const Config& config,
                                       int& iterations_taken) const {
    double du_dt_initial = compute_du_dt(u_old, rho_code, Z, a, config);
    if (du_dt_initial == 0.0) return u_old;

    if (u_old <= u_rad_floor && du_dt_initial < 0.0) {
        iterations_taken = 0;
        return u_rad_floor;
    }

    double u_low, u_high;

    // Negative means cooling, positive means heating
    if (du_dt_initial < 0.0) {
        // Cooling: u_new must be less than u_old
        u_low = u_rad_floor;
        u_high = u_old;

        // Verify the bracket
        double f_low =
            u_low - u_old - dt * compute_du_dt(u_low, rho_code, Z, a, config);

        // If f_low is positive, the root lies below u_rad_floor.
        // This happens when rapid cooling + a large dt push the target
        // below the radiation floor. Just clamp to the floor.
        if (f_low >= 0.0) {
            iterations_taken = 0;
            return u_rad_floor;
        }
    } else {
        // Heating: u_new must be greater than u_old
        u_low = u_old;
        u_high = u_old * 2.0;

        // Find an upper bracket using standard derivative signs
        while ((u_high - u_old -
                dt * compute_du_dt(u_high, rho_code, Z, a, config)) < 0.0) {
            u_high *= 2.0;
            if (u_high > 1e6 * u_old) break;  // Safety escape
        }

        // Explicit fallback if the safety escape triggers
        if ((u_high - u_old -
             dt * compute_du_dt(u_high, rho_code, Z, a, config)) < 0.0) {
            std::cout << "warning: cooling explicit fallback" << std::endl;
            iterations_taken = 1;
            return u_old + (dt * du_dt_initial);
        }
    }

    constexpr double TOLERANCE = 1e-5;
    double u_guess = u_old;
    int iter = 0;

    // Bisection loop
    for (; iter < MAX_ITER; ++iter) {
        u_guess = 0.5 * (u_low + u_high);
        double f_guess = u_guess - u_old -
                         dt * compute_du_dt(u_guess, rho_code, Z, a, config);

        if (std::abs(u_high - u_low) / u_guess < TOLERANCE) {
            break;
        }

        if (f_guess > 0.0) {
            u_high = u_guess;  // Guess was too high
        } else {
            u_low = u_guess;  // Guess was too low
        }
    }

    iterations_taken = iter;
    return std::max(u_guess, u_rad_floor);
}
