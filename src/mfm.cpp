#include "mfm.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>

#include "constants.h"
#include "diagnostics.h"
#include "math_utils.h"

// Enable to use the limiter explained in the Gizmo paper.
// Disable to use the limiter used in the Gizmo code.
// #define THEORETICAL_LIMITER

constexpr double density_floor = 1e-12;
double g_pressure_floor = 0.0;

GasParticleSystem::GasParticleSystem(const Config& config)
    : gas_rho(config.mesh_size), cic_data(config.num_gas_particles) {
    pos_x.reserve(config.num_gas_particles);
    pos_y.reserve(config.num_gas_particles);
    pos_z.reserve(config.num_gas_particles);

    vel_x.reserve(config.num_gas_particles);
    vel_y.reserve(config.num_gas_particles);
    vel_z.reserve(config.num_gas_particles);

    acc_x.reserve(config.num_gas_particles);
    acc_y.reserve(config.num_gas_particles);
    acc_z.reserve(config.num_gas_particles);

    mass.reserve(config.num_gas_particles);

    hydro_acc_x.reserve(config.num_gas_particles);
    hydro_acc_y.reserve(config.num_gas_particles);
    hydro_acc_z.reserve(config.num_gas_particles);

    h.reserve(config.num_gas_particles);
    rho.reserve(config.num_gas_particles);
    pressure.reserve(config.num_gas_particles);
    total_energy.reserve(config.num_gas_particles);
    u.reserve(config.num_gas_particles);
    du_dt.reserve(config.num_gas_particles);
    de_dt.reserve(config.num_gas_particles);
    metal_frac.reserve(config.num_gas_particles);

    entropy.reserve(config.num_gas_particles);
    max_rel_ke.reserve(config.num_gas_particles);
    delta_E_grav.reserve(config.num_gas_particles);

    zeta.reserve(config.num_gas_particles);

    cond_num.reserve(config.num_gas_particles);
    raw_sum_p.reserve(config.num_gas_particles);

    int num_cells = config.mesh_size * config.mesh_size * config.mesh_size;
    pm_cell_list.resize(num_cells, config.num_gas_particles);

    double u_code_K = Cooling::get_internal_energy_from_temp(0.1, 1.0, config);
    pressure_floor = (config.gamma - 1.0) * density_floor * u_code_K;
    g_pressure_floor = pressure_floor;
}

void GasParticleSystem::add_particle(double px, double py, double pz, double vx,
                                     double vy, double vz, double m,
                                     double initial_u, double initial_h,
                                     double z_metal) {
    pos_x.push_back(px);
    pos_y.push_back(py);
    pos_z.push_back(pz);
    vel_x.push_back(vx);
    vel_y.push_back(vy);
    vel_z.push_back(vz);
    acc_x.push_back(0.0);
    acc_y.push_back(0.0);
    acc_z.push_back(0.0);
    mass.push_back(m);

    hydro_acc_x.push_back(0.0);
    hydro_acc_y.push_back(0.0);
    hydro_acc_z.push_back(0.0);

    h.push_back(initial_h);
    rho.push_back(0.0);
    pressure.push_back(0.0);
    u.push_back(initial_u);
    du_dt.push_back(0.0);
    de_dt.push_back(0.0);
    metal_frac.push_back(z_metal);
    B_matrix.push_back(Eigen::Matrix3d::Zero());

    entropy.push_back(0.0);
    max_rel_ke.push_back(0.0);
    delta_E_grav.push_back(0.0);

    zeta.push_back(0.0);

    double initial_ke = 0.5 * (vx * vx + vy * vy + vz * vz);
    total_energy.push_back(initial_u + initial_ke);

    grad_rho.push_back(Eigen::Vector3d::Zero());
    grad_vx.push_back(Eigen::Vector3d::Zero());
    grad_vy.push_back(Eigen::Vector3d::Zero());
    grad_vz.push_back(Eigen::Vector3d::Zero());
    grad_p.push_back(Eigen::Vector3d::Zero());

    cond_num.push_back(0);
    raw_sum_p.push_back(Eigen::Vector3d::Zero());

    num_particles++;
}

void GasParticleSystem::sort_arrays(const std::vector<int>& cell_start,
                                    const std::vector<int>& particle_cell_idx) {
    std::vector<double> new_px(num_particles), new_py(num_particles),
        new_pz(num_particles);
    std::vector<double> new_vx(num_particles), new_vy(num_particles),
        new_vz(num_particles);
    std::vector<double> new_ax(num_particles), new_ay(num_particles),
        new_az(num_particles);
    std::vector<double> new_hax(num_particles), new_hay(num_particles),
        new_haz(num_particles);
    std::vector<double> new_m(num_particles), new_h(num_particles),
        new_rho(num_particles);
    std::vector<double> new_p(num_particles), new_u(num_particles),
        new_dudt(num_particles), new_dedt(num_particles),
        new_total_energy(num_particles);
    std::vector<double> new_metal(num_particles);
    std::vector<CIC_Data> new_cic(num_particles);
    std::vector<Eigen::Matrix3d> new_B(num_particles);
    std::vector<Eigen::Vector3d> new_grad_rho(num_particles),
        new_grad_vx(num_particles);
    std::vector<Eigen::Vector3d> new_grad_vy(num_particles),
        new_grad_vz(num_particles);
    std::vector<Eigen::Vector3d> new_grad_p(num_particles);
    std::vector<double> new_zeta(num_particles);

    std::vector<double> new_entropy(num_particles);
    std::vector<double> new_max_rel_ke(num_particles);
    std::vector<double> new_delta_E_grav(num_particles);

    std::vector<double> new_cond_num(num_particles);
    std::vector<Eigen::Vector3d> new_raw_sum_p(num_particles);

    std::vector<int> write_offset = cell_start;

    for (size_t i = 0; i < num_particles; ++i) {
        int dest = write_offset[particle_cell_idx[i]]++;
        new_px[dest] = pos_x[i];
        new_py[dest] = pos_y[i];
        new_pz[dest] = pos_z[i];
        new_vx[dest] = vel_x[i];
        new_vy[dest] = vel_y[i];
        new_vz[dest] = vel_z[i];
        new_ax[dest] = acc_x[i];
        new_ay[dest] = acc_y[i];
        new_az[dest] = acc_z[i];
        new_m[dest] = mass[i];
        new_hax[dest] = hydro_acc_x[i];
        new_hay[dest] = hydro_acc_y[i];
        new_haz[dest] = hydro_acc_z[i];
        new_h[dest] = h[i];
        new_rho[dest] = rho[i];
        new_p[dest] = pressure[i];
        new_total_energy[dest] = total_energy[i];
        new_u[dest] = u[i];
        new_dudt[dest] = du_dt[i];
        new_dedt[dest] = de_dt[i];
        new_metal[dest] = metal_frac[i];
        new_cic[dest] = cic_data[i];
        new_B[dest] = B_matrix[i];
        new_grad_rho[dest] = grad_rho[i];
        new_grad_vx[dest] = grad_vx[i];
        new_grad_vy[dest] = grad_vy[i];
        new_grad_vz[dest] = grad_vz[i];
        new_grad_p[dest] = grad_p[i];
        new_zeta[dest] = zeta[i];
        new_entropy[dest] = entropy[i];
        new_max_rel_ke[dest] = max_rel_ke[i];
        new_delta_E_grav[dest] = delta_E_grav[i];
        new_cond_num[dest] = cond_num[i];
        new_raw_sum_p[dest] = raw_sum_p[i];
    }

    pos_x = std::move(new_px);
    pos_y = std::move(new_py);
    pos_z = std::move(new_pz);
    vel_x = std::move(new_vx);
    vel_y = std::move(new_vy);
    vel_z = std::move(new_vz);
    acc_x = std::move(new_ax);
    acc_y = std::move(new_ay);
    acc_z = std::move(new_az);
    mass = std::move(new_m);
    hydro_acc_x = std::move(new_hax);
    hydro_acc_y = std::move(new_hay);
    hydro_acc_z = std::move(new_haz);
    h = std::move(new_h);
    rho = std::move(new_rho);
    pressure = std::move(new_p);
    total_energy = std::move(new_total_energy);
    u = std::move(new_u);
    du_dt = std::move(new_dudt);
    de_dt = std::move(new_dedt);
    metal_frac = std::move(new_metal);
    cic_data = std::move(new_cic);
    B_matrix = std::move(new_B);
    grad_rho = std::move(new_grad_rho);
    grad_vx = std::move(new_grad_vx);
    grad_vy = std::move(new_grad_vy);
    grad_vz = std::move(new_grad_vz);
    grad_p = std::move(new_grad_p);
    zeta = std::move(new_zeta);
    entropy = std::move(new_entropy);
    max_rel_ke = std::move(new_max_rel_ke);
    delta_E_grav = std::move(new_delta_E_grav);
    raw_sum_p = std::move(new_raw_sum_p);
    cond_num = std::move(new_cond_num);
}

void GasParticleSystem::build_spatial_hash(double domain_size) {
    if (num_particles == 0) return;

    // Find the maximum smoothing length
    max_h = 0.0;
    for (size_t i = 0; i < num_particles; ++i) {
        if (h[i] > max_h) max_h = h[i];
    }
    if (max_h <= 0.0)
        max_h = domain_size / std::cbrt(num_particles > 0 ? num_particles : 1);

    // Define grid
    hash_cell_size = max_h;
    hash_grid_dim = static_cast<int>(std::floor(domain_size / hash_cell_size));

    if (hash_grid_dim < 1) hash_grid_dim = 1;
    if (hash_grid_dim > 256) hash_grid_dim = 256;

    hash_cell_size = domain_size / hash_grid_dim;
    int num_cells = hash_grid_dim * hash_grid_dim * hash_grid_dim;

    sph_cell_list.resize(num_cells, num_particles);
    std::fill(sph_cell_list.cell_count.begin(), sph_cell_list.cell_count.end(),
              0);
    std::fill(sph_cell_list.cell_start.begin(), sph_cell_list.cell_start.end(),
              0);

    std::vector<int> particle_cell_idx(num_particles);

    // Count particles per cell
    for (size_t i = 0; i < num_particles; ++i) {
        int ix = static_cast<int>(pos_x[i] / hash_cell_size) % hash_grid_dim;
        int iy = static_cast<int>(pos_y[i] / hash_cell_size) % hash_grid_dim;
        int iz = static_cast<int>(pos_z[i] / hash_cell_size) % hash_grid_dim;

        ix = (ix + hash_grid_dim) % hash_grid_dim;
        iy = (iy + hash_grid_dim) % hash_grid_dim;
        iz = (iz + hash_grid_dim) % hash_grid_dim;

        int cell_index =
            iz * hash_grid_dim * hash_grid_dim + iy * hash_grid_dim + ix;
        particle_cell_idx[i] = cell_index;
        sph_cell_list.cell_count[cell_index]++;
    }

    // Prefix sum for offsets
    int current_offset = 0;
    for (int c = 0; c < num_cells; ++c) {
        sph_cell_list.cell_start[c] = current_offset;
        current_offset += sph_cell_list.cell_count[c];
    }

    sort_arrays(sph_cell_list.cell_start, particle_cell_idx);
}

// --------------------------------------------------------------------------------
// Cubic Spline Kernel
// Support radius is 1h. Returns W (value) and dWdh.
// --------------------------------------------------------------------------------
inline void GasParticleSystem::kernel_cubic_spline(double r, double h,
                                                   double& W,
                                                   double& dWdh) const {
    double q = r / h;
    double h3 = h * h * h;
    double norm = 8.0 / (M_PI * h3);  // 3D normalization for 1h support radius
    double dWdr = 0.0;

    if (q < 0.5) {
        // 1 + 6q^2(q - 1) = 1 - 6q^2 + 6q^3
        W = norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
        if (r > 1e-12) {
            // Derivative w.r.t r (using chain rule: dW/dq * dq/dr, where dq/dr
            // = 1/h)
            dWdr = norm * (-12.0 * q + 18.0 * q * q) / h;
        } else {
            dWdr = 0.0;
        }
    } else if (q < 1.0) {
        // 2(1 - q)^3
        double diff = 1.0 - q;
        W = norm * 2.0 * diff * diff * diff;
        // Derivative w.r.t r
        dWdr = norm * -6.0 * diff * diff / h;
    } else {
        // Outside the support radius
        W = 0.0;
        dWdr = 0.0;
    }

    dWdh = -(3.0 / h) * W - q * dWdr;
}

// --------------------------------------------------------------------------------
// Adaptive Gravitational Softening Kernel (GIZMO / Price & Monaghan 2007)
// Returns the modified 1/r^2 force factor and the d(phi)/dh term.
// --------------------------------------------------------------------------------
inline void compute_gravity_kernel_derivatives(double r, double h,
                                               double& dphi_dr,
                                               double& dphi_dh) {
    double q = r / h;
    double q2 = q * q;
    double q3 = q2 * q;
    double q4 = q2 * q2;
    double q5 = q3 * q2;

    double h2 = h * h;
    double h3 = h2 * h;

    if (q < 0.5) {
        // d(phi)/dr : Modified 1/r^2 force
        dphi_dr = (1.0 / h2) *
                  ((10.6666666667) * q - (19.2) * q3 + (10.6666666667) * q4);
        // d(phi)/dh : Softening variation correction
        dphi_dh = (1.0 / h2) * (2.8 - (5.3333333333) * q2 + (11.52) * q4 -
                                (7.1111111111) * q5);
    } else if (q < 1.0) {
        // d(phi)/dr : Modified 1/r^2 force
        dphi_dr = (1.0 / h2) *
                  ((10.6666666667) * q - (19.2) * q3 + (10.6666666667) * q4 -
                   (0.0666666667) / q2 + 0.1 * q5 - (3.2) * q4 +
                   10.6666666667 * q3 - 16.0 * q2 + 11.2 * q - 3.2);
        // d(phi)/dh : Softening variation correction
        dphi_dh = (1.0 / h2) * (3.2 - (10.6666666667) * q + (19.2) * q2 -
                                (11.52) * q3 + (1.0666666667) * q4 - 0.08 / q);
    } else {
        // Pure Newtonian Point Mass outside the softening radius
        dphi_dr = 1.0 / (r * r);
        dphi_dh = 0.0;
    }
}

// --------------------------------------------------------------------------------
// Density and Smoothing Length Iteration
// --------------------------------------------------------------------------------
void GasParticleSystem::compute_density_and_h(const Config& config,
                                              const ParticleSystem& dm) {
    if (num_particles == 0) return;

    build_spatial_hash(config.domain_size);

    double domain_size = config.domain_size;
    double target_N = config.mfm_target_neighbors;
    double tol = config.mfm_neighbor_tolerance;
    int max_iter = config.mfm_max_iterations;
    const double mean_spacing =
        domain_size / std::cbrt(num_particles > 0 ? num_particles : 1);
    const double min_h = 0.05 * mean_spacing;
    //const double max_h = 4.0 * mean_spacing;
    const double max_h = 0.3 * domain_size;

    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];

        double h_low = 0.0;
        double h_high = std::numeric_limits<double>::infinity();
        double h_guess = h[i];  // Start with current h

        int iter = 0;
        double current_n = 0.0;      // Track number density (n_i)
        double current_dn_dh = 0.0;  // Track derivative (dn_i/dh_i)

        while (iter < max_iter) {
            current_n = 0.0;
            current_dn_dh = 0.0;

            // Determine how many hash cells to search based on current h_guess
            double search_radius = h_guess;
            int search_cells =
                static_cast<int>(std::ceil(search_radius / hash_cell_size));
            // Prevent wrapping around and double-counting the domain
            search_cells = std::min(search_cells, hash_grid_dim / 2);

            // Get current particle cell
            int ix = static_cast<int>(p1_x / hash_cell_size) % hash_grid_dim;
            int iy = static_cast<int>(p1_y / hash_cell_size) % hash_grid_dim;
            int iz = static_cast<int>(p1_z / hash_cell_size) % hash_grid_dim;
            ix = (ix + hash_grid_dim) % hash_grid_dim;
            iy = (iy + hash_grid_dim) % hash_grid_dim;
            iz = (iz + hash_grid_dim) % hash_grid_dim;

            // Gather neighbors
            for (int dx = -search_cells; dx <= search_cells; ++dx) {
                for (int dy = -search_cells; dy <= search_cells; ++dy) {
                    for (int dz = -search_cells; dz <= search_cells; ++dz) {
                        // Get cell index
                        int neighbor_ix =
                            (((ix + dx) % hash_grid_dim) + hash_grid_dim) %
                            hash_grid_dim;
                        int neighbor_iy =
                            (((iy + dy) % hash_grid_dim) + hash_grid_dim) %
                            hash_grid_dim;
                        int neighbor_iz =
                            (((iz + dz) % hash_grid_dim) + hash_grid_dim) %
                            hash_grid_dim;

                        int cell_idx =
                            neighbor_iz * hash_grid_dim * hash_grid_dim +
                            neighbor_iy * hash_grid_dim + neighbor_ix;

                        // Loop cell particles
                        int start = sph_cell_list.cell_start[cell_idx];
                        int end = start + sph_cell_list.cell_count[cell_idx];

                        for (int j = start; j < end; ++j) {
                            double dist_x = periodic_displacement(
                                pos_x[j] - p1_x, domain_size);
                            double dist_y = periodic_displacement(
                                pos_y[j] - p1_y, domain_size);
                            double dist_z = periodic_displacement(
                                pos_z[j] - p1_z, domain_size);

                            double r2 = dist_x * dist_x + dist_y * dist_y +
                                        dist_z * dist_z;
                            if (r2 < search_radius * search_radius) {
                                double r = std::sqrt(r2);
                                double W, dWdh;
                                kernel_cubic_spline(r, h_guess, W, dWdh);

                                // MFM relies on particle positions,
                                // summing only the kernel weights
                                current_n += W;
                                current_dn_dh += dWdh;
                            }
                        }
                    }
                }
            }

            // Calculate effective number of neighbors enclosed in the kernel
            // N_enc = (4/3) * pi * h^3 * n_i
            double h3 = h_guess * h_guess * h_guess;
            double N_enc = (4.0 / 3.0) * M_PI * h3 * current_n;

            // Convergence Check
            if (std::abs(N_enc - target_N) < tol) {
                break;
            }

            // Update safety boundaries
            if (N_enc > target_N) {
                h_high = h_guess;
            } else {
                h_low = h_guess;
            }

            // Calculate the derivative of N_enc with respect to h
            // dN_enc/dh = (4/3)*pi * [3*h^2 * n_i + h^3 * dn_i/dh]
            double dN_enc_dh =
                (4.0 / 3.0) * M_PI *
                (3.0 * h_guess * h_guess * current_n + h3 * current_dn_dh);

            // Newton-Raphson Step
            double h_new = h_guess;
            if (dN_enc_dh > 0.0) {
                h_new = h_guess - (N_enc - target_N) / dN_enc_dh;
            }

            // Safe Newton-Raphson: Fall back to bisection if the Newton step
            // overshoots our known boundaries or if the derivative is
            // flat/invalid.
            if (h_new <= h_low || h_new >= h_high || dN_enc_dh <= 0.0) {
                h_guess = std::isinf(h_high) ? (1.26 * h_guess)
                                             : 0.5 * (h_low + h_high);
            } else {
                h_guess = h_new;
            }

            // Clamp maximum smoothing length
            if (h_guess >= max_h) {
                h_guess = max_h;
                if (h_low >= max_h) break;
            }

            if (h_guess < min_h) {
                h_guess = min_h;
                break;
            }

            iter++;
            if (iter == max_iter) {
                std::cout << "Warning: compute_h did not converge" << std::endl;
            }
        }

        // Commit final state for this particle
        h[i] = h_guess;

        // Compute effective volume and density for MFM
        double effective_volume = 1.0 / current_n;
        rho[i] = mass[i] / effective_volume;

        // Adaptive gravity softening corrections (GIZMO App. H2)
        // Equation H10: Omega_a = 1 + (h_a / (n_a * nu)) * (dn_a / dh_a)
        double Omega_i = 1.0 + (h_guess / (current_n * 3.0)) * current_dn_dh;

        // Final gather loop to calculate zeta_i (Equation H9)
        double zeta_sum = 0.0;

        int search_cells =
            static_cast<int>(std::ceil(h_guess / hash_cell_size));
        search_cells = std::min(search_cells, hash_grid_dim / 2);

        int ix = static_cast<int>(p1_x / hash_cell_size) % hash_grid_dim;
        int iy = static_cast<int>(p1_y / hash_cell_size) % hash_grid_dim;
        int iz = static_cast<int>(p1_z / hash_cell_size) % hash_grid_dim;
        ix = (ix + hash_grid_dim) % hash_grid_dim;
        iy = (iy + hash_grid_dim) % hash_grid_dim;
        iz = (iz + hash_grid_dim) % hash_grid_dim;

        for (int dx_c = -search_cells; dx_c <= search_cells; ++dx_c) {
            for (int dy_c = -search_cells; dy_c <= search_cells; ++dy_c) {
                for (int dz_c = -search_cells; dz_c <= search_cells; ++dz_c) {
                    int n_ix = (((ix + dx_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iy = (((iy + dy_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iz = (((iz + dz_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;

                    int cell_idx = n_iz * hash_grid_dim * hash_grid_dim +
                                   n_iy * hash_grid_dim + n_ix;
                    int start = sph_cell_list.cell_start[cell_idx];
                    int end = start + sph_cell_list.cell_count[cell_idx];

                    for (int j = start; j < end; ++j) {
                        double dx =
                            periodic_displacement(pos_x[j] - p1_x, domain_size);
                        double dy =
                            periodic_displacement(pos_y[j] - p1_y, domain_size);
                        double dz =
                            periodic_displacement(pos_z[j] - p1_z, domain_size);
                        double r2 = dx * dx + dy * dy + dz * dz;

                        if (r2 < h_guess * h_guess && r2 > 1e-24) {
                            double r = std::sqrt(r2);
                            double dphi_dr, dphi_dh;
                            compute_gravity_kernel_derivatives(
                                r, h_guess, dphi_dr, dphi_dh);

                            if (use_pm) {
                                double r_scaled = r / (2.0 * r_s);
                                double taper =
                                    std::erfc(r_scaled) +
                                    (r / (std::sqrt(M_PI) * r_s)) *
                                        std::exp(-r_scaled * r_scaled);
                                dphi_dh *= taper;
                            }

                            // Summing m_b * d(phi)/dh
                            zeta_sum += mass[j] * dphi_dh;
                        }
                    }
                }
            }
        }

        // Accumulate zeta from Dark Matter interactions
        // DM uses the PM mesh spacing for neighbor finding
        int dm_dx_start = use_pm ? -static_cast<int>(std::ceil(
                                       config.cutoff_radius / config.cell_size))
                                 : 0;
        int dm_dx_end = use_pm ? static_cast<int>(std::ceil(
                                     config.cutoff_radius / config.cell_size))
                               : config.mesh_size - 1;

        int dm_ix = static_cast<int>(p1_x / config.cell_size);
        int dm_iy = static_cast<int>(p1_y / config.cell_size);
        int dm_iz = static_cast<int>(p1_z / config.cell_size);
        dm_ix = (((dm_ix % config.mesh_size) + config.mesh_size) %
                 config.mesh_size);
        dm_iy = (((dm_iy % config.mesh_size) + config.mesh_size) %
                 config.mesh_size);
        dm_iz = (((dm_iz % config.mesh_size) + config.mesh_size) %
                 config.mesh_size);

        for (int dx_cell = dm_dx_start; dx_cell <= dm_dx_end; ++dx_cell) {
            for (int dy_cell = dm_dx_start; dy_cell <= dm_dx_end; ++dy_cell) {
                for (int dz_cell = dm_dx_start; dz_cell <= dm_dx_end;
                     ++dz_cell) {
                    int n_ix = use_pm
                                   ? (((dm_ix + dx_cell) % config.mesh_size) +
                                      config.mesh_size) %
                                         config.mesh_size
                                   : dx_cell;
                    int n_iy = use_pm
                                   ? (((dm_iy + dy_cell) % config.mesh_size) +
                                      config.mesh_size) %
                                         config.mesh_size
                                   : dy_cell;
                    int n_iz = use_pm
                                   ? (((dm_iz + dz_cell) % config.mesh_size) +
                                      config.mesh_size) %
                                         config.mesh_size
                                   : dz_cell;
                    int cell_idx = n_iz * config.mesh_size * config.mesh_size +
                                   n_iy * config.mesh_size + n_ix;

                    int start = dm.cell_list.cell_start[cell_idx];
                    int end = start + dm.cell_list.cell_count[cell_idx];

                    for (int j = start; j < end; ++j) {
                        double dx = periodic_displacement(dm.pos_x[j] - p1_x,
                                                          domain_size);
                        double dy = periodic_displacement(dm.pos_y[j] - p1_y,
                                                          domain_size);
                        double dz = periodic_displacement(dm.pos_z[j] - p1_z,
                                                          domain_size);
                        double r2 = dx * dx + dy * dy + dz * dz;

                        // DM particles only affect gas zeta if they fall inside
                        // the gas smoothing length
                        if (r2 < h_guess * h_guess && r2 > 1e-24) {
                            if (use_pm && r2 > config.cutoff_radius_squared)
                                continue;

                            double r = std::sqrt(r2);
                            double dphi_dr, dphi_dh;
                            compute_gravity_kernel_derivatives(
                                r, h_guess, dphi_dr, dphi_dh);

                            if (use_pm) {
                                double r_scaled = r / (2.0 * r_s);
                                double taper =
                                    std::erfc(r_scaled) +
                                    (r / (std::sqrt(M_PI) * r_s)) *
                                        std::exp(-r_scaled * r_scaled);
                                dphi_dh *= taper;
                            }

                            zeta_sum += dm.mass[j] * dphi_dh;
                        }
                    }
                }
            }
        }

        // Equation H9: zeta_a = m_a * (h_a / (n_a * nu)) * (1 / Omega_a) *
        // SUM(m_b * dphi/dh)
        zeta[i] = (h_guess / (current_n * 3.0)) * (1.0 / Omega_i) * zeta_sum;
    }

//#define MFM_VOLUME_LIMITER
#ifdef MFM_VOLUME_LIMITER
    // Prevents vacuum particles from artificially borrowing density from shocks
    std::vector<double> limited_rho = rho;
    constexpr double MAX_VOL_RATIO =
        8.0;  // Tunable: Expected between 2.0 and 8.0

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
        double my_vol = mass[i] / rho[i];
        double max_neighbor_vol = 0.0;

        int search_cells = static_cast<int>(std::ceil(h[i] / hash_cell_size));
        search_cells = std::min(search_cells, hash_grid_dim / 2);

        int ix = static_cast<int>(p1_x / hash_cell_size) % hash_grid_dim;
        int iy = static_cast<int>(p1_y / hash_cell_size) % hash_grid_dim;
        int iz = static_cast<int>(p1_z / hash_cell_size) % hash_grid_dim;
        ix = (ix + hash_grid_dim) % hash_grid_dim;
        iy = (iy + hash_grid_dim) % hash_grid_dim;
        iz = (iz + hash_grid_dim) % hash_grid_dim;

        for (int dx_c = -search_cells; dx_c <= search_cells; ++dx_c) {
            for (int dy_c = -search_cells; dy_c <= search_cells; ++dy_c) {
                for (int dz_c = -search_cells; dz_c <= search_cells; ++dz_c) {
                    int n_ix = (((ix + dx_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iy = (((iy + dy_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iz = (((iz + dz_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;

                    int cell_idx = n_iz * hash_grid_dim * hash_grid_dim +
                                   n_iy * hash_grid_dim + n_ix;

                    int start = sph_cell_list.cell_start[cell_idx];
                    int end = start + sph_cell_list.cell_count[cell_idx];

                    for (int j = start; j < end; ++j) {
                        double dx =
                            periodic_displacement(pos_x[j] - p1_x, domain_size);
                        double dy =
                            periodic_displacement(pos_y[j] - p1_y, domain_size);
                        double dz =
                            periodic_displacement(pos_z[j] - p1_z, domain_size);
                        double r2 = dx * dx + dy * dy + dz * dz;

                        if (r2 < h[i] * h[i] && r2 > 1e-24) {
                            double neighbor_vol = mass[j] / rho[j];
                            max_neighbor_vol =
                                std::max(max_neighbor_vol, neighbor_vol);
                        }
                    }
                }
            }
        }

        // Enforce the volume limit
        if (max_neighbor_vol > 0.0 &&
            my_vol > MAX_VOL_RATIO * max_neighbor_vol) {
            double clamped_vol = MAX_VOL_RATIO * max_neighbor_vol;
            limited_rho[i] = mass[i] / clamped_vol;
        }
    }

    rho = std::move(limited_rho);

#endif

    build_spatial_hash(config.domain_size);
}

void GasParticleSystem::bin_and_assign_mass(const Config& config) {
    gas_rho.setZero();
    cic_data.assign(num_particles, {});

    int N = config.mesh_size;
    int num_cells = N * N * N;

    pm_cell_list.resize(num_cells, num_particles);
    std::fill(pm_cell_list.cell_count.begin(), pm_cell_list.cell_count.end(),
              0);
    std::fill(pm_cell_list.cell_start.begin(), pm_cell_list.cell_start.end(),
              0);

    std::vector<int> particle_cell_idx(num_particles);

    // Calculate cells & densities
    for (size_t i = 0; i < num_particles; ++i) {
        double px = pos_x[i], py = pos_y[i], pz = pos_z[i];

        // Must use physical coordinates to align with the PP Gravity solver
        int hash_ix = static_cast<int>(px / config.cell_size) % N;
        int hash_iy = static_cast<int>(py / config.cell_size) % N;
        int hash_iz = static_cast<int>(pz / config.cell_size) % N;

        // Protect against negative bounds
        hash_ix = (hash_ix + N) % N;
        hash_iy = (hash_iy + N) % N;
        hash_iz = (hash_iz + N) % N;
        int cell_index = hash_iz * N * N + hash_iy * N + hash_ix;

        // Cell centered PM grid nodes
        double shifted_x = px - 0.5 * config.cell_size;
        double shifted_y = py - 0.5 * config.cell_size;
        double shifted_z = pz - 0.5 * config.cell_size;

        // Ensure periodic wrap-around bounds
        shifted_x = fmod(shifted_x + config.domain_size, config.domain_size);
        shifted_y = fmod(shifted_y + config.domain_size, config.domain_size);
        shifted_z = fmod(shifted_z + config.domain_size, config.domain_size);

        int ix = static_cast<int>(shifted_x / config.cell_size);
        int iy = static_cast<int>(shifted_y / config.cell_size);
        int iz = static_cast<int>(shifted_z / config.cell_size);

        double frac_x = (shifted_x / config.cell_size) - ix;
        double frac_y = (shifted_y / config.cell_size) - iy;
        double frac_z = (shifted_z / config.cell_size) - iz;

        double w000 = (1 - frac_x) * (1 - frac_y) * (1 - frac_z);
        double w100 = frac_x * (1 - frac_y) * (1 - frac_z);
        double w010 = (1 - frac_x) * frac_y * (1 - frac_z);
        double w110 = frac_x * frac_y * (1 - frac_z);
        double w001 = (1 - frac_x) * (1 - frac_y) * frac_z;
        double w101 = frac_x * (1 - frac_y) * frac_z;
        double w011 = (1 - frac_x) * frac_y * frac_z;
        double w111 = frac_x * frac_y * frac_z;

        cic_data[i] = {ix,   iy,   iz,   w000, w100, w010,
                       w110, w001, w101, w011, w111};

        int ix0 = (ix + N) % N, ix1 = (ix + 1 + N) % N;
        int iy0 = (iy + N) % N, iy1 = (iy + 1 + N) % N;
        int iz0 = (iz + N) % N, iz1 = (iz + 1 + N) % N;

        double m = mass[i];
        gas_rho(ix0, iy0, iz0) += m * w000;
        gas_rho(ix1, iy0, iz0) += m * w100;
        gas_rho(ix0, iy1, iz0) += m * w010;
        gas_rho(ix1, iy1, iz0) += m * w110;
        gas_rho(ix0, iy0, iz1) += m * w001;
        gas_rho(ix1, iy0, iz1) += m * w101;
        gas_rho(ix0, iy1, iz1) += m * w011;
        gas_rho(ix1, iy1, iz1) += m * w111;

        particle_cell_idx[i] = cell_index;
        pm_cell_list.cell_count[cell_index]++;
    }

    gas_rho.data /= config.cell_volume;

    // Prefix sum
    int current_offset = 0;
    for (int c = 0; c < num_cells; ++c) {
        pm_cell_list.cell_start[c] = current_offset;
        current_offset += pm_cell_list.cell_count[c];
    }

    sort_arrays(pm_cell_list.cell_start, particle_cell_idx);
}

void GasParticleSystem::interpolate_cic_forces(const Grid3D& ax_grid,
                                               const Grid3D& ay_grid,
                                               const Grid3D& az_grid,
                                               const Config& config) {
    const int N = config.mesh_size;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        const auto& cd = cic_data[i];

        int ix0 = (cd.ix + N) % N, ix1 = (cd.ix + 1 + N) % N;
        int iy0 = (cd.iy + N) % N, iy1 = (cd.iy + 1 + N) % N;
        int iz0 = (cd.iz + N) % N, iz1 = (cd.iz + 1 + N) % N;

        auto interp = [&](const Grid3D& grid) {
            return grid(ix0, iy0, iz0) * cd.w000 +
                   grid(ix1, iy0, iz0) * cd.w100 +
                   grid(ix0, iy1, iz0) * cd.w010 +
                   grid(ix1, iy1, iz0) * cd.w110 +
                   grid(ix0, iy0, iz1) * cd.w001 +
                   grid(ix1, iy0, iz1) * cd.w101 +
                   grid(ix0, iy1, iz1) * cd.w011 +
                   grid(ix1, iy1, iz1) * cd.w111;
        };

        acc_x[i] = interp(ax_grid);
        acc_y[i] = interp(ay_grid);
        acc_z[i] = interp(az_grid);
    }
}

void GasParticleSystem::apply_cooling(double dt, double a, const Config& config,
                                      Cooling& cooling) {
    if (!config.enable_cooling) return;

    double u_rad_floor = cooling.get_u_rad_floor(a, config);
    double total_radiated = 0.0;
    double total_photoheated = 0.0;
    size_t non_converged_count = 0;
    size_t total_cycles = 0;

    double gamma_minus_1 = config.gamma - 1.0;

#pragma omp parallel for schedule(static)                                 \
    reduction(+ : total_radiated, total_photoheated, non_converged_count, \
                  total_cycles)
    for (size_t i = 0; i < num_particles; ++i) {
        double local_rho = rho[i];
        if (local_rho > 1e-12) {  // Skip vacuum particles
            // In MFM, we directly track the metal mass fraction
            double local_Z_frac = metal_frac[i];
            double u_current = u[i];
            double u_initial = u_current;
            double t_evolved = 0.0;
            int cell_non_converged = 0;

            // Local Particle Subcycling
            while (t_evolved < dt) {
                double du_dt = cooling.compute_du_dt(u_current, local_rho,
                                                     local_Z_frac, a, config);

                double dt_cell = (std::abs(du_dt) > 0.0)
                                     ? 0.1 * (u_current / std::abs(du_dt))
                                     : dt;

                dt_cell = std::min(dt_cell, dt - t_evolved);

                int iters = 0;
                u_current = cooling.solve_cooling_implicit(
                    u_current, local_rho, local_Z_frac, a, dt_cell, u_rad_floor,
                    config, iters);

                if (iters >= Cooling::MAX_ITER) {
                    cell_non_converged++;
                }

                t_evolved += dt_cell;
                total_cycles++;
            }

            if (cell_non_converged > 0) {
                non_converged_count++;
            }

            double delta_u = u_current - u_initial;
            if (std::abs(delta_u) > 0.0) {
                // Update the internal energy
                u[i] = u_current;
                total_energy[i] += delta_u;

                // Sync the entropy so the dual-energy switch
                // doesn't override the update
                entropy[i] =
                    gamma_minus_1 * u[i] / std::pow(rho[i], gamma_minus_1);

                // Track total energy change using the particle mass
                double delta_E = delta_u * mass[i];
                if (delta_u < 0.0) {
                    total_radiated -= delta_E;
                } else {
                    total_photoheated += delta_E;
                }
            }
        }
    }

    this->cooling_total_cycles =
        num_particles > 0 ? (total_cycles / num_particles) : 0;
    this->cooling_failed_cells = non_converged_count;
    this->accumulated_radiated_energy += total_radiated;
    this->accumulated_photoheating_energy += total_photoheated;

    // Resynchronize pressure and dual-energy arrays using the cooled ie
    update_primitive_variables(config, a);
}

// Periodic displacement
inline double mfm_periodic_displacement(double dx, double domain_size) {
    if (dx > 0.5 * domain_size) return dx - domain_size;
    if (dx < -0.5 * domain_size) return dx + domain_size;
    return dx;
}

double GasParticleSystem::get_cfl_timestep(const Config& config) const {
    if (config.hydro_method != HydroMethod::MFM || num_particles == 0) {
        return std::numeric_limits<double>::infinity();
    }

    double min_dt = std::numeric_limits<double>::infinity();
    double gamma = config.gamma;
    double domain_size = config.domain_size;

    // Use dynamic scheduling because the number of neighbors varies per cell
#pragma omp parallel for reduction(min : min_dt) schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
        double v1_x = vel_x[i], v1_y = vel_y[i], v1_z = vel_z[i];

        // Sound speed of particle i
        double c_i =
            (rho[i] > 1e-12) ? std::sqrt(gamma * pressure[i] / rho[i]) : 0.0;
        double h_i = h[i];

        double v_sig_max = 0.0;

        // Hash grid coordinates for particle i
        int ix = static_cast<int>(p1_x / hash_cell_size) % hash_grid_dim;
        int iy = static_cast<int>(p1_y / hash_cell_size) % hash_grid_dim;
        int iz = static_cast<int>(p1_z / hash_cell_size) % hash_grid_dim;
        ix = (ix + hash_grid_dim) % hash_grid_dim;
        iy = (iy + hash_grid_dim) % hash_grid_dim;
        iz = (iz + hash_grid_dim) % hash_grid_dim;

        int search_cells = 1;  // Kernel support spans 1 hash cell in radius

        // Neighbor Search
        for (int dx_c = -search_cells; dx_c <= search_cells; ++dx_c) {
            for (int dy_c = -search_cells; dy_c <= search_cells; ++dy_c) {
                for (int dz_c = -search_cells; dz_c <= search_cells; ++dz_c) {
                    int n_ix = (((ix + dx_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iy = (((iy + dy_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iz = (((iz + dz_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;

                    int cell_idx = n_iz * hash_grid_dim * hash_grid_dim +
                                   n_iy * hash_grid_dim + n_ix;

                    int start = sph_cell_list.cell_start[cell_idx];
                    int end = start + sph_cell_list.cell_count[cell_idx];

                    for (int j = start; j < end; ++j) {
                        if (i == j) continue;  // Skip self

                        double dx = mfm_periodic_displacement(pos_x[j] - p1_x,
                                                              domain_size);
                        double dy = mfm_periodic_displacement(pos_y[j] - p1_y,
                                                              domain_size);
                        double dz = mfm_periodic_displacement(pos_z[j] - p1_z,
                                                              domain_size);
                        double r2 = dx * dx + dy * dy + dz * dz;

                        // Particles interact if they overlap in either support
                        // radius
                        if (r2 < h_i * h_i || r2 < h[j] * h[j]) {
                            if (r2 > 1e-24) {
                                double r = std::sqrt(r2);
                                double c_j =
                                    (rho[j] > 1e-12)
                                        ? std::sqrt(gamma * pressure[j] /
                                                    rho[j])
                                        : 0.0;

                                // Relative velocity: (v_j - v_i)
                                double dvx = vel_x[j] - v1_x;
                                double dvy = vel_y[j] - v1_y;
                                double dvz = vel_z[j] - v1_z;

                                // Dot product: (v_j - v_i) dot (x_j - x_i)
                                // Note: Mathematically identical to (v_i - v_j)
                                // dot (x_i - x_j) from Eq. 25
                                double dv_dot_dx =
                                    dvx * dx + dvy * dy + dvz * dz;

                                // Equation 25: c_{s,i} + c_{s,j}
                                double v_sig_ij = c_i + c_j;

                                // Equation 25: - MIN(0, (dv dot dx) / r)
                                if (dv_dot_dx < 0.0) {
                                    v_sig_ij -= dv_dot_dx / r;
                                }

                                // Equation 25: MAX_j [...]
                                if (v_sig_ij > v_sig_max) {
                                    v_sig_max = v_sig_ij;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Equation 24: dt_i = 2 * C_CFL * (h_i / v_sig_max)
        if (v_sig_max > 1e-9) {
            double dt_i = 2.0 * h_i / v_sig_max;
            if (dt_i < min_dt) {
                min_dt = dt_i;
            }
        }
    }

    return min_dt * config.hydro_courant_factor;
}

double GasParticleSystem::get_cooling_timestep(double a, const Config& config,
                                               Cooling& cooling) const {
    if (!config.enable_cooling || config.hydro_method != HydroMethod::MFM ||
        num_particles == 0) {
        return std::numeric_limits<double>::infinity();
    }

    double min_dt_cool = std::numeric_limits<double>::infinity();
    double u_rad_floor = cooling.get_u_rad_floor(a, config);

#pragma omp parallel for reduction(min : min_dt_cool)
    for (size_t i = 0; i < num_particles; ++i) {
        if (rho[i] > 1e-12) {
            if (u[i] <= u_rad_floor) continue;

            double du_dt_val =
                cooling.compute_du_dt(u[i], rho[i], metal_frac[i], a, config);

            if (std::abs(du_dt_val) > 0.0) {
                // Restrict timestep so internal energy changes at most 10%
                double dt_cool = 0.1 * (u[i] / std::abs(du_dt_val));
                if (dt_cool < min_dt_cool) {
                    min_dt_cool = dt_cool;
                }
            }
        }
    }
    return min_dt_cool;
}

double GasParticleSystem::get_gravity_timestep(const Config& config) const {
    if (num_particles == 0) return std::numeric_limits<double>::infinity();
    double epsilon = std::sqrt(config.softening_squared);
    double a_max = std::sqrt(max_accel_sq);
    double dt_grav = std::sqrt(epsilon / a_max);
    return dt_grav * config.gravity_accuracy_eta;
}

// --------------------------------------------------------------------------------
// Adaptive Gravitational Softening Kernel (GIZMO App. H2 / Price & Monaghan
// 2007) Returns the exact d(phi)/dr potential derivative and d(W)/dr smoothing
// derivative.
// --------------------------------------------------------------------------------
inline void compute_adaptive_gravity_terms(double r, double h, double& dphi_dr,
                                           double& dW_dr) {
    if (r >= h) {
        dphi_dr = 1.0 / (r * r);  // Pure Newtonian point-mass outside h
        dW_dr = 0.0;              // W(r) is zero outside h
        return;
    }

    double q = r / h;
    double q2 = q * q;
    double q3 = q2 * q;
    double q4 = q3 * q;
    double q5 = q4 * q;

    double h2 = h * h;
    double h3 = h2 * h;

    if (q < 0.5) {
        dphi_dr =
            (1.0 / h2) * (10.6666666667 * q - 19.2 * q3 + 10.6666666667 * q4);
    } else {
        dphi_dr =
            (1.0 / h2) * (10.6666666667 * q - 19.2 * q3 + 10.6666666667 * q4 -
                          0.0666666667 / q2 + 0.1 * q5 - 3.2 * q4 +
                          10.6666666667 * q3 - 16.0 * q2 + 11.2 * q - 3.2);
    }

    double norm = 8.0 / (M_PI * h3);
    if (q < 0.5) {
        dW_dr = norm * (-12.0 * q + 18.0 * q2) / h;
    } else {
        double diff = 1.0 - q;
        dW_dr = norm * (-6.0 * diff * diff) / h;
    }
}

void GasParticleSystem::compute_and_add_pp_forces(const Config& config,
                                                  Diagnostics& diag) {
    const int N = config.mesh_size;
    const double cell_size = config.cell_size;
    const double domain_size = config.domain_size;
    const double G = config.G;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * cell_size;
    const bool use_pm = config.use_PM;

    const size_t n_gas = num_particles;
    int dx_start =
        use_pm ? -static_cast<int>(std::ceil(config.cutoff_radius / cell_size))
               : 0;
    int dx_end =
        use_pm ? static_cast<int>(std::ceil(config.cutoff_radius / cell_size))
               : N - 1;
    const int num_cells = N * N * N;

#ifdef USE_GPU
    if (config.enable_GPU) {
        // ========================================================================
        // GPU IMPLEMENTATION (Sorted Array)
        // ========================================================================
        double* d_px = pos_x.data();
        double* d_py = pos_y.data();
        double* d_pz = pos_z.data();
        double* d_m = mass.data();
        double* d_h = h.data();
        double* d_zeta = zeta.data();
        double* d_ax = acc_x.data();
        double* d_ay = acc_y.data();
        double* d_az = acc_z.data();
        int* d_cell_start = pm_cell_list.cell_start.data();
        int* d_cell_count = pm_cell_list.cell_count.data();

        auto start_transfer = std::chrono::high_resolution_clock::now();
#pragma omp target enter data map(                                    \
        to : d_px[0 : n_gas], d_py[0 : n_gas], d_pz[0 : n_gas],       \
            d_m[0 : n_gas], d_h[0 : n_gas], d_zeta[0 : n_gas],        \
            d_cell_start[0 : num_cells], d_cell_count[0 : num_cells], \
            d_ax[0 : n_gas], d_ay[0 : n_gas], d_az[0 : n_gas])
        auto end_transfer = std::chrono::high_resolution_clock::now();

        auto start_compute = std::chrono::high_resolution_clock::now();
#pragma omp target teams distribute parallel for
        for (size_t i = 0; i < n_gas; ++i) {
            double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
            double m_i = d_m[i];
            double h_i = d_h[i];
            double zeta_i = d_zeta[i];

            int ix = static_cast<int>(p1_x / cell_size);
            int iy = static_cast<int>(p1_y / cell_size);
            int iz = static_cast<int>(p1_z / cell_size);
            ix = (((ix % N) + N) % N);
            iy = (((iy % N) + N) % N);
            iz = (((iz % N) + N) % N);

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

            for (int dx_cell = dx_start; dx_cell <= dx_end; ++dx_cell) {
                for (int dy_cell = dx_start; dy_cell <= dx_end; ++dy_cell) {
                    for (int dz_cell = dx_start; dz_cell <= dx_end; ++dz_cell) {
                        int neighbor_ix =
                            use_pm ? (((ix + dx_cell) % N) + N) % N : dx_cell;
                        int neighbor_iy =
                            use_pm ? (((iy + dy_cell) % N) + N) % N : dy_cell;
                        int neighbor_iz =
                            use_pm ? (((iz + dz_cell) % N) + N) % N : dz_cell;
                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        int start_idx = d_cell_start[cell_idx];
                        int end_idx = start_idx + d_cell_count[cell_idx];

                        for (int j = start_idx; j < end_idx; ++j) {
                            if (i == j) continue;

                            double dx = d_px[j] - p1_x;
                            if (dx > 0.5 * domain_size)
                                dx -= domain_size;
                            else if (dx < -0.5 * domain_size)
                                dx += domain_size;

                            double dy = d_py[j] - p1_y;
                            if (dy > 0.5 * domain_size)
                                dy -= domain_size;
                            else if (dy < -0.5 * domain_size)
                                dy += domain_size;

                            double dz = d_pz[j] - p1_z;
                            if (dz > 0.5 * domain_size)
                                dz -= domain_size;
                            else if (dz < -0.5 * domain_size)
                                dz += domain_size;

                            double dist_sq = dx * dx + dy * dy + dz * dz;

                            if (use_pm && dist_sq > cutoff_sq) continue;
                            if (dist_sq < 1e-24) continue;

                            double r = std::sqrt(dist_sq);
                            double m_j = d_m[j];
                            double h_j = d_h[j];
                            double zeta_j = d_zeta[j];

                            // Inline Kernel Derivs for Particle I
                            double dphi_dr_i, dW_dr_i;
                            if (r >= h_i) {
                                dphi_dr_i = 1.0 / (r * r);
                                dW_dr_i = 0.0;
                            } else {
                                double q = r / h_i;
                                double q2 = q * q;
                                double q3 = q2 * q;
                                double q4 = q3 * q;
                                double q5 = q4 * q;
                                double h2 = h_i * h_i;
                                double h3 = h2 * h_i;
                                if (q < 0.5) {
                                    dphi_dr_i = (1.0 / h2) *
                                                (10.6666666667 * q - 19.2 * q3 +
                                                 10.6666666667 * q4);
                                    dW_dr_i = (8.0 / (M_PI * h3)) *
                                              (-12.0 * q + 18.0 * q2) / h_i;
                                } else {
                                    dphi_dr_i = (1.0 / h2) *
                                                (10.6666666667 * q - 19.2 * q3 +
                                                 10.6666666667 * q4 -
                                                 0.0666666667 / q2 + 0.1 * q5 -
                                                 3.2 * q4 + 10.6666666667 * q3 -
                                                 16.0 * q2 + 11.2 * q - 3.2);
                                    double diff = 1.0 - q;
                                    dW_dr_i = (8.0 / (M_PI * h3)) *
                                              (-6.0 * diff * diff) / h_i;
                                }
                            }

                            // Inline Kernel Derivs for Particle J
                            double dphi_dr_j, dW_dr_j;
                            if (r >= h_j) {
                                dphi_dr_j = 1.0 / (r * r);
                                dW_dr_j = 0.0;
                            } else {
                                double q = r / h_j;
                                double q2 = q * q;
                                double q3 = q2 * q;
                                double q4 = q3 * q;
                                double q5 = q4 * q;
                                double h2 = h_j * h_j;
                                double h3 = h2 * h_j;
                                if (q < 0.5) {
                                    dphi_dr_j = (1.0 / h2) *
                                                (10.6666666667 * q - 19.2 * q3 +
                                                 10.6666666667 * q4);
                                    dW_dr_j = (8.0 / (M_PI * h3)) *
                                              (-12.0 * q + 18.0 * q2) / h_j;
                                } else {
                                    dphi_dr_j = (1.0 / h2) *
                                                (10.6666666667 * q - 19.2 * q3 +
                                                 10.6666666667 * q4 -
                                                 0.0666666667 / q2 + 0.1 * q5 -
                                                 3.2 * q4 + 10.6666666667 * q3 -
                                                 16.0 * q2 + 11.2 * q - 3.2);
                                    double diff = 1.0 - q;
                                    dW_dr_j = (8.0 / (M_PI * h3)) *
                                              (-6.0 * diff * diff) / h_j;
                                }
                            }

                            // Symmetrized H8 Equation
                            double force_mag_over_r =
                                (G / 2.0) *
                                ((dphi_dr_i + dphi_dr_j) +
                                 (zeta_i * dW_dr_i) / m_i +
                                 (zeta_j * dW_dr_j) / m_j) /
                                r;

                            if (use_pm) {
                                double r_scaled = r / (2.0 * r_s);
                                double taper =
                                    std::erfc(r_scaled) +
                                    (r / (std::sqrt(M_PI) * r_s)) *
                                        std::exp(-r_scaled * r_scaled);
                                force_mag_over_r *= taper;
                            }

                            local_acc_x += force_mag_over_r * m_j * dx;
                            local_acc_y += force_mag_over_r * m_j * dy;
                            local_acc_z += force_mag_over_r * m_j * dz;
                        }
                    }
                }
            }
            d_ax[i] += local_acc_x;
            d_ay[i] += local_acc_y;
            d_az[i] += local_acc_z;
        }
        auto end_compute = std::chrono::high_resolution_clock::now();

        auto start_return = std::chrono::high_resolution_clock::now();
#pragma omp target exit data map(from : d_ax[0 : n_gas], d_ay[0 : n_gas], \
                                     d_az[0 : n_gas])                     \
    map(delete : d_px[0 : n_gas], d_py[0 : n_gas], d_pz[0 : n_gas],       \
            d_m[0 : n_gas], d_h[0 : n_gas], d_zeta[0 : n_gas],            \
            d_cell_start[0 : num_cells], d_cell_count[0 : num_cells])
        auto end_return = std::chrono::high_resolution_clock::now();

        double diff_transf =
            std::chrono::duration_cast<std::chrono::microseconds>(
                end_transfer - start_transfer)
                .count();
        double diff_comput =
            std::chrono::duration_cast<std::chrono::microseconds>(end_compute -
                                                                  start_compute)
                .count();
        double diff_return =
            std::chrono::duration_cast<std::chrono::microseconds>(end_return -
                                                                  start_return)
                .count();

        diag.add_prof_time(ProfRegion::Transf, diff_transf);
        diag.add_prof_time(ProfRegion::Compute, diff_comput);
        diag.add_prof_time(ProfRegion::Ret, diff_return);
    } else
#endif
    {
        // ========================================================================
        // CPU IMPLEMENTATION (Sorted Array + Chunked Gather Buffer)
        // ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < n_gas; ++i) {
            double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
            double m_i = mass[i];
            double h_i = h[i];
            double zeta_i = zeta[i];

            int ix = static_cast<int>(p1_x / config.cell_size);
            int iy = static_cast<int>(p1_y / config.cell_size);
            int iz = static_cast<int>(p1_z / config.cell_size);
            ix = (((ix % N) + N) % N);
            iy = (((iy % N) + N) % N);
            iz = (((iz % N) + N) % N);

            constexpr int MAX_NEIGHBORS = 1536;
            double n_x[MAX_NEIGHBORS];
            double n_y[MAX_NEIGHBORS];
            double n_z[MAX_NEIGHBORS];
            double n_m[MAX_NEIGHBORS];
            double n_h[MAX_NEIGHBORS];
            double n_zeta[MAX_NEIGHBORS];
            int num_neighbors = 0;

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

            // Lambda to execute the SIMD math and flush the buffer
            auto compute_simd_batch = [&]() {
#pragma omp simd reduction(+ : local_acc_x, local_acc_y, local_acc_z)
                for (int k = 0; k < num_neighbors; ++k) {
                    double dx =
                        mfm_periodic_displacement(n_x[k] - p1_x, domain_size);
                    double dy =
                        mfm_periodic_displacement(n_y[k] - p1_y, domain_size);
                    double dz =
                        mfm_periodic_displacement(n_z[k] - p1_z, domain_size);

                    double dist_sq = dx * dx + dy * dy + dz * dz;

                    if (use_pm && dist_sq > cutoff_sq) continue;
                    if (dist_sq < 1e-24) continue;

                    double r = std::sqrt(dist_sq);
                    double m_j = n_m[k];
                    double h_j = n_h[k];
                    double zeta_j = n_zeta[k];

                    double dphi_dr_i, dW_dr_i;
                    compute_adaptive_gravity_terms(r, h_i, dphi_dr_i, dW_dr_i);

                    double dphi_dr_j, dW_dr_j;
                    compute_adaptive_gravity_terms(r, h_j, dphi_dr_j, dW_dr_j);

                    // Equation H8: Fully Symmetrized Force (divided by r to
                    // multiply by distance vector)
                    double force_mag_over_r =
                        (G / 2.0) *
                        ((dphi_dr_i + dphi_dr_j) + (zeta_i * dW_dr_i) / m_i +
                         (zeta_j * dW_dr_j) / m_j) /
                        r;

                    if (use_pm) {
                        double r_scaled = r / (2.0 * r_s);
                        double taper = std::erfc(r_scaled) +
                                       (r / (std::sqrt(M_PI) * r_s)) *
                                           std::exp(-r_scaled * r_scaled);
                        force_mag_over_r *= taper;
                    }

                    local_acc_x += force_mag_over_r * m_j * dx;
                    local_acc_y += force_mag_over_r * m_j * dy;
                    local_acc_z += force_mag_over_r * m_j * dz;
                }
                num_neighbors = 0;  // Reset the buffer counter
            };

            // Gather: Linearly pull from sorted memory
            for (int dx_cell = dx_start; dx_cell <= dx_end; ++dx_cell) {
                for (int dy_cell = dx_start; dy_cell <= dx_end; ++dy_cell) {
                    for (int dz_cell = dx_start; dz_cell <= dx_end; ++dz_cell) {
                        int neighbor_ix =
                            use_pm ? (((ix + dx_cell) % N) + N) % N : dx_cell;
                        int neighbor_iy =
                            use_pm ? (((iy + dy_cell) % N) + N) % N : dy_cell;
                        int neighbor_iz =
                            use_pm ? (((iz + dz_cell) % N) + N) % N : dz_cell;
                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        int start_idx = pm_cell_list.cell_start[cell_idx];
                        int end_idx =
                            start_idx + pm_cell_list.cell_count[cell_idx];

                        for (int j = start_idx; j < end_idx; ++j) {
                            if (i != j) {
                                n_x[num_neighbors] = pos_x[j];
                                n_y[num_neighbors] = pos_y[j];
                                n_z[num_neighbors] = pos_z[j];
                                n_m[num_neighbors] = mass[j];
                                n_h[num_neighbors] = h[j];
                                n_zeta[num_neighbors] = zeta[j];
                                num_neighbors++;

                                // Flush the buffer if it gets full
                                if (num_neighbors == MAX_NEIGHBORS) {
                                    compute_simd_batch();
                                }
                            }
                        }
                    }
                }
            }

            // Process any remaining particles in the buffer
            if (num_neighbors > 0) {
                compute_simd_batch();
            }

            acc_x[i] += local_acc_x;
            acc_y[i] += local_acc_y;
            acc_z[i] += local_acc_z;
        }
    }
}

void GasParticleSystem::compute_cross_pp_forces(ParticleSystem& dm,
                                                const Config& config,
                                                Diagnostics& diag) {
    const int N = config.mesh_size;
    const double cell_size = config.cell_size;
    const double domain_size = config.domain_size;
    const double G = config.G;
    const double soft_sq = config.softening_squared;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * cell_size;
    const bool use_pm = config.use_PM;

    const size_t n_gas = num_particles;
    const size_t n_dm = dm.num_particles;
    const int num_cells = N * N * N;

    int dx_start =
        use_pm ? -static_cast<int>(std::ceil(config.cutoff_radius / cell_size))
               : 0;
    int dx_end =
        use_pm ? static_cast<int>(std::ceil(config.cutoff_radius / cell_size))
               : N - 1;

#ifdef USE_GPU
    if (config.enable_GPU) {
        // ========================================================================
        // GPU IMPLEMENTATION (Cross-forces)
        // ========================================================================
        double* d_gas_px = pos_x.data();
        double* d_gas_py = pos_y.data();
        double* d_gas_pz = pos_z.data();
        double* d_gas_m = mass.data();
        double* d_gas_h = h.data();
        double* d_gas_zeta = zeta.data();
        double* d_gas_ax = acc_x.data();
        double* d_gas_ay = acc_y.data();
        double* d_gas_az = acc_z.data();

        double* d_dm_px = dm.pos_x.data();
        double* d_dm_py = dm.pos_y.data();
        double* d_dm_pz = dm.pos_z.data();
        double* d_dm_m = dm.mass.data();
        double* d_dm_ax = dm.acc_x.data();
        double* d_dm_ay = dm.acc_y.data();
        double* d_dm_az = dm.acc_z.data();
        int* d_dm_cell_start = dm.cell_list.cell_start.data();
        int* d_dm_cell_count = dm.cell_list.cell_count.data();

        auto start_transfer = std::chrono::high_resolution_clock::now();
#pragma omp target enter data map(                                          \
        to : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], \
            d_gas_m[0 : n_gas], d_gas_h[0 : n_gas], d_gas_zeta[0 : n_gas],  \
            d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas],  \
            d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],        \
            d_dm_m[0 : n_dm], d_dm_cell_start[0 : num_cells],               \
            d_dm_cell_count[0 : num_cells], d_dm_ax[0 : n_dm],              \
            d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm])
        auto end_transfer = std::chrono::high_resolution_clock::now();

        auto start_compute = std::chrono::high_resolution_clock::now();
#pragma omp target teams distribute parallel for
        for (size_t i = 0; i < n_gas; ++i) {
            double p1_x = d_gas_px[i], p1_y = d_gas_py[i], p1_z = d_gas_pz[i];
            int ix = static_cast<int>(p1_x / cell_size);
            int iy = static_cast<int>(p1_y / cell_size);
            int iz = static_cast<int>(p1_z / cell_size);

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
            double m_gas = d_gas_m[i];
            double h_i = d_gas_h[i];
            double zeta_i = d_gas_zeta[i];
            double base_soft = std::sqrt(soft_sq);

            for (int dx_cell = dx_start; dx_cell <= dx_end; ++dx_cell) {
                for (int dy_cell = dx_start; dy_cell <= dx_end; ++dy_cell) {
                    for (int dz_cell = dx_start; dz_cell <= dx_end; ++dz_cell) {
                        int neighbor_ix =
                            use_pm ? (((ix + dx_cell) % N) + N) % N : dx_cell;
                        int neighbor_iy =
                            use_pm ? (((iy + dy_cell) % N) + N) % N : dy_cell;
                        int neighbor_iz =
                            use_pm ? (((iz + dz_cell) % N) + N) % N : dz_cell;
                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        int start = d_dm_cell_start[cell_idx];
                        int end = start + d_dm_cell_count[cell_idx];

                        for (int j = start; j < end; ++j) {
                            double dx = p1_x - d_dm_px[j];
                            if (dx > 0.5 * domain_size)
                                dx -= domain_size;
                            else if (dx < -0.5 * domain_size)
                                dx += domain_size;

                            double dy = p1_y - d_dm_py[j];
                            if (dy > 0.5 * domain_size)
                                dy -= domain_size;
                            else if (dy < -0.5 * domain_size)
                                dy += domain_size;

                            double dz = p1_z - d_dm_pz[j];
                            if (dz > 0.5 * domain_size)
                                dz -= domain_size;
                            else if (dz < -0.5 * domain_size)
                                dz += domain_size;

                            // Flip dx, dy, dz back to (j - i)
                            dx = -dx;
                            dy = -dy;
                            dz = -dz;

                            double dist_sq = dx * dx + dy * dy + dz * dz;

                            if (use_pm && dist_sq > cutoff_sq) continue;

                            double r = std::sqrt(dist_sq + 1e-24);
                            double m_j = d_dm_m[j];

                            // Evaluate kernel derivatives for Gas particle
                            // (uses h_i)
                            double dphi_dr_i, dW_dr_i;
                            if (r >= h_i) {
                                dphi_dr_i = 1.0 / (r * r);
                                dW_dr_i = 0.0;
                            } else {
                                double q = r / h_i;
                                double q2 = q * q;
                                double q3 = q2 * q;
                                double q4 = q3 * q;
                                double q5 = q4 * q;
                                double h2 = h_i * h_i;
                                double h3 = h2 * h_i;

                                if (q < 0.5) {
                                    dphi_dr_i = (1.0 / h2) *
                                                (10.6666666667 * q - 19.2 * q3 +
                                                 10.6666666667 * q4);
                                    dW_dr_i = (8.0 / (M_PI * h3)) *
                                              (-12.0 * q + 18.0 * q2) / h_i;
                                } else {
                                    dphi_dr_i = (1.0 / h2) *
                                                (10.6666666667 * q - 19.2 * q3 +
                                                 10.6666666667 * q4 -
                                                 0.0666666667 / q2 + 0.1 * q5 -
                                                 3.2 * q4 + 10.6666666667 * q3 -
                                                 16.0 * q2 + 11.2 * q - 3.2);
                                    double diff = 1.0 - q;
                                    dW_dr_i = (8.0 / (M_PI * h3)) *
                                              (-6.0 * diff * diff) / h_i;
                                }
                            }

                            // Evaluate kernel derivatives for DM particle (uses
                            // fixed base_soft)
                            double dphi_dr_j;
                            if (r >= base_soft) {
                                dphi_dr_j = 1.0 / (r * r);
                            } else {
                                double q = r / base_soft;
                                double q2 = q * q;
                                double q3 = q2 * q;
                                double q4 = q3 * q;
                                double q5 = q4 * q;
                                double h2 = base_soft * base_soft;

                                if (q < 0.5) {
                                    dphi_dr_j = (1.0 / h2) *
                                                (10.6666666667 * q - 19.2 * q3 +
                                                 10.6666666667 * q4);
                                } else {
                                    dphi_dr_j = (1.0 / h2) *
                                                (10.6666666667 * q - 19.2 * q3 +
                                                 10.6666666667 * q4 -
                                                 0.0666666667 / q2 + 0.1 * q5 -
                                                 3.2 * q4 + 10.6666666667 * q3 -
                                                 16.0 * q2 + 11.2 * q - 3.2);
                                }
                            }

                            // Symmetrized Force Equation (GIZMO App. H8)
                            double force_mag_over_r =
                                (G / 2.0) *
                                ((dphi_dr_i + dphi_dr_j) +
                                 (zeta_i * dW_dr_i) / m_gas) /
                                r;

                            if (use_pm) {
                                double r_scaled = r / (2.0 * r_s);
                                double taper =
                                    std::erfc(r_scaled) +
                                    (r / (std::sqrt(M_PI) * r_s)) *
                                        std::exp(-r_scaled * r_scaled);
                                force_mag_over_r *= taper;
                            }

                            // Pull on MFM Gas
                            local_acc_x += force_mag_over_r * m_j * dx;
                            local_acc_y += force_mag_over_r * m_j * dy;
                            local_acc_z += force_mag_over_r * m_j * dz;

                            // Pull on DM (Equal and opposite)
                            double a_dm = force_mag_over_r * m_gas;
                            double dm_ax = -a_dm * dx;
                            double dm_ay = -a_dm * dy;
                            double dm_az = -a_dm * dz;

#pragma omp atomic
                            d_dm_ax[j] += dm_ax;
#pragma omp atomic
                            d_dm_ay[j] += dm_ay;
#pragma omp atomic
                            d_dm_az[j] += dm_az;
                        }
                    }
                }
            }
            d_gas_ax[i] += local_acc_x;
            d_gas_ay[i] += local_acc_y;
            d_gas_az[i] += local_acc_z;
        }
        auto end_compute = std::chrono::high_resolution_clock::now();

        auto start_return = std::chrono::high_resolution_clock::now();
#pragma omp target exit data map(                                             \
        from : d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas], \
            d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm])          \
    map(delete : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas],                    \
            d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas], d_gas_h[0 : n_gas],      \
            d_gas_zeta[0 : n_gas], d_dm_px[0 : n_dm], d_dm_py[0 : n_dm],      \
            d_dm_pz[0 : n_dm], d_dm_m[0 : n_dm],                              \
            d_dm_cell_start[0 : num_cells], d_dm_cell_count[0 : num_cells])
        auto end_return = std::chrono::high_resolution_clock::now();

        double diff_transf =
            std::chrono::duration_cast<std::chrono::microseconds>(
                end_transfer - start_transfer)
                .count();
        double diff_comput =
            std::chrono::duration_cast<std::chrono::microseconds>(end_compute -
                                                                  start_compute)
                .count();
        double diff_return =
            std::chrono::duration_cast<std::chrono::microseconds>(end_return -
                                                                  start_return)
                .count();

        diag.add_prof_time(ProfRegion::Transf, diff_transf);
        diag.add_prof_time(ProfRegion::Compute, diff_comput);
        diag.add_prof_time(ProfRegion::Ret, diff_return);
    } else
#endif
    {
        // ========================================================================
        // CPU IMPLEMENTATION (Sorted Array + Chunked Gather Buffer)
        // ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < n_gas; ++i) {
            double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
            int ix = static_cast<int>(p1_x / cell_size);
            int iy = static_cast<int>(p1_y / cell_size);
            int iz = static_cast<int>(p1_z / cell_size);
            ix = (((ix % N) + N) % N);
            iy = (((iy % N) + N) % N);
            iz = (((iz % N) + N) % N);

            constexpr int MAX_NEIGHBORS = 1536;
            double n_x[MAX_NEIGHBORS];
            double n_y[MAX_NEIGHBORS];
            double n_z[MAX_NEIGHBORS];
            double n_m[MAX_NEIGHBORS];
            int n_j[MAX_NEIGHBORS];  // Track original DM particle index for
                                     // atomic return force
            int num_neighbors = 0;

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
            double m_gas = mass[i];
            double h_i = h[i];
            double zeta_i = zeta[i];

            // Lambda to execute SIMD math and flush the buffer
            auto compute_simd_batch = [&]() {
                double temp_dm_ax[MAX_NEIGHBORS];
                double temp_dm_ay[MAX_NEIGHBORS];
                double temp_dm_az[MAX_NEIGHBORS];

                double base_soft = std::sqrt(soft_sq);

#pragma omp simd reduction(+ : local_acc_x, local_acc_y, local_acc_z)
                for (int k = 0; k < num_neighbors; ++k) {
                    double dx =
                        mfm_periodic_displacement(n_x[k] - p1_x, domain_size);
                    double dy =
                        mfm_periodic_displacement(n_y[k] - p1_y, domain_size);
                    double dz =
                        mfm_periodic_displacement(n_z[k] - p1_z, domain_size);
                    double dist_sq = dx * dx + dy * dy + dz * dz;

                    if (use_pm && dist_sq > cutoff_sq) {
                        temp_dm_ax[k] = 0.0;
                        temp_dm_ay[k] = 0.0;
                        temp_dm_az[k] = 0.0;
                        continue;
                    }

                    double r = std::sqrt(dist_sq + 1e-24);

                    // Evaluate kernel derivatives for Gas particle (uses h_i)
                    double dphi_dr_i, dW_dr_i;
                    compute_adaptive_gravity_terms(r, h_i, dphi_dr_i, dW_dr_i);

                    // Evaluate kernel derivatives for DM particle (uses fixed
                    // base_soft)
                    double dphi_dr_j, dW_dr_j_dummy;
                    compute_adaptive_gravity_terms(r, base_soft, dphi_dr_j,
                                                   dW_dr_j_dummy);

                    double m_j = n_m[k];

                    // Symmetrized Force Equation (GIZMO App. H8)
                    // Treat DM as a fluid particle with zeta_j = 0
                    double force_mag_over_r =
                        (G / 2.0) *
                        ((dphi_dr_i + dphi_dr_j) + (zeta_i * dW_dr_i) / m_gas) /
                        r;

                    if (use_pm) {
                        double r_scaled = r / (2.0 * r_s);
                        double taper = std::erfc(r_scaled) +
                                       (r / (std::sqrt(M_PI) * r_s)) *
                                           std::exp(-r_scaled * r_scaled);
                        force_mag_over_r *= taper;
                    }

                    local_acc_x += force_mag_over_r * m_j * dx;
                    local_acc_y += force_mag_over_r * m_j * dy;
                    local_acc_z += force_mag_over_r * m_j * dz;

                    double a_dm = force_mag_over_r * m_gas;
                    temp_dm_ax[k] = -a_dm * dx;
                    temp_dm_ay[k] = -a_dm * dy;
                    temp_dm_az[k] = -a_dm * dz;
                }

                // Apply atomics OUTSIDE the SIMD loop to protect vectorization
                for (int k = 0; k < num_neighbors; ++k) {
                    if (temp_dm_ax[k] != 0.0 || temp_dm_ay[k] != 0.0 ||
                        temp_dm_az[k] != 0.0) {
                        int j_orig = n_j[k];
#pragma omp atomic
                        dm.acc_x[j_orig] += temp_dm_ax[k];
#pragma omp atomic
                        dm.acc_y[j_orig] += temp_dm_ay[k];
#pragma omp atomic
                        dm.acc_z[j_orig] += temp_dm_az[k];
                    }
                }

                num_neighbors = 0;  // Reset the buffer counter
            };

            for (int dx_cell = dx_start; dx_cell <= dx_end; ++dx_cell) {
                for (int dy_cell = dx_start; dy_cell <= dx_end; ++dy_cell) {
                    for (int dz_cell = dx_start; dz_cell <= dx_end; ++dz_cell) {
                        int neighbor_ix =
                            use_pm ? (((ix + dx_cell) % N) + N) % N : dx_cell;
                        int neighbor_iy =
                            use_pm ? (((iy + dy_cell) % N) + N) % N : dy_cell;
                        int neighbor_iz =
                            use_pm ? (((iz + dz_cell) % N) + N) % N : dz_cell;
                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        int start = dm.cell_list.cell_start[cell_idx];
                        int end = start + dm.cell_list.cell_count[cell_idx];

                        for (int j = start; j < end; ++j) {
                            n_x[num_neighbors] = dm.pos_x[j];
                            n_y[num_neighbors] = dm.pos_y[j];
                            n_z[num_neighbors] = dm.pos_z[j];
                            n_m[num_neighbors] = dm.mass[j];
                            n_j[num_neighbors] =
                                j;  // Cache the target DM index
                            num_neighbors++;

                            // Flush the buffer if it gets full
                            if (num_neighbors == MAX_NEIGHBORS) {
                                compute_simd_batch();
                            }
                        }
                    }
                }
            }

            // Process any remaining particles in the buffer
            if (num_neighbors > 0) {
                compute_simd_batch();
            }

            acc_x[i] += local_acc_x;
            acc_y[i] += local_acc_y;
            acc_z[i] += local_acc_z;
        }
    }
}

void GasParticleSystem::hydro_step(const Config& config, double a, double dt) {
    update_primitive_variables(config, a);
    compute_gradients(config);
    compute_hydro_forces(config, a, dt);
}

void GasParticleSystem::update_primitive_variables(const Config& config,
                                                   double a) {
    if (num_particles == 0) return;

    static int num_cycles = 0;

    const double gamma_minus_1 = config.gamma - 1.0;

    // Calculate the physical internal energy floor
    const double u_floor =
        Cooling::get_internal_energy_from_temp(config.temp_floor_k, a, config);

    const double alpha_kin = 0.001;
    const double alpha_grav = 0.001;

    double step_floor_heating = 0.0;
    double step_entropy_switch = 0.0;

#pragma omp parallel for simd reduction( \
        + : step_floor_heating, step_entropy_switch) schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        rho[i] = std::max(rho[i], density_floor);

        // Calculate specific kinetic energy
        double ke = 0.5 * (vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] +
                           vel_z[i] * vel_z[i]);

        // Evaluate the Energy-Entropy Switch
        double threshold_kin = alpha_kin * (max_rel_ke[i] + u[i]);
        double threshold_grav = alpha_grav * delta_E_grav[i];

        if (u[i] < threshold_kin || u[i] < threshold_grav) {
            // FALLBACK TRIGGERED: Extreme Mach number detected.
            // Discard the 'u' updated by the Riemann solver and use
            // entropy-based adiabatic evolution. S = (gamma-1) * u /
            // rho^(gamma-1) -> u = S * rho^(gamma-1) / (gamma-1)
            double u_new =
                entropy[i] * std::pow(rho[i], gamma_minus_1) / gamma_minus_1;

            // Track the energy injected (or removed) by the switch
            step_entropy_switch += (u_new - u[i]) * mass[i];
            u[i] = u_new;
        } else {
            // NORMAL REGIME: Trust the integrated internal energy.
            // Re-sync the entropy array to match the shock-heated state.
            entropy[i] = gamma_minus_1 * u[i] / std::pow(rho[i], gamma_minus_1);
        }

        // Enforce thermodynamic floor
        if (u[i] < u_floor) {
            // Track the added energy
            double injected_u = u_floor - u[i];
            step_floor_heating += injected_u * mass[i];

            u[i] = u_floor;
            // Sync the entropy array again if we hit the temperature floor
            entropy[i] = gamma_minus_1 * u[i] / std::pow(rho[i], gamma_minus_1);
        }

        // Calculate pressure based on 'u'
        pressure[i] = gamma_minus_1 * rho[i] * u[i];

        // Passively sync total energy for diagnostic conservation tracking
        // Total energy no longer drives the hydrodynamics
        total_energy[i] = u[i] + ke;

        metal_frac[i] = std::max(0.0, std::min(metal_frac[i], 1.0));
    }
    num_cycles++;

    accumulated_photoheating_energy += step_floor_heating;
    accumulated_entropy_switch_energy += step_entropy_switch;
}

// Cubic spline kernel
inline void compute_kernel(double r, double h, double& W) {
    double q = r / h;
    double norm = 8.0 / (M_PI * h * h * h);
    if (q < 0.5) {
        W = norm * (1.0 - 6.0 * q * q + 6.0 * q * q * q);
    } else if (q < 1.0) {
        double diff = 1.0 - q;
        W = norm * 2.0 * diff * diff * diff;
    } else {
        W = 0.0;
    }
}

#ifdef THEORETICAL_LIMITER
// Helper for the Pairwise Limiter signs
inline int math_sign(double x) { return (x > 0.0) ? 1 : ((x < 0.0) ? -1 : 0); }

// Compute the alpha scaling factor for the Scalar Gradient Limiter
// (Eq. B1 & B2)
inline double compute_gradient_alpha(double d_phi_max, double d_phi_min,
                                     double d_mid_max, double d_mid_min,
                                     double beta) {
    double alpha = 1.0;

    // Check positive deviations
    if (d_mid_max > 0.0) {
        if (d_phi_max <= 0.0)
            alpha = 0.0;
        else
            alpha = std::min(alpha, beta * d_phi_max / d_mid_max);
    }

    // Check negative deviations
    if (d_mid_min < 0.0) {
        if (d_phi_min >= 0.0)
            alpha = 0.0;
        else
            alpha = std::min(alpha, beta * d_phi_min / d_mid_min);
    }
    return alpha;
}

// The Pairwise TVD Limiter applied directly at the face (Eq. B4)
double apply_pairwise_limiter(double phi_L_center, double phi_R_center,
                              double phi_mid_0, double phi_bar) {
    // Tunable constants from the paper
    double psi_1 = 0.5;
    double psi_2 = 0.25;

    double d_phi = std::abs(phi_L_center - phi_R_center);
    if (d_phi < 1e-14) return phi_L_center;

    double phi_min = std::min(phi_L_center, phi_R_center);
    double phi_max = std::max(phi_L_center, phi_R_center);

    double delta_1 = psi_1 * d_phi;
    double delta_2 = psi_2 * d_phi;

    // Calculate bounding limits (phi_minus and phi_plus)
    double phi_minus;
    if (math_sign(phi_min - delta_1) == math_sign(phi_min) || phi_min == 0.0) {
        phi_minus = phi_min - delta_1;
    } else {
        phi_minus =
            phi_min /
            (1.0 + delta_1 / std::abs(phi_min));  // Positivity preserving
    }

    double phi_plus;
    if (math_sign(phi_max + delta_1) == math_sign(phi_max) || phi_max == 0.0) {
        phi_plus = phi_max + delta_1;
    } else {
        phi_plus = phi_max / (1.0 + delta_1 / std::abs(phi_max));
    }

    // Apply the min/max bounds based on the gradient direction
    if (phi_L_center < phi_R_center) {
        return std::max(phi_minus, std::min(phi_bar + delta_2, phi_mid_0));
    } else {
        return std::min(phi_plus, std::max(phi_bar - delta_2, phi_mid_0));
    }
}

#else
// Scalar Gradient Vector Limiter (Translated from GIZMO's local_slopelimiter)
inline void scalar_limiter(Eigen::Vector3d& grad, double valmax, double valmin,
                           double alim, double h, double shoot_tol,
                           bool pos_preserve, double d_max, double val_cen) {
    double d_abs = grad.norm();
    if (d_abs > 0.0) {
        // Inverse change over distance for limiter
        double cfac = 1.0 / (alim * h * d_abs);

        double abs_max = std::abs(valmax);
        double abs_min = std::abs(valmin);

        // Get largest positive/negative deviations, determine smaller in
        // absolute value
        if (abs_max < abs_min) {
            std::swap(abs_max, abs_min);
        }

        // = abs_min for shoot_tol = 0; don't let gradient deviate by more than
        // this in size, slightly larger if 'shoot_tol' allows some overshoot
        // tolerance
        double f_corr_overshoot =
            std::min(abs_min + shoot_tol * abs_max, abs_max);

        // Multiply by the correction factor of interest
        cfac *= f_corr_overshoot;

        // Demand that the limited slope be strictly positivity-preserving over
        // the maximal range to any neighbors
        if (pos_preserve) {
            constexpr double MIN_REAL_NUMBER = 1e-30;

            // Minimum value: smaller of overshoot target or half
            // positive-definite value, but cannot go negative in larger range
            double fmin = std::min(
                val_cen,
                std::max(0.0, std::max(MIN_REAL_NUMBER * val_cen,
                                       std::min(0.5 * (val_cen + valmin),
                                                val_cen - f_corr_overshoot))));

            // Use more conservative limiter, of cfac above or this,
            // over longer range d_max, to restrict here
            cfac = std::min((((val_cen - fmin) / d_max) / d_abs), cfac);
        }

        // Scalar gradient correction
        if (cfac < 1.0) {
            grad *= cfac;
        }
    }
}
#endif

inline double compute_condition_number(const Eigen::Matrix3d& E,
                                       const Eigen::Matrix3d& B) {
    // Eq. C2: The sum of the squared elements of the matrices
    double norm_E_sq = E.squaredNorm();
    double norm_B_sq = B.squaredNorm();

    // Eq. C1: N_cond = (1 / v) * sqrt(||E^-1|| * ||E||) where v = 3
    return (1.0 / 3.0) * std::sqrt(norm_E_sq * norm_B_sq);
}

// MFM Gradient Estimator
ParticleGradients compute_single_particle_gradients(
    const ParticleState& p_i, const std::vector<ParticleState>& neighbors,
    double domain_size) {
    ParticleGradients out;
    out.B_matrix = Eigen::Matrix3d::Zero();
    out.grad_rho = Eigen::Vector3d::Zero();
    out.grad_p = Eigen::Vector3d::Zero();
    out.grad_vx = Eigen::Vector3d::Zero();
    out.grad_vy = Eigen::Vector3d::Zero();
    out.grad_vz = Eigen::Vector3d::Zero();
    out.ill_conditioned = false;

    Eigen::Matrix3d E = Eigen::Matrix3d::Zero();

    // Build the E Matrix
    for (const auto& nj : neighbors) {
        double dx =
            mfm_periodic_displacement(nj.pos.x() - p_i.pos.x(), domain_size);
        double dy =
            mfm_periodic_displacement(nj.pos.y() - p_i.pos.y(), domain_size);
        double dz =
            mfm_periodic_displacement(nj.pos.z() - p_i.pos.z(), domain_size);

        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < p_i.h * p_i.h && r2 > 1e-24) {
            double r = std::sqrt(r2);
            double W;
            compute_kernel(r, p_i.h, W);

            double V_j = nj.mass / nj.rho;
            double weight = W * V_j;

            E(0, 0) += dx * dx * weight;
            E(0, 1) += dx * dy * weight;
            E(0, 2) += dx * dz * weight;
            E(1, 1) += dy * dy * weight;
            E(1, 2) += dy * dz * weight;
            E(2, 2) += dz * dz * weight;
        }
    }

    E(1, 0) = E(0, 1);
    E(2, 0) = E(0, 2);
    E(2, 1) = E(1, 2);

    double det = E.determinant();
    out.ill_conditioned =
        true;  // Assume ill-conditioned until proven otherwise
    out.condition_number = -1.0;

    // We still need a tiny guard so E.inverse() doesn't crash on pure zeroes
    if (std::abs(det) > 1e-30) {
        Eigen::Matrix3d temp_B = E.inverse();

        // Calculate the condition number
        double N_cond = compute_condition_number(E, temp_B);
        out.condition_number = N_cond;

        // Threshold (Hopkins recommends 100 - 1000)
        constexpr double N_cond_crit = 1000.0;

        if (N_cond <= N_cond_crit) {
            out.B_matrix = temp_B;
            out.ill_conditioned = false;
        }
    }

    if (!out.ill_conditioned) {
        out.B_matrix = E.inverse();
        Eigen::Vector3d sum_rho = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_p = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vx = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vy = Eigen::Vector3d::Zero();
        Eigen::Vector3d sum_vz = Eigen::Vector3d::Zero();

        // Track the min/max differences for the face-extrapolated limiter
        double d_rho_max = 0.0, d_rho_min = 0.0;
        double d_p_max = 0.0, d_p_min = 0.0;
        double d_vx_max = 0.0, d_vx_min = 0.0;
        double d_vy_max = 0.0, d_vy_min = 0.0;
        double d_vz_max = 0.0, d_vz_min = 0.0;

        double r_max = 1e-12;

        // Calculate sums and record extrema
        for (const auto& nj : neighbors) {
            double dx = mfm_periodic_displacement(nj.pos.x() - p_i.pos.x(),
                                                  domain_size);
            double dy = mfm_periodic_displacement(nj.pos.y() - p_i.pos.y(),
                                                  domain_size);
            double dz = mfm_periodic_displacement(nj.pos.z() - p_i.pos.z(),
                                                  domain_size);

            double r2 = dx * dx + dy * dy + dz * dz;
            // Limiter bounds: Record extrema
            if ((r2 < p_i.h * p_i.h || r2 < nj.h * nj.h) && r2 > 1e-24) {
                double r = std::sqrt(r2);
                r_max = std::max(r_max, r);

                // Update min/max neighbor differences
                double d_rho = nj.rho - p_i.rho;
                d_rho_max = std::max(d_rho_max, d_rho);
                d_rho_min = std::min(d_rho_min, d_rho);

                double d_p = nj.pressure - p_i.pressure;
                d_p_max = std::max(d_p_max, d_p);
                d_p_min = std::min(d_p_min, d_p);

                double d_vx = nj.vel.x() - p_i.vel.x();
                d_vx_max = std::max(d_vx_max, d_vx);
                d_vx_min = std::min(d_vx_min, d_vx);

                double d_vy = nj.vel.y() - p_i.vel.y();
                d_vy_max = std::max(d_vy_max, d_vy);
                d_vy_min = std::min(d_vy_min, d_vy);

                double d_vz = nj.vel.z() - p_i.vel.z();
                d_vz_max = std::max(d_vz_max, d_vz);
                d_vz_min = std::min(d_vz_min, d_vz);
            }

            // Gradient sums
            if (r2 < p_i.h * p_i.h && r2 > 1e-24) {
                double r = std::sqrt(r2);
                double W;
                compute_kernel(r, p_i.h, W);

                double V_j = nj.mass / nj.rho;
                double weight = W * V_j;
                Eigen::Vector3d dx_vec(dx, dy, dz);

                sum_rho += (nj.rho - p_i.rho) * dx_vec * weight;
                sum_p += (nj.pressure - p_i.pressure) * dx_vec * weight;
                sum_vx += (nj.vel.x() - p_i.vel.x()) * dx_vec * weight;
                sum_vy += (nj.vel.y() - p_i.vel.y()) * dx_vec * weight;
                sum_vz += (nj.vel.z() - p_i.vel.z()) * dx_vec * weight;
            }
        }

        // Apply the B Matrix
        out.grad_rho = out.B_matrix * sum_rho;
        out.grad_p = out.B_matrix * sum_p;
        out.grad_vx = out.B_matrix * sum_vx;
        out.grad_vy = out.B_matrix * sum_vy;
        out.grad_vz = out.B_matrix * sum_vz;
        out.raw_sum_p = sum_p;

#ifdef THEORETICAL_LIMITER
        // Scalar Gradient Limiter
        double phi_mid_max_rho = 0.0, phi_mid_min_rho = 0.0;
        double phi_mid_max_p = 0.0, phi_mid_min_p = 0.0;
        double phi_mid_max_vx = 0.0, phi_mid_min_vx = 0.0;
        double phi_mid_max_vy = 0.0, phi_mid_min_vy = 0.0;
        double phi_mid_max_vz = 0.0, phi_mid_min_vz = 0.0;

        // Second pass: Find the extrapolated mid-point extrema
        for (const auto& nj : neighbors) {
            double dx = mfm_periodic_displacement(nj.pos.x() - p_i.pos.x(),
                                                  domain_size);
            double dy = mfm_periodic_displacement(nj.pos.y() - p_i.pos.y(),
                                                  domain_size);
            double dz = mfm_periodic_displacement(nj.pos.z() - p_i.pos.z(),
                                                  domain_size);

            double r2 = dx * dx + dy * dy + dz * dz;
            if ((r2 < p_i.h * p_i.h || r2 < nj.h * nj.h) && r2 > 1e-24) {
                // Find where the face actually is
                double fraction_i = p_i.h / (p_i.h + nj.h);
                Eigen::Vector3d dx_face_i =
                    fraction_i * Eigen::Vector3d(dx, dy, dz);

                // Extrapolated deltas (phi_mid - phi_i)
                double d_mid_rho = out.grad_rho.dot(dx_face_i);
                phi_mid_max_rho = std::max(phi_mid_max_rho, d_mid_rho);
                phi_mid_min_rho = std::min(phi_mid_min_rho, d_mid_rho);

                double d_mid_p = out.grad_p.dot(dx_face_i);
                phi_mid_max_p = std::max(phi_mid_max_p, d_mid_p);
                phi_mid_min_p = std::min(phi_mid_min_p, d_mid_p);

                double d_mid_vx = out.grad_vx.dot(dx_face_i);
                phi_mid_max_vx = std::max(phi_mid_max_vx, d_mid_vx);
                phi_mid_min_vx = std::min(phi_mid_min_vx, d_mid_vx);

                double d_mid_vy = out.grad_vy.dot(dx_face_i);
                phi_mid_max_vy = std::max(phi_mid_max_vy, d_mid_vy);
                phi_mid_min_vy = std::min(phi_mid_min_vy, d_mid_vy);

                double d_mid_vz = out.grad_vz.dot(dx_face_i);
                phi_mid_max_vz = std::max(phi_mid_max_vz, d_mid_vz);
                phi_mid_min_vz = std::min(phi_mid_min_vz, d_mid_vz);
            }
        }

        // Tunable aggressiveness parameter beta (between 1.0 and 2.0)
        double beta = 1.5;

        // Compute alpha scalars and apply them to the gradients
        out.grad_rho *= compute_gradient_alpha(
            d_rho_max, d_rho_min, phi_mid_max_rho, phi_mid_min_rho, beta);
        out.grad_p *= compute_gradient_alpha(d_p_max, d_p_min, phi_mid_max_p,
                                             phi_mid_min_p, beta);
        out.grad_vx *= compute_gradient_alpha(
            d_vx_max, d_vx_min, phi_mid_max_vx, phi_mid_min_vx, beta);
        out.grad_vy *= compute_gradient_alpha(
            d_vy_max, d_vy_min, phi_mid_max_vy, phi_mid_min_vy, beta);
        out.grad_vz *= compute_gradient_alpha(
            d_vz_max, d_vz_min, phi_mid_max_vz, phi_mid_min_vz, beta);
#else
        // Apply the Scalar Limiter
        double alim = 0.5;  // 0.5 is standard for aggressive limiting
        double h_lim = std::max(p_i.h, r_max);
        double d_max = h_lim;
        double stol = 0.1;  // overshoot tolerance for pressure/velocity

        // Density: no overshoot tolerance (0.0), positivity preserving (true)
        scalar_limiter(out.grad_rho, d_rho_max, d_rho_min, alim, h_lim, 0.0,
                       true, d_max, p_i.rho);

        // Pressure: standard overshoot tolerance (stol), positivity preserving
        // (true)
        scalar_limiter(out.grad_p, d_p_max, d_p_min, alim, h_lim, stol, true,
                       d_max, p_i.pressure);

        // Velocity: standard overshoot tolerance (stol), NOT positivity
        // preserving (false)
        scalar_limiter(out.grad_vx, d_vx_max, d_vx_min, alim, h_lim, stol,
                       false, d_max, p_i.vel.x());
        scalar_limiter(out.grad_vy, d_vy_max, d_vy_min, alim, h_lim, stol,
                       false, d_max, p_i.vel.y());
        scalar_limiter(out.grad_vz, d_vz_max, d_vz_min, alim, h_lim, stol,
                       false, d_max, p_i.vel.z());
#endif
    } else {
        // Matrix is ill-conditioned (pathological alignment).
        // Fall back to standard SPH gradient estimator (Hopkins 2015, Eq. C4).
        out.B_matrix = Eigen::Matrix3d::Identity();  // Dummy valid matrix
        out.ill_conditioned = true;

        // Use a dimensionally correct isotropic average for the dummy B-matrix
        // to prevent Face Area explosions in the Riemann solver.
        /*double trace_E = E(0, 0) + E(1, 1) + E(2, 2);
        if (trace_E > 1e-24) {
            out.B_matrix = (3.0 / trace_E) * Eigen::Matrix3d::Identity();
        } else {
            // Absolute fallback if particle is completely isolated
            //double h2 = p_i.h * p_i.h;
            out.B_matrix = (1.0 / h2) * Eigen::Matrix3d::Identity();
        }*/

        for (const auto& nj : neighbors) {
            double dx = mfm_periodic_displacement(nj.pos.x() - p_i.pos.x(),
                                                  domain_size);
            double dy = mfm_periodic_displacement(nj.pos.y() - p_i.pos.y(),
                                                  domain_size);
            double dz = mfm_periodic_displacement(nj.pos.z() - p_i.pos.z(),
                                                  domain_size);
            double r2 = dx * dx + dy * dy + dz * dz;

            // SPH gradient sums
            if (r2 < p_i.h * p_i.h && r2 > 1e-24) {
                double r = std::sqrt(r2);

                // Fetch the derivative of the kernel respect to 'r' using the
                // existing gravity helper
                double dphi_dr_dummy, dW_dr;
                compute_adaptive_gravity_terms(r, p_i.h, dphi_dr_dummy, dW_dr);

                double V_j = nj.mass / nj.rho;

                // grad_i W = - (dW/dr) * (x_j - x_i) / r
                Eigen::Vector3d dx_vec(dx, dy, dz);
                Eigen::Vector3d grad_W_i = -dW_dr * dx_vec / r;

                // SPH gradient estimator: grad f_i = sum_j V_j (f_j - f_i)
                // grad_i W_ij
                out.grad_rho += V_j * (nj.rho - p_i.rho) * grad_W_i;
                out.grad_p += V_j * (nj.pressure - p_i.pressure) * grad_W_i;
                out.grad_vx += V_j * (nj.vel.x() - p_i.vel.x()) * grad_W_i;
                out.grad_vy += V_j * (nj.vel.y() - p_i.vel.y()) * grad_W_i;
                out.grad_vz += V_j * (nj.vel.z() - p_i.vel.z()) * grad_W_i;
            }
        }
    }

    return out;
}

void GasParticleSystem::compute_gradients(const Config& config) {
    if (num_particles == 0) return;
    double domain_size = config.domain_size;

    size_t step_ill_conditioned = 0;

#pragma omp parallel for reduction(+ : step_ill_conditioned) \
    schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        // Construct the isolated state for particle i
        ParticleState p_i;
        p_i.pos = Eigen::Vector3d(pos_x[i], pos_y[i], pos_z[i]);
        p_i.vel = Eigen::Vector3d(vel_x[i], vel_y[i], vel_z[i]);
        p_i.mass = mass[i];
        p_i.rho = rho[i];
        p_i.pressure = pressure[i];
        p_i.h = h[i];

        double local_max_rel_ke = 0.0;

        // Gather neighbors into a standard vector
        std::vector<ParticleState> neighbors;

        int search_cells = 1;

        int ix = static_cast<int>(p_i.pos.x() / hash_cell_size) % hash_grid_dim;
        int iy = static_cast<int>(p_i.pos.y() / hash_cell_size) % hash_grid_dim;
        int iz = static_cast<int>(p_i.pos.z() / hash_cell_size) % hash_grid_dim;
        ix = (ix + hash_grid_dim) % hash_grid_dim;
        iy = (iy + hash_grid_dim) % hash_grid_dim;
        iz = (iz + hash_grid_dim) % hash_grid_dim;

        for (int dx_c = -search_cells; dx_c <= search_cells; ++dx_c) {
            for (int dy_c = -search_cells; dy_c <= search_cells; ++dy_c) {
                for (int dz_c = -search_cells; dz_c <= search_cells; ++dz_c) {
                    int n_ix = (((ix + dx_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iy = (((iy + dy_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iz = (((iz + dz_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;

                    int cell_idx = n_iz * hash_grid_dim * hash_grid_dim +
                                   n_iy * hash_grid_dim + n_ix;
                    int start = sph_cell_list.cell_start[cell_idx];
                    int end = start + sph_cell_list.cell_count[cell_idx];

                    for (int j = start; j < end; ++j) {
                        ParticleState nj;
                        nj.pos = Eigen::Vector3d(pos_x[j], pos_y[j], pos_z[j]);
                        nj.vel = Eigen::Vector3d(vel_x[j], vel_y[j], vel_z[j]);
                        nj.mass = mass[j];
                        nj.rho = rho[j];
                        nj.pressure = pressure[j];
                        nj.h = h[j];
                        neighbors.push_back(nj);

                        // Calculate periodic distance
                        double dx = mfm_periodic_displacement(
                            pos_x[j] - p_i.pos.x(), domain_size);
                        double dy = mfm_periodic_displacement(
                            pos_y[j] - p_i.pos.y(), domain_size);
                        double dz = mfm_periodic_displacement(
                            pos_z[j] - p_i.pos.z(), domain_size);
                        double r2 = dx * dx + dy * dy + dz * dz;

                        // evaluate if they are actually interacting neighbors
                        if ((r2 < p_i.h * p_i.h || r2 < h[j] * h[j]) &&
                            r2 > 1e-24) {
                            double rel_vx = vel_x[j] - p_i.vel.x();
                            double rel_vy = vel_y[j] - p_i.vel.y();
                            double rel_vz = vel_z[j] - p_i.vel.z();
                            double rel_v2 = rel_vx * rel_vx + rel_vy * rel_vy +
                                            rel_vz * rel_vz;

                            local_max_rel_ke =
                                std::max(local_max_rel_ke, 0.5 * rel_v2);
                        }
                    }
                }
            }
        }

        max_rel_ke[i] = local_max_rel_ke;

        // Calculate delta E_grav = |a_grav| * h
        double a_grav_mag = std::sqrt(
            acc_x[i] * acc_x[i] + acc_y[i] * acc_y[i] + acc_z[i] * acc_z[i]);
        delta_E_grav[i] = a_grav_mag * p_i.h;

        ParticleGradients grads =
            compute_single_particle_gradients(p_i, neighbors, domain_size);

        // Map results back to SOA
        B_matrix[i] = grads.B_matrix;
        grad_rho[i] = grads.grad_rho;
        grad_p[i] = grads.grad_p;
        grad_vx[i] = grads.grad_vx;
        grad_vy[i] = grads.grad_vy;
        grad_vz[i] = grads.grad_vz;
        cond_num[i] = grads.condition_number;
        raw_sum_p[i] = grads.raw_sum_p;

        if (grads.ill_conditioned) {
            step_ill_conditioned++;
        }
    }

    ill_conditioned_cases += step_ill_conditioned;
}

// MFM Face Reconstruction
ReconstructedFace compute_face_reconstruction(const ParticleState& p_i,
                                              const ParticleGradients& grad_i,
                                              const ParticleState& p_j,
                                              const ParticleGradients& grad_j,
                                              double domain_size) {
    ReconstructedFace face;
    face.is_valid = false;

    double dx =
        mfm_periodic_displacement(p_j.pos.x() - p_i.pos.x(), domain_size);
    double dy =
        mfm_periodic_displacement(p_j.pos.y() - p_i.pos.y(), domain_size);
    double dz =
        mfm_periodic_displacement(p_j.pos.z() - p_i.pos.z(), domain_size);

    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < 1e-24) return face;

    face.r = std::sqrt(r2);
    Eigen::Vector3d dx_vec(dx, dy, dz);
    face.n = dx_vec / face.r;

    // Volume-weighted quadrature fraction
    double fraction_i = p_i.h / (p_i.h + p_j.h);
    double fraction_j = 1.0 - fraction_i;

    Eigen::Vector3d dx_face_i = fraction_i * dx_vec;
    Eigen::Vector3d dx_face_j = -fraction_j * dx_vec;

// #define DISABLE_LIMITER
#ifndef DISABLE_LIMITER
#ifdef THEORETICAL_LIMITER
    // Linearly interpolated "bar" values at the face
    double rho_bar = p_i.rho + fraction_i * (p_j.rho - p_i.rho);
    double p_bar = p_i.pressure + fraction_i * (p_j.pressure - p_i.pressure);
    double vx_bar = p_i.vel.x() + fraction_i * (p_j.vel.x() - p_i.vel.x());
    double vy_bar = p_i.vel.y() + fraction_i * (p_j.vel.y() - p_i.vel.y());
    double vz_bar = p_i.vel.z() + fraction_i * (p_j.vel.z() - p_i.vel.z());

    // Pairwise TVD Limiter Applied to the Reconstructed States
    face.rho_L = apply_pairwise_limiter(
        p_i.rho, p_j.rho, p_i.rho + grad_i.grad_rho.dot(dx_face_i), rho_bar);
    face.rho_R = apply_pairwise_limiter(
        p_j.rho, p_i.rho, p_j.rho + grad_j.grad_rho.dot(dx_face_j), rho_bar);

    face.p_L = apply_pairwise_limiter(
        p_i.pressure, p_j.pressure, p_i.pressure + grad_i.grad_p.dot(dx_face_i),
        p_bar);
    face.p_R = apply_pairwise_limiter(
        p_j.pressure, p_i.pressure, p_j.pressure + grad_j.grad_p.dot(dx_face_j),
        p_bar);

    face.v_L.x() = apply_pairwise_limiter(
        p_i.vel.x(), p_j.vel.x(), p_i.vel.x() + grad_i.grad_vx.dot(dx_face_i),
        vx_bar);
    face.v_R.x() = apply_pairwise_limiter(
        p_j.vel.x(), p_i.vel.x(), p_j.vel.x() + grad_j.grad_vx.dot(dx_face_j),
        vx_bar);

    face.v_L.y() = apply_pairwise_limiter(
        p_i.vel.y(), p_j.vel.y(), p_i.vel.y() + grad_i.grad_vy.dot(dx_face_i),
        vy_bar);
    face.v_R.y() = apply_pairwise_limiter(
        p_j.vel.y(), p_i.vel.y(), p_j.vel.y() + grad_j.grad_vy.dot(dx_face_j),
        vy_bar);

    face.v_L.z() = apply_pairwise_limiter(
        p_i.vel.z(), p_j.vel.z(), p_i.vel.z() + grad_i.grad_vz.dot(dx_face_i),
        vz_bar);
    face.v_R.z() = apply_pairwise_limiter(
        p_j.vel.z(), p_i.vel.z(), p_j.vel.z() + grad_j.grad_vz.dot(dx_face_j),
        vz_bar);

    // Enforce thermodynamic floors after the limiter
    face.rho_L = std::max(face.rho_L, density_floor);
    face.rho_R = std::max(face.rho_R, density_floor);
    face.p_L = std::max(face.p_L, g_pressure_floor);
    face.p_R = std::max(face.p_R, g_pressure_floor);

#else
    // 2nd-order reconstruction using the scalar-limited gradients
    face.rho_L =
        std::max(p_i.rho + grad_i.grad_rho.dot(dx_face_i), density_floor);
    face.rho_R =
        std::max(p_j.rho + grad_j.grad_rho.dot(dx_face_j), density_floor);

    face.p_L =
        std::max(p_i.pressure + grad_i.grad_p.dot(dx_face_i), g_pressure_floor);
    face.p_R =
        std::max(p_j.pressure + grad_j.grad_p.dot(dx_face_j), g_pressure_floor);

    face.v_L.x() = p_i.vel.x() + grad_i.grad_vx.dot(dx_face_i);
    face.v_L.y() = p_i.vel.y() + grad_i.grad_vy.dot(dx_face_i);
    face.v_L.z() = p_i.vel.z() + grad_i.grad_vz.dot(dx_face_i);

    face.v_R.x() = p_j.vel.x() + grad_j.grad_vx.dot(dx_face_j);
    face.v_R.y() = p_j.vel.y() + grad_j.grad_vy.dot(dx_face_j);
    face.v_R.z() = p_j.vel.z() + grad_j.grad_vz.dot(dx_face_j);
#endif
#else
    // TEMPORARY DIAGNOSTIC: PURE 2ND-ORDER RECONSTRUCTION (NO LIMITERS)

    // Density (with absolute physical floor to prevent NaN)
    face.rho_L =
        std::max(p_i.rho + grad_i.grad_rho.dot(dx_face_i), density_floor);
    face.rho_R =
        std::max(p_j.rho + grad_j.grad_rho.dot(dx_face_j), density_floor);

    // Pressure (with absolute physical floor to prevent NaN)
    face.p_L =
        std::max(p_i.pressure + grad_i.grad_p.dot(dx_face_i), g_pressure_floor);
    face.p_R =
        std::max(p_j.pressure + grad_j.grad_p.dot(dx_face_j), g_pressure_floor);

    // Velocity (completely unbounded)
    face.v_L.x() = p_i.vel.x() + grad_i.grad_vx.dot(dx_face_i);
    face.v_L.y() = p_i.vel.y() + grad_i.grad_vy.dot(dx_face_i);
    face.v_L.z() = p_i.vel.z() + grad_i.grad_vz.dot(dx_face_i);

    face.v_R.x() = p_j.vel.x() + grad_j.grad_vx.dot(dx_face_j);
    face.v_R.y() = p_j.vel.y() + grad_j.grad_vy.dot(dx_face_j);
    face.v_R.z() = p_j.vel.z() + grad_j.grad_vz.dot(dx_face_j);
#endif

    face.is_valid = true;
    return face;
}

// MFM Riemann Solver (Frame Boosted)
MFMFaceFlux solve_mfm_riemann(const ReconstructedFace& face,
                              const Eigen::Vector3d& v_frame, double gamma) {
    // Boost and Project Left & Right States to the face frame
    double vn_L = (face.v_L - v_frame).dot(face.n);
    double vn_R = (face.v_R - v_frame).dot(face.n);

    // Compute Sound Speeds
    double cs_L =
        (face.rho_L > 1e-12) ? std::sqrt(gamma * face.p_L / face.rho_L) : 0.0;
    double cs_R =
        (face.rho_R > 1e-12) ? std::sqrt(gamma * face.p_R / face.rho_R) : 0.0;

    // Wave Speed Estimates
    double S_L = std::min(vn_L - cs_L, vn_R - cs_R);
    double S_R = std::max(vn_L + cs_L, vn_R + cs_R);

    double P_star = 0.0;
    double S_star = 0.0;

    // Solve the HLLC Riemann Problem for the Star State
    if (S_L >= 0.0) {
        P_star = face.p_L;
        S_star = vn_L;
    } else if (S_R <= 0.0) {
        P_star = face.p_R;
        S_star = vn_R;
    } else {
        // Subsonic Contact Wave
        double den_star = face.rho_L * (S_L - vn_L) - face.rho_R * (S_R - vn_R);
        if (std::abs(den_star) < 1e-15)
            den_star = (den_star > 0 ? 1e-15 : -1e-15);

        S_star = (face.p_R - face.p_L + face.rho_L * vn_L * (S_L - vn_L) -
                  face.rho_R * vn_R * (S_R - vn_R)) /
                 den_star;
        P_star = face.p_L + face.rho_L * (S_L - vn_L) * (S_star - vn_L);

        if (P_star < 0.0) P_star = 0.0;  // Floor to physical values
    }

    // Assemble MFM Fluxes (Mass flux is zero)
    MFMFaceFlux out;
    out.flux_mom = P_star * face.n;
    out.P_star = P_star;
    out.S_star = S_star;  // Return the relative face speed

    return out;
}

Eigen::Vector3d compute_mfm_face_area_vector(const ParticleState& p_i,
                                             const Eigen::Matrix3d& B_i,
                                             const ParticleState& p_j,
                                             const Eigen::Matrix3d& B_j,
                                             const ReconstructedFace& face) {
    double V_i = p_i.mass / p_i.rho;
    double V_j = p_j.mass / p_j.rho;

    double W_i, W_j;
    compute_kernel(face.r, p_i.h, W_i);
    compute_kernel(face.r, p_j.h, W_j);

    Eigen::Vector3d dx_vec = face.n * face.r;
    Eigen::Vector3d Area_vec =
        (V_i * V_j * W_i * (B_i * dx_vec)) + (V_i * V_j * W_j * (B_j * dx_vec));

    return Area_vec;
}

void GasParticleSystem::compute_hydro_forces(const Config& config, double a,
                                             double dt) {
    if (num_particles == 0) return;

    double domain_size = config.domain_size;

    // Reset accumulators
    std::fill(hydro_acc_x.begin(), hydro_acc_x.end(), 0.0);
    std::fill(hydro_acc_y.begin(), hydro_acc_y.end(), 0.0);
    std::fill(hydro_acc_z.begin(), hydro_acc_z.end(), 0.0);
    std::fill(du_dt.begin(), du_dt.end(), 0.0);
    std::fill(de_dt.begin(), de_dt.end(), 0.0);

    // Search radius is 1 cell (max_h) to catch all pairwise
    // interactions
    int search_cells = 1;

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        // Build the isolated state proxy for particle i
        ParticleState p_i;
        p_i.pos = Eigen::Vector3d(pos_x[i], pos_y[i], pos_z[i]);
        p_i.vel = Eigen::Vector3d(vel_x[i], vel_y[i], vel_z[i]);
        p_i.mass = mass[i];
        p_i.rho = rho[i];
        p_i.pressure = pressure[i];
        p_i.h = h[i];

        ParticleGradients grad_i;
        grad_i.B_matrix = B_matrix[i];
        grad_i.grad_rho = grad_rho[i];
        grad_i.grad_p = grad_p[i];
        grad_i.grad_vx = grad_vx[i];
        grad_i.grad_vy = grad_vy[i];
        grad_i.grad_vz = grad_vz[i];

        int ix = static_cast<int>(p_i.pos.x() / hash_cell_size) % hash_grid_dim;
        int iy = static_cast<int>(p_i.pos.y() / hash_cell_size) % hash_grid_dim;
        int iz = static_cast<int>(p_i.pos.z() / hash_cell_size) % hash_grid_dim;
        ix = (ix + hash_grid_dim) % hash_grid_dim;
        iy = (iy + hash_grid_dim) % hash_grid_dim;
        iz = (iz + hash_grid_dim) % hash_grid_dim;

        for (int dx_c = -search_cells; dx_c <= search_cells; ++dx_c) {
            for (int dy_c = -search_cells; dy_c <= search_cells; ++dy_c) {
                for (int dz_c = -search_cells; dz_c <= search_cells; ++dz_c) {
                    int n_ix = (((ix + dx_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iy = (((iy + dy_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;
                    int n_iz = (((iz + dz_c) % hash_grid_dim) + hash_grid_dim) %
                               hash_grid_dim;

                    int cell_idx = n_iz * hash_grid_dim * hash_grid_dim +
                                   n_iy * hash_grid_dim + n_ix;
                    int start = sph_cell_list.cell_start[cell_idx];
                    int end = start + sph_cell_list.cell_count[cell_idx];

                    for (int j = start; j < end; ++j) {
                        // Process each pair just once
                        if (j <= i) continue;

                        double dx = mfm_periodic_displacement(
                            pos_x[j] - p_i.pos.x(), domain_size);
                        double dy = mfm_periodic_displacement(
                            pos_y[j] - p_i.pos.y(), domain_size);
                        double dz = mfm_periodic_displacement(
                            pos_z[j] - p_i.pos.z(), domain_size);

                        double r2 = dx * dx + dy * dy + dz * dz;
                        double hj = h[j];

                        // Skip overlaps to prevent division by zero
                        if (r2 < 1e-24) continue;

                        // Particles interact if they fall within either
                        // smoothing length
                        if (r2 < p_i.h * p_i.h || r2 < hj * hj) {
                            // Build the isolated state proxy for particle j
                            ParticleState p_j;
                            p_j.pos =
                                Eigen::Vector3d(pos_x[j], pos_y[j], pos_z[j]);
                            p_j.vel =
                                Eigen::Vector3d(vel_x[j], vel_y[j], vel_z[j]);
                            p_j.mass = mass[j];
                            p_j.rho = rho[j];
                            p_j.pressure = pressure[j];
                            p_j.h = hj;

                            ParticleGradients grad_j;
                            grad_j.B_matrix = B_matrix[j];
                            grad_j.grad_rho = grad_rho[j];
                            grad_j.grad_p = grad_p[j];
                            grad_j.grad_vx = grad_vx[j];
                            grad_j.grad_vy = grad_vy[j];
                            grad_j.grad_vz = grad_vz[j];

                            // SOLVER PIPELINE

                            // Face Reconstruction (Extrapolates states to
                            // the midpoint)
                            ReconstructedFace face =
                                compute_face_reconstruction(
                                    p_i, grad_i, p_j, grad_j, domain_size);
                            if (!face.is_valid) continue;

                            // Compute the Area Vector
                            Eigen::Vector3d Area_vec =
                                compute_mfm_face_area_vector(
                                    p_i, grad_i.B_matrix, p_j, grad_j.B_matrix,
                                    face);

                            // Check if the area is zero to prevent division
                            // by zero
                            double A_mag = Area_vec.norm();
                            if (A_mag < 1e-20) continue;

                            // Re-align the face normal to
                            // the Area Vector. The Riemann problem MUST be
                            // solved along the A_ij axis.
                            face.n = Area_vec / A_mag;

                            // Frame Boosting & Riemann Solver
                            double fraction_i = p_i.h / (p_i.h + p_j.h);
                            double fraction_j = 1.0 - fraction_i;
                            Eigen::Vector3d v_frame =
                                p_i.vel + fraction_i * (p_j.vel - p_i.vel);

                            // Because face.n is Area_vec.normalized(), vn_L
                            // and vn_R inside this function will be
                            // correctly projected onto the true face.
                            MFMFaceFlux flux_1d =
                                solve_mfm_riemann(face, v_frame, config.gamma);

                            // Apply Vector Fluxes
                            Eigen::Vector3d Force_mom =
                                flux_1d.P_star * Area_vec;
                            Eigen::Vector3d v_star_lab =
                                v_frame + (flux_1d.S_star * face.n);

                            double work_i =
                                flux_1d.P_star *
                                (v_star_lab - p_i.vel).dot(Area_vec);
                            double work_j =
                                flux_1d.P_star *
                                (v_star_lab - p_j.vel).dot(Area_vec);

                            double du_dt_i = -work_i / p_i.mass;
                            double du_dt_j = work_j / p_j.mass;

                            double Rate_energy =
                                flux_1d.P_star * v_star_lab.dot(Area_vec);
                            double de_dt_i = -Rate_energy / p_i.mass;
                            double de_dt_j = Rate_energy / p_j.mass;

                            if (config.disable_hydro_forces) {
                                continue;
                            }

#pragma omp atomic
                            hydro_acc_x[i] -= Force_mom.x() / p_i.mass;
#pragma omp atomic
                            hydro_acc_y[i] -= Force_mom.y() / p_i.mass;
#pragma omp atomic
                            hydro_acc_z[i] -= Force_mom.z() / p_i.mass;
#pragma omp atomic
                            du_dt[i] += du_dt_i;
#pragma omp atomic
                            de_dt[i] += de_dt_i;

#pragma omp atomic
                            hydro_acc_x[j] += Force_mom.x() / p_j.mass;
#pragma omp atomic
                            hydro_acc_y[j] += Force_mom.y() / p_j.mass;
#pragma omp atomic
                            hydro_acc_z[j] += Force_mom.z() / p_j.mass;
#pragma omp atomic
                            du_dt[j] += du_dt_j;
#pragma omp atomic
                            de_dt[j] += de_dt_j;
                        }
                    }
                }
            }
        }
    }
}