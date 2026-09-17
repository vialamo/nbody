#include "mfm.h"

#include <chrono>

#include "cic.h"
#include "constants.h"
#include "diagnostics.h"
#include "kernels.h"
#include "math_utils.h"

// #define UNCORRECTED_GRAVITY

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

    grad_rho.reserve(config.num_gas_particles);
    grad_vx.reserve(config.num_gas_particles);
    grad_vy.reserve(config.num_gas_particles);
    grad_vz.reserve(config.num_gas_particles);
    grad_p.reserve(config.num_gas_particles);
    B_matrix.reserve(config.num_gas_particles);

    entropy.reserve(config.num_gas_particles);
    max_rel_ke.reserve(config.num_gas_particles);
    delta_E_grav.reserve(config.num_gas_particles);

    zeta.reserve(config.num_gas_particles);

    cond_num.reserve(config.num_gas_particles);
    raw_sum_p.reserve(config.num_gas_particles);
    n_enc_final.reserve(config.num_gas_particles);

    size_t num_nodes = 2 * config.num_gas_particles - 1;
    if (config.num_dm_particles > 0) {
        morton_codes.reserve(config.num_gas_particles);
        sorted_indices.reserve(config.num_gas_particles);
        bvh_nodes.reserve(num_nodes);
    }

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
    n_enc_final.push_back(0);

    num_particles++;
}

void GasParticleSystem::sort_arrays(const std::vector<int>& sorted_indices) {
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
    std::vector<double> new_n_enc_final(num_particles);

    for (size_t i = 0; i < num_particles; ++i) {
        int src = sorted_indices[i];
        new_px[i] = pos_x[src];
        new_py[i] = pos_y[src];
        new_pz[i] = pos_z[src];
        new_vx[i] = vel_x[src];
        new_vy[i] = vel_y[src];
        new_vz[i] = vel_z[src];
        new_ax[i] = acc_x[src];
        new_ay[i] = acc_y[src];
        new_az[i] = acc_z[src];
        new_m[i] = mass[src];
        new_hax[i] = hydro_acc_x[src];
        new_hay[i] = hydro_acc_y[src];
        new_haz[i] = hydro_acc_z[src];
        new_h[i] = h[src];
        new_rho[i] = rho[src];
        new_p[i] = pressure[src];
        new_total_energy[i] = total_energy[src];
        new_u[i] = u[src];
        new_dudt[i] = du_dt[src];
        new_dedt[i] = de_dt[src];
        new_metal[i] = metal_frac[src];
        new_cic[i] = cic_data[src];
        new_B[i] = B_matrix[src];
        new_grad_rho[i] = grad_rho[src];
        new_grad_vx[i] = grad_vx[src];
        new_grad_vy[i] = grad_vy[src];
        new_grad_vz[i] = grad_vz[src];
        new_grad_p[i] = grad_p[src];
        new_zeta[i] = zeta[src];
        new_entropy[i] = entropy[src];
        new_max_rel_ke[i] = max_rel_ke[src];
        new_delta_E_grav[i] = delta_E_grav[src];
        new_cond_num[i] = cond_num[src];
        new_raw_sum_p[i] = raw_sum_p[src];
        new_n_enc_final[i] = n_enc_final[src];
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
    n_enc_final = std::move(new_n_enc_final);
}

void GasParticleSystem::build_lbvh(const Config& config) {
    if (num_particles == 0) return;

    // Generate the Morton codes and the index map
    LBVH::compute_morton_and_sort_indices(num_particles, config.domain_size,
                                          pos_x, pos_y, pos_z, morton_codes,
                                          sorted_indices);

    // Shuffle arrays
    sort_arrays(sorted_indices);

    // Build the tree topology and compute Centers of Mass / Bounding Boxes
    LBVH::build_topology_and_aggregate(num_particles, pos_x, pos_y, pos_z, mass,
                                       &h, morton_codes, bvh_nodes);
}

void GasParticleSystem::evaluate_density_sum(size_t particle_idx,
                                             double h_guess, double domain_size,
                                             double& out_n,
                                             double& out_dn_dh) const {
    out_n = 0.0;
    out_dn_dh = 0.0;
    double search_sq = h_guess * h_guess;
    double p1_x = pos_x[particle_idx];
    double p1_y = pos_y[particle_idx];
    double p1_z = pos_z[particle_idx];

    int stack[128];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    while (stack_ptr > 0) {
        int node_idx = stack[--stack_ptr];
        const BVHNode& node = bvh_nodes[node_idx];

        double dist_sq = min_periodic_dist_sq(p1_x, node.bbox.min_x,
                                              node.bbox.max_x, domain_size) +
                         min_periodic_dist_sq(p1_y, node.bbox.min_y,
                                              node.bbox.max_y, domain_size) +
                         min_periodic_dist_sq(p1_z, node.bbox.min_z,
                                              node.bbox.max_z, domain_size);

        if (dist_sq > search_sq) continue;

        if (node.particle_idx != -1) {
            int j = node.particle_idx;
            double dx = periodic_displacement(pos_x[j] - p1_x, domain_size);
            double dy = periodic_displacement(pos_y[j] - p1_y, domain_size);
            double dz = periodic_displacement(pos_z[j] - p1_z, domain_size);
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 < search_sq) {
                double r = std::sqrt(r2);
                double W, dWdh;
                Kernels::cubic_spline(r, h_guess, W, dWdh);
                out_n += W;
                out_dn_dh += dWdh;
            }
        } else {
            stack[stack_ptr++] = node.left_child;
            stack[stack_ptr++] = node.right_child;
        }
    }
}

double GasParticleSystem::compute_zeta_contribution(
    double p1_x, double p1_y, double p1_z, double h_i,
    const std::vector<double>& target_x, const std::vector<double>& target_y,
    const std::vector<double>& target_z, const std::vector<double>& target_mass,
    const std::vector<BVHNode>& target_bvh, const Config& config) const {
    double zeta_sum = 0.0;
    double search_sq_zeta = h_i * h_i;
    double domain_size = config.domain_size;

    bool use_pm = config.use_PM;
    double r_s = config.PM_smoothing_cells * config.cell_size;

    int stack[128];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    while (stack_ptr > 0) {
        int node_idx = stack[--stack_ptr];
        const BVHNode& node = target_bvh[node_idx];

        double dist_sq = min_periodic_dist_sq(p1_x, node.bbox.min_x,
                                              node.bbox.max_x, domain_size) +
                         min_periodic_dist_sq(p1_y, node.bbox.min_y,
                                              node.bbox.max_y, domain_size) +
                         min_periodic_dist_sq(p1_z, node.bbox.min_z,
                                              node.bbox.max_z, domain_size);

        if (dist_sq > search_sq_zeta) continue;

        if (node.particle_idx != -1) {
            int j = node.particle_idx;
            double dx = periodic_displacement(target_x[j] - p1_x, domain_size);
            double dy = periodic_displacement(target_y[j] - p1_y, domain_size);
            double dz = periodic_displacement(target_z[j] - p1_z, domain_size);
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 < search_sq_zeta) {
                if (use_pm && r2 > config.cutoff_radius_squared) continue;

                double r = std::sqrt(r2);
                double dphi_dr, dphi_dh;
                Kernels::gravity_derivatives(r, h_i, dphi_dr, dphi_dh);

                if (use_pm) {
                    double r_scaled = r / (2.0 * r_s);
                    double taper = std::erfc(r_scaled) +
                                   (r / (std::sqrt(M_PI) * r_s)) *
                                       std::exp(-r_scaled * r_scaled);
                    dphi_dh *= taper;
                }
                zeta_sum += target_mass[j] * dphi_dh;
            }
        } else {
            stack[stack_ptr++] = node.left_child;
            stack[stack_ptr++] = node.right_child;
        }
    }
    return zeta_sum;
}

// --------------------------------------------------------------------------------
// Density and Smoothing Length Iteration
// --------------------------------------------------------------------------------
void GasParticleSystem::compute_density_and_h(const Config& config,
                                              const ParticleSystem& dm) {
    if (num_particles == 0) return;

    build_lbvh(config);

    double domain_size = config.domain_size;
    double target_N = config.mfm_target_neighbors;
    double tol = config.mfm_neighbor_tolerance;
    int max_iter = config.mfm_max_iterations;
    const double mean_spacing =
        domain_size / std::cbrt(num_particles > 0 ? num_particles : 1);
    const double min_h_cap = 0.05 * mean_spacing;
    const double max_h_cap = 0.5 * domain_size;
    constexpr double MAX_H_GROWTH = 8.0;

    size_t num_h_clamped = 0;

#pragma omp parallel for schedule(dynamic, 64) reduction(+ : num_h_clamped)
    for (size_t i = 0; i < num_particles; ++i) {
        double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
        double h_low = 0.0;
        double h_high = std::numeric_limits<double>::infinity();
        double h_guess = h[i];
        double step_max_h = std::min(max_h_cap, h[i] * MAX_H_GROWTH);

        int iter = 0;
        double current_n = 0.0;
        double current_dn_dh = 0.0;
        bool is_converged = false;
        bool h_clamped = false;

        // Newton-Raphson Solver for h_i
        while (iter < max_iter) {
            evaluate_density_sum(i, h_guess, domain_size, current_n,
                                 current_dn_dh);

            double h3 = h_guess * h_guess * h_guess;
            double N_enc = (4.0 / 3.0) * M_PI * h3 * current_n;
            n_enc_final[i] = N_enc;

            if (std::abs(N_enc - target_N) < tol) {
                is_converged = true;
                break;
            }

            if (N_enc > target_N)
                h_high = h_guess;
            else
                h_low = h_guess;

            double dN_enc_dh =
                (4.0 / 3.0) * M_PI *
                (3.0 * h_guess * h_guess * current_n + h3 * current_dn_dh);
            double h_new = h_guess;

            if (dN_enc_dh > 0.0)
                h_new = h_guess - (N_enc - target_N) / dN_enc_dh;

            if (h_new <= h_low || h_new >= h_high || dN_enc_dh <= 0.0) {
                h_guess = std::isinf(h_high) ? (1.26 * h_guess)
                                             : 0.5 * (h_low + h_high);
            } else {
                h_guess = h_new;
            }

            if (h_guess >= step_max_h) {
                h_guess = step_max_h;
                h_clamped = true;
                break;
            }
            if (h_guess < min_h_cap) {
                h_guess = min_h_cap;
                h_clamped = true;
                break;
            }
            iter++;
        }

        if (h_clamped) {
            num_h_clamped++;
        }

        // IMPORTANT: If we exited the loop due to max_iter OR clamping,
        // h_guess has been updated but current_n and current_dn_dh are stale.
        // We must re-evaluate
        if (!is_converged) {
            // Sync the properties using the finalized, clamped h_guess
            evaluate_density_sum(i, h_guess, domain_size, current_n,
                                 current_dn_dh);
            double h3 = h_guess * h_guess * h_guess;
            n_enc_final[i] = (4.0 / 3.0) * M_PI * h3 * current_n;
        }

        // Commit finalized state
        h[i] = h_guess;
        rho[i] = mass[i] * current_n;

#ifndef UNCORRECTED_GRAVITY
        // Adaptive gravity correction (Zeta)
        double Omega_i = std::max(
            1.0 + (h_guess / (current_n * 3.0)) * current_dn_dh, 1e-12);

        // Sum Gas contribution
        double zeta_sum =
            compute_zeta_contribution(p1_x, p1_y, p1_z, h_guess, pos_x, pos_y,
                                      pos_z, mass, bvh_nodes, config);

        // Sum DM contribution
        if (dm.num_particles > 0 && !dm.bvh_nodes.empty()) {
            zeta_sum += compute_zeta_contribution(
                p1_x, p1_y, p1_z, h_guess, dm.pos_x, dm.pos_y, dm.pos_z,
                dm.mass, dm.bvh_nodes, config);
        }

        zeta[i] = (h_guess / (current_n * 3.0)) * (1.0 / Omega_i) * zeta_sum;
#endif
    }

    clamped_h_cases += num_h_clamped;
}

void GasParticleSystem::bin_and_assign_mass(const Config& config) {
    CIC::bin_and_assign_mass(config, num_particles, pos_x, pos_y, pos_z, mass,
                             cic_data, gas_rho);
}

void GasParticleSystem::interpolate_cic_forces(const Grid3D& ax_grid,
                                               const Grid3D& ay_grid,
                                               const Grid3D& az_grid,
                                               const Config& config) {
    CIC::interpolate_forces(config, num_particles, cic_data, ax_grid, ay_grid,
                            az_grid, acc_x, acc_y, acc_z);
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
            // Track the metal mass fraction
            double local_Z_frac = metal_frac[i];
            double u_current = u[i];
            double u_initial = u_current;
            double t_evolved = 0.0;
            int cell_non_converged = 0;

            // Local Particle Subcycling
            while (t_evolved < dt) {
                double du_dt = cooling.compute_du_dt(u_current, local_rho,
                                                     local_Z_frac, a, config);

                double dt_cell;
                if (u_current <= u_rad_floor && du_dt < 0.0) {
                    // The particle is at the temperature floor and trying to
                    // cool. It is in thermal equilibrium. Consume the rest of
                    // the step
                    dt_cell = dt - t_evolved;
                } else {
                    dt_cell = (std::abs(du_dt) > 0.0)
                                  ? 0.1 * (u_current / std::abs(du_dt))
                                  : dt;
                }

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

    // Resynchronize pressure and dual-energy arrays using the ie
    update_primitive_variables(config, a);
}

double GasParticleSystem::get_cfl_timestep(double a,
                                           const Config& config) const {
    if (config.hydro_method != HydroMethod::MFM || num_particles == 0) {
        return std::numeric_limits<double>::infinity();
    }

    double min_dt = std::numeric_limits<double>::infinity();
    double gamma = config.gamma;
    double domain_size = config.domain_size;

    // Precompute cosmology conversion factors
    double a_inv = 1.0 / a;
    double a_inv3 = a_inv * a_inv * a_inv;

#pragma omp parallel for reduction(min : min_dt) schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];

        // Evaluate particle i in PHYSICAL units
        double rho_phys_i = rho[i] * a_inv3;
        double p_phys_i = pressure[i] * a_inv;
        double c_phys_i = (rho_phys_i > 1e-12)
                              ? std::sqrt(gamma * p_phys_i / rho_phys_i)
                              : 0.0;

        double v1_x_phys = vel_x[i] * a;
        double v1_y_phys = vel_y[i] * a;
        double v1_z_phys = vel_z[i] * a;

        double h_i = h[i];
        double v_sig_max_phys = 0.0;

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;  // Push root

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = bvh_nodes[node_idx];

            double dist_sq =
                min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                     domain_size) +
                min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                     domain_size) +
                min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                     domain_size);

            double eff_h = std::max(h_i, node.max_h);
            if (dist_sq > eff_h * eff_h) continue;

            if (node.particle_idx != -1) {
                int j = node.particle_idx;
                if (i == j) continue;

                double dx = periodic_displacement(pos_x[j] - p1_x, domain_size);
                double dy = periodic_displacement(pos_y[j] - p1_y, domain_size);
                double dz = periodic_displacement(pos_z[j] - p1_z, domain_size);
                double r2 = dx * dx + dy * dy + dz * dz;

                if (r2 < h_i * h_i || r2 < h[j] * h[j]) {
                    if (r2 > 1e-24) {
                        double r = std::sqrt(r2);

                        // Evaluate particle j in PHYSICAL units
                        double rho_phys_j = rho[j] * a_inv3;
                        double p_phys_j = pressure[j] * a_inv;
                        double c_phys_j =
                            (rho_phys_j > 1e-12)
                                ? std::sqrt(gamma * p_phys_j / rho_phys_j)
                                : 0.0;

                        // Physical relative velocity
                        double dvx_phys = (vel_x[j] * a) - v1_x_phys;
                        double dvy_phys = (vel_y[j] * a) - v1_y_phys;
                        double dvz_phys = (vel_z[j] * a) - v1_z_phys;

                        double dv_dot_dx_phys =
                            dvx_phys * dx + dvy_phys * dy + dvz_phys * dz;
                        double v_sig_ij_phys = c_phys_i + c_phys_j;

                        if (dv_dot_dx_phys < 0.0) {
                            v_sig_ij_phys -= dv_dot_dx_phys / r;
                        }

                        if (v_sig_ij_phys > v_sig_max_phys) {
                            v_sig_max_phys = v_sig_ij_phys;
                        }
                    }
                }
            } else {
                stack[stack_ptr++] = node.left_child;
                stack[stack_ptr++] = node.right_child;
            }
        }

        if (v_sig_max_phys > 1e-9) {
            double dt_i = (a * h_i) / v_sig_max_phys;
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

void GasParticleSystem::compute_and_add_pp_forces(double a,
                                                  const Config& config,
                                                  Diagnostics& diag) {
    if (num_particles == 0) return;

    const double domain_size = config.domain_size;
    // Pre-scale G so all output forces are comoving accelerations
    const double G = config.G / (a * a * a);
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t n_gas = num_particles;
    const size_t num_nodes = 2 * n_gas - 1;

#ifdef UNCORRECTED_GRAVITY
    const double soft_sq = config.softening_squared;
    const double base_soft = std::sqrt(soft_sq);
#endif

    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();

    // Extract pointers
    double* d_px = pos_x.data();
    double* d_py = pos_y.data();
    double* d_pz = pos_z.data();
    double* d_m = mass.data();
    double* d_h = h.data();
    double* d_zeta = zeta.data();
    double* d_ax = acc_x.data();
    double* d_ay = acc_y.data();
    double* d_az = acc_z.data();
    BVHNode* d_bvh_nodes = bvh_nodes.data();

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                        \
        to : d_px[0 : n_gas], d_py[0 : n_gas], d_pz[0 : n_gas],           \
            d_m[0 : n_gas], d_h[0 : n_gas], d_zeta[0 : n_gas],            \
            d_bvh_nodes[0 : num_nodes], d_ax[0 : n_gas], d_ay[0 : n_gas], \
            d_az[0 : n_gas])

    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

// Swap the OpenMP execution pragma based on the compile-time target
#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t i = 0; i < n_gas; ++i) {
        double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
        double m_i = d_m[i];
        double h_i = d_h[i];
        double zeta_i = d_zeta[i];

        double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;  // Push root node

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = d_bvh_nodes[node_idx];

            // AABB Culling against cutoff radius
            double aabb_dist_sq =
                min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                     domain_size) +
                min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                     domain_size) +
                min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                     domain_size);

            if (aabb_dist_sq > search_sq) continue;

            if (node.particle_idx != -1) {
                int j = node.particle_idx;
                if (static_cast<size_t>(j) == i) continue;

                double dx = p1_x - d_px[j];
                if (dx > 0.5 * domain_size)
                    dx -= domain_size;
                else if (dx < -0.5 * domain_size)
                    dx += domain_size;

                double dy = p1_y - d_py[j];
                if (dy > 0.5 * domain_size)
                    dy -= domain_size;
                else if (dy < -0.5 * domain_size)
                    dy += domain_size;

                double dz = p1_z - d_pz[j];
                if (dz > 0.5 * domain_size)
                    dz -= domain_size;
                else if (dz < -0.5 * domain_size)
                    dz += domain_size;

                dx = -dx;
                dy = -dy;
                dz = -dz;

                double dist_sq = dx * dx + dy * dy + dz * dz;

                if (use_pm && dist_sq > cutoff_sq) continue;
                if (dist_sq < 1e-24) continue;

                double r = std::sqrt(dist_sq);
                double m_j = d_m[j];

#ifdef UNCORRECTED_GRAVITY
                // Pure Plummer Softening (Matches DM exactly, bypasses Kernels)
                double pp_dist_sq = dist_sq + soft_sq;
                double pp_dist = std::sqrt(pp_dist_sq);

                // We define force_mag_over_r so that when it is later
                // multiplied by (m_j * dx), it equals the DM equation: (G * m_j
                // * dx) / (pp_dist^3)
                double force_mag_over_r = G / (pp_dist_sq * pp_dist);

                // Override 'r' to match the DM's PM taper behavior perfectly
                // r = pp_dist;
#else
                double h_j = d_h[j];
                double zeta_j = d_zeta[j];

                double dphi_dr_i, dW_dr_i, dphi_dr_j, dW_dr_j;
                Kernels::adaptive_gravity_terms(r, h_i, dphi_dr_i, dW_dr_i);
                Kernels::adaptive_gravity_terms(r, h_j, dphi_dr_j, dW_dr_j);

                double force_mag_over_r =
                    (G / 2.0) *
                    ((dphi_dr_i + dphi_dr_j) + (zeta_i * dW_dr_i) / m_i +
                     (zeta_j * dW_dr_j) / m_j) /
                    r;
#endif

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
            } else {
                stack[stack_ptr++] = node.left_child;
                stack[stack_ptr++] = node.right_child;
            }
        }

        d_ax[i] += local_acc_x;
        d_ay[i] += local_acc_y;
        d_az[i] += local_acc_z;
    }

#ifdef USE_GPU
    auto end_compute = std::chrono::high_resolution_clock::now();
    auto start_return = std::chrono::high_resolution_clock::now();

#pragma omp target exit data map(from : d_ax[0 : n_gas], d_ay[0 : n_gas], \
                                     d_az[0 : n_gas])                     \
    map(delete : d_px[0 : n_gas], d_py[0 : n_gas], d_pz[0 : n_gas],       \
            d_m[0 : n_gas], d_h[0 : n_gas], d_zeta[0 : n_gas],            \
            d_bvh_nodes[0 : num_nodes])

    auto end_return = std::chrono::high_resolution_clock::now();

    diag.add_prof_time(ProfRegion::Transf,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_transfer - start_transfer)
                           .count());
    diag.add_prof_time(ProfRegion::Compute,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_compute - start_compute)
                           .count());
    diag.add_prof_time(ProfRegion::Ret,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_return - start_return)
                           .count());
#endif
}

void GasParticleSystem::compute_cross_pp_forces(double a, ParticleSystem& dm,
                                                const Config& config,
                                                Diagnostics& diag) {
    if (num_particles == 0 || dm.num_particles == 0) return;

    const double domain_size = config.domain_size;
    // Pre-scale G so all output forces are comoving accelerations
    const double G = config.G / (a * a * a);
    const double soft_sq = config.softening_squared;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t n_gas = num_particles;
    const size_t n_dm = dm.num_particles;
    const size_t dm_num_nodes = 2 * n_dm - 1;

    // Search radius is the PM cutoff. If not using PM, it's infinite (N^2
    // tree traversal)
    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();
    const double base_soft = std::sqrt(soft_sq);
    const double spline_equivalent_h = 2.8 * base_soft;

    // Extract pointers
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

    BVHNode* d_dm_bvh = dm.bvh_nodes.data();

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                           \
        to : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas],  \
            d_gas_m[0 : n_gas], d_gas_h[0 : n_gas], d_gas_zeta[0 : n_gas],   \
            d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas],   \
            d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],         \
            d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes], d_dm_ax[0 : n_dm], \
            d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm])

    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

// Swap the OpenMP execution pragma based on the compile-time target
#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t i = 0; i < n_gas; ++i) {
        double p1_x = d_gas_px[i], p1_y = d_gas_py[i], p1_z = d_gas_pz[i];
        double m_gas = d_gas_m[i];
        double h_i = d_gas_h[i];
        double zeta_i = d_gas_zeta[i];

        double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;  // Push root of DM tree

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = d_dm_bvh[node_idx];

            // Spatial Culling: Does the search sphere intersect this node's
            // AABB?
            double aabb_dist_sq =
                min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                     domain_size) +
                min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                     domain_size) +
                min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                     domain_size);

            if (aabb_dist_sq > search_sq) continue;  // Prune branch completely

            if (node.particle_idx != -1) {
                // Exact P-P interaction at the leaf
                int j = node.particle_idx;

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

                // Flip back to (j - i) vector direction
                dx = -dx;
                dy = -dy;
                dz = -dz;

                double dist_sq = dx * dx + dy * dy + dz * dz;

                if (use_pm && dist_sq > cutoff_sq) continue;
                double r = std::sqrt(dist_sq + 1e-24);
                double m_j = d_dm_m[j];

#ifdef UNCORRECTED_GRAVITY
                // Pure Plummer Softening (Matches DM exactly, bypasses Kernels)
                double pp_dist_sq = dist_sq + soft_sq;
                double pp_dist = std::sqrt(pp_dist_sq);

                // We define force_mag_over_r so that when it is later
                // multiplied by (m_j * dx), it equals the DM equation: (G * m_j
                // * dx) / (pp_dist^3)
                double force_mag_over_r = G / (pp_dist_sq * pp_dist);

                // Override 'r' to match the DM's PM taper behavior perfectly
                // r = pp_dist;
#else
                // Gas kernel derivatives
                double dphi_dr_i, dW_dr_i, dphi_dr_j, dummy_dphi_dh;
                Kernels::adaptive_gravity_terms(r, h_i, dphi_dr_i, dW_dr_i);

                // DM kernel derivatives
                Kernels::gravity_derivatives(r, spline_equivalent_h, dphi_dr_j,
                                             dummy_dphi_dh);

                // Force magnitude
                double force_mag_over_r =
                    (G / 2.0) *
                    ((dphi_dr_i + dphi_dr_j) + (zeta_i * dW_dr_i) / m_gas) / r;
#endif

                if (use_pm) {
                    double r_scaled = r / (2.0 * r_s);
                    double taper = std::erfc(r_scaled) +
                                   (r / (std::sqrt(M_PI) * r_s)) *
                                       std::exp(-r_scaled * r_scaled);
                    force_mag_over_r *= taper;
                }

                // Pull on MFM Gas
                local_acc_x += force_mag_over_r * m_j * dx;
                local_acc_y += force_mag_over_r * m_j * dy;
                local_acc_z += force_mag_over_r * m_j * dz;

                // Pull on DM (Equal and opposite, N3L)
                double a_dm = force_mag_over_r * m_gas;
#pragma omp atomic
                d_dm_ax[j] -= a_dm * dx;
#pragma omp atomic
                d_dm_ay[j] -= a_dm * dy;
#pragma omp atomic
                d_dm_az[j] -= a_dm * dz;

            } else {
                // Push children
                stack[stack_ptr++] = node.left_child;
                stack[stack_ptr++] = node.right_child;
            }
        }

        d_gas_ax[i] += local_acc_x;
        d_gas_ay[i] += local_acc_y;
        d_gas_az[i] += local_acc_z;
    }

#ifdef USE_GPU
    auto end_compute = std::chrono::high_resolution_clock::now();
    auto start_return = std::chrono::high_resolution_clock::now();

#pragma omp target exit data map(                                             \
        from : d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas], \
            d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm])          \
    map(delete : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas],                    \
            d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas], d_gas_h[0 : n_gas],      \
            d_gas_zeta[0 : n_gas], d_dm_px[0 : n_dm], d_dm_py[0 : n_dm],      \
            d_dm_pz[0 : n_dm], d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes])

    auto end_return = std::chrono::high_resolution_clock::now();

    diag.add_prof_time(ProfRegion::Transf,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_transfer - start_transfer)
                           .count());
    diag.add_prof_time(ProfRegion::Compute,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_compute - start_compute)
                           .count());
    diag.add_prof_time(ProfRegion::Ret,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_return - start_return)
                           .count());
#endif
}

void GasParticleSystem::hydro_step(const Config& config, double a, double H,
                                   double dt) {
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
            // Re-sync the entropy array to match the shock-heated state
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
        total_energy[i] = u[i] + ke;

        metal_frac[i] = std::max(0.0, std::min(metal_frac[i], 1.0));
    }
    num_cycles++;

    accumulated_photoheating_energy += step_floor_heating;
    accumulated_entropy_switch_energy += step_entropy_switch;
}

void GasParticleSystem::compute_gradients(const Config& config) {
    if (num_particles == 0) return;
    double domain_size = config.domain_size;
    size_t step_ill_conditioned = 0;

#pragma omp parallel for reduction(+ : step_ill_conditioned) \
    schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        // Construct the isolated state for particle i
        Reconstruction::ParticleState p_i;
        p_i.pos = Eigen::Vector3d(pos_x[i], pos_y[i], pos_z[i]);
        p_i.vel = Eigen::Vector3d(vel_x[i], vel_y[i], vel_z[i]);
        p_i.mass = mass[i];
        p_i.rho = rho[i];
        p_i.pressure = pressure[i];
        p_i.h = h[i];

        double local_max_rel_ke = 0.0;
        std::vector<Reconstruction::ParticleState> neighbors;

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;  // Push root

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = bvh_nodes[node_idx];

            double dist_sq =
                min_periodic_dist_sq(p_i.pos.x(), node.bbox.min_x,
                                     node.bbox.max_x, domain_size) +
                min_periodic_dist_sq(p_i.pos.y(), node.bbox.min_y,
                                     node.bbox.max_y, domain_size) +
                min_periodic_dist_sq(p_i.pos.z(), node.bbox.min_z,
                                     node.bbox.max_z, domain_size);

            double eff_h = std::max(p_i.h, node.max_h);
            if (dist_sq > eff_h * eff_h) continue;

            if (node.particle_idx != -1) {
                int j = node.particle_idx;

                double dx =
                    periodic_displacement(pos_x[j] - p_i.pos.x(), domain_size);
                double dy =
                    periodic_displacement(pos_y[j] - p_i.pos.y(), domain_size);
                double dz =
                    periodic_displacement(pos_z[j] - p_i.pos.z(), domain_size);
                double r2 = dx * dx + dy * dy + dz * dz;

                if ((r2 < p_i.h * p_i.h || r2 < h[j] * h[j]) && r2 > 1e-24) {
                    Reconstruction::ParticleState nj;
                    nj.pos = Eigen::Vector3d(pos_x[j], pos_y[j], pos_z[j]);
                    nj.vel = Eigen::Vector3d(vel_x[j], vel_y[j], vel_z[j]);
                    nj.mass = mass[j];
                    nj.rho = rho[j];
                    nj.pressure = pressure[j];
                    nj.h = h[j];
                    neighbors.push_back(nj);

                    double rel_vx = vel_x[j] - p_i.vel.x();
                    double rel_vy = vel_y[j] - p_i.vel.y();
                    double rel_vz = vel_z[j] - p_i.vel.z();
                    double rel_v2 =
                        rel_vx * rel_vx + rel_vy * rel_vy + rel_vz * rel_vz;

                    local_max_rel_ke = std::max(local_max_rel_ke, 0.5 * rel_v2);
                }
            } else {
                stack[stack_ptr++] = node.left_child;
                stack[stack_ptr++] = node.right_child;
            }
        }

        max_rel_ke[i] = local_max_rel_ke;

        // Calculate delta E_grav = |a_grav| * h
        double a_grav_mag = std::sqrt(
            acc_x[i] * acc_x[i] + acc_y[i] * acc_y[i] + acc_z[i] * acc_z[i]);
        delta_E_grav[i] = a_grav_mag * p_i.h;

        Reconstruction::ParticleGradients grads =
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

// MFM Riemann Solver (Frame Boosted)
MFMFaceFlux solve_mfm_riemann(const Reconstruction::ReconstructedFace& face,
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

        // Enforce physical bounds on the Riemann fan.
        // Prevents S_star from exploding to +/- infinity if den_star is
        // corrupted
        S_star = std::max(S_L, std::min(S_star, S_R));

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

Eigen::Vector3d compute_mfm_face_area_vector(
    const Reconstruction::ParticleState& p_i, const Eigen::Matrix3d& B_i,
    const Reconstruction::ParticleState& p_j, const Eigen::Matrix3d& B_j,
    const Reconstruction::ReconstructedFace& face) {
    double V_i = p_i.mass / p_i.rho;
    double V_j = p_j.mass / p_j.rho;

    double W_i, W_j;
    Kernels::cubic_spline_value(face.r, p_i.h, W_i);
    Kernels::cubic_spline_value(face.r, p_j.h, W_j);

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

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        // Build the isolated state proxy for particle i
        Reconstruction::ParticleState p_i;
        p_i.pos = Eigen::Vector3d(pos_x[i], pos_y[i], pos_z[i]);
        p_i.vel = Eigen::Vector3d(vel_x[i], vel_y[i], vel_z[i]);
        p_i.mass = mass[i];
        p_i.rho = rho[i];
        p_i.pressure = pressure[i];
        p_i.h = h[i];

        Reconstruction::ParticleGradients grad_i;
        grad_i.B_matrix = B_matrix[i];
        grad_i.grad_rho = grad_rho[i];
        grad_i.grad_p = grad_p[i];
        grad_i.grad_vx = grad_vx[i];
        grad_i.grad_vy = grad_vy[i];
        grad_i.grad_vz = grad_vz[i];

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;  // Push root

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = bvh_nodes[node_idx];

            double dist_sq =
                min_periodic_dist_sq(p_i.pos.x(), node.bbox.min_x,
                                     node.bbox.max_x, domain_size) +
                min_periodic_dist_sq(p_i.pos.y(), node.bbox.min_y,
                                     node.bbox.max_y, domain_size) +
                min_periodic_dist_sq(p_i.pos.z(), node.bbox.min_z,
                                     node.bbox.max_z, domain_size);

            double eff_h = std::max(p_i.h, node.max_h);
            if (dist_sq > eff_h * eff_h) continue;

            if (node.particle_idx != -1) {
                int j = node.particle_idx;

                // Process each pair just once to maintain N3L conservation
                if (j <= static_cast<int>(i)) continue;

                double dx =
                    periodic_displacement(pos_x[j] - p_i.pos.x(), domain_size);
                double dy =
                    periodic_displacement(pos_y[j] - p_i.pos.y(), domain_size);
                double dz =
                    periodic_displacement(pos_z[j] - p_i.pos.z(), domain_size);
                double r2 = dx * dx + dy * dy + dz * dz;

                if (r2 < 1e-24) continue;

                double hj = h[j];

                // Particles interact if they fall within either smoothing
                // length
                if (r2 < p_i.h * p_i.h || r2 < hj * hj) {
                    Reconstruction::ParticleState p_j;
                    p_j.pos = Eigen::Vector3d(pos_x[j], pos_y[j], pos_z[j]);
                    p_j.vel = Eigen::Vector3d(vel_x[j], vel_y[j], vel_z[j]);
                    p_j.mass = mass[j];
                    p_j.rho = rho[j];
                    p_j.pressure = pressure[j];
                    p_j.h = hj;

                    Reconstruction::ParticleGradients grad_j;
                    grad_j.B_matrix = B_matrix[j];
                    grad_j.grad_rho = grad_rho[j];
                    grad_j.grad_p = grad_p[j];
                    grad_j.grad_vx = grad_vx[j];
                    grad_j.grad_vy = grad_vy[j];
                    grad_j.grad_vz = grad_vz[j];

                    // SOLVER PIPELINE
                    Reconstruction::ReconstructedFace face =
                        compute_face_reconstruction(p_i, grad_i, p_j, grad_j,
                                                    domain_size, density_floor,
                                                    this->pressure_floor);
                    if (!face.is_valid) continue;

                    Eigen::Vector3d Area_vec = compute_mfm_face_area_vector(
                        p_i, grad_i.B_matrix, p_j, grad_j.B_matrix, face);
                    double A_mag = Area_vec.norm();
                    if (A_mag < 1e-20) continue;

                    face.n = Area_vec / A_mag;

                    double fraction_i = p_i.h / (p_i.h + p_j.h);
                    double fraction_j = 1.0 - fraction_i;
                    Eigen::Vector3d v_frame =
                        p_i.vel + fraction_i * (p_j.vel - p_i.vel);

                    // Convert to physical units
                    double a_inv = 1.0 / a;
                    double a_inv3 = a_inv * a_inv * a_inv;

                    Reconstruction::ReconstructedFace face_phys = face;
                    face_phys.rho_L = face.rho_L * a_inv3;
                    face_phys.rho_R = face.rho_R * a_inv3;
                    face_phys.p_L = face.p_L * a_inv;
                    face_phys.p_R = face.p_R * a_inv;
                    face_phys.v_L = face.v_L * a;
                    face_phys.v_R = face.v_R * a;

                    Eigen::Vector3d v_frame_phys = v_frame * a;

                    MFMFaceFlux flux_1d_phys = solve_mfm_riemann(
                        face_phys, v_frame_phys, config.gamma);

                    // Back to code units
                    double P_star_com = flux_1d_phys.P_star * a;
                    double S_star_com = flux_1d_phys.S_star * a_inv;

                    Eigen::Vector3d Force_mom = P_star_com * Area_vec;
                    Eigen::Vector3d v_star_lab =
                        v_frame + (S_star_com * face.n);

                    double work_i =
                        P_star_com * (v_star_lab - p_i.vel).dot(Area_vec);
                    double work_j =
                        P_star_com * (v_star_lab - p_j.vel).dot(Area_vec);

                    double du_dt_i = -work_i / p_i.mass;
                    double du_dt_j = work_j / p_j.mass;

                    double Rate_energy = P_star_com * v_star_lab.dot(Area_vec);
                    double de_dt_i = -Rate_energy / p_i.mass;
                    double de_dt_j = Rate_energy / p_j.mass;

                    if (config.disable_hydro_forces) continue;

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
            } else {
                stack[stack_ptr++] = node.left_child;
                stack[stack_ptr++] = node.right_child;
            }
        }
    }
}