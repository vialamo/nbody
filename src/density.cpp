#include "density.h"

#include <omp.h>

#include <cmath>
#include <iostream>

#include "kernels.h"
#include "math_utils.h"

namespace DensitySolver {

namespace {  // Anonymous namespace for private helper functions

// Evaluates the condition number of a geometric E matrix. Returns -1.0 if
// singular
inline double evaluate_matrix_condition_number(const Eigen::Matrix3d& E,
                                               const Eigen::Matrix3d& B) {
    double det = E.determinant();
    if (std::abs(det) < 1e-30) return -1.0;
    return (1.0 / 3.0) * std::sqrt(E.squaredNorm() * B.squaredNorm());
}

void evaluate_density_sum(const GasParticleSystem& gas, size_t particle_idx,
                          double h_guess, const Vec3& domain_size,
                          double& out_n, double& out_dn_dh,
                          int& out_num_neighb) {
    out_n = 0.0;
    out_dn_dh = 0.0;
    out_num_neighb = 0;
    double search_sq = h_guess * h_guess;
    double p1_x = gas.pos_x[particle_idx];
    double p1_y = gas.pos_y[particle_idx];
    double p1_z = gas.pos_z[particle_idx];

    int stack[128];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    while (stack_ptr > 0) {
        int node_idx = stack[--stack_ptr];
        const BVHNode& node = gas.bvh_nodes[node_idx];

        double dist_sq = min_periodic_dist_sq(p1_x, node.bbox.min_x,
                                              node.bbox.max_x, domain_size.x) +
                         min_periodic_dist_sq(p1_y, node.bbox.min_y,
                                              node.bbox.max_y, domain_size.y) +
                         min_periodic_dist_sq(p1_z, node.bbox.min_z,
                                              node.bbox.max_z, domain_size.z);

        if (dist_sq > search_sq) continue;

        if (node.particle_idx != -1) {
            int j = node.particle_idx;
            double dx =
                periodic_displacement(gas.pos_x[j] - p1_x, domain_size.x);
            double dy =
                periodic_displacement(gas.pos_y[j] - p1_y, domain_size.y);
            double dz =
                periodic_displacement(gas.pos_z[j] - p1_z, domain_size.z);
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 < search_sq) {
                double r = std::sqrt(r2);
                double W, dWdh;
                Kernels::cubic_spline(r, h_guess, W, dWdh);
                out_n += W;
                out_dn_dh += dWdh;
                out_num_neighb++;
            }
        } else {
            stack[stack_ptr++] = node.left_child;
            stack[stack_ptr++] = node.right_child;
        }
    }
}

double build_matrix_and_compute_condition(const GasParticleSystem& gas,
                                          size_t particle_idx, double h_guess,
                                          const Vec3& domain_size,
                                          Eigen::Matrix3d& B,
                                          bool volume_weighted) {
    Eigen::Matrix3d E = Eigen::Matrix3d::Zero();
    double p1_x = gas.pos_x[particle_idx];
    double p1_y = gas.pos_y[particle_idx];
    double p1_z = gas.pos_z[particle_idx];
    double search_sq = h_guess * h_guess;

    int stack[128];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    while (stack_ptr > 0) {
        int node_idx = stack[--stack_ptr];
        const BVHNode& node = gas.bvh_nodes[node_idx];

        double dist_sq = min_periodic_dist_sq(p1_x, node.bbox.min_x,
                                              node.bbox.max_x, domain_size.x) +
                         min_periodic_dist_sq(p1_y, node.bbox.min_y,
                                              node.bbox.max_y, domain_size.y) +
                         min_periodic_dist_sq(p1_z, node.bbox.min_z,
                                              node.bbox.max_z, domain_size.z);

        if (dist_sq > search_sq) continue;

        if (node.particle_idx != -1) {
            int j = node.particle_idx;
            if (static_cast<size_t>(j) == particle_idx) continue;

            double dx =
                periodic_displacement(gas.pos_x[j] - p1_x, domain_size.x);
            double dy =
                periodic_displacement(gas.pos_y[j] - p1_y, domain_size.y);
            double dz =
                periodic_displacement(gas.pos_z[j] - p1_z, domain_size.z);
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 < search_sq && r2 > 1e-24) {
                double r = std::sqrt(r2);
                double W;
                Kernels::cubic_spline_value(r, h_guess, W);
                double weight = W;
                if (volume_weighted) {
                    // Use a safe density floor (1e-12) to avoid division by
                    // zero
                    double V_j = gas.mass[j] / std::max(gas.rho[j], 1e-12);
                    weight = W * V_j;
                }
                E(0, 0) += dx * dx * weight;
                E(0, 1) += dx * dy * weight;
                E(0, 2) += dx * dz * weight;
                E(1, 1) += dy * dy * weight;
                E(1, 2) += dy * dz * weight;
                E(2, 2) += dz * dz * weight;
            }
        } else {
            stack[stack_ptr++] = node.left_child;
            stack[stack_ptr++] = node.right_child;
        }
    }

    E(1, 0) = E(0, 1);
    E(2, 0) = E(0, 2);
    E(2, 1) = E(1, 2);

    B = E.inverse();
    return evaluate_matrix_condition_number(E, B);
}

#ifdef USE_ADAPTIVE_SOFTENING
double compute_zeta_contribution(
    double p1_x, double p1_y, double p1_z, double h_i,
    const std::vector<double>& target_x, const std::vector<double>& target_y,
    const std::vector<double>& target_z, const std::vector<double>& target_mass,
    const std::vector<BVHNode>& target_bvh, const Config& config) {
    double zeta_sum = 0.0;
    double search_sq_zeta = h_i * h_i;
    const Vec3& domain_size = config.domain_size;
    bool use_pm = config.use_PM;
    double r_s = config.PM_smoothing_cells * config.cell_size;

    int stack[128];
    int stack_ptr = 0;
    stack[stack_ptr++] = 0;

    while (stack_ptr > 0) {
        int node_idx = stack[--stack_ptr];
        const BVHNode& node = target_bvh[node_idx];

        double dist_sq = min_periodic_dist_sq(p1_x, node.bbox.min_x,
                                              node.bbox.max_x, domain_size.x) +
                         min_periodic_dist_sq(p1_y, node.bbox.min_y,
                                              node.bbox.max_y, domain_size.y) +
                         min_periodic_dist_sq(p1_z, node.bbox.min_z,
                                              node.bbox.max_z, domain_size.z);

        if (dist_sq > search_sq_zeta) continue;

        if (node.particle_idx != -1) {
            int j = node.particle_idx;
            double dx =
                periodic_displacement(target_x[j] - p1_x, domain_size.x);
            double dy =
                periodic_displacement(target_y[j] - p1_y, domain_size.y);
            double dz =
                periodic_displacement(target_z[j] - p1_z, domain_size.z);
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
#endif

}  // end anonymous namespace

// --------------------------------------------------------------------------------
// Density and Smoothing Length Iteration
// --------------------------------------------------------------------------------
void compute_density_and_h(GasParticleSystem& gas, const Config& config,
                           const ParticleSystem& dm) {
    if (gas.num_particles == 0) return;

    const Vec3& domain_size = config.domain_size;
    double target_N = config.mfm_target_neighbors;
    double tol = config.mfm_neighbor_tolerance;
    int max_iter = config.mfm_max_iterations;

    const double min_h_cap =
#ifdef USE_ADAPTIVE_SOFTENING
        2.8 *
#endif
        config.mfm_min_hsml_fraction * config.softening;

    const double max_h_cap = 0.5 * config.domain_size_1d;

    size_t num_h_clamped = 0;
    size_t num_non_converged = 0;
    size_t step_neighbor_increased_cases = 0;

    if (gas.num_active == 0) return;

#pragma omp parallel for schedule(dynamic, 64) reduction( \
        + : num_h_clamped, num_non_converged, step_neighbor_increased_cases)
    for (size_t k = 0; k < gas.num_active; ++k) {
        size_t i = gas.active_indices[k];

        double p1_x = gas.pos_x[i], p1_y = gas.pos_y[i], p1_z = gas.pos_z[i];
        double h_low = 0.0;
        double h_high = std::numeric_limits<double>::infinity();
        double h_guess = gas.h[i];
        double current_n_enc = gas.n_enc_final[i];

        double current_target_N = target_N;
        double current_n = 0.0;
        double current_dn_dh = 0.0;
        int current_num_neighb = 0;

        bool matrix_converged = false;
        int matrix_iters = 0;
        bool h_clamped = false;
        bool is_converged = false;

        // OUTER LOOP: Matrix Regularization via Target Expansion
        while (!matrix_converged && matrix_iters < 10 && h_guess < max_h_cap) {
            int nr_iter = 0;
            is_converged = false;
            h_clamped = false;

            // INNER LOOP: Standard Newton-Raphson Solver
            while (nr_iter < max_iter) {
                evaluate_density_sum(gas, i, h_guess, domain_size, current_n,
                                     current_dn_dh, current_num_neighb);

                double h3 = h_guess * h_guess * h_guess;
                double N_enc = (4.0 / 3.0) * M_PI * h3 * current_n;
                current_n_enc = N_enc;

                if (std::abs(N_enc - current_target_N) < tol) {
                    is_converged = true;
                    break;
                }

                // Floor Termination
                if (h_guess <= min_h_cap && N_enc > current_target_N) {
                    h_clamped = true;
                    is_converged = true;
                    break;
                }

                // Ceiling Termination
                if (h_guess >= max_h_cap && N_enc < current_target_N) {
                    h_clamped = true;
                    is_converged = true;
                    break;
                }

                // Update bounds
                if (N_enc > current_target_N)
                    h_high = h_guess;
                else
                    h_low = h_guess;

                // Calculate next step
                double dN_enc_dh =
                    (4.0 / 3.0) * M_PI *
                    (3.0 * h_guess * h_guess * current_n + h3 * current_dn_dh);
                double h_new = h_guess;

                if (dN_enc_dh > 0.0)
                    h_new = h_guess - (N_enc - current_target_N) / dN_enc_dh;

                // Bisection fallback
                if (h_new <= h_low || h_new >= h_high || dN_enc_dh <= 0.0) {
                    h_guess = std::isinf(h_high) ? (1.26 * h_guess)
                                                 : std::sqrt(h_low * h_high);
                } else {
                    h_guess = h_new;
                }

                // Safe Clamping
                if (h_guess > max_h_cap) h_guess = max_h_cap;
                if (h_guess < min_h_cap) h_guess = min_h_cap;

                nr_iter++;
            }

            // Sync properties if we exited without convergence
            if (h_clamped || !is_converged) {
                evaluate_density_sum(gas, i, h_guess, domain_size, current_n,
                                     current_dn_dh, current_num_neighb);
                double h3 = h_guess * h_guess * h_guess;
                current_n_enc = (4.0 / 3.0) * M_PI * h3 * current_n;
            }

            // If clamped at the ceiling, further expansion is impossible
            if (h_clamped && h_guess >= max_h_cap) break;

            // Evaluate matrix condition
            double cond = build_matrix_and_compute_condition(
                gas, i, h_guess, domain_size, gas.B_matrix[i], false);

            constexpr int min_num_neighb = 20;
            if ((cond < 0.0 || cond > N_cond_crit ||
                 current_num_neighb < min_num_neighb) &&
                current_target_N < 3.0 * target_N) {
                // Find more neighbors
                current_target_N *= 1.2;

                // Reset the bounds because the new root is larger
                h_low = h_guess;
                h_high = std::numeric_limits<double>::infinity();
                matrix_iters++;
                step_neighbor_increased_cases++;
            } else {
                matrix_converged = true;
            }
        }

        // Commit finalized state
        gas.cond_num[i] = build_matrix_and_compute_condition(
            gas, i, h_guess, domain_size, gas.B_matrix[i], false);
        gas.h[i] = h_guess;
        gas.rho[i] = gas.mass[i] * current_n;
        gas.n_enc_final[i] = current_n_enc;

        if (h_clamped) num_h_clamped++;
        if (!is_converged) num_non_converged++;
    }

// #define VOLUME_RENORMALIZATION
#ifdef VOLUME_RENORMALIZATION
    // Only valid when all particles are active
    if (gas.num_active == gas.num_particles) {
        double alpha = 1.0;
        double volume_error = 1.0;
        int vol_iter = 0;
        double box_volume = domain_size.x * domain_size.y * domain_size.z;

        // Save the baseline topological lengths
        std::vector<double> h_base = gas.h;

        while (std::abs(volume_error) > 1e-5 && vol_iter < 20) {
            double total_sim_volume = 0.0;

#pragma omp parallel for reduction(+ : total_sim_volume) schedule(dynamic, 64)
            for (size_t k = 0; k < gas.num_active; ++k) {
                size_t i = gas.active_indices[k];
                gas.h[i] = h_base[i] * alpha;

                // Re-evaluate from the tree. h and rho remain synced
                double current_n = 0.0, dummy_dn_dh = 0.0;
                int current_num_neighb = 0;
                evaluate_density_sum(gas, i, gas.h[i], domain_size, current_n,
                                     dummy_dn_dh, current_num_neighb);

                gas.rho[i] = gas.mass[i] * current_n;
                total_sim_volume += gas.mass[i] / std::max(gas.rho[i], 1e-12);
            }

            // V_sim scales with alpha^3. Use cbrt for rapid convergence
            volume_error = (total_sim_volume / box_volume) - 1.0;
            alpha *= std::cbrt(box_volume / total_sim_volume);
            vol_iter++;
        }
    }
#endif

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t k = 0; k < gas.num_active; ++k) {
        size_t i = gas.active_indices[k];

        // Sync final geometric N_enc for diagnostics (will slightly deviate
        // from Target_N)
        double current_n = gas.rho[i] / gas.mass[i];
        gas.n_enc_final[i] =
            (4.0 / 3.0) * M_PI * std::pow(gas.h[i], 3) * current_n;

        gas.cond_num[i] = build_matrix_and_compute_condition(
            gas, i, gas.h[i], domain_size, gas.B_matrix[i], true);

#ifdef USE_ADAPTIVE_SOFTENING
        if (config.use_PP) {
            double current_dn_dh = 0.0, dummy_n = 0.0;
            int current_num_neighb = 0;
            evaluate_density_sum(gas, i, gas.h[i], domain_size, dummy_n,
                                 current_dn_dh, current_num_neighb);

            // Adaptive gravity correction (Zeta)
            double Omega_i = 1.0;
            double h_over_3n = 0.0;
            if (current_n > 1e-12) {
                h_over_3n = gas.h[i] / (current_n * 3.0);
                Omega_i = std::max(1.0 + h_over_3n * current_dn_dh, 1e-12);
            }

            // Sum Gas contribution
            double p1_x = gas.pos_x[i], p1_y = gas.pos_y[i],
                   p1_z = gas.pos_z[i];
            double zeta_sum = compute_zeta_contribution(
                p1_x, p1_y, p1_z, gas.h[i], gas.pos_x, gas.pos_y, gas.pos_z,
                gas.mass, gas.bvh_nodes, config);

            // Sum DM contribution
            if (dm.num_particles > 0 && !dm.bvh_nodes.empty()) {
                zeta_sum += compute_zeta_contribution(
                    p1_x, p1_y, p1_z, gas.h[i], dm.pos_x, dm.pos_y, dm.pos_z,
                    dm.mass, dm.bvh_nodes, config);
            }

            gas.zeta[i] = h_over_3n * (1.0 / Omega_i) * zeta_sum;
        }
#endif
    }

    gas.clamped_h_cases += num_h_clamped;
    gas.non_converged_h_cases += num_non_converged;
    gas.neighbor_increased_cases += step_neighbor_increased_cases;

    // Refresh Bounding Boxes using the final converged smoothing lengths
    LBVH::update_h_and_bboxes(gas.num_particles, gas.pos_x, gas.pos_y,
                              gas.pos_z, gas.h, gas.bvh_nodes);
}

}  // namespace DensitySolver
