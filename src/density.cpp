#include "density.h"

#include <omp.h>

#include "config.h"
#include "debug_utils.h"
#include "kernels.h"
#include "math_utils.h"
#include "mfm.h"

#define USE_NBG_LIST
// #define USE_ADAPTIVE_CONDITION_NUMBER

namespace DensitySolver {

namespace {  // Anonymous namespace for private helper functions

// Evaluates the condition number of a geometric E matrix. Returns -1.0 if
// singular
static inline double evaluate_matrix_condition_number(
    const Eigen::Matrix3d& E, const Eigen::Matrix3d& B) {
    double det = E.determinant();
    if (std::abs(det) < 1e-30) return -1.0;
    return (1.0 / 3.0) * std::sqrt(E.squaredNorm() * B.squaredNorm());
}

#ifdef USE_ADAPTIVE_CONDITION_NUMBER
static void evaluate_density_sum(const GasParticleSystem& gas,
                                 size_t particle_idx, double h,
                                 const Vec3& domain_size, double& out_n,
                                 double& out_dn_dh, int& out_num_neighb) {
    out_n = 0.0;
    out_dn_dh = 0.0;
    out_num_neighb = 0;
    double search_sq = h * h;
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
                Kernels::cubic_spline(r, h, W, dWdh);
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
#endif

#ifdef USE_NBG_LIST
static void gather_neighbors(const GasParticleSystem& gas, size_t particle_idx,
                             double h, const Vec3& domain_size,
                             std::vector<int>& out_indices,
                             std::vector<double>& out_r2) {
    out_indices.clear();
    out_r2.clear();
    double search_sq = h * h;
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

            if (r2 <= search_sq) {
                out_indices.push_back(j);
                out_r2.push_back(r2);
            }
        } else {
            stack[stack_ptr++] = node.left_child;
            stack[stack_ptr++] = node.right_child;
        }
    }
}

static void evaluate_density_sum_from_cache(
    double h, const std::vector<int>& cached_indices,
    const std::vector<double>& cached_r2, double& out_n, double& out_dn_dh,
    int& out_num_neighb) {
    out_n = 0.0;
    out_dn_dh = 0.0;
    out_num_neighb = 0;
    double search_sq = h * h;

    for (size_t k = 0; k < cached_indices.size(); ++k) {
        double r2 = cached_r2[k];
        if (r2 < search_sq) {
            double r = std::sqrt(r2);
            double W, dWdh;
            Kernels::cubic_spline(r, h, W, dWdh);
            out_n += W;
            out_dn_dh += dWdh;
            out_num_neighb++;
        }
    }
}
#endif

static double build_matrix_and_compute_condition(const GasParticleSystem& gas,
                                                 size_t particle_idx, double h,
                                                 const Vec3& domain_size,
                                                 Eigen::Matrix3d& B) {
    Eigen::Matrix3d E = Eigen::Matrix3d::Zero();

    double p1_x = gas.pos_x[particle_idx];
    double p1_y = gas.pos_y[particle_idx];
    double p1_z = gas.pos_z[particle_idx];
    double search_sq = h * h;

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
                Kernels::cubic_spline_value(r, h, W);

#ifdef USE_GEOMETRIC_MATRIX
                double weight = W;
#else
                double V_j = gas.mass[j] / std::max(gas.rho[j], 1e-12);
                double weight = W * V_j;
#endif

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

// #define USE_PSEUDOINVERSE
#ifndef USE_PSEUDOINVERSE
    B = E.inverse();
    return evaluate_matrix_condition_number(E, B);
#else
    // Compute Moore-Penrose pseudoinverse via Eigendecomposition.
    // (Mathematically identical to SVD for symmetric positive semi-definite
    // matrices)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(E);
    Eigen::Vector3d eigenvalues = solver.eigenvalues();
    double max_eigenvalue = eigenvalues.maxCoeff();

    // Scale-invariant absolute floor. Tied to the smoothing length (h^2)
    // to prevent matrix conditioning from breaking if simulation units change
    double absolute_tol = 1e-14 * h * h;

    // Relative geometric floor. A high threshold (e.g., 1e-1 to 1e-2)
    // aggressively collapses the solver to zeroth-order in highly deformed
    // spatial distributions (like 1D strings or 2D sheets), preventing
    // unphysical thermal explosions in voids
    constexpr double relative_tolerance_factor = 1e-4;
    double relative_tol = relative_tolerance_factor * max_eigenvalue;

    double threshold = std::max(absolute_tol, relative_tol);

    Eigen::Vector3d inv_eigenvalues = Eigen::Vector3d::Zero();

#define USE_TIKHONOV_REGULARIZATION
#ifdef TIKHONOV_REGULARIZATION
    // Tikhonov Regularization (Damped Pseudoinverse)
    // Smoothly suppresses degenerate dimensions to prevent abrupt mathematical
    // discontinuities as particles enter and exit poorly-sampled configurations
    double threshold_sq = threshold * threshold;
    for (int i = 0; i < 3; ++i) {
        inv_eigenvalues(i) =
            eigenvalues(i) / (eigenvalues(i) * eigenvalues(i) + threshold_sq);
    }
#else
    // Truncated Pseudoinverse (Hard Cutoff)
    // Completely zeroes out the gradient and flux along severely
    // ill-conditioned axes
    for (int i = 0; i < 3; ++i) {
        if (eigenvalues(i) > threshold) {
            inv_eigenvalues(i) = 1.0 / eigenvalues(i);
        }
    }
#endif

    // Reconstruct the regularized shape matrix (B)
    B = solver.eigenvectors() * inv_eigenvalues.asDiagonal() *
        solver.eigenvectors().transpose();

    // Return dummy condition number since the pseudoinverse is unconditionally
    // stable
    return 1.0;

#endif
}

}  // end anonymous namespace

// --------------------------------------------------------------------------------
// Density and Smoothing Length Iteration
// --------------------------------------------------------------------------------
void compute_density_and_h(GasParticleSystem& gas, const Config& config) {
    if (gas.num_particles == 0) return;

    const Vec3& domain_size = config.domain_size;
    double target_N = config.mfm_target_neighbors;
    double tol = config.mfm_neighbor_tolerance;
    int max_iter = config.mfm_max_iterations;

    const double min_h_cap = config.mfm_min_hsml_fraction * config.softening;
    const double max_h_cap = 0.5 * config.domain_size_1d;

    constexpr double sphere_volume_factor = (4.0 / 3.0) * M_PI;

    size_t num_h_clamped = 0;
    size_t num_non_converged = 0;
    size_t step_neighbor_increased_cases = 0;
    size_t num_not_enough_neighbors = 0;

    if (gas.num_active == 0) return;

    // Cache to pass dn_dh to the adaptive softening loop without re-evaluating
    // the tree
    std::vector<double> cached_dn_dh(gas.num_particles, 0.0);

    // Cache the individual target_N so Pass 2 starts from Pass 1's expanded
    // target
    std::vector<double> particle_target_N(gas.num_particles, target_N);

    bool diagnostic_printed = false;

// ========================================================================
// PASS 1: Baseline Density and h Computation (Enforces Minimum Neighbors)
// ========================================================================
#ifdef USE_NBG_LIST
#pragma omp parallel reduction(+ : num_h_clamped, num_non_converged, \
                                   step_neighbor_increased_cases)
#endif
    {
#ifdef USE_NBG_LIST
        // Thread-local caches to prevent constant reallocation
        std::vector<int> cached_indices;
        std::vector<double> cached_r2;
        cached_indices.reserve(256);
        cached_r2.reserve(256);
#pragma omp for schedule(dynamic, 64)
#else
#pragma omp parallel for schedule(dynamic, 64) reduction( \
        + : num_h_clamped, num_non_converged, step_neighbor_increased_cases)
#endif
        for (size_t k = 0; k < gas.num_active; ++k) {
            size_t i = gas.active_indices[k];

            double h_guess = gas.h[i];
            if (h_guess < min_h_cap) h_guess = min_h_cap;
            if (h_guess > max_h_cap) h_guess = max_h_cap;

            double current_target_N = target_N;

            double number_density = 0.0;
            double dn_dh = 0.0;
            int num_neighb = 0;
            double effective_neighbors = gas.n_enc_final[i];

            bool h_clamped = false;
            bool is_converged = false;

            bool neighbors_increased = false;

            double h_low = 0.0;
            double h_high = std::numeric_limits<double>::infinity();
            int nr_iter = 0;
            double cached_h = -1.0;  // Track the current cache radius

            is_converged = false;
            h_clamped = false;

            // Newton-Raphson Solver with Bisection fallback
            while (nr_iter < max_iter) {
#ifdef USE_NBG_LIST
                // Only rebuild the cache if the new guess exceeds our
                // gathered radius
                if (h_guess > cached_h) {
                    cached_h = h_guess;
                    gather_neighbors(gas, i, cached_h, domain_size,
                                     cached_indices, cached_r2);
                }

                evaluate_density_sum_from_cache(h_guess, cached_indices,
                                                cached_r2, number_density,
                                                dn_dh, num_neighb);
#else
                evaluate_density_sum(gas, i, h_guess, domain_size,
                                     number_density, dn_dh, num_neighb);
#endif

                double h3 = h_guess * h_guess * h_guess;
                effective_neighbors =
                    sphere_volume_factor * h3 * number_density;

                if (std::abs(effective_neighbors - current_target_N) < tol) {
                    is_converged = true;
                    break;
                }

                // Floor / Ceiling Termination
                if (h_guess <= min_h_cap &&
                    effective_neighbors > current_target_N) {
                    h_clamped = true;
                    is_converged = true;
                    break;
                }
                if (h_guess >= max_h_cap &&
                    effective_neighbors < current_target_N) {
                    h_clamped = true;
                    is_converged = true;
                    break;
                }

                // Update bounds
                if (effective_neighbors > current_target_N)
                    h_high = h_guess;
                else
                    h_low = h_guess;

                // Calculate next step
                double dN_enc_dh =
                    sphere_volume_factor *
                    (3.0 * h_guess * h_guess * number_density + h3 * dn_dh);
                double h_new = h_guess;

                if (dN_enc_dh > 0.0)
                    h_new = h_guess - (effective_neighbors - current_target_N) /
                                          dN_enc_dh;

                // Bisection fallback
                if (h_new <= h_low || h_new >= h_high || dN_enc_dh <= 0.0) {
                    h_new = std::isinf(h_high) ? (1.26 * h_guess)
                                               : std::sqrt(h_low * h_high);
                }

                // Safe Clamping
                if (h_new > max_h_cap) h_new = max_h_cap;
                if (h_new < min_h_cap) h_new = min_h_cap;

                nr_iter++;

                if (nr_iter < max_iter) {
                    h_guess = h_new;
                }
            }

            // Commit finalized state for Pass 1
            gas.h[i] = h_guess;
            gas.rho[i] = gas.mass[i] * number_density;
            gas.n_enc_final[i] = effective_neighbors;
            cached_dn_dh[i] = dn_dh;
            particle_target_N[i] = current_target_N;  // Save for Pass 2

            if (h_clamped) num_h_clamped++;
            if (!is_converged) {
                num_non_converged++;
                if (!h_clamped) {
                    DebugUtils::print_non_converged_diagnostic(
                        gas, i, h_guess, h_low, h_high, effective_neighbors,
                        current_target_N, target_N, tol, num_neighb, nr_iter,
                        max_iter, h_clamped, diagnostic_printed);
                }
            }
            if (neighbors_increased) step_neighbor_increased_cases++;
        }
    }

#ifdef USE_ADAPTIVE_CONDITION_NUMBER
// ========================================================================
// PASS 2: Matrix Evaluation & Target Expansion (using synced densities)
// ========================================================================
#pragma omp parallel for schedule(dynamic, 64) \
    reduction(+ : step_neighbor_increased_cases, num_not_enough_neighbors)
    for (size_t k = 0; k < gas.num_active; ++k) {
        size_t i = gas.active_indices[k];

        double h_guess = gas.h[i];
        double current_target_N =
            particle_target_N[i];  // Start from Pass 1's target

        int anisotropy_retries = 0;
        bool matrix_converged = false;
        bool neighbors_increased = false;
        double cond = 0.0;

        while (!matrix_converged && anisotropy_retries < max_iter &&
               h_guess < max_h_cap) {
            Eigen::Matrix3d temp_B;
            cond = build_matrix_and_compute_condition(gas, i, h_guess,
                                                      domain_size, temp_B);

            // Matrix regularization expansion loop
            if ((cond < 0.0 || cond > N_cond_crit) &&
                current_target_N < 1.2 * target_N) {
                // Find more neighbors
                current_target_N *= 1.2;
                anisotropy_retries++;
                neighbors_increased = true;

                double h_low = h_guess;
                double h_high = std::numeric_limits<double>::infinity();
                int nr_iter = 0;
                double effective_neighbors = 0;
                bool h_clamped = false;

                double number_density = 0.0;
                double dn_dh = 0.0;
                int num_neighb = 0;

                // Re-solve with Newton-Raphson
                while (nr_iter < max_iter) {
                    evaluate_density_sum(gas, i, h_guess, domain_size,
                                         number_density, dn_dh, num_neighb);

                    double h3 = h_guess * h_guess * h_guess;
                    effective_neighbors =
                        sphere_volume_factor * h3 * number_density;

                    if (std::abs(effective_neighbors - current_target_N) < tol)
                        break;

                    if (h_guess <= min_h_cap &&
                        effective_neighbors > current_target_N) {
                        h_clamped = true;
                        break;
                    }
                    if (h_guess >= max_h_cap &&
                        effective_neighbors < current_target_N) {
                        h_clamped = true;
                        break;
                    }

                    if (effective_neighbors > current_target_N)
                        h_high = h_guess;
                    else
                        h_low = h_guess;

                    double dN_enc_dh =
                        sphere_volume_factor *
                        (3.0 * h_guess * h_guess * number_density + h3 * dn_dh);
                    double h_new = h_guess;

                    if (dN_enc_dh > 0.0)
                        h_new =
                            h_guess - (effective_neighbors - current_target_N) /
                                          dN_enc_dh;

                    if (h_new <= h_low || h_new >= h_high || dN_enc_dh <= 0.0) {
                        h_new = std::isinf(h_high) ? (1.26 * h_guess)
                                                   : std::sqrt(h_low * h_high);
                    }

                    if (h_new > max_h_cap) h_new = max_h_cap;
                    if (h_new < min_h_cap) h_new = min_h_cap;

                    nr_iter++;
                    if (nr_iter < max_iter) h_guess = h_new;
                }

                gas.h[i] = h_guess;
                gas.rho[i] = gas.mass[i] * number_density;
                gas.n_enc_final[i] = effective_neighbors;
                cached_dn_dh[i] = dn_dh;

                if (h_clamped && h_guess >= max_h_cap) break;
            } else {
                matrix_converged = true;
            }
        }

        // Count neighbors locally from the last evaluate_density_sum
        double final_number_density = 0.0;
        double final_dn_dh = 0.0;
        int final_num_neighb = 0;
        evaluate_density_sum(gas, i, gas.h[i], domain_size,
                             final_number_density, final_dn_dh,
                             final_num_neighb);
        if (final_num_neighb < min_num_neighb) num_not_enough_neighbors++;

        if (neighbors_increased) step_neighbor_increased_cases++;
    }
#endif

    // ========================================================================
    // PASS 3: Final Matrix Sync & Adaptive Gravity
    // ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
    for (size_t k = 0; k < gas.num_active; ++k) {
        size_t i = gas.active_indices[k];

#ifdef USE_ADAPTIVE_CONDITION_NUMBER
        // Sync final geometric N_enc for diagnostics
        double number_density = gas.rho[i] / gas.mass[i];
        gas.n_enc_final[i] = sphere_volume_factor * gas.h[i] * gas.h[i] *
                             gas.h[i] * number_density;
#endif

        // Assign the final volume-weighted B matrix and condition
        gas.cond_num[i] = build_matrix_and_compute_condition(
            gas, i, gas.h[i], domain_size, gas.B_matrix[i]);

        if (gas.n_enc_final[i] > 3.0 * config.mfm_target_neighbors) {
            gas.cond_num[i] = 2 * N_cond_crit;
        }
    }

    gas.clamped_h_cases += num_h_clamped;
    gas.non_converged_h_cases += num_non_converged;
    gas.neighbor_increased_cases += step_neighbor_increased_cases;
    gas.not_enough_neighbors_cases += num_not_enough_neighbors;

    // Refresh Bounding Boxes using the final smoothing lengths
    LBVH::update_h_and_bboxes(gas.num_particles, gas.pos_x, gas.pos_y,
                              gas.pos_z, gas.h, gas.bvh_nodes);
}

}  // namespace DensitySolver
