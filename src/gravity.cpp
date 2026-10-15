#include "gravity.h"

#include <chrono>
#include <cmath>

#include "kernels.h"
#include "math_utils.h"
#include "pm.h"

namespace GravitySolver {

void update_softening(SimState& state, Config& config) {
    double a = state.scale_factor;
    double cap_a = config.physical_softening_cap_a;
    double current_comoving_softening = config.base_comoving_softening;

    if (a > cap_a && config.expanding_universe) {
        // Shrink comoving softening so physical softening stays constant
        current_comoving_softening =
            config.base_comoving_softening * (cap_a / a);
    }
    config.softening = current_comoving_softening;
}

void compute_and_add_gas_pp_forces(GasParticleSystem& gas, double a,
                                   const Config& config, Diagnostics& diag) {
    if (gas.num_particles == 0) return;

    const Vec3& domain_size = config.domain_size;
    const double G = config.G / (a * a * a);
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t n_gas = gas.num_particles;
    const size_t num_nodes = 2 * n_gas - 1;
#ifndef USE_ADAPTIVE_SOFTENING
    double epsilon = config.softening;
#endif

    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();

    double* d_px = gas.pos_x.data();
    double* d_py = gas.pos_y.data();
    double* d_pz = gas.pos_z.data();
    double* d_m = gas.mass.data();
    double* d_h = gas.h.data();
#ifdef USE_ADAPTIVE_SOFTENING
    double* d_zeta = gas.zeta.data();
#endif
    double* d_ax = gas.acc_x.data();
    double* d_ay = gas.acc_y.data();
    double* d_az = gas.acc_z.data();
    BVHNode* d_bvh_nodes = gas.bvh_nodes.data();
    int* d_active_idx = gas.active_indices.data();
    size_t n_active = gas.num_active;

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();
#ifdef USE_ADAPTIVE_SOFTENING
#pragma omp target enter data map(                                        \
        to : d_px[0 : n_gas], d_py[0 : n_gas], d_pz[0 : n_gas],           \
            d_m[0 : n_gas], d_h[0 : n_gas], d_zeta[0 : n_gas],            \
            d_bvh_nodes[0 : num_nodes], d_ax[0 : n_gas], d_ay[0 : n_gas], \
            d_az[0 : n_gas], d_active_idx[0 : n_active])
#else
#pragma omp target enter data map(                                      \
        to : d_px[0 : n_gas], d_py[0 : n_gas], d_pz[0 : n_gas],         \
            d_m[0 : n_gas], d_h[0 : n_gas], d_bvh_nodes[0 : num_nodes], \
            d_ax[0 : n_gas], d_ay[0 : n_gas], d_az[0 : n_gas],          \
            d_active_idx[0 : n_active])
#endif
    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t k = 0; k < n_active; ++k) {
        size_t i = d_active_idx[k];

        double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
#ifdef USE_ADAPTIVE_SOFTENING
        double m_i = d_m[i];
        double h_i = d_h[i];
        double zeta_i = d_zeta[i];
#endif
        double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = d_bvh_nodes[node_idx];

            // AABB Culling against cutoff radius
            double aabb_dist_sq =
                min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                     domain_size.x) +
                min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                     domain_size.y) +
                min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                     domain_size.z);

            if (aabb_dist_sq > search_sq) continue;

            if (node.particle_idx != -1) {
                int j = node.particle_idx;
                if (static_cast<size_t>(j) == i) continue;

                double dx = p1_x - d_px[j];
                if (dx > 0.5 * domain_size.x)
                    dx -= domain_size.x;
                else if (dx < -0.5 * domain_size.x)
                    dx += domain_size.x;

                double dy = p1_y - d_py[j];
                if (dy > 0.5 * domain_size.y)
                    dy -= domain_size.y;
                else if (dy < -0.5 * domain_size.y)
                    dy += domain_size.y;

                double dz = p1_z - d_pz[j];
                if (dz > 0.5 * domain_size.z)
                    dz -= domain_size.z;
                else if (dz < -0.5 * domain_size.z)
                    dz += domain_size.z;

                dx = -dx;
                dy = -dy;
                dz = -dz;
                double dist_sq = dx * dx + dy * dy + dz * dz;

                if (use_pm && dist_sq > cutoff_sq) continue;
                if (dist_sq < 1e-24) continue;

                double r = std::sqrt(dist_sq);
                double m_j = d_m[j];

#ifdef USE_ADAPTIVE_SOFTENING
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
#else
                // TreePM PP Force
                double dphi_dr, dummy_dphi_dh;
                Kernels::gravity_derivatives(r, epsilon, dphi_dr,
                                             dummy_dphi_dh);
                double force_mag_over_r = G * dphi_dr / r;
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
#ifdef USE_ADAPTIVE_SOFTENING
#pragma omp target exit data map(from : d_ax[0 : n_gas], d_ay[0 : n_gas], \
                                     d_az[0 : n_gas])                     \
    map(delete : d_px[0 : n_gas], d_py[0 : n_gas], d_pz[0 : n_gas],       \
            d_m[0 : n_gas], d_h[0 : n_gas], d_zeta[0 : n_gas],            \
            d_bvh_nodes[0 : num_nodes], d_active_idx[0 : n_active])
#else
#pragma omp target exit data map(from : d_ax[0 : n_gas], d_ay[0 : n_gas], \
                                     d_az[0 : n_gas])                     \
    map(delete : d_px[0 : n_gas], d_py[0 : n_gas], d_pz[0 : n_gas],       \
            d_m[0 : n_gas], d_h[0 : n_gas], d_bvh_nodes[0 : num_nodes],   \
            d_active_idx[0 : n_active])
#endif
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

static void compute_gas_cross_pp_forces(GasParticleSystem& gas, double a,
                                        ParticleSystem& dm,
                                        const Config& config,
                                        Diagnostics& diag) {
    if (gas.num_particles == 0 || dm.num_particles == 0) return;

    const Vec3& domain_size = config.domain_size;
    // Pre-scale G so all output forces are comoving accelerations
    const double G = config.G / (a * a * a);
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t n_gas = gas.num_particles;
    const size_t n_dm = dm.num_particles;
    const size_t dm_num_nodes = 2 * n_dm - 1;

    // Search radius is the PM cutoff. If not using PM, it's infinite (N^2
    // tree traversal)
    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();
    const double base_soft = config.softening;

    double* d_gas_px = gas.pos_x.data();
    double* d_gas_py = gas.pos_y.data();
    double* d_gas_pz = gas.pos_z.data();
    double* d_gas_m = gas.mass.data();
    double* d_gas_h = gas.h.data();
#ifdef USE_ADAPTIVE_SOFTENING
    double* d_gas_zeta = gas.zeta.data();
#endif
    double* d_gas_ax = gas.acc_x.data();
    double* d_gas_ay = gas.acc_y.data();
    double* d_gas_az = gas.acc_z.data();

    double* d_dm_px = dm.pos_x.data();
    double* d_dm_py = dm.pos_y.data();
    double* d_dm_pz = dm.pos_z.data();
    double* d_dm_m = dm.mass.data();
    double* d_dm_ax = dm.acc_x.data();
    double* d_dm_ay = dm.acc_y.data();
    double* d_dm_az = dm.acc_z.data();

    BVHNode* d_dm_bvh = dm.bvh_nodes.data();
    int* d_active_idx = gas.active_indices.data();
    size_t n_active = gas.num_active;

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();
#ifdef USE_ADAPTIVE_SOFTENING
#pragma omp target enter data map(                                            \
        to : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas],   \
            d_gas_m[0 : n_gas], d_gas_h[0 : n_gas], d_gas_zeta[0 : n_gas],    \
            d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas],    \
            d_active_idx[0 : n_active], d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], \
            d_dm_pz[0 : n_dm], d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes])
#else
#pragma omp target enter data map(                                            \
        to : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas],   \
            d_gas_m[0 : n_gas], d_gas_h[0 : n_gas], d_gas_ax[0 : n_gas],      \
            d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas],                         \
            d_active_idx[0 : n_active], d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], \
            d_dm_pz[0 : n_dm], d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes])
#endif
    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t k = 0; k < n_active; ++k) {
        size_t i = d_active_idx[k];

        double p1_x = d_gas_px[i], p1_y = d_gas_py[i], p1_z = d_gas_pz[i];
#ifdef USE_ADAPTIVE_SOFTENING
        double m_gas = d_gas_m[i];
        double h_i = d_gas_h[i];
        double zeta_i = d_gas_zeta[i];
#endif
        double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = d_dm_bvh[node_idx];

            // Spatial Culling
            double aabb_dist_sq =
                min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                     domain_size.x) +
                min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                     domain_size.y) +
                min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                     domain_size.z);

            if (aabb_dist_sq > search_sq) continue;

            if (node.particle_idx != -1) {
                // PP interaction at the leaf
                int j = node.particle_idx;

                double dx = p1_x - d_dm_px[j];
                if (dx > 0.5 * domain_size.x)
                    dx -= domain_size.x;
                else if (dx < -0.5 * domain_size.x)
                    dx += domain_size.x;

                double dy = p1_y - d_dm_py[j];
                if (dy > 0.5 * domain_size.y)
                    dy -= domain_size.y;
                else if (dy < -0.5 * domain_size.y)
                    dy += domain_size.y;

                double dz = p1_z - d_dm_pz[j];
                if (dz > 0.5 * domain_size.z)
                    dz -= domain_size.z;
                else if (dz < -0.5 * domain_size.z)
                    dz += domain_size.z;

                // Flip back to (j - i) vector direction
                dx = -dx;
                dy = -dy;
                dz = -dz;
                double dist_sq = dx * dx + dy * dy + dz * dz;

                if (use_pm && dist_sq > cutoff_sq) continue;
                double r = std::sqrt(dist_sq + 1e-24);
                double m_j = d_dm_m[j];

#ifdef USE_ADAPTIVE_SOFTENING
                // Gas kernel derivatives
                double dphi_dr_i, dW_dr_i, dphi_dr_j, dummy_dphi_dh;
                Kernels::adaptive_gravity_terms(r, h_i, dphi_dr_i, dW_dr_i);
                // DM kernel derivatives
                Kernels::gravity_derivatives(r, base_soft, dphi_dr_j,
                                             dummy_dphi_dh);
                // Force magnitude
                double force_mag_over_r =
                    (G / 2.0) *
                    ((dphi_dr_i + dphi_dr_j) + (zeta_i * dW_dr_i) / m_gas) / r;
#else
                // TreePM PP Force
                double dphi_dr, dummy_dphi_dh;
                Kernels::gravity_derivatives(r, base_soft, dphi_dr,
                                             dummy_dphi_dh);
                double force_mag_over_r = G * dphi_dr / r;
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
#ifdef USE_ADAPTIVE_SOFTENING
#pragma omp target exit data map(from : d_gas_ax[0 : n_gas],                   \
                                     d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas]) \
    map(delete : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas],                     \
            d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas], d_gas_h[0 : n_gas],       \
            d_gas_zeta[0 : n_gas], d_active_idx[0 : n_active],                 \
            d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],           \
            d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes])
#else
#pragma omp target exit data map(from : d_gas_ax[0 : n_gas],                   \
                                     d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas]) \
    map(delete : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas],                     \
            d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas], d_gas_h[0 : n_gas],       \
            d_active_idx[0 : n_active], d_dm_px[0 : n_dm], d_dm_py[0 : n_dm],  \
            d_dm_pz[0 : n_dm], d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes])
#endif
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

void compute_and_add_dm_pp_forces(ParticleSystem& dm, double a,
                                  const Config& config, Diagnostics& diag) {
    if (dm.num_particles == 0) return;

    const size_t n_parts = dm.num_particles;
    const Vec3& domain_size = config.domain_size;
    // Pre-scale G so output forces are comoving accelerations
    const double G = config.G / (a * a * a);
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t num_nodes = 2 * n_parts - 1;
    const double base_soft = config.softening;

    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();

    // Extract pointers
    double* d_px = dm.pos_x.data();
    double* d_py = dm.pos_y.data();
    double* d_pz = dm.pos_z.data();
    double* d_m = dm.mass.data();
    double* d_ax = dm.acc_x.data();
    double* d_ay = dm.acc_y.data();
    double* d_az = dm.acc_z.data();
    const BVHNode* d_bvh_nodes = dm.bvh_nodes.data();
    int* d_active_idx = dm.active_indices.data();
    size_t n_active = dm.num_active;

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                           \
        to : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts],        \
            d_m[0 : n_parts], d_bvh_nodes[0 : num_nodes], d_ax[0 : n_parts], \
            d_ay[0 : n_parts], d_az[0 : n_parts], d_active_idx[0 : n_active])

    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t k = 0; k < n_active; ++k) {
        size_t i = d_active_idx[k];

        double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
        double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = d_bvh_nodes[node_idx];

            double aabb_dist_sq =
                min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                     domain_size.x) +
                min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                     domain_size.y) +
                min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                     domain_size.z);

            if (aabb_dist_sq > search_sq) continue;

            if (node.particle_idx != -1) {
                int j = node.particle_idx;
                if (static_cast<size_t>(j) == i) continue;

                double dx = p1_x - d_px[j];
                if (dx > 0.5 * domain_size.x)
                    dx -= domain_size.x;
                else if (dx < -0.5 * domain_size.x)
                    dx += domain_size.x;

                double dy = p1_y - d_py[j];
                if (dy > 0.5 * domain_size.y)
                    dy -= domain_size.y;
                else if (dy < -0.5 * domain_size.y)
                    dy += domain_size.y;

                double dz = p1_z - d_pz[j];
                if (dz > 0.5 * domain_size.z)
                    dz -= domain_size.z;
                else if (dz < -0.5 * domain_size.z)
                    dz += domain_size.z;

                dx = -dx;
                dy = -dy;
                dz = -dz;

                double dist_sq = dx * dx + dy * dy + dz * dz;

                if (use_pm && dist_sq > cutoff_sq) continue;

                double r = std::sqrt(dist_sq + 1e-24);
                double force_mag_over_r = 0.0;

                // Cubic Spline Gravity
                double dphi_dr, dummy_dphi_dh;
#ifdef USE_ADAPTIVE_SOFTENING
                Kernels::gravity_derivatives(r, base_soft, dphi_dr,
                                             dummy_dphi_dh);
#else
                Kernels::gravity_derivatives(r, base_soft, dphi_dr,
                                             dummy_dphi_dh);
#endif

                // DM-DM has identical kernels and no zeta, so (dphi_i +
                // dphi_j)/2 = dphi_dr
                force_mag_over_r = G * dphi_dr / r;

                if (use_pm) {
                    // Taper evaluated safely using true geometric distance 'r'
                    double r_scaled = r / (2.0 * r_s);
                    double taper = std::erfc(r_scaled) +
                                   (r / (std::sqrt(M_PI) * r_s)) *
                                       std::exp(-r_scaled * r_scaled);
                    force_mag_over_r *= taper;
                }

                local_acc_x += force_mag_over_r * node.mass * dx;
                local_acc_y += force_mag_over_r * node.mass * dy;
                local_acc_z += force_mag_over_r * node.mass * dz;
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

#pragma omp target exit data map(from : d_ax[0 : n_parts], d_ay[0 : n_parts], \
                                     d_az[0 : n_parts])                       \
    map(delete : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts],     \
            d_m[0 : n_parts], d_bvh_nodes[0 : num_nodes],                     \
            d_active_idx[0 : n_active])

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

static void compute_dm_cross_pp_forces(ParticleSystem& dm,
                                       const GasParticleSystem& gas, double a,
                                       const Config& config,
                                       Diagnostics& diag) {
    if (dm.num_particles == 0 || gas.num_particles == 0) return;

    const Vec3& domain_size = config.domain_size;

    // Pre-scale G so all output forces are comoving accelerations
    const double G = config.G / (a * a * a);
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;

    const size_t n_dm = dm.num_particles;
    const size_t n_gas = gas.num_particles;
    const size_t gas_num_nodes = 2 * n_gas - 1;

    // Search radius is the PM cutoff. If not using PM, it's infinite
    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();
    const double base_soft = config.softening;

    // Extract DM pointers
    double* d_dm_px = dm.pos_x.data();
    double* d_dm_py = dm.pos_y.data();
    double* d_dm_pz = dm.pos_z.data();
    uint8_t* d_dm_is_active = dm.is_active.data();
    double* d_dm_ax = dm.acc_x.data();
    double* d_dm_ay = dm.acc_y.data();
    double* d_dm_az = dm.acc_z.data();

    // Extract Gas pointers (read-only)
    const double* d_gas_px = gas.pos_x.data();
    const double* d_gas_py = gas.pos_y.data();
    const double* d_gas_pz = gas.pos_z.data();
    const double* d_gas_m = gas.mass.data();
    const double* d_gas_h = gas.h.data();
#ifdef USE_ADAPTIVE_SOFTENING
    const double* d_gas_zeta = gas.zeta.data();
#endif
    const BVHNode* d_gas_bvh = gas.bvh_nodes.data();
    int* d_active_idx = dm.active_indices.data();  // Extract indirection array
    size_t n_active = dm.num_active;

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();
#ifdef USE_ADAPTIVE_SOFTENING
#pragma omp target enter data map(                                            \
        to : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],         \
            d_active_idx[0 : n_active], d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], \
            d_dm_az[0 : n_dm], d_gas_px[0 : n_gas], d_gas_py[0 : n_gas],      \
            d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas], d_gas_h[0 : n_gas],      \
            d_gas_zeta[0 : n_gas], d_gas_bvh[0 : gas_num_nodes])
#else
#pragma omp target enter data map(                                            \
        to : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],         \
            d_active_idx[0 : n_active], d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], \
            d_dm_az[0 : n_dm], d_gas_px[0 : n_gas], d_gas_py[0 : n_gas],      \
            d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas], d_gas_h[0 : n_gas],      \
            d_gas_bvh[0 : gas_num_nodes])
#endif
    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

// Swap the OpenMP execution pragma based on the compile-time target
#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t k = 0; k < n_active; ++k) {
        size_t i = d_active_idx[k];

        double p1_x = d_dm_px[i], p1_y = d_dm_py[i], p1_z = d_dm_pz[i];
        double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;  // Push root of Gas tree

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = d_gas_bvh[node_idx];

            // Spatial Culling against the gas node's AABB
            double aabb_dist_sq =
                min_periodic_dist_sq(p1_x, node.bbox.min_x, node.bbox.max_x,
                                     domain_size.x) +
                min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                     domain_size.y) +
                min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                     domain_size.z);

            if (aabb_dist_sq > search_sq) continue;  // Prune branch

            if (node.particle_idx != -1) {
                // Exact P-P interaction at the leaf
                int j = node.particle_idx;

                double dx = p1_x - d_gas_px[j];
                if (dx > 0.5 * domain_size.x)
                    dx -= domain_size.x;
                else if (dx < -0.5 * domain_size.x)
                    dx += domain_size.x;

                double dy = p1_y - d_gas_py[j];
                if (dy > 0.5 * domain_size.y)
                    dy -= domain_size.y;
                else if (dy < -0.5 * domain_size.y)
                    dy += domain_size.y;

                double dz = p1_z - d_gas_pz[j];
                if (dz > 0.5 * domain_size.z)
                    dz -= domain_size.z;
                else if (dz < -0.5 * domain_size.z)
                    dz += domain_size.z;

                // Flip to (gas - dm) direction so force pulls DM toward Gas
                dx = -dx;
                dy = -dy;
                dz = -dz;

                double dist_sq = dx * dx + dy * dy + dz * dz;

                if (use_pm && dist_sq > cutoff_sq) continue;

                double r = std::sqrt(dist_sq + 1e-24);
                double m_gas = d_gas_m[j];
#ifdef USE_ADAPTIVE_SOFTENING
                double h_gas = d_gas_h[j];
                double zeta_gas = d_gas_zeta[j];

                double dphi_dr_i, dummy_dphi_dh, dphi_dr_j, dW_dr_j;

                // DM kernel derivatives (particle i)
                Kernels::gravity_derivatives(r, base_soft, dphi_dr_i,
                                             dummy_dphi_dh);

                // Gas kernel derivatives (particle j)
                Kernels::adaptive_gravity_terms(r, h_gas, dphi_dr_j, dW_dr_j);

                // Force magnitude calculation (only Gas contributes zeta, DM
                // softening is fixed)
                double force_mag_over_r =
                    (G / 2.0) *
                    ((dphi_dr_i + dphi_dr_j) + (zeta_gas * dW_dr_j) / m_gas) /
                    r;

#else
                // TreePM PP Force
                double dphi_dr, dummy_dphi_dh;
                Kernels::gravity_derivatives(r, base_soft, dphi_dr,
                                             dummy_dphi_dh);
                double force_mag_over_r = G * dphi_dr / r;
#endif

                if (use_pm) {
                    double r_scaled = r / (2.0 * r_s);
                    double taper = std::erfc(r_scaled) +
                                   (r / (std::sqrt(M_PI) * r_s)) *
                                       std::exp(-r_scaled * r_scaled);
                    force_mag_over_r *= taper;
                }

                // Pull on the Dark Matter particle
                local_acc_x += force_mag_over_r * m_gas * dx;
                local_acc_y += force_mag_over_r * m_gas * dy;
                local_acc_z += force_mag_over_r * m_gas * dz;

            } else {
                // Push children
                stack[stack_ptr++] = node.left_child;
                stack[stack_ptr++] = node.right_child;
            }
        }

        // Apply accumulated local force to the DM arrays
        d_dm_ax[i] += local_acc_x;
        d_dm_ay[i] += local_acc_y;
        d_dm_az[i] += local_acc_z;
    }

#ifdef USE_GPU
    auto end_compute = std::chrono::high_resolution_clock::now();
    auto start_return = std::chrono::high_resolution_clock::now();
#ifdef USE_ADAPTIVE_SOFTENING
#pragma omp target exit data map(from : d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], \
                                     d_dm_az[0 : n_dm])                       \
    map(delete : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],     \
            d_active_idx[0 : n_active], d_gas_px[0 : n_gas],                  \
            d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas],     \
            d_gas_h[0 : n_gas], d_gas_zeta[0 : n_gas],                        \
            d_gas_bvh[0 : gas_num_nodes])
#else
#pragma omp target exit data map(from : d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], \
                                     d_dm_az[0 : n_dm])                       \
    map(delete : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],     \
            d_active_idx[0 : n_active], d_gas_px[0 : n_gas],                  \
            d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas],     \
            d_gas_h[0 : n_gas], d_gas_bvh[0 : gas_num_nodes])
#endif
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

void compute_gas_dm_cross_pp_forces(ParticleSystem& dm, GasParticleSystem& gas,
                                    double a, const Config& config,
                                    Diagnostics& diag) {
    // Gas-DM and DM-Gas
    compute_gas_cross_pp_forces(gas, a, dm, config, diag);

    // Gas pulling on DM
    compute_dm_cross_pp_forces(dm, gas, a, config, diag);
}

double get_gravity_timestep(const GasParticleSystem& gas,
                            const Config& config) {
    if (gas.num_particles == 0) return std::numeric_limits<double>::infinity();
    if (config.individual_particle_timesteps) {
        return std::numeric_limits<double>::infinity();
    } else {
        double epsilon = config.softening;
        double a_max = std::sqrt(gas.max_accel_sq);
        double dt_grav = std::sqrt(epsilon / a_max);
        return dt_grav * config.gravity_accuracy_eta;
    }
}

double get_gravity_timestep(const ParticleSystem& dm, const Config& config) {
    if (dm.num_particles == 0) return std::numeric_limits<double>::infinity();
    if (config.individual_particle_timesteps) {
        // Bootstrap: On the first step, request safe step to initialize
        // forces
        if (dm.global_time == 0.0) {
            return config.fixed_dt / pow(2, 16);
        }

        // Return the time until the NEXT active DM particle finishes its step
        double next_time = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < dm.num_particles; ++i) {
            // Only look at future bins ignoring float precision drift matches
            if (dm.t_end[i] > dm.global_time + 1e-10) {
                next_time = std::min(next_time, dm.t_end[i]);
            }
        }

        double dt_next = next_time - dm.global_time;
        return (dt_next > 1e-10) ? dt_next
                                 : std::numeric_limits<double>::infinity();
    } else {
        double epsilon = config.softening;
        double a_max = std::sqrt(dm.max_accel_sq);
        double dt_grav = std::sqrt(epsilon / a_max);

        return dt_grav * config.gravity_accuracy_eta;
    }
}

}  // namespace GravitySolver
