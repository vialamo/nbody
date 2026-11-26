#include "gravity.h"

#include <chrono>
#include <cmath>

#include "diagnostics.h"
#include "kernels.h"
#include "math_utils.h"
#include "state.h"

namespace {
// Microscopic floor for squared distances to prevent 1/r singularities
// when particles perfectly overlap
constexpr double MIN_R2_EPSILON = 1e-24;

#ifdef USE_GPU
#pragma omp declare target
#endif
inline double compute_pm_taper(double r, double inv_two_r_s,
                               double taper_coeff) {
    double r_scaled = r * inv_two_r_s;
    return std::erfc(r_scaled) +
           (r * taper_coeff) * std::exp(-r_scaled * r_scaled);
}
#ifdef USE_GPU
#pragma omp end declare target
#endif
}  // namespace

#ifdef GRAVITY_GPU_MEMORY
void GravitySolver::prepare_arrays(ParticleSystem& dm, GasParticleSystem& gas,
                                   const Config& config, Diagnostics& diag) {
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas.num_particles;
    size_t n_total = n_dm + n_gas;

    if (n_total == 0) return;

    auto t_start = std::chrono::high_resolution_clock::now();

    // GPU Allocation
    if (n_total > current_gpu_capacity) {
        int dev = omp_get_default_device();
        // Free old GPU memory if capacity is expanding
        if (current_gpu_capacity > 0) {
            omp_target_free(d_u_px, dev);
            omp_target_free(d_t_px, dev);
            omp_target_free(d_u_py, dev);
            omp_target_free(d_t_py, dev);
            omp_target_free(d_u_pz, dev);
            omp_target_free(d_t_pz, dev);
            omp_target_free(d_u_m, dev);
            omp_target_free(d_t_m, dev);
            omp_target_free(d_u_ax, dev);
            omp_target_free(d_u_ay, dev);
            omp_target_free(d_u_az, dev);
            omp_target_free(d_u_orig_idx, dev);
            omp_target_free(d_t_orig_idx, dev);
            omp_target_free(d_u_type, dev);
            omp_target_free(d_t_type, dev);
            omp_target_free(d_t_is_active, dev);
            omp_target_free(d_u_morton, dev);
            omp_target_free(d_u_sorted_idx, dev);
            omp_target_free(d_u_eps, dev);
            omp_target_free(d_t_eps, dev);
            omp_target_free(d_u_gradh_fac, dev);
            omp_target_free(d_t_gradh_fac, dev);
            omp_target_free(d_u_xi, dev);

            // Only use map(delete) for the host-backed active array
            uint8_t* old_is_active = u_is_active.data();
#pragma omp target exit data map( \
        delete : old_is_active[0 : current_gpu_capacity])
        }

        u_is_active.resize(n_total);
        uint8_t* d_is_act = u_is_active.data();
#pragma omp target enter data map(alloc : d_is_act[0 : n_total])

        // Allocate on the GPU
        d_u_px = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_t_px = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_py = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_t_py = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_pz = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_t_pz = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_m = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_t_m = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_ax = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_ay = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_az = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_orig_idx = (int*)omp_target_alloc(n_total * sizeof(int), dev);
        d_t_orig_idx = (int*)omp_target_alloc(n_total * sizeof(int), dev);
        d_u_type = (uint8_t*)omp_target_alloc(n_total * sizeof(uint8_t), dev);
        d_t_type = (uint8_t*)omp_target_alloc(n_total * sizeof(uint8_t), dev);
        d_t_is_active =
            (uint8_t*)omp_target_alloc(n_total * sizeof(uint8_t), dev);
        d_u_morton =
            (uint64_t*)omp_target_alloc(n_total * sizeof(uint64_t), dev);
        d_u_sorted_idx = (int*)omp_target_alloc(n_total * sizeof(int), dev);

        d_u_eps = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_t_eps = (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_gradh_fac =
            (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_t_gradh_fac =
            (double*)omp_target_alloc(n_total * sizeof(double), dev);
        d_u_xi = (double*)omp_target_alloc(n_total * sizeof(double), dev);

        current_gpu_capacity = n_total;
    }

    // Extract all class members to local stack variables for OpenMP
    double *loc_u_px = d_u_px, *loc_u_py = d_u_py, *loc_u_pz = d_u_pz,
           *loc_u_m = d_u_m;
    double *loc_u_ax = d_u_ax, *loc_u_ay = d_u_ay, *loc_u_az = d_u_az;
    int *loc_u_orig_idx = d_u_orig_idx, *loc_u_sorted_idx = d_u_sorted_idx;
    uint8_t *loc_u_type = d_u_type, *loc_u_is_active = u_is_active.data();

    double *loc_t_px = d_t_px, *loc_t_py = d_t_py, *loc_t_pz = d_t_pz,
           *loc_t_m = d_t_m;
    int* loc_t_orig_idx = d_t_orig_idx;
    uint8_t *loc_t_type = d_t_type, *loc_t_is_active = d_t_is_active;

    double *loc_u_eps = d_u_eps, *loc_t_eps = d_t_eps;
    double *loc_u_gradh_fac = d_u_gradh_fac, *loc_t_gradh_fac = d_t_gradh_fac;

    // Zero out the accelerations on the GPU for the current step
#pragma omp target teams distribute parallel for is_device_ptr( \
        loc_u_ax, loc_u_ay, loc_u_az)
    for (size_t i = 0; i < n_total; ++i) {
        loc_u_ax[i] = 0.0;
        loc_u_ay[i] = 0.0;
        loc_u_az[i] = 0.0;
    }

    const bool use_adaptive = config.use_adaptive_softening;
    const double eps_min = softening_length;
    const double eps_max = config.use_PM
                               ? std::sqrt(config.cutoff_radius_squared)
                               : std::numeric_limits<double>::infinity();

    // On-Device Concatenation
    // Map the CPU particle arrays to the device temporarily to build the
    // unified arrays
    double* d_dm_px = dm.pos_x.data();
    double* d_dm_py = dm.pos_y.data();
    double* d_dm_pz = dm.pos_z.data();
    double* d_dm_m = dm.mass.data();
    uint8_t* d_dm_is_active = dm.is_active.data();

#pragma omp target teams distribute parallel for is_device_ptr(            \
        loc_u_px, loc_u_py, loc_u_pz, loc_u_m, loc_u_orig_idx, loc_u_type, \
            loc_u_eps, loc_u_gradh_fac)                                    \
    map(to : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],      \
            d_dm_m[0 : n_dm], d_dm_is_active[0 : n_dm])
    for (size_t i = 0; i < n_dm; ++i) {
        loc_u_px[i] = d_dm_px[i];
        loc_u_py[i] = d_dm_py[i];
        loc_u_pz[i] = d_dm_pz[i];
        loc_u_m[i] = d_dm_m[i];
        loc_u_orig_idx[i] = i;
        loc_u_type[i] = 0;  // Type 0: DM
        loc_u_is_active[i] = d_dm_is_active[i];
        loc_u_eps[i] = eps_min;
        loc_u_gradh_fac[i] = 0.0;
    }

    if (n_gas > 0) {
        double* d_gas_px = gas.pos_x.data();
        double* d_gas_py = gas.pos_y.data();
        double* d_gas_pz = gas.pos_z.data();
        double* d_gas_m = gas.mass.data();
        uint8_t* d_gas_is_active = gas.is_active.data();
        double* d_gas_h = gas.h.data();
        double* d_gas_rho = gas.rho.data();
        double* d_gas_omega = gas.omega.data();

#pragma omp target teams distribute parallel for is_device_ptr(             \
        loc_u_px, loc_u_py, loc_u_pz, loc_u_m, loc_u_orig_idx, loc_u_type,  \
            loc_u_eps, loc_u_gradh_fac)                                     \
    map(to : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], \
            d_gas_m[0 : n_gas], d_gas_is_active[0 : n_gas],                 \
            d_gas_h[0 : n_gas], d_gas_rho[0 : n_gas], d_gas_omega[0 : n_gas])
        for (size_t i = 0; i < n_gas; ++i) {
            size_t offset = n_dm + i;
            loc_u_px[offset] = d_gas_px[i];
            loc_u_py[offset] = d_gas_py[i];
            loc_u_pz[offset] = d_gas_pz[i];
            loc_u_m[offset] = d_gas_m[i];
            loc_u_orig_idx[offset] = i;
            loc_u_type[offset] = 1;  // Type 1: Gas
            loc_u_is_active[offset] = d_gas_is_active[i];
            double eps_val = eps_min;
            double gradh_fac = 0.0;
            if (use_adaptive) {
                double h_i = d_gas_h[i];
                if (h_i <= eps_min) {
                    eps_val = eps_min;
                } else if (h_i >= eps_max) {
                    eps_val = eps_max;
                } else {
                    eps_val = h_i;
                    gradh_fac = (d_gas_m[i] * h_i) /
                                (3.0 * d_gas_rho[i] * d_gas_omega[i]);
                }
            }
            loc_u_eps[offset] = eps_val;
            loc_u_gradh_fac[offset] = gradh_fac;
        }
    }

    auto t_transfer = std::chrono::high_resolution_clock::now();

    // Run Morton generation
    LBVH::compute_morton_and_sort_indices_gpu(
        lbvh_ws, n_total, config.domain_size, loc_u_px, loc_u_py, loc_u_pz,
        d_u_morton, d_u_sorted_idx);

// Permute into temporary device arrays
#pragma omp target teams distribute parallel for is_device_ptr(                \
        loc_u_px, loc_u_py, loc_u_pz, loc_u_m, loc_u_orig_idx, loc_u_type,     \
            loc_u_sorted_idx, loc_t_px, loc_t_py, loc_t_pz, loc_t_m,           \
            loc_t_orig_idx, loc_t_type, loc_t_is_active, loc_u_eps, loc_t_eps, \
            loc_u_gradh_fac, loc_t_gradh_fac)
    for (size_t i = 0; i < n_total; ++i) {
        int idx = loc_u_sorted_idx[i];
        loc_t_px[i] = loc_u_px[idx];
        loc_t_py[i] = loc_u_py[idx];
        loc_t_pz[i] = loc_u_pz[idx];
        loc_t_m[i] = loc_u_m[idx];
        loc_t_orig_idx[i] = loc_u_orig_idx[idx];
        loc_t_type[i] = loc_u_type[idx];
        loc_t_is_active[i] = loc_u_is_active[idx];
        loc_t_eps[i] = loc_u_eps[idx];
        loc_t_gradh_fac[i] = loc_u_gradh_fac[idx];
    }

// Copy back to the persistent device arrays
#pragma omp target teams distribute parallel for is_device_ptr(                \
        loc_u_px, loc_u_py, loc_u_pz, loc_u_m, loc_u_orig_idx, loc_u_type,     \
            loc_t_px, loc_t_py, loc_t_pz, loc_t_m, loc_t_orig_idx, loc_t_type, \
            loc_t_is_active, loc_u_eps, loc_t_eps, loc_u_gradh_fac,            \
            loc_t_gradh_fac)
    for (size_t i = 0; i < n_total; ++i) {
        loc_u_px[i] = loc_t_px[i];
        loc_u_py[i] = loc_t_py[i];
        loc_u_pz[i] = loc_t_pz[i];
        loc_u_m[i] = loc_t_m[i];
        loc_u_orig_idx[i] = loc_t_orig_idx[i];
        loc_u_type[i] = loc_t_type[i];
        loc_u_is_active[i] = loc_t_is_active[i];
        loc_u_eps[i] = loc_t_eps[i];
        loc_u_gradh_fac[i] = loc_t_gradh_fac[i];
    }

    // Persistent GPU allocation for the tree
    size_t required_nodes = (n_total > 0) ? (2 * n_total - 1) : 0;

    if (required_nodes > current_bvh_capacity) {
        int dev = omp_get_default_device();
        if (current_bvh_capacity > 0) {
            omp_target_free(d_u_bvh, dev);
        }
        d_u_bvh =
            (BVHNode*)omp_target_alloc(required_nodes * sizeof(BVHNode), dev);
        current_bvh_capacity = required_nodes;
    }

    // Call the device-side builder using the GPU-resident arrays
    LBVH::build_topology_and_aggregate_gpu(lbvh_ws, n_total, loc_u_px, loc_u_py,
                                           loc_u_pz, loc_u_m, d_u_morton,
                                           d_u_bvh);

    uint8_t* d_u_is_active_end = u_is_active.data();
#pragma omp target update from(d_u_is_active_end[0 : n_total])

    // Rebuild the dense array sequentially to preserve Morton spatial
    // locality
    u_active_indices.clear();
    for (size_t i = 0; i < n_total; ++i) {
        if (u_is_active[i]) u_active_indices.push_back(i);
    }

    auto t_tree = std::chrono::high_resolution_clock::now();

    // Log the times to diagnostics
    diag.add_prof_time(GPUProfRegion::Transf,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           t_transfer - t_start)
                           .count());
    diag.add_prof_time(GPUProfRegion::Tree,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           t_tree - t_transfer)
                           .count());
}

#endif

// Map accelerations back to native structures
void GravitySolver::map_results(ParticleSystem& dm, GasParticleSystem& gas,
                                Diagnostics& diag) {
#ifdef GRAVITY_GPU_MEMORY
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas.num_particles;
    size_t n_active = u_active_indices.size();

    if (n_active == 0) return;

    int* d_active_idx = u_active_indices.data();

    double *loc_u_ax = d_u_ax, *loc_u_ay = d_u_ay, *loc_u_az = d_u_az;
    int* loc_u_orig_idx = d_u_orig_idx;
    uint8_t* loc_u_type = d_u_type;

    auto start_return = std::chrono::high_resolution_clock::now();
#pragma omp target data map(to : d_active_idx[0 : n_active])
    {
        if (n_dm > 0) {
            double* d_dm_ax = dm.acc_x.data();
            double* d_dm_ay = dm.acc_y.data();
            double* d_dm_az = dm.acc_z.data();

#pragma omp target teams distribute parallel for map(                     \
        tofrom : d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm]) \
    is_device_ptr(loc_u_ax, loc_u_ay, loc_u_az, loc_u_orig_idx, loc_u_type)
            for (size_t k = 0; k < n_active; ++k) {
                size_t i = d_active_idx[k];
                if (loc_u_type[i] == 0) {
                    int orig = loc_u_orig_idx[i];
                    d_dm_ax[orig] += loc_u_ax[i];
                    d_dm_ay[orig] += loc_u_ay[i];
                    d_dm_az[orig] += loc_u_az[i];
                }
            }
        }

        if (n_gas > 0) {
            double* d_gas_ax = gas.acc_x.data();
            double* d_gas_ay = gas.acc_y.data();
            double* d_gas_az = gas.acc_z.data();

#pragma omp target teams distribute parallel for map(      \
        tofrom : d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas], \
            d_gas_az[0 : n_gas])                           \
    is_device_ptr(loc_u_ax, loc_u_ay, loc_u_az, loc_u_orig_idx, loc_u_type)
            for (size_t k = 0; k < n_active; ++k) {
                size_t i = d_active_idx[k];
                if (loc_u_type[i] == 1) {
                    int orig = loc_u_orig_idx[i];
                    d_gas_ax[orig] += loc_u_ax[i];
                    d_gas_ay[orig] += loc_u_ay[i];
                    d_gas_az[orig] += loc_u_az[i];
                }
            }
        }
    }

    auto end_return = std::chrono::high_resolution_clock::now();
    diag.add_prof_time(GPUProfRegion::Ret,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_return - start_return)
                           .count());
#endif
}

#ifdef GRAVITY_GPU_MEMORY
void GravitySolver::add_pp_accel(ParticleSystem& dm, GasParticleSystem& gas,
                                 double a, const Config& config,
                                 Diagnostics& diag) {
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas.num_particles;
    size_t n_total = n_dm + n_gas;

    if (n_total == 0) return;

    prepare_arrays(dm, gas, config, diag);

    size_t n_active = u_active_indices.size();
    if (n_active == 0) return;

    // Shared Force Constants
    const Vec3& domain_size = config.domain_size;
    const double G = config.G / (a * a * a);
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const double epsilon = softening_length;
    const double search_sq = use_pm ? config.cutoff_radius_squared
                                    : std::numeric_limits<double>::infinity();
    const double inv_two_r_s = 1.0 / (2.0 * r_s);
    const double taper_coeff = 1.0 / (std::sqrt(M_PI) * r_s);
    const bool use_adaptive = config.use_adaptive_softening;

    // Raw pointers for GPU Offload
    int* d_active_idx = u_active_indices.data();
    double* d_px = d_u_px;
    double* d_py = d_u_py;
    double* d_pz = d_u_pz;
    double* d_m = d_u_m;
    double* d_ax = d_u_ax;
    double* d_ay = d_u_ay;
    double* d_az = d_u_az;
    BVHNode* d_bvh = d_u_bvh;
    double* d_eps = d_u_eps;
    double* d_gradh_fac = d_u_gradh_fac;
    double* d_xi = d_u_xi;
    uint8_t* d_type = d_u_type;

    auto start_compute = std::chrono::high_resolution_clock::now();

    if (use_adaptive) {
#pragma omp target teams distribute parallel for is_device_ptr( \
        d_px, d_py, d_pz, d_m, d_eps, d_gradh_fac, d_xi, d_bvh)
        for (size_t i = 0; i < n_total; ++i) {
            double fac = d_gradh_fac[i];
            if (fac == 0.0) {
                d_xi[i] = 0.0;
                continue;
            }

            double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
            double eps_i = d_eps[i];
            double eps_sq = eps_i * eps_i;
            double sum_dphi_dh = 0.0;

            LBVH::traverse(
                d_bvh, d_px, d_py, d_pz, p1_x, p1_y, p1_z, eps_sq, domain_size,
                [=, &sum_dphi_dh](int j, double dx, double dy, double dz,
                                  double r2) {
                    if (r2 > 0.0 && r2 < eps_sq) {
                        double r = std::sqrt(r2 + MIN_R2_EPSILON);
                        double term =
                            d_m[j] * Kernels::gravity_dphi_dh(r, eps_i);
                        if (use_pm) {
                            term *= std::erfc(r * inv_two_r_s);
                        }
                        sum_dphi_dh += term;
                    }
                });

            d_xi[i] = fac * sum_dphi_dh;
        }
    }

#pragma omp target teams distribute parallel for map(                          \
        to : d_active_idx[0 : n_active])                                       \
    is_device_ptr(d_px, d_py, d_pz, d_m, d_eps, d_bvh, d_ax, d_ay, d_az, d_xi, \
                      d_type)
    for (size_t k = 0; k < n_active; ++k) {
        size_t i = d_active_idx[k];

        double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
        double eps_i = d_eps[i];
        double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
        double m_i = d_m[i];
        double xi_i = d_xi[i];
        uint8_t type_i = d_type[i];

        LBVH::traverse(
            d_bvh, d_px, d_py, d_pz, p1_x, p1_y, p1_z, search_sq, domain_size,
            [=, &local_acc_x, &local_acc_y, &local_acc_z](
                int j, double dx, double dy, double dz, double r2) {
                if (r2 < search_sq) {
                    double r = std::sqrt(r2 + MIN_R2_EPSILON);
                    double m_j = d_m[j];
                    double eps_j = d_eps[j];

                    double dphi_dr;
                    if (use_adaptive) {
                        dphi_dr =
                            0.5 * (Kernels::gravity_derivatives(r, eps_i) +
                                   Kernels::gravity_derivatives(r, eps_j));
                    } else {
                        dphi_dr = Kernels::gravity_derivatives(r, epsilon);
                    }

                    double force_mag_over_r = G * dphi_dr / r;

                    if (use_pm) {
                        force_mag_over_r *=
                            compute_pm_taper(r, inv_two_r_s, taper_coeff);
                    }

                    local_acc_x += force_mag_over_r * m_j * dx;
                    local_acc_y += force_mag_over_r * m_j * dy;
                    local_acc_z += force_mag_over_r * m_j * dz;

                    // Price & Monaghan (2007) adaptive softening correction
                    // (gas-gas pairs only; dW/dr <= 0 so leading minus makes it
                    // attractive)
                    if (use_adaptive && type_i == 1 && d_type[j] == 1) {
                        double xi_j = d_xi[j];
                        if (xi_i != 0.0 || xi_j != 0.0) {
                            double dW_dr_i =
                                (xi_i != 0.0)
                                    ? Kernels::cubic_spline_gradient(r, eps_i)
                                    : 0.0;
                            double dW_dr_j =
                                (xi_j != 0.0)
                                    ? Kernels::cubic_spline_gradient(r, eps_j)
                                    : 0.0;

                            double corr_over_r =
                                -0.5 * (G / m_i) *
                                (m_i * xi_i * dW_dr_i + m_j * xi_j * dW_dr_j) /
                                r;

                            local_acc_x += corr_over_r * dx;
                            local_acc_y += corr_over_r * dy;
                            local_acc_z += corr_over_r * dz;
                        }
                    }
                }
            });
        d_ax[i] += local_acc_x;
        d_ay[i] += local_acc_y;
        d_az[i] += local_acc_z;
    }

    auto end_compute = std::chrono::high_resolution_clock::now();
    diag.add_prof_time(GPUProfRegion::Compute,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_compute - start_compute)
                           .count());
}

#else
void GravitySolver::add_pp_accel(ParticleSystem& dm, GasParticleSystem& gas,
                                 double a, const Config& config,
                                 Diagnostics& diag) {
    if (dm.num_particles == 0 && gas.num_particles == 0) return;

    const Vec3& domain_size = config.domain_size;
    const double G = config.G / (a * a * a);
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const double epsilon = softening_length;
    const bool use_adaptive = config.use_adaptive_softening;
    const double eps_min = softening_length;
    const double eps_max = use_pm ? std::sqrt(config.cutoff_radius_squared)
                                  : std::numeric_limits<double>::infinity();
    const double inv_two_r_s = 1.0 / (2.0 * r_s);
    const double taper_coeff = 1.0 / (std::sqrt(M_PI) * r_s);

    const double search_sq = use_pm ? config.cutoff_radius_squared
                                    : std::numeric_limits<double>::infinity();

    const size_t n_dm = dm.num_particles;
    const size_t dm_num_nodes = (n_dm > 0) ? (2 * n_dm - 1) : 0;
    const size_t n_dm_active = dm.get_num_active();

    const size_t n_gas = gas.num_particles;
    const size_t gas_num_nodes = (n_gas > 0) ? (2 * n_gas - 1) : 0;
    const size_t n_gas_active = gas.get_num_active();

    // Extract DM pointers
    double* d_dm_px = dm.pos_x.data();
    double* d_dm_py = dm.pos_y.data();
    double* d_dm_pz = dm.pos_z.data();
    double* d_dm_m = dm.mass.data();
    double* d_dm_ax = dm.acc_x.data();
    double* d_dm_ay = dm.acc_y.data();
    double* d_dm_az = dm.acc_z.data();
    BVHNode* d_dm_bvh = dm.bvh_nodes.data();
    int* d_dm_active_idx = dm.active_indices.data();

    // Extract Gas pointers
    double* d_gas_px = gas.pos_x.data();
    double* d_gas_py = gas.pos_y.data();
    double* d_gas_pz = gas.pos_z.data();
    double* d_gas_m = gas.mass.data();
    double* d_gas_ax = gas.acc_x.data();
    double* d_gas_ay = gas.acc_y.data();
    double* d_gas_az = gas.acc_z.data();
    BVHNode* d_gas_bvh = gas.bvh_nodes.data();
    int* d_gas_active_idx = gas.active_indices.data();
    double* d_gas_h = gas.h.data();

    std::vector<double> gas_xi(n_gas, 0.0);
    double* d_gas_xi = gas_xi.data();
    double* d_gas_rho = gas.rho.data();
    double* d_gas_omega = gas.omega.data();

    // Accumulate xi_i for all gas particles (gas + DM neighbors)
    if (use_adaptive && n_gas > 0) {
#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < n_gas; ++i) {
            double h_i = d_gas_h[i];
            if (h_i <= eps_min || h_i >= eps_max) {
                d_gas_xi[i] = 0.0;
                continue;
            }

            double p1_x = d_gas_px[i], p1_y = d_gas_py[i], p1_z = d_gas_pz[i];
            double eps_sq = h_i * h_i;
            double sum_dphi_dh = 0.0;

            // Sum over Gas neighbors
            LBVH::traverse(d_gas_bvh, d_gas_px, d_gas_py, d_gas_pz, p1_x, p1_y,
                           p1_z, eps_sq, domain_size,
                           [=, &sum_dphi_dh](int j, double dx, double dy,
                                             double dz, double r2) {
                               if (r2 > 0.0 && r2 < eps_sq) {
                                   double r = std::sqrt(r2 + MIN_R2_EPSILON);
                                   double term =
                                       d_gas_m[j] *
                                       Kernels::gravity_dphi_dh(r, h_i);
                                   if (use_pm) {
                                       term *= std::erfc(r * inv_two_r_s);
                                   }
                                   sum_dphi_dh += term;
                               }
                           });

            // Sum over DM neighbors
            if (n_dm > 0) {
                LBVH::traverse(
                    d_dm_bvh, d_dm_px, d_dm_py, d_dm_pz, p1_x, p1_y, p1_z,
                    eps_sq, domain_size,
                    [=, &sum_dphi_dh](int j, double dx, double dy, double dz,
                                      double r2) {
                        if (r2 > 0.0 && r2 < eps_sq) {
                            double r = std::sqrt(r2 + MIN_R2_EPSILON);
                            double term =
                                d_dm_m[j] * Kernels::gravity_dphi_dh(r, h_i);
                            if (use_pm) {
                                term *= std::erfc(r * inv_two_r_s);
                            }
                            sum_dphi_dh += term;
                        }
                    });
            }

            double fac =
                (d_gas_m[i] * h_i) / (3.0 * d_gas_rho[i] * d_gas_omega[i]);
            d_gas_xi[i] = fac * sum_dphi_dh;
        }
    }

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                            \
        to : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],         \
            d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes],                     \
            d_dm_active_idx[0 : n_dm_active], d_dm_ax[0 : n_dm],              \
            d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm], d_gas_px[0 : n_gas],        \
            d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas],     \
            d_gas_bvh[0 : gas_num_nodes], d_gas_active_idx[0 : n_gas_active], \
            d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas],    \
            d_gas_h[0 : n_gas], d_gas_xi[0 : n_gas])

    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

    // STAGE 1: Target = DM, Trees = DM + Gas
    if (n_dm > 0) {
#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (size_t k = 0; k < n_dm_active; ++k) {
            size_t i = d_dm_active_idx[k];
            double p1_x = d_dm_px[i], p1_y = d_dm_py[i], p1_z = d_dm_pz[i];
            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
            LBVH::traverse(
                d_dm_bvh, d_dm_px, d_dm_py, d_dm_pz, p1_x, p1_y, p1_z,
                search_sq, domain_size,
                [=, &local_acc_x, &local_acc_y, &local_acc_z](
                    int j, double dx, double dy, double dz, double r2) {
                    if (r2 > search_sq) return;

                    double r = std::sqrt(r2 + MIN_R2_EPSILON);
                    double m_j = d_dm_m[j];

                    double dphi_dr = Kernels::gravity_derivatives(r, epsilon);
                    double force_mag_over_r = G * dphi_dr / r;

                    if (use_pm) {
                        force_mag_over_r *=
                            compute_pm_taper(r, inv_two_r_s, taper_coeff);
                    }

                    local_acc_x += force_mag_over_r * m_j * dx;
                    local_acc_y += force_mag_over_r * m_j * dy;
                    local_acc_z += force_mag_over_r * m_j * dz;
                });

            // Evaluate against Gas Tree
            if (n_gas > 0) {
                LBVH::traverse(
                    d_gas_bvh, d_gas_px, d_gas_py, d_gas_pz, p1_x, p1_y, p1_z,
                    search_sq, domain_size,
                    [=, &local_acc_x, &local_acc_y, &local_acc_z](
                        int j, double dx, double dy, double dz, double r2) {
                        if (r2 > search_sq) return;

                        double r = std::sqrt(r2 + MIN_R2_EPSILON);
                        double m_j = d_gas_m[j];

                        double dphi_dr;
                        if (use_adaptive) {
                            double eps_j = d_gas_h[j];
                            if (eps_j < eps_min) eps_j = eps_min;
                            if (eps_j > eps_max) eps_j = eps_max;
                            dphi_dr =
                                0.5 *
                                (Kernels::gravity_derivatives(r, epsilon) +
                                 Kernels::gravity_derivatives(r, eps_j));
                        } else {
                            dphi_dr = Kernels::gravity_derivatives(r, epsilon);
                        }
                        double force_mag_over_r = G * dphi_dr / r;

                        if (use_pm) {
                            force_mag_over_r *=
                                compute_pm_taper(r, inv_two_r_s, taper_coeff);
                        }

                        local_acc_x += force_mag_over_r * m_j * dx;
                        local_acc_y += force_mag_over_r * m_j * dy;
                        local_acc_z += force_mag_over_r * m_j * dz;
                    });
            }

            d_dm_ax[i] += local_acc_x;
            d_dm_ay[i] += local_acc_y;
            d_dm_az[i] += local_acc_z;
        }
    }

    // STAGE 2: Target = Gas, Trees = Gas + DM
    if (n_gas > 0) {
#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (size_t k = 0; k < n_gas_active; ++k) {
            size_t i = d_gas_active_idx[k];
            double p1_x = d_gas_px[i], p1_y = d_gas_py[i], p1_z = d_gas_pz[i];
            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
            double eps_i = epsilon;
            if (use_adaptive) {
                eps_i = d_gas_h[i];
                if (eps_i < eps_min) eps_i = eps_min;
                if (eps_i > eps_max) eps_i = eps_max;
            }
            double m_i = d_gas_m[i];
            double xi_i = d_gas_xi[i];

            LBVH::traverse(
                d_gas_bvh, d_gas_px, d_gas_py, d_gas_pz, p1_x, p1_y, p1_z,
                search_sq, domain_size,
                [=, &local_acc_x, &local_acc_y, &local_acc_z](
                    int j, double dx, double dy, double dz, double r2) {
                    if (r2 > search_sq) return;

                    double r = std::sqrt(r2 + MIN_R2_EPSILON);
                    double m_j = d_gas_m[j];
                    double eps_j = d_gas_h[j];

                    double dphi_dr;
                    if (use_adaptive) {
                        if (eps_j < eps_min) eps_j = eps_min;
                        if (eps_j > eps_max) eps_j = eps_max;
                        dphi_dr =
                            0.5 * (Kernels::gravity_derivatives(r, eps_i) +
                                   Kernels::gravity_derivatives(r, eps_j));
                    } else {
                        dphi_dr = Kernels::gravity_derivatives(r, epsilon);
                    }
                    double force_mag_over_r = G * dphi_dr / r;

                    if (use_pm) {
                        force_mag_over_r *=
                            compute_pm_taper(r, inv_two_r_s, taper_coeff);
                    }

                    local_acc_x += force_mag_over_r * m_j * dx;
                    local_acc_y += force_mag_over_r * m_j * dy;
                    local_acc_z += force_mag_over_r * m_j * dz;

                    // Price & Monaghan (2007) adaptive softening correction
                    if (use_adaptive) {
                        double xi_j = d_gas_xi[j];
                        if (xi_i != 0.0 || xi_j != 0.0) {
                            double dW_dr_i =
                                (xi_i != 0.0)
                                    ? Kernels::cubic_spline_gradient(r, eps_i)
                                    : 0.0;
                            double dW_dr_j =
                                (xi_j != 0.0)
                                    ? Kernels::cubic_spline_gradient(r, eps_j)
                                    : 0.0;

                            double corr_over_r =
                                -0.5 * (G / m_i) *
                                (m_i * xi_i * dW_dr_i + m_j * xi_j * dW_dr_j) /
                                r;

                            local_acc_x += corr_over_r * dx;
                            local_acc_y += corr_over_r * dy;
                            local_acc_z += corr_over_r * dz;
                        }
                    }
                });

            // Evaluate against DM Tree
            if (n_dm > 0) {
                LBVH::traverse(
                    d_dm_bvh, d_dm_px, d_dm_py, d_dm_pz, p1_x, p1_y, p1_z,
                    search_sq, domain_size,
                    [=, &local_acc_x, &local_acc_y, &local_acc_z](
                        int j, double dx, double dy, double dz, double r2) {
                        if (r2 > search_sq) return;

                        double r = std::sqrt(r2 + MIN_R2_EPSILON);
                        double m_j = d_dm_m[j];

                        double dphi_dr =
                            use_adaptive
                                ? 0.5 *
                                      (Kernels::gravity_derivatives(r, eps_i) +
                                       Kernels::gravity_derivatives(r, epsilon))
                                : Kernels::gravity_derivatives(r, epsilon);
                        double force_mag_over_r = G * dphi_dr / r;

                        if (use_pm) {
                            force_mag_over_r *=
                                compute_pm_taper(r, inv_two_r_s, taper_coeff);
                        }

                        local_acc_x += force_mag_over_r * m_j * dx;
                        local_acc_y += force_mag_over_r * m_j * dy;
                        local_acc_z += force_mag_over_r * m_j * dz;
                    });
            }

            d_gas_ax[i] += local_acc_x;
            d_gas_ay[i] += local_acc_y;
            d_gas_az[i] += local_acc_z;
        }
    }

#ifdef USE_GPU
    auto end_compute = std::chrono::high_resolution_clock::now();
    auto start_return = std::chrono::high_resolution_clock::now();

#pragma omp target exit data map(from : d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm],  \
                                     d_dm_az[0 : n_dm], d_gas_ax[0 : n_gas],   \
                                     d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas]) \
    map(delete : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],      \
            d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes],                      \
            d_dm_active_idx[0 : n_dm_active], d_gas_px[0 : n_gas],             \
            d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas],      \
            d_gas_bvh[0 : gas_num_nodes], d_gas_active_idx[0 : n_gas_active],  \
            d_gas_h[0 : n_gas], d_gas_xi[0 : n_gas])

    auto end_return = std::chrono::high_resolution_clock::now();

    diag.add_prof_time(GPUProfRegion::Transf,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_transfer - start_transfer)
                           .count());
    diag.add_prof_time(GPUProfRegion::Compute,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_compute - start_compute)
                           .count());
    diag.add_prof_time(GPUProfRegion::Ret,
                       std::chrono::duration_cast<std::chrono::microseconds>(
                           end_return - start_return)
                           .count());
#endif
}
#endif

namespace {
double calculate_max_acceleration(const ParticleSystem& dm) {
    double local_max = 1e-9;
    const size_t n = dm.num_particles;
    const double* ax = dm.acc_x.data();
    const double* ay = dm.acc_y.data();
    const double* az = dm.acc_z.data();

#pragma omp parallel for reduction(max : local_max) schedule(static)
    for (size_t i = 0; i < n; ++i) {
        double accel_sq = ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i];
        local_max = (accel_sq > local_max) ? accel_sq : local_max;
    }

    return std::sqrt(local_max);
}

double calculate_max_acceleration(const GasParticleSystem& gas) {
    double local_max = 1e-9;
    const size_t n = gas.num_particles;
    const double* ax = gas.acc_x.data();
    const double* ay = gas.acc_y.data();
    const double* az = gas.acc_z.data();

#pragma omp parallel for reduction(max : local_max) schedule(static)
    for (size_t i = 0; i < n; ++i) {
        double accel_sq = ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i];
        local_max = (accel_sq > local_max) ? accel_sq : local_max;
    }
    return std::sqrt(local_max);
}
}  // namespace

double GravitySolver::get_timestep(const ParticleSystem& dm,
                                   const GasParticleSystem& gas,
                                   const Config& config) const {
    if (dm.num_particles == 0 && gas.num_particles == 0)
        return std::numeric_limits<double>::infinity();

    const double eps_min = softening_length;
    const double eps_max = config.use_PM
                               ? std::sqrt(config.cutoff_radius_squared)
                               : std::numeric_limits<double>::infinity();

    double dt_min = std::numeric_limits<double>::infinity();

    // DM: fixed softening
    if (dm.num_particles > 0) {
        double a_dm_max = calculate_max_acceleration(dm);
        double dt_dm =
            std::sqrt((2.0 * config.gravity_accuracy_eta * eps_min) / a_dm_max);
        dt_min = std::min(dt_min, dt_dm);
    }

    // Gas: per-particle adaptive softening (or fixed if disabled)
    if (gas.num_particles > 0) {
        if (!config.use_adaptive_softening) {
            double a_gas_max = calculate_max_acceleration(gas);
            double dt_gas = std::sqrt(
                (2.0 * config.gravity_accuracy_eta * eps_min) / a_gas_max);
            dt_min = std::min(dt_min, dt_gas);
        } else {
            double max_a_over_eps = 1e-9;
            const size_t n = gas.num_particles;
            const double* ax = gas.acc_x.data();
            const double* ay = gas.acc_y.data();
            const double* az = gas.acc_z.data();
            const double* h = gas.h.data();

#pragma omp parallel for reduction(max : max_a_over_eps) schedule(static)
            for (size_t i = 0; i < n; ++i) {
                double eps_i = h[i];
                if (eps_i < eps_min) eps_i = eps_min;
                if (eps_i > eps_max) eps_i = eps_max;

                double a_mag =
                    std::sqrt(ax[i] * ax[i] + ay[i] * ay[i] + az[i] * az[i]);
                double ratio = a_mag / eps_i;
                if (ratio > max_a_over_eps) max_a_over_eps = ratio;
            }

            double dt_gas =
                std::sqrt((2.0 * config.gravity_accuracy_eta) / max_a_over_eps);
            dt_min = std::min(dt_min, dt_gas);
        }
    }

    return dt_min;
}

void GravitySolver::update_softening(double a, const Config& config) {
    double cap_a = config.softening_cap_a;
    double current_comoving_softening = config.base_comoving_softening;

    if (a > cap_a && config.expanding_universe) {
        // Shrink comoving softening so physical softening stays constant
        current_comoving_softening =
            config.base_comoving_softening * (cap_a / a);
    }
    softening_length = current_comoving_softening;
}

void GravitySolver::reset_accels(SimState& state, const Config& config) const {
    std::fill(state.dm.acc_x.begin(), state.dm.acc_x.end(), 0.0);
    std::fill(state.dm.acc_y.begin(), state.dm.acc_y.end(), 0.0);
    std::fill(state.dm.acc_z.begin(), state.dm.acc_z.end(), 0.0);
    std::fill(state.gas.acc_x.begin(), state.gas.acc_x.end(), 0.0);
    std::fill(state.gas.acc_y.begin(), state.gas.acc_y.end(), 0.0);
    std::fill(state.gas.acc_z.begin(), state.gas.acc_z.end(), 0.0);
}

GravitySolver::GravitySolver() {}

GravitySolver::~GravitySolver() {
#ifdef GRAVITY_GPU_MEMORY
    int dev = omp_get_default_device();

    if (current_gpu_capacity > 0) {
        omp_target_free(d_u_px, dev);
        omp_target_free(d_t_px, dev);
        omp_target_free(d_u_py, dev);
        omp_target_free(d_t_py, dev);
        omp_target_free(d_u_pz, dev);
        omp_target_free(d_t_pz, dev);
        omp_target_free(d_u_m, dev);
        omp_target_free(d_t_m, dev);
        omp_target_free(d_u_ax, dev);
        omp_target_free(d_u_ay, dev);
        omp_target_free(d_u_az, dev);
        omp_target_free(d_u_orig_idx, dev);
        omp_target_free(d_t_orig_idx, dev);
        omp_target_free(d_u_type, dev);
        omp_target_free(d_t_type, dev);
        omp_target_free(d_t_is_active, dev);
        omp_target_free(d_u_morton, dev);
        omp_target_free(d_u_sorted_idx, dev);

        omp_target_free(d_u_eps, dev);
        omp_target_free(d_t_eps, dev);
        omp_target_free(d_u_gradh_fac, dev);
        omp_target_free(d_t_gradh_fac, dev);
        omp_target_free(d_u_xi, dev);

        uint8_t* old_is_active = u_is_active.data();
#pragma omp target exit data map( \
        delete : old_is_active[0 : current_gpu_capacity])

        current_gpu_capacity = 0;
    }

    if (current_bvh_capacity > 0) {
        omp_target_free(d_u_bvh, dev);
        current_bvh_capacity = 0;
    }
#endif
}
