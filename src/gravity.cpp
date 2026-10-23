#include "gravity.h"

#include <chrono>
#include <cmath>

#include "diagnostics.h"
#include "kernels.h"
#include "math_utils.h"
#include "state.h"

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

#ifdef GRAVITY_GPU_MEMORY
static double *d_u_px = nullptr, *d_u_py = nullptr, *d_u_pz = nullptr,
              *d_u_m = nullptr;
static double *d_u_ax = nullptr, *d_u_ay = nullptr, *d_u_az = nullptr;
static double *d_t_px = nullptr, *d_t_py = nullptr, *d_t_pz = nullptr,
              *d_t_m = nullptr;
static int *d_u_orig_idx = nullptr, *d_t_orig_idx = nullptr,
           *d_u_sorted_idx = nullptr;
static uint8_t *d_u_type = nullptr, *d_t_type = nullptr,
               *d_t_is_active = nullptr;
static uint64_t* d_u_morton = nullptr;
static BVHNode* d_u_bvh = nullptr;

// Arrays requiring CPU interaction
static std::vector<uint8_t> u_is_active;
static std::vector<int> u_active_indices;

static void prepare_arrays(ParticleSystem& dm, GasParticleSystem* gas,
                           const Config& config) {
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas ? gas->num_particles : 0;
    size_t n_total = n_dm + n_gas;

    if (n_total == 0) return;

    static size_t current_gpu_capacity = 0;

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

            // Only use map(delete) for the host-backed active array
            uint8_t* old_is_active = u_is_active.data();
#pragma omp target exit data map( \
        delete : old_is_active[0 : current_gpu_capacity])
        }

        u_is_active.resize(n_total);
        uint8_t* d_is_act = u_is_active.data();
#pragma omp target enter data map(alloc : d_is_act[0 : n_total])

        // Allocate directly on the GPU without copying uninitialized host data
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

        current_gpu_capacity = n_total;
    }

    // Zero out the accelerations on the GPU for the current step
#pragma omp target teams distribute parallel for is_device_ptr(d_u_ax, d_u_ay, \
                                                                   d_u_az)
    for (size_t i = 0; i < n_total; ++i) {
        d_u_ax[i] = 0.0;
        d_u_ay[i] = 0.0;
        d_u_az[i] = 0.0;
    }

    // Pointers for Device Mapping
    uint8_t* d_u_is_active = u_is_active.data();

    // On-Device Concatenation
    // Map the CPU particle arrays to the device temporarily to build the
    // unified arrays
    double* d_dm_px = dm.pos_x.data();
    double* d_dm_py = dm.pos_y.data();
    double* d_dm_pz = dm.pos_z.data();
    double* d_dm_m = dm.mass.data();
    uint8_t* d_dm_is_active = dm.is_active.data();

#pragma omp target teams distribute parallel for is_device_ptr(       \
        d_u_px, d_u_py, d_u_pz, d_u_m, d_u_orig_idx, d_u_type)        \
    map(to : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm], \
            d_dm_m[0 : n_dm], d_dm_is_active[0 : n_dm])               \
    map(alloc : d_u_is_active[0 : n_total])
    for (size_t i = 0; i < n_dm; ++i) {
        d_u_px[i] = d_dm_px[i];
        d_u_py[i] = d_dm_py[i];
        d_u_pz[i] = d_dm_pz[i];
        d_u_m[i] = d_dm_m[i];
        d_u_orig_idx[i] = i;
        d_u_type[i] = 0;  // Type 0: DM
        d_u_is_active[i] = d_dm_is_active[i];
    }

    if (gas && n_gas > 0) {
        double* d_gas_px = gas->pos_x.data();
        double* d_gas_py = gas->pos_y.data();
        double* d_gas_pz = gas->pos_z.data();
        double* d_gas_m = gas->mass.data();
        uint8_t* d_gas_is_active = gas->is_active.data();

#pragma omp target teams distribute parallel for is_device_ptr(             \
        d_u_px, d_u_py, d_u_pz, d_u_m, d_u_orig_idx, d_u_type)              \
    map(to : d_gas_px[0 : n_gas], d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], \
            d_gas_m[0 : n_gas], d_gas_is_active[0 : n_gas])                 \
    map(alloc : d_u_is_active[0 : n_total])
        for (size_t i = 0; i < n_gas; ++i) {
            size_t offset = n_dm + i;
            d_u_px[offset] = d_gas_px[i];
            d_u_py[offset] = d_gas_py[i];
            d_u_pz[offset] = d_gas_pz[i];
            d_u_m[offset] = d_gas_m[i];
            d_u_orig_idx[offset] = i;
            d_u_type[offset] = 1;  // Type 1: Gas
            d_u_is_active[offset] = d_gas_is_active[i];
        }
    }

    // Run Morton generation
    LBVH::compute_morton_and_sort_indices_gpu(n_total, config.domain_size,
                                              d_u_px, d_u_py, d_u_pz,
                                              d_u_morton, d_u_sorted_idx);

// Permute into temporary device arrays
#pragma omp target teams distribute parallel for is_device_ptr(                \
        d_u_px, d_u_py, d_u_pz, d_u_m, d_u_orig_idx, d_u_type, d_u_sorted_idx, \
            d_t_px, d_t_py, d_t_pz, d_t_m, d_t_orig_idx, d_t_type,             \
            d_t_is_active) map(alloc : d_u_is_active[0 : n_total])
    for (size_t i = 0; i < n_total; ++i) {
        int idx = d_u_sorted_idx[i];
        d_t_px[i] = d_u_px[idx];
        d_t_py[i] = d_u_py[idx];
        d_t_pz[i] = d_u_pz[idx];
        d_t_m[i] = d_u_m[idx];
        d_t_orig_idx[i] = d_u_orig_idx[idx];
        d_t_type[i] = d_u_type[idx];
        d_t_is_active[i] = d_u_is_active[idx];
    }

// Copy back to the persistent device arrays
#pragma omp target teams distribute parallel for is_device_ptr(                \
        d_u_px, d_u_py, d_u_pz, d_u_m, d_u_orig_idx, d_u_type, d_t_px, d_t_py, \
            d_t_pz, d_t_m, d_t_orig_idx, d_t_type, d_t_is_active)              \
    map(alloc : d_u_is_active[0 : n_total])
    for (size_t i = 0; i < n_total; ++i) {
        d_u_px[i] = d_t_px[i];
        d_u_py[i] = d_t_py[i];
        d_u_pz[i] = d_t_pz[i];
        d_u_m[i] = d_t_m[i];
        d_u_orig_idx[i] = d_t_orig_idx[i];
        d_u_type[i] = d_t_type[i];
        d_u_is_active[i] = d_t_is_active[i];
    }

    // Persistent GPU allocation for the tree
    static size_t current_bvh_capacity = 0;
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
    LBVH::build_topology_and_aggregate_gpu(n_total, d_u_px, d_u_py, d_u_pz,
                                           d_u_m, d_u_morton, d_u_bvh);

    uint8_t* d_u_is_active_end = u_is_active.data();
#pragma omp target update from(d_u_is_active_end[0 : n_total])

    // Rebuild the dense array sequentially to preserve Morton spatial
    // locality
    u_active_indices.clear();
    for (size_t i = 0; i < n_total; ++i) {
        if (u_is_active[i]) u_active_indices.push_back(i);
    }
}

// Map accelerations back to native structures
static void map_results(ParticleSystem& dm, GasParticleSystem* gas) {
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas ? gas->num_particles : 0;
    size_t n_total = n_dm + n_gas;
    size_t n_active = u_active_indices.size();

    if (n_active == 0) return;

    int* d_active_idx = u_active_indices.data();

#pragma omp target data map(to : d_active_idx[0 : n_active])
    {
        if (n_dm > 0) {
            double* d_dm_ax = dm.acc_x.data();
            double* d_dm_ay = dm.acc_y.data();
            double* d_dm_az = dm.acc_z.data();

#pragma omp target teams distribute parallel for map(                     \
        tofrom : d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm]) \
    map(alloc : d_active_idx[0 : n_active])                               \
    is_device_ptr(d_u_ax, d_u_ay, d_u_az, d_u_orig_idx, d_u_type)
            for (size_t k = 0; k < n_active; ++k) {
                size_t i = d_active_idx[k];
                if (d_u_type[i] == 0) {
                    int orig = d_u_orig_idx[i];
                    d_dm_ax[orig] += d_u_ax[i];
                    d_dm_ay[orig] += d_u_ay[i];
                    d_dm_az[orig] += d_u_az[i];
                }
            }
        }

        if (n_gas > 0) {
            double* d_gas_ax = gas->acc_x.data();
            double* d_gas_ay = gas->acc_y.data();
            double* d_gas_az = gas->acc_z.data();

#pragma omp target teams distribute parallel for map(                    \
        tofrom : d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas],               \
            d_gas_az[0 : n_gas]) map(alloc : d_active_idx[0 : n_active]) \
    is_device_ptr(d_u_ax, d_u_ay, d_u_az, d_u_orig_idx, d_u_type)
            for (size_t k = 0; k < n_active; ++k) {
                size_t i = d_active_idx[k];
                if (d_u_type[i] == 1) {
                    int orig = d_u_orig_idx[i];
                    d_gas_ax[orig] += d_u_ax[i];
                    d_gas_ay[orig] += d_u_ay[i];
                    d_gas_az[orig] += d_u_az[i];
                }
            }
        }
    }
}

void add_pp_accel(ParticleSystem& dm, GasParticleSystem* gas, double a,
                  const Config& config, Diagnostics& diag) {
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas ? gas->num_particles : 0;
    size_t n_total = n_dm + n_gas;

    if (n_total == 0) return;

    prepare_arrays(dm, gas, config);

    size_t n_active = u_active_indices.size();
    if (n_active == 0) return;

    // Shared Force Constants
    const Vec3 domain_size = config.domain_size;
    const double G = config.G / (a * a * a);
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const double epsilon = config.softening;
    const double search_sq = use_pm ? config.cutoff_radius_squared
                                    : std::numeric_limits<double>::infinity();
    const size_t num_nodes = 2 * n_total - 1;

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

    auto start_transfer = std::chrono::high_resolution_clock::now();
    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();

#pragma omp target teams distribute parallel for map( \
        to : d_active_idx[0 : n_active])              \
    is_device_ptr(d_px, d_py, d_pz, d_m, d_bvh, d_ax, d_ay, d_az)
    for (size_t k = 0; k < n_active; ++k) {
        size_t i = d_active_idx[k];

        double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
        double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
        int stack[128];
        int stack_ptr = 0;
        stack[stack_ptr++] = 0;

        while (stack_ptr > 0) {
            int node_idx = stack[--stack_ptr];
            const BVHNode& node = d_bvh[node_idx];

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

                if (dist_sq > search_sq) continue;

                double r = std::sqrt(dist_sq + 1e-24);
                double m_j = d_m[j];

                double dphi_dr = Kernels::gravity_derivatives(r, epsilon);
                double force_mag_over_r = G * dphi_dr / r;

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

    auto end_compute = std::chrono::high_resolution_clock::now();
    auto start_return = std::chrono::high_resolution_clock::now();
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

    map_results(dm, gas);
}

#else
void add_pp_accel(ParticleSystem& dm, GasParticleSystem* gas, double a,
                  const Config& config, Diagnostics& diag) {
    if (dm.num_particles == 0 &&
        (gas == nullptr || (gas && gas->num_particles == 0)))
        return;

    const Vec3& domain_size = config.domain_size;
    const double G = config.G / (a * a * a);
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const double epsilon = config.softening;

    const double search_sq = use_pm ? config.cutoff_radius_squared
                                    : std::numeric_limits<double>::infinity();

    const size_t n_dm = dm.num_particles;
    const size_t dm_num_nodes = (n_dm > 0) ? (2 * n_dm - 1) : 0;
    const size_t n_dm_active = dm.num_active;

    const size_t n_gas = gas ? gas->num_particles : 0;
    const size_t gas_num_nodes = (n_gas > 0) ? (2 * n_gas - 1) : 0;
    const size_t n_gas_active = gas ? gas->num_active : 0;

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
    double* d_gas_px = gas ? gas->pos_x.data() : nullptr;
    double* d_gas_py = gas ? gas->pos_y.data() : nullptr;
    double* d_gas_pz = gas ? gas->pos_z.data() : nullptr;
    double* d_gas_m = gas ? gas->mass.data() : nullptr;
    double* d_gas_ax = gas ? gas->acc_x.data() : nullptr;
    double* d_gas_ay = gas ? gas->acc_y.data() : nullptr;
    double* d_gas_az = gas ? gas->acc_z.data() : nullptr;
    BVHNode* d_gas_bvh = gas ? gas->bvh_nodes.data() : nullptr;
    int* d_gas_active_idx = gas ? gas->active_indices.data() : nullptr;

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                            \
        to : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],         \
            d_dm_m[0 : n_dm], d_dm_bvh[0 : dm_num_nodes],                     \
            d_dm_active_idx[0 : n_dm_active], d_dm_ax[0 : n_dm],              \
            d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm], d_gas_px[0 : n_gas],        \
            d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas],     \
            d_gas_bvh[0 : gas_num_nodes], d_gas_active_idx[0 : n_gas_active], \
            d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas], d_gas_az[0 : n_gas])

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
            int stack[128];

            // Evaluate against DM Tree
            int stack_ptr = 0;
            stack[stack_ptr++] = 0;
            while (stack_ptr > 0) {
                int node_idx = stack[--stack_ptr];
                const BVHNode& node = d_dm_bvh[node_idx];

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

                    dx = -dx;
                    dy = -dy;
                    dz = -dz;
                    double dist_sq = dx * dx + dy * dy + dz * dz;

                    if (dist_sq > search_sq) continue;

                    double r = std::sqrt(dist_sq + 1e-24);
                    double m_j = d_dm_m[j];

                    double dphi_dr = Kernels::gravity_derivatives(r, epsilon);
                    double force_mag_over_r = G * dphi_dr / r;

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

            // Evaluate against Gas Tree
            if (n_gas > 0) {
                stack_ptr = 0;
                stack[stack_ptr++] = 0;
                while (stack_ptr > 0) {
                    int node_idx = stack[--stack_ptr];
                    const BVHNode& node = d_gas_bvh[node_idx];

                    double aabb_dist_sq =
                        min_periodic_dist_sq(p1_x, node.bbox.min_x,
                                             node.bbox.max_x, domain_size.x) +
                        min_periodic_dist_sq(p1_y, node.bbox.min_y,
                                             node.bbox.max_y, domain_size.y) +
                        min_periodic_dist_sq(p1_z, node.bbox.min_z,
                                             node.bbox.max_z, domain_size.z);

                    if (aabb_dist_sq > search_sq) continue;

                    if (node.particle_idx != -1) {
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

                        dx = -dx;
                        dy = -dy;
                        dz = -dz;
                        double dist_sq = dx * dx + dy * dy + dz * dz;

                        if (dist_sq > search_sq) continue;

                        double r = std::sqrt(dist_sq + 1e-24);
                        double m_j = d_gas_m[j];

                        double dphi_dr =
                            Kernels::gravity_derivatives(r, epsilon);
                        double force_mag_over_r = G * dphi_dr / r;

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
            int stack[128];

            // Evaluate against Gas Tree
            int stack_ptr = 0;
            stack[stack_ptr++] = 0;
            while (stack_ptr > 0) {
                int node_idx = stack[--stack_ptr];
                const BVHNode& node = d_gas_bvh[node_idx];

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

                    dx = -dx;
                    dy = -dy;
                    dz = -dz;
                    double dist_sq = dx * dx + dy * dy + dz * dz;

                    if (dist_sq > search_sq) continue;

                    double r = std::sqrt(dist_sq + 1e-24);
                    double m_j = d_gas_m[j];

                    double dphi_dr = Kernels::gravity_derivatives(r, epsilon);
                    double force_mag_over_r = G * dphi_dr / r;

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

            // Evaluate against DM Tree
            if (n_dm > 0) {
                stack_ptr = 0;
                stack[stack_ptr++] = 0;
                while (stack_ptr > 0) {
                    int node_idx = stack[--stack_ptr];
                    const BVHNode& node = d_dm_bvh[node_idx];

                    double aabb_dist_sq =
                        min_periodic_dist_sq(p1_x, node.bbox.min_x,
                                             node.bbox.max_x, domain_size.x) +
                        min_periodic_dist_sq(p1_y, node.bbox.min_y,
                                             node.bbox.max_y, domain_size.y) +
                        min_periodic_dist_sq(p1_z, node.bbox.min_z,
                                             node.bbox.max_z, domain_size.z);

                    if (aabb_dist_sq > search_sq) continue;

                    if (node.particle_idx != -1) {
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

                        dx = -dx;
                        dy = -dy;
                        dz = -dz;
                        double dist_sq = dx * dx + dy * dy + dz * dz;

                        if (dist_sq > search_sq) continue;

                        double r = std::sqrt(dist_sq + 1e-24);
                        double m_j = d_dm_m[j];

                        double dphi_dr =
                            Kernels::gravity_derivatives(r, epsilon);
                        double force_mag_over_r = G * dphi_dr / r;

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
            d_gas_bvh[0 : gas_num_nodes], d_gas_active_idx[0 : n_gas_active])

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
#endif

static double calculate_max_acceleration(const ParticleSystem& dm) {
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

static double calculate_max_acceleration(const GasParticleSystem& gas) {
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

double get_timestep(const GasParticleSystem& gas, const Config& config) {
    if (gas.num_particles == 0) return std::numeric_limits<double>::infinity();

    double epsilon = config.softening;
    double a_max = calculate_max_acceleration(gas);
    double dt_grav = std::sqrt(epsilon / a_max);
    return dt_grav * config.gravity_accuracy_eta;
}

double get_timestep(const ParticleSystem& dm, const Config& config) {
    if (dm.num_particles == 0) return std::numeric_limits<double>::infinity();

    double epsilon = config.softening;
    double a_max = calculate_max_acceleration(dm);
    double dt_grav = std::sqrt(epsilon / a_max);
    return dt_grav * config.gravity_accuracy_eta;
}

}  // namespace GravitySolver
