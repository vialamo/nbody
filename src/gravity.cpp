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
static std::vector<double> u_px, u_py, u_pz, u_m, u_ax, u_ay, u_az;
static std::vector<double> t_px, t_py, t_pz, t_m;
static std::vector<int> u_original_idx, t_original_idx;
static std::vector<uint8_t> u_type, t_type;
static std::vector<uint8_t> u_is_active, t_is_active;
static std::vector<uint64_t> u_morton;
static std::vector<int> u_sorted_indices;
static std::vector<int> u_active_indices;
static std::vector<BVHNode> u_bvh;

static void prepare_arrays(ParticleSystem& dm, GasParticleSystem* gas,
                           const Config& config) {
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas ? gas->num_particles : 0;
    size_t n_total = n_dm + n_gas;

    if (n_total == 0) return;

    static size_t current_gpu_capacity = 0;

    // GPU Allocation
    if (n_total > current_gpu_capacity) {
        // Free old GPU memory if capacity is expanding
        if (current_gpu_capacity > 0) {
            double* old_px = u_px.data();
            double* old_t_px = t_px.data();
            double* old_py = u_py.data();
            double* old_t_py = t_py.data();
            double* old_pz = u_pz.data();
            double* old_t_pz = t_pz.data();
            double* old_m = u_m.data();
            double* old_t_m = t_m.data();
            double* old_ax = u_ax.data();
            double* old_ay = u_ay.data();
            double* old_az = u_az.data();
            int* old_orig_idx = u_original_idx.data();
            int* old_t_orig_idx = t_original_idx.data();
            uint8_t* old_type = u_type.data();
            uint8_t* old_t_type = t_type.data();
            uint8_t* old_is_active = u_is_active.data();
            uint8_t* old_t_is_active = t_is_active.data();
            uint64_t* old_u_morton = u_morton.data();
            int* old_u_sorted_idx = u_sorted_indices.data();

#pragma omp target exit data map(                                              \
        delete : old_px[0 : current_gpu_capacity],                             \
            old_py[0 : current_gpu_capacity],                                  \
            old_pz[0 : current_gpu_capacity], old_m[0 : current_gpu_capacity], \
            old_ax[0 : current_gpu_capacity],                                  \
            old_ay[0 : current_gpu_capacity],                                  \
            old_az[0 : current_gpu_capacity],                                  \
            old_orig_idx[0 : current_gpu_capacity],                            \
            old_type[0 : current_gpu_capacity],                                \
            old_is_active[0 : current_gpu_capacity],                           \
            old_u_morton[0 : current_gpu_capacity],                            \
            old_u_sorted_idx[0 : current_gpu_capacity],                        \
            old_t_px[0 : current_gpu_capacity],                                \
            old_t_py[0 : current_gpu_capacity],                                \
            old_t_pz[0 : current_gpu_capacity],                                \
            old_t_m[0 : current_gpu_capacity],                                 \
            old_t_orig_idx[0 : current_gpu_capacity],                          \
            old_t_type[0 : current_gpu_capacity],                              \
            old_t_is_active[0 : current_gpu_capacity])
        }

        // Resize host vectors to get new valid pointers
        u_px.resize(n_total);
        u_py.resize(n_total);
        u_pz.resize(n_total);
        u_m.resize(n_total);
        u_original_idx.resize(n_total);
        u_type.resize(n_total);
        u_is_active.resize(n_total);
        u_morton.resize(n_total);
        u_sorted_indices.resize(n_total);
        u_ax.assign(n_total, 0.0);
        u_ay.assign(n_total, 0.0);
        u_az.assign(n_total, 0.0);

        t_px.resize(n_total);
        t_py.resize(n_total);
        t_pz.resize(n_total);
        t_m.resize(n_total);
        t_original_idx.resize(n_total);
        t_type.resize(n_total);
        t_is_active.resize(n_total);

        // Allocate directly on the GPU without copying uninitialized host data
        double* d_u_px = u_px.data();
        double* d_t_px = t_px.data();
        double* d_u_py = u_py.data();
        double* d_t_py = t_py.data();
        double* d_u_pz = u_pz.data();
        double* d_t_pz = t_pz.data();
        double* d_u_m = u_m.data();
        double* d_t_m = t_m.data();
        double* d_u_ax = u_ax.data();
        double* d_u_ay = u_ay.data();
        double* d_u_az = u_az.data();
        int* d_u_orig_idx = u_original_idx.data();
        int* d_t_orig_idx = t_original_idx.data();
        uint8_t* d_u_type = u_type.data();
        uint8_t* d_t_type = t_type.data();
        uint8_t* d_u_is_active = u_is_active.data();
        uint8_t* d_t_is_active = t_is_active.data();
        uint64_t* d_u_morton = u_morton.data();
        int* d_u_sorted_idx = u_sorted_indices.data();

#pragma omp target enter data map(                                             \
        alloc : d_u_px[0 : n_total], d_u_py[0 : n_total], d_u_pz[0 : n_total], \
            d_u_m[0 : n_total], d_u_ax[0 : n_total], d_u_ay[0 : n_total],      \
            d_u_az[0 : n_total], d_u_orig_idx[0 : n_total],                    \
            d_u_type[0 : n_total], d_u_is_active[0 : n_total],                 \
            d_u_morton[0 : n_total], d_u_sorted_idx[0 : n_total],              \
            d_t_px[0 : n_total], d_t_py[0 : n_total], d_t_pz[0 : n_total],     \
            d_t_m[0 : n_total], d_t_orig_idx[0 : n_total],                     \
            d_t_type[0 : n_total], d_t_is_active[0 : n_total])

        current_gpu_capacity = n_total;
    }

    // Zero out the accelerations on the GPU for the current step
    double* d_u_ax = u_ax.data();
    double* d_u_ay = u_ay.data();
    double* d_u_az = u_az.data();

#pragma omp target teams distribute parallel for map(       \
        present : d_u_ax[0 : n_total], d_u_ay[0 : n_total], \
            d_u_az[0 : n_total])
    for (size_t i = 0; i < n_total; ++i) {
        d_u_ax[i] = 0.0;
        d_u_ay[i] = 0.0;
        d_u_az[i] = 0.0;
    }

    // Pointers for Device Mapping
    double* d_u_px = u_px.data();
    double* d_u_py = u_py.data();
    double* d_u_pz = u_pz.data();
    double* d_u_m = u_m.data();
    int* d_u_orig_idx = u_original_idx.data();
    uint8_t* d_u_type = u_type.data();
    uint8_t* d_u_is_active = u_is_active.data();

    // On-Device Concatenation
    // Map the CPU particle arrays to the device temporarily to build the
    // unified arrays
    double* d_dm_px = dm.pos_x.data();
    double* d_dm_py = dm.pos_y.data();
    double* d_dm_pz = dm.pos_z.data();
    double* d_dm_m = dm.mass.data();
    uint8_t* d_dm_is_active = dm.is_active.data();

#pragma omp target teams distribute parallel for map(                         \
        present : d_u_px[0 : n_total], d_u_py[0 : n_total],                   \
            d_u_pz[0 : n_total], d_u_m[0 : n_total],                          \
            d_u_orig_idx[0 : n_total], d_u_type[0 : n_total],                 \
            d_u_is_active[0 : n_total], d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], \
            d_dm_pz[0 : n_dm], d_dm_m[0 : n_dm], d_dm_is_active[0 : n_dm])
    for (size_t i = 0; i < n_dm; ++i) {
        d_u_px[i] = d_dm_px[i];
        d_u_py[i] = d_dm_py[i];
        d_u_pz[i] = d_dm_pz[i];
        d_u_m[i] = d_dm_m[i];
        d_u_orig_idx[i] = i;
        d_u_type[i] = 0;  // Type 0: DM
        d_u_is_active[i] = d_dm_is_active[i];
    }

#pragma omp target exit data map(                                         \
        delete : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm], \
            d_dm_m[0 : n_dm], d_dm_is_active[0 : n_dm])

    if (gas && n_gas > 0) {
        double* d_gas_px = gas->pos_x.data();
        double* d_gas_py = gas->pos_y.data();
        double* d_gas_pz = gas->pos_z.data();
        double* d_gas_m = gas->mass.data();
        uint8_t* d_gas_is_active = gas->is_active.data();

#pragma omp target teams distribute parallel for map(                     \
        present : d_u_px[0 : n_total], d_u_py[0 : n_total],               \
            d_u_pz[0 : n_total], d_u_m[0 : n_total],                      \
            d_u_orig_idx[0 : n_total], d_u_type[0 : n_total],             \
            d_u_is_active[0 : n_total], d_gas_px[0 : n_gas],              \
            d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas], \
            d_gas_is_active[0 : n_gas])
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

#pragma omp target exit data map(delete : d_gas_px[0 : n_gas],                 \
                                     d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], \
                                     d_gas_m[0 : n_gas],                       \
                                     d_gas_is_active[0 : n_gas])
    }

    // Run Morton generation
    LBVH::compute_morton_and_sort_indices_gpu(
        n_total, config.domain_size, d_u_px, d_u_py, d_u_pz, u_morton.data(),
        u_sorted_indices.data());

    // Allocate Temp Arrays on Host and Device
    double* d_t_px = t_px.data();
    double* d_t_py = t_py.data();
    double* d_t_pz = t_pz.data();
    double* d_t_m = t_m.data();
    int* d_t_orig_idx = t_original_idx.data();
    uint8_t* d_t_type = t_type.data();
    int* d_sorted_idx = u_sorted_indices.data();
    uint8_t* d_t_is_active = t_is_active.data();

// Permute into temporary device arrays
#pragma omp target teams distribute parallel for map(                      \
        present : d_u_px[0 : n_total], d_u_py[0 : n_total],                \
            d_u_pz[0 : n_total], d_u_m[0 : n_total],                       \
            d_u_orig_idx[0 : n_total], d_u_type[0 : n_total],              \
            d_u_is_active[0 : n_total], d_sorted_idx[0 : n_total],         \
            d_t_px[0 : n_total], d_t_py[0 : n_total], d_t_pz[0 : n_total], \
            d_t_m[0 : n_total], d_t_orig_idx[0 : n_total],                 \
            d_t_type[0 : n_total], d_t_is_active[0 : n_total])
    for (size_t i = 0; i < n_total; ++i) {
        int idx = d_sorted_idx[i];
        d_t_px[i] = d_u_px[idx];
        d_t_py[i] = d_u_py[idx];
        d_t_pz[i] = d_u_pz[idx];
        d_t_m[i] = d_u_m[idx];
        d_t_orig_idx[i] = d_u_orig_idx[idx];
        d_t_type[i] = d_u_type[idx];
        d_t_is_active[i] = d_u_is_active[idx];
    }

// Copy back to the persistent device arrays
#pragma omp target teams distribute parallel for map(                     \
        present : d_u_px[0 : n_total], d_u_py[0 : n_total],               \
            d_u_pz[0 : n_total], d_u_m[0 : n_total],                      \
            d_u_orig_idx[0 : n_total], d_u_type[0 : n_total],             \
            d_u_is_active[0 : n_total], d_t_px[0 : n_total],              \
            d_t_py[0 : n_total], d_t_pz[0 : n_total], d_t_m[0 : n_total], \
            d_t_orig_idx[0 : n_total], d_t_type[0 : n_total],             \
            d_t_is_active[0 : n_total])
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
        if (current_bvh_capacity > 0) {
            BVHNode* old_bvh = u_bvh.data();
#pragma omp target exit data map(delete : old_bvh[0 : current_bvh_capacity])
        }
        u_bvh.resize(required_nodes);
        BVHNode* d_u_bvh = u_bvh.data();
#pragma omp target enter data map(alloc : d_u_bvh[0 : required_nodes])
        current_bvh_capacity = required_nodes;
    }

    // Call the device-side builder using the GPU-resident arrays
    LBVH::build_topology_and_aggregate_gpu(n_total, u_px.data(), u_py.data(),
                                           u_pz.data(), u_m.data(),
                                           u_morton.data(), u_bvh.data());

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

    double* d_u_ax = u_ax.data();
    double* d_u_ay = u_ay.data();
    double* d_u_az = u_az.data();
    int* d_u_orig_idx = u_original_idx.data();
    uint8_t* d_u_type = u_type.data();
    int* d_active_idx = u_active_indices.data();

#pragma omp target data map(to : d_active_idx[0 : n_active])
    {
        if (n_dm > 0) {
            double* d_dm_ax = dm.acc_x.data();
            double* d_dm_ay = dm.acc_y.data();
            double* d_dm_az = dm.acc_z.data();

#pragma omp target teams distribute parallel for map(                     \
        tofrom : d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], d_dm_az[0 : n_dm]) \
    map(present : d_u_ax[0 : n_total], d_u_ay[0 : n_total],               \
            d_u_az[0 : n_total], d_u_orig_idx[0 : n_total],               \
            d_u_type[0 : n_total], d_active_idx[0 : n_active])
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

#pragma omp target teams distribute parallel for map(       \
        tofrom : d_gas_ax[0 : n_gas], d_gas_ay[0 : n_gas],  \
            d_gas_az[0 : n_gas])                            \
    map(present : d_u_ax[0 : n_total], d_u_ay[0 : n_total], \
            d_u_az[0 : n_total], d_u_orig_idx[0 : n_total], \
            d_u_type[0 : n_total], d_active_idx[0 : n_active])
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
    size_t n_active = dm.num_active + (gas ? gas->num_active : 0);

    if (n_active == 0) return;

    prepare_arrays(dm, gas, config);

    // Shared Force Constants
    const Vec3& domain_size = config.domain_size;
    const double G = config.G / (a * a * a);
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const double epsilon = config.softening;
    const double search_sq = use_pm ? config.cutoff_radius_squared
                                    : std::numeric_limits<double>::infinity();
    const size_t num_nodes = 2 * n_total - 1;

    // Raw pointers for GPU Offload
    double* d_px = u_px.data();
    double* d_py = u_py.data();
    double* d_pz = u_pz.data();
    double* d_m = u_m.data();
    double* d_ax = u_ax.data();
    double* d_ay = u_ay.data();
    double* d_az = u_az.data();
    BVHNode* d_bvh = u_bvh.data();
    int* d_active_idx = u_active_indices.data();

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();
#pragma omp target enter data map(to : d_active_idx[0 : n_active])
    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();

#pragma omp target teams distribute parallel for map(                      \
        present : d_px[0 : n_total], d_py[0 : n_total], d_pz[0 : n_total], \
            d_m[0 : n_total], d_bvh[0 : num_nodes], d_ax[0 : n_total],     \
            d_ay[0 : n_total], d_az[0 : n_total], d_active_idx[0 : n_active])
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    // Unified Kernel Execution
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

#ifdef USE_GPU
    auto end_compute = std::chrono::high_resolution_clock::now();
    auto start_return = std::chrono::high_resolution_clock::now();

    // Copy the computed accelerations to the CPU
#pragma omp target update from(d_ax[0 : n_total], d_ay[0 : n_total], \
                                   d_az[0 : n_total])
#pragma omp target exit data map(delete : d_active_idx[0 : n_active])

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
