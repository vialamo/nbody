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

#ifdef GRAVITY_OWN_MEMORY
static std::vector<double> u_px, u_py, u_pz, u_m, u_ax, u_ay, u_az;
static std::vector<double> t_px, t_py, t_pz, t_m;
static std::vector<int> u_original_idx, t_original_idx;
static std::vector<uint8_t> u_type, t_type;
static std::vector<uint64_t> u_morton;
static std::vector<int> u_sorted_indices;
static std::vector<int> u_active_indices;
static std::vector<BVHNode> u_bvh;
#endif

void prepare_arrays(ParticleSystem& dm, GasParticleSystem* gas,
                    const Config& config) {
#ifdef GRAVITY_OWN_MEMORY
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas ? gas->num_particles : 0;
    size_t n_total = n_dm + n_gas;

    if (n_total == 0) return;

    // Resize global arrays
    u_px.resize(n_total);
    u_py.resize(n_total);
    u_pz.resize(n_total);
    u_m.resize(n_total);
    u_original_idx.resize(n_total);
    u_type.resize(n_total);
    u_ax.assign(n_total, 0.0);
    u_ay.assign(n_total, 0.0);
    u_az.assign(n_total, 0.0);

// Data Concatenation
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_dm; ++i) {
        u_px[i] = dm.pos_x[i];
        u_py[i] = dm.pos_y[i];
        u_pz[i] = dm.pos_z[i];
        u_m[i] = dm.mass[i];
        u_original_idx[i] = i;
        u_type[i] = 0;  // Type 0: DM
    }

    if (gas) {
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < n_gas; ++i) {
            size_t offset = n_dm + i;
            u_px[offset] = gas->pos_x[i];
            u_py[offset] = gas->pos_y[i];
            u_pz[offset] = gas->pos_z[i];
            u_m[offset] = gas->mass[i];
            u_original_idx[offset] = i;
            u_type[offset] = 1;  // Type 1: Gas
        }
    }

    // Morton Code Generation and Array Sorting
    LBVH::compute_morton_and_sort_indices(n_total, config.domain_size, u_px,
                                          u_py, u_pz, u_morton,
                                          u_sorted_indices);

    // Reorder data inline using temp arrays
    t_px.resize(n_total);
    t_py.resize(n_total);
    t_pz.resize(n_total);
    t_m.resize(n_total);
    t_original_idx.resize(n_total);
    t_type.resize(n_total);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_total; ++i) {
        int idx = u_sorted_indices[i];
        t_px[i] = u_px[idx];
        t_py[i] = u_py[idx];
        t_pz[i] = u_pz[idx];
        t_m[i] = u_m[idx];
        t_original_idx[i] = u_original_idx[idx];
        t_type[i] = u_type[idx];
    }

    u_px.swap(t_px);
    u_py.swap(t_py);
    u_pz.swap(t_pz);
    u_m.swap(t_m);
    u_original_idx.swap(t_original_idx);
    u_type.swap(t_type);

    // Build Tree & Gather Active Indices
    LBVH::build_topology_and_aggregate(n_total, u_px, u_py, u_pz, u_m, nullptr,
                                       u_morton, u_bvh);

    u_active_indices.clear();
    for (size_t i = 0; i < n_total; ++i) {
        int orig = u_original_idx[i];
        bool is_active =
            (u_type[i] == 0) ? dm.is_active[orig] : gas->is_active[orig];
        if (is_active) u_active_indices.push_back(i);
    }
#endif
}

// Map accelerations back to native structures
void map_results(ParticleSystem& dm, GasParticleSystem* gas) {
#ifdef GRAVITY_OWN_MEMORY
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas ? gas->num_particles : 0;
    size_t n_total = n_dm + n_gas;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_total; ++i) {
        int orig = u_original_idx[i];
        if (u_type[i] == 0) {
            dm.acc_x[orig] += u_ax[i];
            dm.acc_y[orig] += u_ay[i];
            dm.acc_z[orig] += u_az[i];
        } else {
            gas->acc_x[orig] += u_ax[i];
            gas->acc_y[orig] += u_ay[i];
            gas->acc_z[orig] += u_az[i];
        }
    }
#endif
}

#ifdef GRAVITY_OWN_MEMORY
void add_pp_accel(ParticleSystem& dm, GasParticleSystem* gas, double a,
                  const Config& config, Diagnostics& diag) {
    size_t n_dm = dm.num_particles;
    size_t n_gas = gas ? gas->num_particles : 0;
    size_t n_total = n_dm + n_gas;

    if (n_total == 0) return;

    size_t n_active = u_active_indices.size();
    if (n_active == 0) return;

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

#pragma omp target enter data map(                                            \
        to : d_px[0 : n_total], d_py[0 : n_total], d_pz[0 : n_total],         \
            d_m[0 : n_total], d_bvh[0 : num_nodes],                           \
            d_active_idx[0 : n_active], d_ax[0 : n_total], d_ay[0 : n_total], \
            d_az[0 : n_total])

    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();

#pragma omp target teams distribute parallel for
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

#pragma omp target exit data map(from : d_ax[0 : n_total], d_ay[0 : n_total], \
                                     d_az[0 : n_total])                       \
    map(delete : d_px[0 : n_total], d_py[0 : n_total], d_pz[0 : n_total],     \
            d_m[0 : n_total], d_bvh[0 : num_nodes],                           \
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
