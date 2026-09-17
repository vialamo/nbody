#include "particles.h"

#include <omp.h>

#include "cic.h"
#include "diagnostics.h"
#include "math_utils.h"

ParticleSystem::ParticleSystem(const Config& config)
    : dm_rho(config.mesh_size),
      cic_data(config.num_dm_particles),
      max_accel_sq(0.0) {
    pos_x.reserve(config.num_dm_particles);
    pos_y.reserve(config.num_dm_particles);
    pos_z.reserve(config.num_dm_particles);

    vel_x.reserve(config.num_dm_particles);
    vel_y.reserve(config.num_dm_particles);
    vel_z.reserve(config.num_dm_particles);

    acc_x.reserve(config.num_dm_particles);
    acc_y.reserve(config.num_dm_particles);
    acc_z.reserve(config.num_dm_particles);

    mass.reserve(config.num_dm_particles);

    size_t num_nodes = 2 * config.num_dm_particles - 1;
    if (config.num_dm_particles > 0) {
        morton_codes.reserve(config.num_dm_particles);
        sorted_indices.reserve(config.num_dm_particles);
        bvh_nodes.reserve(num_nodes);
    }
}

void ParticleSystem::add_particle(double px, double py, double pz, double vx,
                                  double vy, double vz, double m) {
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
    num_particles++;
}

void ParticleSystem::bin_and_assign_mass(const Config& config) {
    CIC::bin_and_assign_mass(config, num_particles, pos_x, pos_y, pos_z, mass,
                             cic_data, dm_rho);
}

void ParticleSystem::interpolate_cic_forces(const Grid3D& ax_grid,
                                            const Grid3D& ay_grid,
                                            const Grid3D& az_grid,
                                            const Config& config) {
    CIC::interpolate_forces(config, num_particles, cic_data, ax_grid, ay_grid,
                            az_grid, acc_x, acc_y, acc_z);
}

void ParticleSystem::sort_arrays(const std::vector<int>& sorted_indices) {
    if (num_particles == 0) return;

    std::vector<double> new_px(num_particles), new_py(num_particles),
        new_pz(num_particles);
    std::vector<double> new_vx(num_particles), new_vy(num_particles),
        new_vz(num_particles);
    std::vector<double> new_ax(num_particles), new_ay(num_particles),
        new_az(num_particles);
    std::vector<double> new_m(num_particles);

    // CIC_Data is only needed for PM gravity.
    std::vector<CIC_Data> new_cic;
    bool has_cic = !cic_data.empty();
    if (has_cic) {
        new_cic.resize(num_particles);
    }

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

        if (has_cic) new_cic[i] = cic_data[src];
    }

    // Move sorted data back
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

    if (has_cic) cic_data = std::move(new_cic);
}

void ParticleSystem::build_lbvh(const Config& config) {
    if (num_particles == 0) return;

    // Generate Morton codes
    LBVH::compute_morton_and_sort_indices(num_particles, config.domain_size,
                                          pos_x, pos_y, pos_z, morton_codes,
                                          sorted_indices);

    // Shuffle the DM arrays
    sort_arrays(sorted_indices);

    // Build topology (Passing nullptr for 'h' since DM doesn't use it)
    LBVH::build_topology_and_aggregate(num_particles, pos_x, pos_y, pos_z, mass,
                                       nullptr, morton_codes, bvh_nodes);
}

void ParticleSystem::compute_and_add_pp_forces(double a, const Config& config,
                                               Diagnostics& diag) {
    if (num_particles == 0) return;

    const size_t n_parts = num_particles;
    const double domain_size = config.domain_size;
    // Pre-scale G so all output forces are comoving accelerations
    const double G = config.G / (a * a * a);
    const double soft_sq = config.softening_squared;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t num_nodes = 2 * n_parts - 1;

    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();

    // Extract pointers
    double* d_px = pos_x.data();
    double* d_py = pos_y.data();
    double* d_pz = pos_z.data();
    double* d_m = mass.data();
    double* d_ax = acc_x.data();
    double* d_ay = acc_y.data();
    double* d_az = acc_z.data();
    const BVHNode* d_bvh_nodes = bvh_nodes.data();

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                           \
        to : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts],        \
            d_m[0 : n_parts], d_bvh_nodes[0 : num_nodes], d_ax[0 : n_parts], \
            d_ay[0 : n_parts], d_az[0 : n_parts])

    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t i = 0; i < n_parts; ++i) {
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

                double pp_dist_sq = dist_sq + soft_sq;
                double pp_dist = std::sqrt(pp_dist_sq);
                double a_pp = G * node.mass / pp_dist_sq;

                if (use_pm) {
                    double r = std::sqrt(dist_sq);
                    double r_scaled = r / (2.0 * r_s);
                    a_pp *= (std::erfc(r_scaled) +
                             (r / (std::sqrt(M_PI) * r_s)) *
                                 std::exp(-r_scaled * r_scaled));
                }

                local_acc_x += a_pp * dx / pp_dist;
                local_acc_y += a_pp * dy / pp_dist;
                local_acc_z += a_pp * dz / pp_dist;
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
            d_m[0 : n_parts], d_bvh_nodes[0 : num_nodes])

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

double ParticleSystem::get_gravity_timestep(const Config& config) const {
    if (num_particles == 0) return std::numeric_limits<double>::infinity();

    double epsilon = std::sqrt(config.softening_squared);
    double a_max = std::sqrt(max_accel_sq);
    double dt_grav = std::sqrt(epsilon / a_max);

    return dt_grav * config.gravity_accuracy_eta;
}