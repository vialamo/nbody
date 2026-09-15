#include "particles.h"

#include <omp.h>

#include "cic.h"
#include "diagnostics.h"
#include "gas.h"
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

void ParticleSystem::compute_and_add_pp_forces(const Config& config,
                                               Diagnostics& diag) {
    if (num_particles == 0) return;

    const size_t n_parts = num_particles;
    const double domain_size = config.domain_size;
    const double G = config.G;
    const double soft_sq = config.softening_squared;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t num_nodes = 2 * n_parts - 1;

    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();

    // Extract pointers for OpenMP GPU map clauses
    double* d_px = pos_x.data();
    double* d_py = pos_y.data();
    double* d_pz = pos_z.data();
    double* d_m = mass.data();
    double* d_ax = acc_x.data();
    double* d_ay = acc_y.data();
    double* d_az = acc_z.data();
    const BVHNode* d_bvh_nodes = bvh_nodes.data();

#ifdef USE_GPU
    if (config.enable_GPU) {
        auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                           \
        to : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts],        \
            d_m[0 : n_parts], d_bvh_nodes[0 : num_nodes], d_ax[0 : n_parts], \
            d_ay[0 : n_parts], d_az[0 : n_parts])

        auto end_transfer = std::chrono::high_resolution_clock::now();
        auto start_compute = std::chrono::high_resolution_clock::now();

#pragma omp target teams distribute parallel for
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
                        double r_scaled = pp_dist / (2.0 * r_s);
                        a_pp *= (std::erfc(r_scaled) +
                                 (pp_dist / (std::sqrt(M_PI) * r_s)) *
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

        auto end_compute = std::chrono::high_resolution_clock::now();
        auto start_return = std::chrono::high_resolution_clock::now();

#pragma omp target exit data map(from : d_ax[0 : n_parts], d_ay[0 : n_parts], \
                                     d_az[0 : n_parts])                       \
    map(delete : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts],     \
            d_m[0 : n_parts], d_bvh_nodes[0 : num_nodes])

        auto end_return = std::chrono::high_resolution_clock::now();

        diag.add_prof_time(
            ProfRegion::Transf,
            std::chrono::duration_cast<std::chrono::microseconds>(
                end_transfer - start_transfer)
                .count());
        diag.add_prof_time(
            ProfRegion::Compute,
            std::chrono::duration_cast<std::chrono::microseconds>(end_compute -
                                                                  start_compute)
                .count());
        diag.add_prof_time(
            ProfRegion::Ret,
            std::chrono::duration_cast<std::chrono::microseconds>(end_return -
                                                                  start_return)
                .count());

    } else
#endif
    {
        // ========================================================================
        // CPU IMPLEMENTATION
        // ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];
            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;

            int stack[128];
            int stack_ptr = 0;
            stack[stack_ptr++] = 0;

            while (stack_ptr > 0) {
                int node_idx = stack[--stack_ptr];
                const BVHNode& node = bvh_nodes[node_idx];

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

                    double dx = p1_x - pos_x[j];
                    if (dx > 0.5 * domain_size)
                        dx -= domain_size;
                    else if (dx < -0.5 * domain_size)
                        dx += domain_size;

                    double dy = p1_y - pos_y[j];
                    if (dy > 0.5 * domain_size)
                        dy -= domain_size;
                    else if (dy < -0.5 * domain_size)
                        dy += domain_size;

                    double dz = p1_z - pos_z[j];
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
                        double r_scaled = pp_dist / (2.0 * r_s);
                        a_pp *= (std::erfc(r_scaled) +
                                 (pp_dist / (std::sqrt(M_PI) * r_s)) *
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

            acc_x[i] += local_acc_x;
            acc_y[i] += local_acc_y;
            acc_z[i] += local_acc_z;
        }
    }
}

void ParticleSystem::compute_gas_dm_pp_forces(const GasGrid& gas,
                                              Grid3D& grav_x, Grid3D& grav_y,
                                              Grid3D& grav_z,
                                              const Config& config,
                                              Diagnostics& diag) {
    if (config.hydro_method != HydroMethod::Eulerian) return;

    // Shared setup variables
    int search_radius =
        config.use_PM
            ? static_cast<int>(ceil(config.cutoff_radius / config.cell_size))
            : config.mesh_size / 2;
    const int N = config.mesh_size;
    const double domain_size = config.domain_size;
    const double G = config.G;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t n_parts = num_particles;

    const double cell_vol = config.cell_volume;
    const double cell_size = config.cell_size;

#ifdef USE_GPU
    if (config.enable_GPU) {
        // ========================================================================
        // GPU IMPLEMENTATION
        // ========================================================================
        double* d_px = pos_x.data();
        double* d_py = pos_y.data();
        double* d_pz = pos_z.data();
        double* d_m = mass.data();
        double* d_ax = acc_x.data();
        double* d_ay = acc_y.data();
        double* d_az = acc_z.data();

        const double* d_gas_rho = gas.get_density().raw_data();
        double* d_grav_x = grav_x.data.data();
        double* d_grav_y = grav_y.data.data();
        double* d_grav_z = grav_z.data.data();

#pragma omp target enter data map(                                    \
        to : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts], \
            d_m[0 : n_parts], d_ax[0 : n_parts], d_ay[0 : n_parts],   \
            d_az[0 : n_parts], d_gas_rho[0 : N * N * N],              \
            d_grav_x[0 : N * N * N], d_grav_y[0 : N * N * N],         \
            d_grav_z[0 : N * N * N])

#pragma omp target teams distribute parallel for
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = d_px[i], p1_y = d_py[i], p1_z = d_pz[i];
            int ix = static_cast<int>(p1_x / cell_size);
            int iy = static_cast<int>(p1_y / cell_size);
            int iz = static_cast<int>(p1_z / cell_size);

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
            double dm_mass = d_m[i];

            for (int dx_cell = -search_radius; dx_cell <= search_radius;
                 ++dx_cell) {
                for (int dy_cell = -search_radius; dy_cell <= search_radius;
                     ++dy_cell) {
                    for (int dz_cell = -search_radius; dz_cell <= search_radius;
                         ++dz_cell) {
                        int neighbor_ix = (((ix + dx_cell) % N) + N) % N;
                        int neighbor_iy = (((iy + dy_cell) % N) + N) % N;
                        int neighbor_iz = (((iz + dz_cell) % N) + N) % N;

                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        double cell_px = (neighbor_ix + 0.5) * cell_size;
                        double cell_py = (neighbor_iy + 0.5) * cell_size;
                        double cell_pz = (neighbor_iz + 0.5) * cell_size;
                        double gas_mass = d_gas_rho[cell_idx] * cell_vol;

                        double dx =
                            periodic_displacement(cell_px - p1_x, domain_size);
                        double dy =
                            periodic_displacement(cell_py - p1_y, domain_size);
                        double dz =
                            periodic_displacement(cell_pz - p1_z, domain_size);

                        double dist_sq = dx * dx + dy * dy + dz * dz;
                        if (use_pm && dist_sq > cutoff_sq) continue;
                        if (dist_sq == 0.0)
                            continue;  // Prevent NaN if particle is dead center

                        double r = std::sqrt(dist_sq);
                        double q = r / cell_size;

                        // Monaghan 1992 Cubic Spline Weighting (Models internal
                        // fluid Shell Theorem)
                        double spline_weight = 1.0;
                        if (q < 1.0) {
                            double q3 = q * q * q;
                            if (q < 0.5) {
                                double q5 = q3 * q * q;
                                double q6 = q5 * q;
                                spline_weight = (32.0 / 3.0) * q3 -
                                                (192.0 / 5.0) * q5 + 32.0 * q6;
                            } else {
                                double q4 = q3 * q;
                                double q5 = q4 * q;
                                double q6 = q5 * q;
                                spline_weight = -1.0 / 15.0 +
                                                (64.0 / 3.0) * q3 - 48.0 * q4 +
                                                (192.0 / 5.0) * q5 -
                                                (32.0 / 3.0) * q6;
                            }
                        }

                        // True Newtonian magnitude, scaled smoothly down to 0
                        // inside the cell
                        double a_pp_mag =
                            (G * gas_mass / dist_sq) * spline_weight;
                        double a_gas_mag =
                            (G * dm_mass / dist_sq) * spline_weight;

                        if (use_pm) {
                            double r_scaled = r / (2.0 * r_s);
                            double pm_filter =
                                (std::erfc(r_scaled) +
                                 (r / (std::sqrt(M_PI) * r_s)) *
                                     std::exp(-r_scaled * r_scaled));
                            a_pp_mag *= pm_filter;
                            a_gas_mag *= pm_filter;
                        }

                        // Accumulate locally for DM particle
                        local_acc_x += a_pp_mag * dx / r;
                        local_acc_y += a_pp_mag * dy / r;
                        local_acc_z += a_pp_mag * dz / r;

                        // Write equal and opposite reaction to Gas grid
                        double gas_ax = -a_gas_mag * dx / r;
                        double gas_ay = -a_gas_mag * dy / r;
                        double gas_az = -a_gas_mag * dz / r;

#pragma omp atomic
                        d_grav_x[cell_idx] += gas_ax;
#pragma omp atomic
                        d_grav_y[cell_idx] += gas_ay;
#pragma omp atomic
                        d_grav_z[cell_idx] += gas_az;
                    }
                }
            }

            d_ax[i] += local_acc_x;
            d_ay[i] += local_acc_y;
            d_az[i] += local_acc_z;
        }
#pragma omp target exit data map(                                         \
        from : d_ax[0 : n_parts], d_ay[0 : n_parts], d_az[0 : n_parts],   \
            d_grav_x[0 : N * N * N], d_grav_y[0 : N * N * N],             \
            d_grav_z[0 : N * N * N])                                      \
    map(delete : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts], \
            d_m[0 : n_parts])
    } else
#endif
    {
        // ========================================================================
        // CPU IMPLEMENTATION
        // ========================================================================
#pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < n_parts; ++i) {
            double p1_x = pos_x[i], p1_y = pos_y[i], p1_z = pos_z[i];

            int ix = static_cast<int>(p1_x / cell_size);
            int iy = static_cast<int>(p1_y / cell_size);
            int iz = static_cast<int>(p1_z / cell_size);

            double local_acc_x = 0.0, local_acc_y = 0.0, local_acc_z = 0.0;
            double dm_mass = mass[i];

            for (int dx_cell = -search_radius; dx_cell <= search_radius;
                 ++dx_cell) {
                for (int dy_cell = -search_radius; dy_cell <= search_radius;
                     ++dy_cell) {
                    for (int dz_cell = -search_radius; dz_cell <= search_radius;
                         ++dz_cell) {
                        int neighbor_ix = (((ix + dx_cell) % N) + N) % N;
                        int neighbor_iy = (((iy + dy_cell) % N) + N) % N;
                        int neighbor_iz = (((iz + dz_cell) % N) + N) % N;

                        int cell_idx =
                            neighbor_iz * N * N + neighbor_iy * N + neighbor_ix;

                        double cell_px = (neighbor_ix + 0.5) * cell_size;
                        double cell_py = (neighbor_iy + 0.5) * cell_size;
                        double cell_pz = (neighbor_iz + 0.5) * cell_size;
                        double gas_mass =
                            gas.get_density().data[cell_idx] * cell_vol;

                        double dx =
                            periodic_displacement(cell_px - p1_x, domain_size);
                        double dy =
                            periodic_displacement(cell_py - p1_y, domain_size);
                        double dz =
                            periodic_displacement(cell_pz - p1_z, domain_size);

                        double dist_sq = dx * dx + dy * dy + dz * dz;
                        if (use_pm && dist_sq > cutoff_sq) continue;
                        if (dist_sq == 0.0)
                            continue;  // Prevent NaN if particle is dead center

                        double r = std::sqrt(dist_sq);
                        double q = r / cell_size;

                        // Monaghan 1992 Cubic Spline Weighting
                        double spline_weight = 1.0;
                        if (q < 1.0) {
                            double q3 = q * q * q;
                            if (q < 0.5) {
                                double q5 = q3 * q * q;
                                double q6 = q5 * q;
                                spline_weight = (32.0 / 3.0) * q3 -
                                                (192.0 / 5.0) * q5 + 32.0 * q6;
                            } else {
                                double q4 = q3 * q;
                                double q5 = q4 * q;
                                double q6 = q5 * q;
                                spline_weight = -1.0 / 15.0 +
                                                (64.0 / 3.0) * q3 - 48.0 * q4 +
                                                (192.0 / 5.0) * q5 -
                                                (32.0 / 3.0) * q6;
                            }
                        }

                        // True Newtonian magnitude, scaled smoothly down to 0
                        // inside the cell
                        double a_pp_mag =
                            (G * gas_mass / dist_sq) * spline_weight;
                        double a_gas_mag =
                            (G * dm_mass / dist_sq) * spline_weight;

                        if (use_pm) {
                            double r_scaled = r / (2.0 * r_s);
                            double pm_filter =
                                (std::erfc(r_scaled) +
                                 (r / (std::sqrt(M_PI) * r_s)) *
                                     std::exp(-r_scaled * r_scaled));
                            a_pp_mag *= pm_filter;
                            a_gas_mag *= pm_filter;
                        }

                        // Accumulate locally for DM particle
                        local_acc_x += a_pp_mag * dx / r;
                        local_acc_y += a_pp_mag * dy / r;
                        local_acc_z += a_pp_mag * dz / r;

                        // Write equal and opposite reaction to Gas grid
                        double gas_ax = -a_gas_mag * dx / r;
                        double gas_ay = -a_gas_mag * dy / r;
                        double gas_az = -a_gas_mag * dz / r;

#pragma omp atomic
                        grav_x.data[cell_idx] += gas_ax;
#pragma omp atomic
                        grav_y.data[cell_idx] += gas_ay;
#pragma omp atomic
                        grav_z.data[cell_idx] += gas_az;
                    }
                }
            }

            acc_x[i] += local_acc_x;
            acc_y[i] += local_acc_y;
            acc_z[i] += local_acc_z;
        }
    }
}

double ParticleSystem::get_gravity_timestep(const Config& config) const {
    if (num_particles == 0) return std::numeric_limits<double>::infinity();

    double epsilon = std::sqrt(config.softening_squared);
    double a_max = std::sqrt(max_accel_sq);
    double dt_grav = std::sqrt(epsilon / a_max);

    return dt_grav * config.gravity_accuracy_eta;
}