#include "particles.h"

#include <omp.h>

#include "cic.h"
#include "diagnostics.h"
#include "kernels.h"
#include "math_utils.h"
#include "mfm.h"

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

    is_active.reserve(config.num_dm_particles);
    time_bin.reserve(config.num_dm_particles);
    dt_step.reserve(config.num_dm_particles);
    t_current.reserve(config.num_dm_particles);
    t_end.reserve(config.num_dm_particles);
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

    time_bin.push_back(0);
    dt_step.push_back(0.0);
    t_current.push_back(0.0);
    t_end.push_back(0.0);
    is_active.push_back(1);

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

    std::vector<int> new_time_bin(num_particles);
    std::vector<double> new_dt_step(num_particles),
        new_t_current(num_particles), new_t_end(num_particles);
    std::vector<uint8_t> new_is_active(num_particles);

    // CIC_Data is only needed for PM gravity
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

        new_time_bin[i] = time_bin[src];
        new_dt_step[i] = dt_step[src];
        new_t_current[i] = t_current[src];
        new_t_end[i] = t_end[src];
        new_is_active[i] = is_active[src];

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

    time_bin = std::move(new_time_bin);
    dt_step = std::move(new_dt_step);
    t_current = std::move(new_t_current);
    t_end = std::move(new_t_end);
    is_active = std::move(new_is_active);

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
    // Pre-scale G so output forces are comoving accelerations
    const double G = config.G / (a * a * a);
    const double soft_sq = config.softening_squared;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;
    const size_t num_nodes = 2 * n_parts - 1;
    const double base_soft = std::sqrt(soft_sq);
    const double spline_h = 2.8 * base_soft;

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
    uint8_t* d_is_active = is_active.data();

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();

#pragma omp target enter data map(                                           \
        to : d_px[0 : n_parts], d_py[0 : n_parts], d_pz[0 : n_parts],        \
            d_m[0 : n_parts], d_bvh_nodes[0 : num_nodes], d_ax[0 : n_parts], \
            d_ay[0 : n_parts], d_az[0 : n_parts], d_is_active[0 : n_parts])

    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t i = 0; i < n_parts; ++i) {
        if (!d_is_active[i]) continue;

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

                double r = std::sqrt(dist_sq + 1e-24);
                double force_mag_over_r = 0.0;

                // Cubic Spline Gravity
                double dphi_dr, dummy_dphi_dh;
                Kernels::gravity_derivatives(r, spline_h, dphi_dr,
                                             dummy_dphi_dh);

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
            d_is_active[0 : n_parts])

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

void ParticleSystem::compute_cross_pp_forces(double a,
                                             const GasParticleSystem& gas,
                                             const Config& config,
                                             Diagnostics& diag) {
    if (num_particles == 0 || gas.num_particles == 0) return;

    const double domain_size = config.domain_size;

    // Pre-scale G so all output forces are comoving accelerations
    const double G = config.G / (a * a * a);
    const double soft_sq = config.softening_squared;
    const double cutoff_sq = config.cutoff_radius_squared;
    const double r_s = config.PM_smoothing_cells * config.cell_size;
    const bool use_pm = config.use_PM;

    const size_t n_dm = num_particles;
    const size_t n_gas = gas.num_particles;
    const size_t gas_num_nodes = 2 * n_gas - 1;

    // Search radius is the PM cutoff. If not using PM, it's infinite
    const double search_sq =
        use_pm ? cutoff_sq : std::numeric_limits<double>::infinity();
    const double base_soft = std::sqrt(soft_sq);
    const double spline_equivalent_h = 2.8 * base_soft;

    // Extract DM pointers
    double* d_dm_px = pos_x.data();
    double* d_dm_py = pos_y.data();
    double* d_dm_pz = pos_z.data();
    uint8_t* d_dm_is_active = is_active.data();
    double* d_dm_ax = acc_x.data();
    double* d_dm_ay = acc_y.data();
    double* d_dm_az = acc_z.data();

    // Extract Gas pointers (read-only)
    const double* d_gas_px = gas.pos_x.data();
    const double* d_gas_py = gas.pos_y.data();
    const double* d_gas_pz = gas.pos_z.data();
    const double* d_gas_m = gas.mass.data();
    const double* d_gas_h = gas.h.data();
    const double* d_gas_zeta = gas.zeta.data();
    const BVHNode* d_gas_bvh = gas.bvh_nodes.data();

#ifdef USE_GPU
    auto start_transfer = std::chrono::high_resolution_clock::now();
#pragma omp target enter data map(                                          \
        to : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],       \
            d_dm_is_active[0 : n_dm], d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], \
            d_dm_az[0 : n_dm], d_gas_px[0 : n_gas], d_gas_py[0 : n_gas],    \
            d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas], d_gas_h[0 : n_gas],    \
            d_gas_zeta[0 : n_gas], d_gas_bvh[0 : gas_num_nodes])
    auto end_transfer = std::chrono::high_resolution_clock::now();
    auto start_compute = std::chrono::high_resolution_clock::now();
#endif

// Swap the OpenMP execution pragma based on the compile-time target
#ifdef USE_GPU
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (size_t i = 0; i < n_dm; ++i) {
        if (!d_dm_is_active[i]) continue;  // Skip sleeping DM particles

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
                                     domain_size) +
                min_periodic_dist_sq(p1_y, node.bbox.min_y, node.bbox.max_y,
                                     domain_size) +
                min_periodic_dist_sq(p1_z, node.bbox.min_z, node.bbox.max_z,
                                     domain_size);

            if (aabb_dist_sq > search_sq) continue;  // Prune branch completely

            if (node.particle_idx != -1) {
                // Exact P-P interaction at the leaf
                int j = node.particle_idx;

                double dx = p1_x - d_gas_px[j];
                if (dx > 0.5 * domain_size)
                    dx -= domain_size;
                else if (dx < -0.5 * domain_size)
                    dx += domain_size;

                double dy = p1_y - d_gas_py[j];
                if (dy > 0.5 * domain_size)
                    dy -= domain_size;
                else if (dy < -0.5 * domain_size)
                    dy += domain_size;

                double dz = p1_z - d_gas_pz[j];
                if (dz > 0.5 * domain_size)
                    dz -= domain_size;
                else if (dz < -0.5 * domain_size)
                    dz += domain_size;

                // Flip to (gas - dm) direction so force pulls DM toward Gas
                dx = -dx;
                dy = -dy;
                dz = -dz;

                double dist_sq = dx * dx + dy * dy + dz * dz;

                if (use_pm && dist_sq > cutoff_sq) continue;

                double r = std::sqrt(dist_sq + 1e-24);
                double m_gas = d_gas_m[j];
                double h_gas = d_gas_h[j];
                double zeta_gas = d_gas_zeta[j];

                double dphi_dr_i, dummy_dphi_dh, dphi_dr_j, dW_dr_j;

                // DM kernel derivatives (particle i)
                Kernels::gravity_derivatives(r, spline_equivalent_h, dphi_dr_i,
                                             dummy_dphi_dh);

                // Gas kernel derivatives (particle j)
                Kernels::adaptive_gravity_terms(r, h_gas, dphi_dr_j, dW_dr_j);

                // Force magnitude calculation (only Gas contributes zeta, DM
                // softening is fixed)
                double force_mag_over_r =
                    (G / 2.0) *
                    ((dphi_dr_i + dphi_dr_j) + (zeta_gas * dW_dr_j) / m_gas) /
                    r;

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
#pragma omp target exit data map(from : d_dm_ax[0 : n_dm], d_dm_ay[0 : n_dm], \
                                     d_dm_az[0 : n_dm])                       \
    map(delete : d_dm_px[0 : n_dm], d_dm_py[0 : n_dm], d_dm_pz[0 : n_dm],     \
            d_dm_is_active[0 : n_dm], d_gas_px[0 : n_gas],                    \
            d_gas_py[0 : n_gas], d_gas_pz[0 : n_gas], d_gas_m[0 : n_gas],     \
            d_gas_h[0 : n_gas], d_gas_zeta[0 : n_gas],                        \
            d_gas_bvh[0 : gas_num_nodes])
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
    if (config.enable_individual_timesteps) {
        // Bootstrap: On the first step, request safe step to initialize
        // forces
        if (global_time == 0.0) {
            double dt_boot = config.fixed_dt;
            while (dt_boot > 1e-6) dt_boot *= 0.5;
            return dt_boot;
        }

        // Return the time until the NEXT active DM particle finishes its step
        double next_time = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < num_particles; ++i) {
            // Only look at future bins ignoring float precision drift matches
            if (t_end[i] > global_time + 1e-10) {
                next_time = std::min(next_time, t_end[i]);
            }
        }

        double dt_next = next_time - global_time;
        return (dt_next > 1e-10) ? dt_next
                                 : std::numeric_limits<double>::infinity();
    } else {
        double epsilon = std::sqrt(config.softening_squared);
        double a_max = std::sqrt(max_accel_sq);
        double dt_grav = std::sqrt(epsilon / a_max);

        return dt_grav * config.gravity_accuracy_eta;
    }
}

void ParticleSystem::sync_and_activate(double dt, const Config& config) {
    if (!config.enable_individual_timesteps) {
        global_time += dt;
        for (size_t i = 0; i < num_particles; ++i) {
            is_active[i] = 1;
            dt_step[i] = dt;
            t_end[i] = global_time;  // Keep timeline synced
        }
        return;
    }

    double new_global_time = global_time + dt;

    // Advance global clock
    global_time = new_global_time;

    // Flag active particles
    for (size_t i = 0; i < num_particles; ++i) {
        // A particle is active if its block step ends at the new global time
        is_active[i] = (std::abs(t_end[i] - global_time) < 1e-10) ? 1 : 0;

        // BOOTSTRAP: If it's the very first step (t=0), everyone is active
        if (global_time <= dt + 1e-12 && t_end[i] == 0.0) {
            is_active[i] = 1;
            dt_step[i] =
                dt;  // Give them a tiny initial dt so kicks work properly
            t_end[i] = global_time;
        }
    }
}

void ParticleSystem::update_particle_timesteps(double dt_max,
                                               const Config& config) {
    if (num_particles == 0) return;

    if (!config.enable_individual_timesteps) {
        // Under global timesteps, all particles just take the macro step
        for (size_t i = 0; i < num_particles; ++i) {
            dt_step[i] = dt_max;
            time_bin[i] = 0;
            // t_end is already synced in sync_and_activate for the global
            // scheme
        }
        return;
    }

    double epsilon = std::sqrt(config.softening_squared);

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        // ONLY update the clock for particles that are currently active
        if (!is_active[i]) continue;

        double dt_ideal = dt_max;

        // Gravity Timestep
        double a_mag = std::sqrt(acc_x[i] * acc_x[i] + acc_y[i] * acc_y[i] +
                                 acc_z[i] * acc_z[i]);

        if (a_mag > 0.0) {
            double dt_grav =
                std::sqrt(epsilon / a_mag) * config.gravity_accuracy_eta;
            dt_ideal = std::min(dt_ideal, dt_grav);
        }

        // POWER-OF-TWO BINNING
        int n = 0;
        double dt_bin = dt_max;

        // Halve the timestep until it safely bounds the ideal timestep
        // constraint
        while (dt_bin > dt_ideal) {
            dt_bin *= 0.5;
            n++;
        }

        // Commit the new hierarchical timestep configuration
        time_bin[i] = n;
        dt_step[i] = dt_bin;

        // Project the end time of the new step based on the current global
        // clock
        t_end[i] = global_time + dt_bin;
    }
}