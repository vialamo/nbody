#include "particles.h"

#include <omp.h>

#include <cmath>

#include "config.h"

ParticleSystem::ParticleSystem(const Config& config) {
    cic_data.reserve(config.num_dm_particles);

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

    h.reserve(config.num_dm_particles);
    rho.reserve(config.num_dm_particles);
    omega.reserve(config.num_dm_particles);

    size_t num_nodes = 2 * config.num_dm_particles - 1;
    if (config.num_dm_particles > 0) {
        morton_codes.reserve(config.num_dm_particles);
        sorted_indices.reserve(config.num_dm_particles);
        bvh_nodes.reserve(num_nodes);
    }

    is_active.reserve(config.num_dm_particles);
    active_indices.reserve(config.num_dm_particles);

    ti_step.reserve(config.num_dm_particles);
    ti_current.reserve(config.num_dm_particles);
    ti_end.reserve(config.num_dm_particles);
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

    h.push_back(0.05);
    rho.push_back(0.0);
    omega.push_back(1.0);

    ti_step.push_back(0);
    ti_current.push_back(0);
    ti_end.push_back(0);

    cic_data.push_back({});

    is_active.push_back(1);
    active_indices.push_back(num_particles);

    num_particles++;
}

void ParticleSystem::sort_arrays(const std::vector<int>& sorted_indices) {
    if (num_particles == 0) return;

    // Doubles
    permute_array(pos_x, scratch_double, sorted_indices);
    permute_array(pos_y, scratch_double, sorted_indices);
    permute_array(pos_z, scratch_double, sorted_indices);
    permute_array(vel_x, scratch_double, sorted_indices);
    permute_array(vel_y, scratch_double, sorted_indices);
    permute_array(vel_z, scratch_double, sorted_indices);
    permute_array(acc_x, scratch_double, sorted_indices);
    permute_array(acc_y, scratch_double, sorted_indices);
    permute_array(acc_z, scratch_double, sorted_indices);
    permute_array(mass, scratch_double, sorted_indices);
    permute_array(h, scratch_double, sorted_indices);
    permute_array(rho, scratch_double, sorted_indices);
    permute_array(omega, scratch_double, sorted_indices);

    // Integers & flags
    permute_array(ti_step, scratch_uint64, sorted_indices);
    permute_array(ti_current, scratch_uint64, sorted_indices);
    permute_array(ti_end, scratch_uint64, sorted_indices);
    permute_array(is_active, scratch_uint8, sorted_indices);

    // Others
    permute_array(cic_data, scratch_cic, sorted_indices);

    active_indices.clear();
    for (size_t i = 0; i < num_particles; ++i) {
        if (is_active[i]) {
            active_indices.push_back(i);
        }
    }
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

void ParticleSystem::sync_and_activate(uint64_t ti_dt, const Config& config) {
    if (!config.individual_particle_timesteps) {
        ti_global_time += ti_dt;
        for (size_t i = 0; i < num_particles; ++i) {
            is_active[i] = 1;
            ti_step[i] = ti_dt;
            ti_end[i] = ti_global_time;
        }
        return;
    }

    uint64_t new_ti_global_time = ti_global_time + ti_dt;

    // Advance global clock
    ti_global_time = new_ti_global_time;

    // Flag active particles
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        // A particle is active if its block step ends at the new global time
        is_active[i] = (ti_end[i] == ti_global_time) ? 1 : 0;

        // BOOTSTRAP: If it's the very first step (t=0), everyone is active
        if (ti_global_time == ti_dt && ti_end[i] == 0) {
            is_active[i] = 1;
            ti_step[i] = ti_dt;
            ti_end[i] = ti_global_time;
        }
    }

    active_indices.clear();
    for (size_t i = 0; i < num_particles; ++i) {
        if (is_active[i]) {
            active_indices.push_back(i);
        }
    }
}

void ParticleSystem::update_particle_timesteps(const Config& config,
                                               double softening_length,
                                               double hubble_param) {
    if (num_particles == 0) return;

    double dt_max = config.fixed_dt;

    // Apply cosmological expansion limiter
    if (config.expanding_universe && hubble_param > 0.0) {
        dt_max = std::min(dt_max, config.max_expansion_step / hubble_param);
    }

    if (!config.individual_particle_timesteps) {
        // Under global timesteps, calculate the integer step from dt_max
        uint64_t ti_dt =
            static_cast<uint64_t>(dt_max / config.timebase_to_float);
        if (ti_dt == 0) ti_dt = 1;

        for (size_t i = 0; i < num_particles; ++i) {
            ti_step[i] = ti_dt;
            ti_end[i] = ti_global_time + ti_dt;
        }
        return;
    }

    double epsilon = softening_length;

    // Prevent OpenMP from deadlocking on 0-trip dynamic loops
    size_t num_active = get_num_active();
    if (num_active == 0) return;

    const double inv_timebase = 1.0 / config.timebase_to_float;
    const uint64_t max_aligned_ti =
        (ti_global_time > 0) ? (ti_global_time & (~ti_global_time + 1))
                             : TimeConstants::TIMEBASE;

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t k = 0; k < num_active; ++k) {
        size_t i = active_indices[k];

        double dt_ideal = dt_max;

        // Gravity Timestep
        double a_mag = std::sqrt(acc_x[i] * acc_x[i] + acc_y[i] * acc_y[i] +
                                 acc_z[i] * acc_z[i]);

        if (a_mag > 0.0) {
            double eps_i = epsilon;
            if (config.use_adaptive_softening_dm) {
                double eps_max = config.use_PM
                                     ? std::sqrt(config.cutoff_radius_squared)
                                     : std::numeric_limits<double>::infinity();
                eps_i = std::clamp(h[i], softening_length, eps_max);
            }
            double dt_grav =
                std::sqrt((2.0 * config.gravity_accuracy_eta * eps_i) / a_mag);
            dt_ideal = std::min(dt_ideal, dt_grav);
        }

        // Power-of-two binning
        uint64_t ti_bin;
        if (dt_ideal >= config.fixed_dt) {
            ti_bin = TimeConstants::TIMEBASE;
        } else {
            uint64_t ti_ideal = static_cast<uint64_t>(dt_ideal * inv_timebase);
            // Largest power of 2 <= ti_ideal (clamped to 1)
            ti_bin = (ti_ideal <= 1)
                         ? 1ULL
                         : (1ULL << (63 - __builtin_clzll(ti_ideal)));
        }

        // A particle can only adopt a timestep bin if the current global clock
        // is aligned to an integer multiple of that bin's boundary
        ti_bin = std::min(ti_bin, max_aligned_ti);

        // Update
        ti_step[i] = ti_bin;
        ti_end[i] = ti_global_time + ti_bin;
    }
}
