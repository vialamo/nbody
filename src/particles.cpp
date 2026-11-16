#include "particles.h"

#include <omp.h>

#include <cmath>

#include "config.h"

ParticleSystem::ParticleSystem(const Config& config)
    : cic_data(config.num_dm_particles) {
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

    ti_step.push_back(0);
    ti_current.push_back(0);
    ti_end.push_back(0);

    is_active.push_back(1);
    active_indices.push_back(num_particles);

    num_particles++;
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

    std::vector<uint8_t> new_is_active(num_particles);

    std::vector<uint64_t> new_ti_step(num_particles);
    std::vector<uint64_t> new_ti_current(num_particles);
    std::vector<uint64_t> new_ti_end(num_particles);

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

        new_ti_step[i] = ti_step[src];
        new_ti_current[i] = ti_current[src];
        new_ti_end[i] = ti_end[src];

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

    is_active = std::move(new_is_active);

    ti_step = std::move(new_ti_step);
    ti_current = std::move(new_ti_current);
    ti_end = std::move(new_ti_end);

    if (has_cic) cic_data = std::move(new_cic);

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
                                               double softening_length) {
    if (num_particles == 0) return;

    double dt_max = config.fixed_dt;

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
    // during micro-cycles where all particles are asleep
    size_t num_active = get_num_active();
    if (num_active == 0) return;

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t k = 0; k < num_active; ++k) {
        size_t i = active_indices[k];

        double dt_ideal = dt_max;

        // Gravity Timestep
        double a_mag = std::sqrt(acc_x[i] * acc_x[i] + acc_y[i] * acc_y[i] +
                                 acc_z[i] * acc_z[i]);

        if (a_mag > 0.0) {
            double dt_grav = std::sqrt(
                (2.0 * config.gravity_accuracy_eta * epsilon) / a_mag);
            dt_ideal = std::min(dt_ideal, dt_grav);
        }

        // Integer Power-of-two binning
        int n = 0;
        uint64_t ti_bin = TimeConstants::TIMEBASE;
        double dt_bin = config.fixed_dt;  // Always anchor to the fixed_dt

        // Halve the timestep until it bounds the ideal timestep
        while (dt_bin > dt_ideal && n < TimeConstants::MAX_TIMESTEP_BIN_DEPTH) {
            ti_bin >>= 1;  // Fast integer divide by 2
            dt_bin *= 0.5;
            n++;
        }

        // A particle can ONLY adopt a timestep bin if the current global clock
        // is aligned to an integer multiple of that bin's boundary
        while ((ti_global_time % ti_bin) != 0 &&
               n < TimeConstants::MAX_TIMESTEP_BIN_DEPTH) {
            ti_bin >>= 1;
            dt_bin *= 0.5;
            n++;
        }

        // Update the Integer Source of Truth
        ti_step[i] = ti_bin;
        ti_end[i] = ti_global_time + ti_bin;
    }
}
