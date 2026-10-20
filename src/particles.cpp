#include "particles.h"

#include <omp.h>

#include "config.h"

ParticleSystem::ParticleSystem(const Config& config)
    : grid_rho(config.mesh_size), cic_data(config.num_dm_particles) {
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
    active_indices.reserve(config.num_dm_particles);
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

    active_indices.push_back(num_particles);
    num_active++;

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

    active_indices.clear();
    for (size_t i = 0; i < num_particles; ++i) {
        if (is_active[i]) {
            active_indices.push_back(i);
        }
    }
    num_active = active_indices.size();
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

void ParticleSystem::sync_and_activate(double dt, const Config& config) {
    if (!config.individual_particle_timesteps) {
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

    active_indices.clear();
    for (size_t i = 0; i < num_particles; ++i) {
        if (is_active[i]) {
            active_indices.push_back(i);
        }
    }
    num_active = active_indices.size();
}

void ParticleSystem::update_particle_timesteps(double dt_max,
                                               const Config& config) {
    if (num_particles == 0) return;

    if (!config.individual_particle_timesteps) {
        // Under global timesteps, all particles just take the macro step
        for (size_t i = 0; i < num_particles; ++i) {
            dt_step[i] = dt_max;
            time_bin[i] = 0;
            // t_end is already synced in sync_and_activate for the global
            // scheme
        }
        return;
    }

    double epsilon = config.softening;

    // Prevent OpenMP from deadlocking on 0-trip dynamic loops
    // during micro-cycles where all particles are asleep
    if (num_active == 0) return;

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t k = 0; k < num_active; ++k) {
        size_t i = active_indices[k];

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

        // A particle can ONLY adopt a timestep bin if the current global clock
        // is aligned to an integer multiple of that bin's boundary.
        // If not aligned, it must drop to a smaller bin that IS aligned
        double ratio = global_time / dt_bin;
        while (std::abs(ratio - std::round(ratio)) > 1e-5 && n < 30) {
            dt_bin *= 0.5;
            n++;
            ratio = global_time / dt_bin;
        }

        time_bin[i] = n;
        dt_step[i] = dt_bin;

        // Project the end time of the new step
        t_end[i] = global_time + dt_bin;
    }
}
