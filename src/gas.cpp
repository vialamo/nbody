#include "gas.h"

#include "config.h"
#include "hydro.h"
#include "integrator.h"

GasParticleSystem::GasParticleSystem(const Config& config) {
    cic_data.reserve(config.num_gas_particles);

    pos_x.reserve(config.num_gas_particles);
    pos_y.reserve(config.num_gas_particles);
    pos_z.reserve(config.num_gas_particles);

    vel_x.reserve(config.num_gas_particles);
    vel_y.reserve(config.num_gas_particles);
    vel_z.reserve(config.num_gas_particles);

    acc_x.reserve(config.num_gas_particles);
    acc_y.reserve(config.num_gas_particles);
    acc_z.reserve(config.num_gas_particles);

    mass.reserve(config.num_gas_particles);

    hydro_acc_x.reserve(config.num_gas_particles);
    hydro_acc_y.reserve(config.num_gas_particles);
    hydro_acc_z.reserve(config.num_gas_particles);

    h.reserve(config.num_gas_particles);
    rho.reserve(config.num_gas_particles);
    pressure.reserve(config.num_gas_particles);
    total_energy.reserve(config.num_gas_particles);
    u.reserve(config.num_gas_particles);
    du_dt.reserve(config.num_gas_particles);
    de_dt.reserve(config.num_gas_particles);
    metal_frac.reserve(config.num_gas_particles);

    grad_rho.reserve(config.num_gas_particles);
    grad_vx.reserve(config.num_gas_particles);
    grad_vy.reserve(config.num_gas_particles);
    grad_vz.reserve(config.num_gas_particles);
    grad_p.reserve(config.num_gas_particles);
    B_matrix.reserve(config.num_gas_particles);

    entropy.reserve(config.num_gas_particles);
    max_rel_ke.reserve(config.num_gas_particles);

    vel_x_pred.reserve(config.num_gas_particles);
    vel_y_pred.reserve(config.num_gas_particles);
    vel_z_pred.reserve(config.num_gas_particles);
    pressure_pred.reserve(config.num_gas_particles);

    v_sig_max.reserve(config.num_gas_particles);
    cond_num.reserve(config.num_gas_particles);
    n_enc_final.reserve(config.num_gas_particles);

    ti_step.reserve(config.num_gas_particles);
    ti_current.reserve(config.num_gas_particles);
    ti_end.reserve(config.num_gas_particles);

    is_active.reserve(config.num_gas_particles);
    needs_wakeup.reserve(config.num_gas_particles);
    active_indices.reserve(config.num_gas_particles);

    size_t num_nodes =
        (config.num_gas_particles > 0) ? (2 * config.num_gas_particles - 1) : 0;
    if (config.num_gas_particles > 0) {
        morton_codes.reserve(config.num_gas_particles);
        sorted_indices.reserve(config.num_gas_particles);
        bvh_nodes.reserve(num_nodes);
    }
}

void GasParticleSystem::add_particle(double px, double py, double pz, double vx,
                                     double vy, double vz, double m,
                                     double initial_u, double initial_h,
                                     double z_metal) {
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

    hydro_acc_x.push_back(0.0);
    hydro_acc_y.push_back(0.0);
    hydro_acc_z.push_back(0.0);

    h.push_back(initial_h);
    rho.push_back(0.0);
    pressure.push_back(0.0);
    u.push_back(initial_u);
    du_dt.push_back(0.0);
    de_dt.push_back(0.0);
    metal_frac.push_back(z_metal);
    B_matrix.push_back(Eigen::Matrix3d::Zero());

    entropy.push_back(0.0);
    max_rel_ke.push_back(0.0);

    vel_x_pred.push_back(vx);
    vel_y_pred.push_back(vy);
    vel_z_pred.push_back(vz);
    pressure_pred.push_back(0.0);

    v_sig_max.push_back(0.0);

    double initial_ke = 0.5 * (vx * vx + vy * vy + vz * vz);
    total_energy.push_back(initial_u + initial_ke);

    grad_rho.push_back(Eigen::Vector3d::Zero());
    grad_vx.push_back(Eigen::Vector3d::Zero());
    grad_vy.push_back(Eigen::Vector3d::Zero());
    grad_vz.push_back(Eigen::Vector3d::Zero());
    grad_p.push_back(Eigen::Vector3d::Zero());

    cond_num.push_back(0);
    n_enc_final.push_back(0);

    cic_data.push_back({});

    ti_step.push_back(0);
    ti_current.push_back(0);
    ti_end.push_back(0);

    is_active.push_back(1);  // All particles start active at t=0
    needs_wakeup.push_back(0);

    active_indices.push_back(num_particles);

    num_particles++;
}

void GasParticleSystem::sort_arrays(const std::vector<int>& sorted_indices) {
    std::vector<double> new_px(num_particles), new_py(num_particles),
        new_pz(num_particles);
    std::vector<double> new_vx(num_particles), new_vy(num_particles),
        new_vz(num_particles);
    std::vector<double> new_ax(num_particles), new_ay(num_particles),
        new_az(num_particles);
    std::vector<double> new_hax(num_particles), new_hay(num_particles),
        new_haz(num_particles);
    std::vector<double> new_m(num_particles), new_h(num_particles),
        new_rho(num_particles);
    std::vector<double> new_p(num_particles), new_u(num_particles),
        new_dudt(num_particles), new_dedt(num_particles),
        new_total_energy(num_particles);
    std::vector<double> new_metal(num_particles);
    std::vector<CIC_Data> new_cic(num_particles);
    std::vector<Eigen::Matrix3d> new_B(num_particles);
    std::vector<Eigen::Vector3d> new_grad_rho(num_particles),
        new_grad_vx(num_particles);
    std::vector<Eigen::Vector3d> new_grad_vy(num_particles),
        new_grad_vz(num_particles);
    std::vector<Eigen::Vector3d> new_grad_p(num_particles);

    std::vector<double> new_vxp(num_particles), new_vyp(num_particles),
        new_vzp(num_particles), new_ppred(num_particles);

    std::vector<double> new_entropy(num_particles);
    std::vector<double> new_max_rel_ke(num_particles);

    std::vector<double> new_cond_num(num_particles);
    std::vector<double> new_n_enc_final(num_particles);

    std::vector<uint8_t> new_is_active(num_particles);
    std::vector<uint8_t> new_needs_wakeup(num_particles);

    std::vector<uint64_t> new_ti_step(num_particles);
    std::vector<uint64_t> new_ti_current(num_particles);
    std::vector<uint64_t> new_ti_end(num_particles);

    std::vector<double> new_v_sig_max(num_particles);

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
        new_hax[i] = hydro_acc_x[src];
        new_hay[i] = hydro_acc_y[src];
        new_haz[i] = hydro_acc_z[src];
        new_h[i] = h[src];
        new_rho[i] = rho[src];
        new_p[i] = pressure[src];
        new_total_energy[i] = total_energy[src];
        new_u[i] = u[src];
        new_dudt[i] = du_dt[src];
        new_dedt[i] = de_dt[src];
        new_metal[i] = metal_frac[src];
        new_cic[i] = cic_data[src];
        new_B[i] = B_matrix[src];
        new_grad_rho[i] = grad_rho[src];
        new_grad_vx[i] = grad_vx[src];
        new_grad_vy[i] = grad_vy[src];
        new_grad_vz[i] = grad_vz[src];
        new_grad_p[i] = grad_p[src];

        new_vxp[i] = vel_x_pred[src];
        new_vyp[i] = vel_y_pred[src];
        new_vzp[i] = vel_z_pred[src];
        new_ppred[i] = pressure_pred[src];

        new_entropy[i] = entropy[src];
        new_max_rel_ke[i] = max_rel_ke[src];
        new_cond_num[i] = cond_num[src];
        new_n_enc_final[i] = n_enc_final[src];

        new_ti_step[i] = ti_step[src];
        new_ti_current[i] = ti_current[src];
        new_ti_end[i] = ti_end[src];

        new_is_active[i] = is_active[src];
        new_needs_wakeup[i] = needs_wakeup[src];
        new_v_sig_max[i] = v_sig_max[src];
    }

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
    hydro_acc_x = std::move(new_hax);
    hydro_acc_y = std::move(new_hay);
    hydro_acc_z = std::move(new_haz);
    h = std::move(new_h);
    rho = std::move(new_rho);
    pressure = std::move(new_p);
    total_energy = std::move(new_total_energy);
    u = std::move(new_u);
    du_dt = std::move(new_dudt);
    de_dt = std::move(new_dedt);
    metal_frac = std::move(new_metal);
    cic_data = std::move(new_cic);
    B_matrix = std::move(new_B);
    grad_rho = std::move(new_grad_rho);
    grad_vx = std::move(new_grad_vx);
    grad_vy = std::move(new_grad_vy);
    grad_vz = std::move(new_grad_vz);
    grad_p = std::move(new_grad_p);

    vel_x_pred = std::move(new_vxp);
    vel_y_pred = std::move(new_vyp);
    vel_z_pred = std::move(new_vzp);
    pressure_pred = std::move(new_ppred);

    entropy = std::move(new_entropy);
    max_rel_ke = std::move(new_max_rel_ke);
    cond_num = std::move(new_cond_num);
    n_enc_final = std::move(new_n_enc_final);

    ti_step = std::move(new_ti_step);
    ti_current = std::move(new_ti_current);
    ti_end = std::move(new_ti_end);

    is_active = std::move(new_is_active);
    needs_wakeup = std::move(new_needs_wakeup);
    v_sig_max = std::move(new_v_sig_max);

    active_indices.clear();
    for (size_t i = 0; i < num_particles; ++i) {
        if (is_active[i]) {
            active_indices.push_back(i);
        }
    }
}

void GasParticleSystem::build_lbvh(const Config& config) {
    if (num_particles == 0) return;

    // Generate the Morton codes and the index map
    LBVH::compute_morton_and_sort_indices(num_particles, config.domain_size,
                                          pos_x, pos_y, pos_z, morton_codes,
                                          sorted_indices);

    sort_arrays(sorted_indices);

    // Build the tree topology and compute Centers of Mass / Bounding Boxes
    LBVH::build_topology_and_aggregate(num_particles, pos_x, pos_y, pos_z, mass,
                                       &h, morton_codes, bvh_nodes);
}

void GasParticleSystem::sync_and_activate(Integrator& integrator,
                                          uint64_t ti_dt, double a, double H,
                                          const Config& config) {
    double timebase_to_float = config.timebase_to_float;
    if (!config.individual_particle_timesteps) {
        ti_global_time += ti_dt;
        for (size_t i = 0; i < num_particles; ++i) {
            is_active[i] = 1;
            ti_step[i] = ti_dt;
            ti_end[i] = ti_global_time;
            needs_wakeup[i] = 0;
        }
        return;
    }

    uint64_t new_ti_global_time = ti_global_time + ti_dt;

    // Process wakeups from previous cycle
    for (size_t i = 0; i < num_particles; ++i) {
        if (needs_wakeup[i]) {
            uint64_t ti_start_old = ti_end[i] - ti_step[i];
            uint64_t ti_kicked_to = ti_start_old + (ti_step[i] / 2);

            // If the particle was kicked past the current global time, rewind
            // it
            if (ti_kicked_to > ti_global_time) {
                uint64_t ti_excess = ti_kicked_to - ti_global_time;
                double dt_excess = ti_excess * timebase_to_float;

                Integrator::KickWork work = integrator.kick_single_gas_particle(
                    *this, i, -dt_excess, a, H, config);

                budget.gravitational_work += work.grav_work;
                budget.expansion_work += work.exp_work;
            }

            ti_step[i] = ti_dt;
            ti_end[i] = new_ti_global_time;
            needs_wakeup[i] = 0;
        }
    }

    // Advance global clock
    ti_global_time = new_ti_global_time;

    // Flag active particles using
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        is_active[i] = (ti_end[i] == ti_global_time) ? 1 : 0;

        // BOOTSTRAP: First step assignment
        if (ti_global_time == ti_dt && ti_end[i] == 0) {
            is_active[i] = 1;
            ti_step[i] = ti_dt;
            ti_end[i] = ti_global_time;
        }
    }

    active_indices.clear();
    for (size_t i = 0; i < num_particles; ++i) {
        if (is_active[i]) active_indices.push_back(i);
    }

    diag.active_particles_fraction +=
        num_particles > 0
            ? static_cast<double>(get_num_active()) / num_particles
            : 0;
    diag.active_particles_num_cycles++;
}

void GasParticleSystem::update_particle_timesteps(const Config& config,
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

        // CFL Timestep (Using cached signal velocity from compute_hydro_forces)
        if (v_sig_max[i] > HydroSolver::MIN_SIGNAL_VELOCITY) {
            double dt_cfl = (h[i] / v_sig_max[i]) * config.hydro_courant_factor;
            dt_ideal = std::min(dt_ideal, dt_cfl);
        }

        // Integer Power-of-two binning
        int n = 0;
        uint64_t ti_bin = TimeConstants::TIMEBASE;
        double dt_bin = config.fixed_dt;

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

        // Diagnostic tracker
        if (n == TimeConstants::MAX_TIMESTEP_BIN_DEPTH && ti_global_time > 0) {
#pragma omp atomic
            diag.particles_at_min_dt++;
        }

        // Update the Integer
        ti_step[i] = ti_bin;
        ti_end[i] = ti_global_time + ti_bin;
    }
}

void GasParticleSystem::reset_diagnostic_counters() {
    diag.active_particles_fraction = 0.0;
    diag.active_particles_num_cycles = 0;
    diag.ill_conditioned = 0;
    diag.clamped_h = 0;
    diag.neighbor_increased = 0;
}
