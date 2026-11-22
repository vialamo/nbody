#include "gas.h"

#include "config.h"
#include "hydro.h"
#include "integrator.h"

namespace {
inline void atomic_update_min_uint64(uint64_t* address, uint64_t val) {
    uint64_t old_val = *address;
    while (old_val > val) {
        uint64_t actual_old =
            __sync_val_compare_and_swap(address, old_val, val);
        if (actual_old == old_val) break;
        old_val = actual_old;
    }
}
}  // namespace

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
    changed_timestep.reserve(config.num_gas_particles);

    is_active.reserve(config.num_gas_particles);
    needs_wakeup.reserve(config.num_gas_particles);
    active_indices.reserve(config.num_gas_particles);

    hydro_dv_x.reserve(config.num_gas_particles);
    hydro_dv_y.reserve(config.num_gas_particles);
    hydro_dv_z.reserve(config.num_gas_particles);
    hydro_du.reserve(config.num_gas_particles);
    hydro_de.reserve(config.num_gas_particles);

    snap_vel_x.reserve(config.num_gas_particles);
    snap_vel_y.reserve(config.num_gas_particles);
    snap_vel_z.reserve(config.num_gas_particles);
    snap_u.reserve(config.num_gas_particles);
    snap_total_energy.reserve(config.num_gas_particles);
    snap_pressure.reserve(config.num_gas_particles);
    snap_entropy.reserve(config.num_gas_particles);

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
    changed_timestep.push_back(0);

    is_active.push_back(1);  // All particles start active at t=0
    needs_wakeup.push_back(0);

    active_indices.push_back(num_particles);

    hydro_dv_x.push_back(0.0);
    hydro_dv_y.push_back(0.0);
    hydro_dv_z.push_back(0.0);
    hydro_du.push_back(0.0);
    hydro_de.push_back(0.0);

    snap_vel_x.push_back(0.0);
    snap_vel_y.push_back(0.0);
    snap_vel_z.push_back(0.0);
    snap_u.push_back(0.0);
    snap_total_energy.push_back(0.0);
    snap_pressure.push_back(0.0);
    snap_entropy.push_back(0.0);

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
    std::vector<uint8_t> new_changed_timestep(num_particles);

    std::vector<double> new_v_sig_max(num_particles);

    std::vector<double> new_hydro_dv_x(num_particles);
    std::vector<double> new_hydro_dv_y(num_particles);
    std::vector<double> new_hydro_dv_z(num_particles);
    std::vector<double> new_hydro_du(num_particles);
    std::vector<double> new_hydro_de(num_particles);

    std::vector<double> new_snap_vel_x(num_particles);
    std::vector<double> new_snap_vel_y(num_particles);
    std::vector<double> new_snap_vel_z(num_particles);
    std::vector<double> new_snap_u(num_particles);
    std::vector<double> new_snap_total_energy(num_particles);
    std::vector<double> new_snap_pressure(num_particles);
    std::vector<double> new_snap_entropy(num_particles);

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
        new_changed_timestep[i] = changed_timestep[src];

        new_is_active[i] = is_active[src];
        new_needs_wakeup[i] = needs_wakeup[src];
        new_v_sig_max[i] = v_sig_max[src];

        new_hydro_dv_x[i] = hydro_dv_x[src];
        new_hydro_dv_y[i] = hydro_dv_y[src];
        new_hydro_dv_z[i] = hydro_dv_z[src];
        new_hydro_du[i] = hydro_du[src];
        new_hydro_de[i] = hydro_de[src];

        new_snap_vel_x[i] = snap_vel_x[src];
        new_snap_vel_y[i] = snap_vel_y[src];
        new_snap_vel_z[i] = snap_vel_z[src];
        new_snap_u[i] = snap_u[src];
        new_snap_total_energy[i] = snap_total_energy[src];
        new_snap_pressure[i] = snap_pressure[src];
        new_snap_entropy[i] = snap_entropy[src];
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
    changed_timestep = std::move(new_changed_timestep);

    is_active = std::move(new_is_active);
    needs_wakeup = std::move(new_needs_wakeup);
    v_sig_max = std::move(new_v_sig_max);

    hydro_dv_x = std::move(new_hydro_dv_x);
    hydro_dv_y = std::move(new_hydro_dv_y);
    hydro_dv_z = std::move(new_hydro_dv_z);
    hydro_du = std::move(new_hydro_du);
    hydro_de = std::move(new_hydro_de);

    snap_vel_x = std::move(new_snap_vel_x);
    snap_vel_y = std::move(new_snap_vel_y);
    snap_vel_z = std::move(new_snap_vel_z);
    snap_u = std::move(new_snap_u);
    snap_total_energy = std::move(new_snap_total_energy);
    snap_pressure = std::move(new_snap_pressure);
    snap_entropy = std::move(new_snap_entropy);

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

    enforce_timestep_restrictions(ti_dt, config);

    uint64_t new_ti_global_time = ti_global_time + ti_dt;

    // Process wakeups from previous cycle
    for (size_t i = 0; i < num_particles; ++i) {
        if (needs_wakeup[i]) {
#define DYNAMIC_WAKEUP
#ifdef DYNAMIC_WAKEUP
            uint64_t ti_start_old = ti_end[i] - ti_step[i];

            // Restore state to the beginning of sleep (ti_start_old)
            vel_x[i] = snap_vel_x[i];
            vel_y[i] = snap_vel_y[i];
            vel_z[i] = snap_vel_z[i];
            u[i] = snap_u[i];
            total_energy[i] = snap_total_energy[i];
            pressure[i] = snap_pressure[i];
            entropy[i] = snap_entropy[i];

            // Calculate exact time asleep to catch up to current clock
            uint64_t ti_actual_asleep = ti_global_time - ti_start_old;
            double dt_actual = ti_actual_asleep * config.timebase_to_float;

            // Discard the accumulated impulses
            hydro_dv_x[i] = 0.0;
            hydro_dv_y[i] = 0.0;
            hydro_dv_z[i] = 0.0;
            hydro_du[i] = 0.0;
            hydro_de[i] = 0.0;

            integrator.kick_wakeup_catchup(*this, i, dt_actual, a, H, config);

            changed_timestep[i] = (ti_step[i] == ti_dt ? 0 : 1);
            ti_step[i] = ti_dt;
            ti_end[i] = new_ti_global_time;

            diag.forced_wakeup++;
#endif
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
            changed_timestep[i] = 1;
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

        // A particle can only adopt a timestep bin if the current global clock
        // is aligned to an integer multiple of that bin's boundary
        while ((ti_global_time % ti_bin) != 0 &&
               n < TimeConstants::MAX_TIMESTEP_BIN_DEPTH) {
            ti_bin >>= 1;
            n++;
        }

        // Diagnostic tracker
        if (n == TimeConstants::MAX_TIMESTEP_BIN_DEPTH && ti_global_time > 0) {
#pragma omp atomic
            diag.particles_at_min_dt++;
        }

        // Update
        changed_timestep[i] = (ti_step[i] == ti_bin ? 0 : 1);
        ti_step[i] = ti_bin;
        ti_end[i] = ti_global_time + ti_bin;
    }
}

void GasParticleSystem::enforce_timestep_restrictions(uint64_t ti_dt,
                                                      const Config& config) {
    if (!config.individual_particle_timesteps || ti_dt == 0) return;

    // Enforce the timestep limiter on neighbors
    // Saitoh & Makino (2009) Timestep Limiter
    const std::vector<uint64_t> calc_ti_step = ti_step;

    // The absolute maximum ti_step any particle can possibly have
    const uint64_t max_system_ti_step =
        static_cast<uint64_t>(config.fixed_dt / config.timebase_to_float);

#pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < num_particles; ++i) {
        // Process only particles that will be active and changed their
        // timestep
        if (ti_global_time > 0 &&
            (!changed_timestep[i] || ti_end[i] > ti_global_time + ti_dt)) {
            continue;
        }

        changed_timestep[i] = 0;

        double search_sq = h[i] * h[i];
        double px = pos_x[i];
        double py = pos_y[i];
        double pz = pos_z[i];

        LBVH::traverse(
            bvh_nodes.data(), pos_x.data(), pos_y.data(), pos_z.data(), px, py,
            pz, search_sq, config.domain_size,
            [&](int j, double dx, double dy, double dz, double r2) {
                // Exclude self and check interacting distance
                if (j != static_cast<int>(i) &&
                    (r2 < search_sq || r2 < h[j] * h[j])) {
                    uint64_t max_tj_allowed =
                        calc_ti_step[i] * TimeConstants::MAX_NGB_TIME_FACTOR;
                    uint64_t max_ti_allowed =
                        calc_ti_step[j] * TimeConstants::MAX_NGB_TIME_FACTOR;

                    // Helper lambda to apply the restriction, align to
                    // power-of-two, and update
                    auto apply_timestep_limit = [&](int target_idx,
                                                    uint64_t current_t,
                                                    uint64_t max_allowed) {
                        uint64_t new_t = max_allowed;
                        int depth = 0;
                        while ((ti_global_time % new_t) != 0 &&
                               depth < TimeConstants::MAX_TIMESTEP_BIN_DEPTH) {
                            new_t >>= 1;
                            depth++;
                        }

                        if (is_active[target_idx]) {
                            if (new_t < current_t) {
#pragma omp atomic write
                                changed_timestep[target_idx] = 1;
                                atomic_update_min_uint64(&ti_step[target_idx],
                                                         new_t);
                                atomic_update_min_uint64(
                                    &ti_end[target_idx],
                                    ti_global_time + new_t);
                            }
                        } else {
#pragma omp atomic write
                            needs_wakeup[target_idx] = 1;
                        }
                    };

                    // Enforce dt_j
                    if (max_tj_allowed < max_system_ti_step) {
                        uint64_t current_tj;
#pragma omp atomic read
                        current_tj = ti_step[j];
                        if (current_tj > max_tj_allowed) {
                            apply_timestep_limit(j, current_tj, max_tj_allowed);
                        }
                    }

                    // Enforce dt_i
                    if (max_ti_allowed < max_system_ti_step) {
                        uint64_t current_ti;
#pragma omp atomic read
                        current_ti = ti_step[i];
                        if (current_ti > max_ti_allowed) {
                            apply_timestep_limit(i, current_ti, max_ti_allowed);
                        }
                    }
                }
            });
    }
}

void GasParticleSystem::reset_diagnostic_counters() {
    diag.active_particles_fraction = 0.0;
    diag.active_particles_num_cycles = 0;
    diag.ill_conditioned = 0;
    diag.clamped_h = 0;
    diag.neighbor_increased = 0;
    diag.forced_wakeup = 0;
}
