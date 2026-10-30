#include "mfm.h"

#include "config.h"
#include "density.h"
#include "hydro.h"
#include "integrator.h"
#include "math_utils.h"
#include "thermodynamics.h"

GasParticleSystem::GasParticleSystem(const Config& config)
    : cic_data(config.num_gas_particles) {
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
    delta_E_grav.reserve(config.num_gas_particles);

    vel_x_pred.reserve(config.num_gas_particles);
    vel_y_pred.reserve(config.num_gas_particles);
    vel_z_pred.reserve(config.num_gas_particles);
    pressure_pred.reserve(config.num_gas_particles);

    v_sig_max.reserve(config.num_gas_particles);
    cond_num.reserve(config.num_gas_particles);
    n_enc_final.reserve(config.num_gas_particles);

    time_bin.reserve(config.num_gas_particles);
    dt_step.reserve(config.num_gas_particles);
    t_current.reserve(config.num_gas_particles);
    t_end.reserve(config.num_gas_particles);
    is_active.reserve(config.num_gas_particles);
    needs_wakeup.reserve(config.num_gas_particles);
    active_indices.reserve(config.num_gas_particles);

    size_t num_nodes = 2 * config.num_gas_particles - 1;
    if (config.num_dm_particles > 0) {
        morton_codes.reserve(config.num_gas_particles);
        sorted_indices.reserve(config.num_gas_particles);
        bvh_nodes.reserve(num_nodes);
    }

    double u_code_K =
        Thermodynamics::get_internal_energy_from_temp(0.1, 1.0, config);
    pressure_floor = (config.gamma - 1.0) * density_floor * u_code_K;
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
    delta_E_grav.push_back(0.0);

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

    time_bin.push_back(0);
    dt_step.push_back(0.0);
    t_current.push_back(0.0);
    t_end.push_back(0.0);
    is_active.push_back(1);  // All particles start active at t=0
    needs_wakeup.push_back(0);

    active_indices.push_back(num_particles);
    num_active++;

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
    std::vector<double> new_delta_E_grav(num_particles);

    std::vector<double> new_cond_num(num_particles);
    std::vector<double> new_n_enc_final(num_particles);

    std::vector<int> new_time_bin(num_particles);
    std::vector<double> new_dt_step(num_particles),
        new_t_current(num_particles), new_t_end(num_particles);
    std::vector<uint8_t> new_is_active(num_particles);
    std::vector<uint8_t> new_needs_wakeup(num_particles);
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
        new_delta_E_grav[i] = delta_E_grav[src];
        new_cond_num[i] = cond_num[src];
        new_n_enc_final[i] = n_enc_final[src];
        new_time_bin[i] = time_bin[src];
        new_dt_step[i] = dt_step[src];
        new_t_current[i] = t_current[src];
        new_t_end[i] = t_end[src];
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
    delta_E_grav = std::move(new_delta_E_grav);
    cond_num = std::move(new_cond_num);
    n_enc_final = std::move(new_n_enc_final);
    time_bin = std::move(new_time_bin);
    dt_step = std::move(new_dt_step);
    t_current = std::move(new_t_current);
    t_end = std::move(new_t_end);
    is_active = std::move(new_is_active);
    needs_wakeup = std::move(new_needs_wakeup);
    v_sig_max = std::move(new_v_sig_max);

    active_indices.clear();
    for (size_t i = 0; i < num_particles; ++i) {
        if (is_active[i]) {
            active_indices.push_back(i);
        }
    }
    num_active = active_indices.size();
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

void GasParticleSystem::sync_and_activate(double dt, double a, double H,
                                          const Config& config) {
    if (!config.individual_particle_timesteps) {
        global_time += dt;
        for (size_t i = 0; i < num_particles; ++i) {
            is_active[i] = 1;
            dt_step[i] = dt;
            t_end[i] = global_time;  // Keep timeline synced
            needs_wakeup[i] = 0;
        }
        return;
    }

    double new_global_time = global_time + dt;

    // Process wakeups from the previous cycle
    for (size_t i = 0; i < num_particles; ++i) {
        if (needs_wakeup[i]) {
            double dt_old = dt_step[i];
            double dt_new = dt;

            // The particle previously received Kick 1 for dt_old / 2.
            // This projected its velocity to t_kicked_to
            double t_start_old = t_end[i] - dt_old;
            double t_kicked_to = t_start_old + (dt_old / 2.0);

            // We must rewind the velocity back to the CURRENT
            // global_time. Then the engine will natively apply Kick 1 for
            // dt_new / 2
            double dt_excess = t_kicked_to - global_time;

            // Rewind kicks using negative time
            // Note that we should better compute the a and H which
            // originally kicked the particle
            Integrator::KickWork w_grav = Integrator::kick_single_gas_particle(
                *this, i, -dt_excess, a, H, config);

            // Adjust global accumulators (adding negative work subtracts the
            // phantom work)
            accumulated_gravitational_work += w_grav.grav_work;
            accumulated_expansion_work += w_grav.exp_work;

            // GIZMO Wake-up logic:
            // Snap the particle to the current active micro-step size.
            // It will drift the full distance since it last moved (using
            // t_current), but its kicks and cooling will use this tiny
            // micro-step
            dt_step[i] = dt_new;
            t_end[i] = new_global_time;  // Force sync to the current clock
            needs_wakeup[i] = 0;
        }
    }

    // Advance global clock
    global_time = new_global_time;

    // Flag active particles
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        // A particle is active if its block step ends at the new global time
        is_active[i] = (std::abs(t_end[i] - global_time) < 1e-10) ? 1 : 0;

        // BOOTSTRAP: If it's the first step (t=0), everyone is active
        if (global_time <= dt + 1e-12 && t_end[i] == 0.0) {
            is_active[i] = 1;
            dt_step[i] = dt;  // Give them a tiny initial dt so Kick 2 works
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

    active_particles_fraction +=
        static_cast<double>(num_active) / static_cast<double>(num_particles);
    active_particles_num_cycles++;
}

double GasParticleSystem::get_cfl_timestep(double a,
                                           const Config& config) const {
    if (config.hydro_method != HydroMethod::MFM || num_particles == 0 ||
        config.individual_particle_timesteps) {
        return std::numeric_limits<double>::infinity();
    }

    double min_dt = std::numeric_limits<double>::infinity();

#pragma omp parallel for reduction(min : min_dt) schedule(static)
    for (size_t i = 0; i < num_particles; ++i) {
        // v_sig_max is cached in physical units during compute_hydro_forces
        if (v_sig_max[i] > 1e-9) {
            double dt_i = (a * h[i]) / v_sig_max[i];
            if (dt_i < min_dt) {
                min_dt = dt_i;
            }
        }
    }

    // BOOTSTRAP: On the very first step (t = 0), forces haven't been computed
    // yet, so v_sig_max is 0. Provide a tiny, safe initial timestep to get the
    // engine started
    if (min_dt == std::numeric_limits<double>::infinity()) {
        return config.fixed_dt / std::pow(2.0, 10.0);
    }

    return min_dt * config.hydro_courant_factor;
}

void GasParticleSystem::update_particle_timesteps(double dt_max, double a,
                                                  const Config& config,
                                                  Cooling& cooling) {
    if (num_particles == 0) return;

    if (!config.individual_particle_timesteps) {
        // Under global timesteps, the engine already calculated dt_max as the
        // global minimum
        for (size_t i = 0; i < num_particles; ++i) {
            dt_step[i] = dt_max;
            time_bin[i] = 0;
        }
        return;
    }

    double u_floor =
        config.enable_cooling ? Thermodynamics::get_u_floor(a, config) : 0.0;
    double epsilon = config.softening;

    if (num_active == 0) return;

    double current_min_dt = dt_max;
    for (size_t j = 0; j < num_particles; ++j) {
        if (dt_step[j] > 0.0) {
            current_min_dt = std::min(current_min_dt, dt_step[j]);
        }
    }
    double cooling_clamp_floor = std::min(current_min_dt, dt_max / 16.0);

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

        // CFL Timestep (Using cached signal velocity from compute_hydro_forces)
        double v_sig_max_phys = v_sig_max[i];
        if (v_sig_max_phys > 1e-9) {
            double dt_cfl =
                ((a * h[i]) / v_sig_max_phys) * config.hydro_courant_factor;
            dt_ideal = std::min(dt_ideal, dt_cfl);
        }

        // Power-of-two binning
        int n = 0;
        double dt_bin = dt_max;

        // Halve the timestep until it bounds the ideal timestep
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

void GasParticleSystem::hydro_step(const Config& config, double a) {
    Thermodynamics::predict_primitive_variables(*this, config, a);
    HydroSolver::compute_gradients(*this, config);
    HydroSolver::compute_forces(*this, config, a);
}

void GasParticleSystem::sync_thermo(const Config& config, double a,
                                    bool first_call) {
    if (first_call) {
        for (size_t i = 0; i < num_particles; ++i) {
            double ke = 0.5 * (vel_x[i] * vel_x[i] + vel_y[i] * vel_y[i] +
                               vel_z[i] * vel_z[i]);
            total_energy[i] = u[i] + ke;
        }
    }
    Thermodynamics::update_primitive_variables(*this, config, a);
}
