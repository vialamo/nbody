#include "integrator.h"

#include <omp.h>

#include "config.h"
#include "cosmology.h"
#include "density.h"
#include "diagnostics.h"
#include "gravity.h"
#include "particles.h"
#include "pm.h"
#include "state.h"
#include "types.h"

void apply_mesh_gas_gravity_kick(GasGrid& gas, const Grid3D& grav_x,
                                 const Grid3D& grav_y, const Grid3D& grav_z,
                                 double dt, double a, double H,
                                 const Config& config) {
    if (config.hydro_method != HydroMethod::Eulerian) return;

    gas.update_primitive_variables(a);

    Grid3D total_ax_gas(config.mesh_size), total_ay_gas(config.mesh_size),
        total_az_gas(config.mesh_size);
    total_ax_gas.data = grav_x.array() - (2 * H * gas.get_velocity_x().array());
    total_ay_gas.data = grav_y.array() - (2 * H * gas.get_velocity_y().array());
    total_az_gas.data = grav_z.array() - (2 * H * gas.get_velocity_z().array());

    Grid3D g_mom_x_source(config.mesh_size), g_mom_y_source(config.mesh_size),
        g_mom_z_source(config.mesh_size);
    g_mom_x_source.data = gas.get_density().array() * total_ax_gas.array();
    g_mom_y_source.data = gas.get_density().array() * total_ay_gas.array();
    g_mom_z_source.data = gas.get_density().array() * total_az_gas.array();

    // Calculate the change in momentum (dp)
    auto dp_x = g_mom_x_source.array() * dt;
    auto dp_y = g_mom_y_source.array() * dt;
    auto dp_z = g_mom_z_source.array() * dt;

    // Fetch current momentum
    const auto& p_x = gas.get_momentum_x().array();
    const auto& p_y = gas.get_momentum_y().array();
    const auto& p_z = gas.get_momentum_z().array();

    // Kinetic Energy work done by gravity: (P_new^2 - P_old^2) / (2*rho)
    Eigen::ArrayXd gravity_work_density =
        ((p_x + dp_x).square() + (p_y + dp_y).square() + (p_z + dp_z).square() -
         (p_x.square() + p_y.square() + p_z.square())) /
        (2.0 * gas.get_density().array());

    // Apply updates
    gas.momentum_x.array() += dp_x;
    gas.momentum_y.array() += dp_y;
    gas.momentum_z.array() += dp_z;
    gas.energy.array() += gravity_work_density;

    // Adiabatic Expansion Cooling (PdV Work)
    //
    // In physical units, adiabatic expansion causes the thermal energy to decay
    // at a rate of -3(gamma - 1)H.
    //
    // However, we track comoving velocity as v = dx/dt, meaning
    // physical kinetic energy scales as a^2 * v_code^2. To conserve total
    // energy when adding kinetic and internal arrays together, our comoving
    // internal energy must also absorb this scaling: E_phys = a^2 * E_code.
    //
    // Applying the chain rule (d/dt) to this a^2 transformation introduces an
    // extra -2H mathematical decay just to keep the arrays synchronized.
    // The total required numerical decay is -(3*gamma - 1)H.
    double expansion_factor = (3.0 * config.gamma - 1.0) * H * dt;
    Eigen::ArrayXd expansion_cooling =
        expansion_factor * gas.get_internal_energy().array();

    gas.energy.array() -= expansion_cooling;
    gas.internal_energy.array() -= expansion_cooling;

    // Accumulate global work for diagnostics
    // sum() adds up the energy density of all cells. Multiply by volume to get
    // total code energy.
    double step_grav_work = gravity_work_density.sum() * config.cell_volume;
    double step_exp_work = expansion_cooling.sum() * config.cell_volume;

    gas.add_gravitational_work(step_grav_work);
    gas.add_expansion_work(step_exp_work);
}

static void apply_dm_kick(ParticleSystem& dm, double dt, double a, double H,
                          const Config& config, bool is_kick1) {
    double step_grav_work = 0.0;
    double step_exp_work = 0.0;
    const size_t n = dm.num_particles;

    // dt passed to this function is half of the macro step
    double dt_macro = 2.0 * dt;

#pragma omp parallel for reduction(+ : step_grav_work, step_exp_work) \
    schedule(static)
    for (size_t i = 0; i < n; i++) {
        // Gate the kick based on the individual timeline boundaries
        if (config.individual_particle_timesteps) {
            double t_start = dm.t_end[i] - dm.dt_step[i];
            bool do_kick = false;

            if (is_kick1) {
                // Kick 1 only applies if the particle's block step just started
                if (std::abs(dm.global_time - dt_macro - t_start) < 1e-10)
                    do_kick = true;
            } else {
                // Kick 2 only applies if the particle's block step just ended
                if (std::abs(dm.global_time - dm.t_end[i]) < 1e-10)
                    do_kick = true;
            }
            if (!do_kick) continue;
        } else {
            if (!dm.is_active[i]) continue;
        }

        // Use the exact half-step of the individual particle
        double half_dt = dm.dt_step[i] / 2.0;

        // Calculate cosmological damping factor tailored to this specific step
        double drag_factor = 1.0;
        if (H > 0.0 && a > 0.0) {
            double a_next = a + a * H * half_dt;
            drag_factor =
                (a / a_next) * (a / a_next);  // (1 / (1 + H*half_dt))^2
        }

        double m = dm.mass[i];
        double vx_old = dm.vel_x[i];
        double vy_old = dm.vel_y[i];
        double vz_old = dm.vel_z[i];

        // Comoving Gravitational Acceleration
        double gx = dm.acc_x[i];
        double gy = dm.acc_y[i];
        double gz = dm.acc_z[i];

        // Intermediate velocity after pure gravity kick
        double vx_g = vx_old + gx * half_dt;
        double vy_g = vy_old + gy * half_dt;
        double vz_g = vz_old + gz * half_dt;

        // Gravitational Work (dW_grav = dK_grav)
        double ke_old =
            0.5 * m * (vx_old * vx_old + vy_old * vy_old + vz_old * vz_old);
        double ke_g = 0.5 * m * (vx_g * vx_g + vy_g * vy_g + vz_g * vz_g);
        step_grav_work += (ke_g - ke_old);

        // Final velocity after Hubble Drag decay
        double vx_new = vx_g * drag_factor;
        double vy_new = vy_g * drag_factor;
        double vz_new = vz_g * drag_factor;

        // Expansion Loss Work (dW_exp = dK_drag)
        double ke_new =
            0.5 * m * (vx_new * vx_new + vy_new * vy_new + vz_new * vz_new);
        step_exp_work += (ke_g - ke_new);

        dm.vel_x[i] = vx_new;
        dm.vel_y[i] = vy_new;
        dm.vel_z[i] = vz_new;
    }

    dm.accumulated_gravitational_work += step_grav_work;
    dm.accumulated_expansion_work += step_exp_work;
}

static void apply_dm_drift(ParticleSystem& dm, double dt,
                           const Vec3& domain_size) {
    double inv_domain_x = 1.0 / domain_size.x;
    double inv_domain_y = 1.0 / domain_size.y;
    double inv_domain_z = 1.0 / domain_size.z;
    const size_t n = dm.num_particles;

#pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < n; i++) {
        double drift_dt = dm.global_time - dm.t_current[i];
        if (drift_dt <= 0.0) continue;

        double nx = dm.pos_x[i] + dm.vel_x[i] * drift_dt;
        double ny = dm.pos_y[i] + dm.vel_y[i] * drift_dt;
        double nz = dm.pos_z[i] + dm.vel_z[i] * drift_dt;

        // SIMD-friendly periodic boundary using independent axes
        dm.pos_x[i] = nx - domain_size.x * std::floor(nx * inv_domain_x);
        dm.pos_y[i] = ny - domain_size.y * std::floor(ny * inv_domain_y);
        dm.pos_z[i] = nz - domain_size.z * std::floor(nz * inv_domain_z);

        dm.t_current[i] = dm.global_time;
    }
}

KickWork kick_single_gas_particle_gravity(GasParticleSystem& gas, size_t i,
                                          double dt, double a, double H,
                                          const Config& config) {
    KickWork work;

    // Cosmological factors parameterized by the local dt
    double drag_factor = 1.0;
    if (H > 0.0 && a > 0.0) {
        double a_next = a + a * H * dt;
        drag_factor = (a / a_next) * (a / a_next);
    }
    double expansion_factor = (3.0 * config.gamma - 1.0) * H * dt;

    double m = gas.mass[i];
    double vx_old = gas.vel_x[i];
    double vy_old = gas.vel_y[i];
    double vz_old = gas.vel_z[i];

    // Intermediate velocity after gravity kick
    double vx_g = vx_old + gas.acc_x[i] * dt;
    double vy_g = vy_old + gas.acc_y[i] * dt;
    double vz_g = vz_old + gas.acc_z[i] * dt;

    // Gravitational Work (dW_grav = dK_grav)
    double ke_old =
        0.5 * m * (vx_old * vx_old + vy_old * vy_old + vz_old * vz_old);
    double ke_g = 0.5 * m * (vx_g * vx_g + vy_g * vy_g + vz_g * vz_g);
    work.grav_work = (ke_g - ke_old);

    // Final velocity after Hubble Drag decay
    double vx_new = vx_g * drag_factor;
    double vy_new = vy_g * drag_factor;
    double vz_new = vz_g * drag_factor;

    // Apply Cosmological Adiabatic Cooling (Thermal energy lost to PdV work)
    double u_old = gas.u[i];
    double u_cooling = u_old * expansion_factor;

    // Expansion Work Accumulation
    double ke_new =
        0.5 * m * (vx_new * vx_new + vy_new * vy_new + vz_new * vz_new);
    work.exp_work = (ke_g - ke_new) + (u_cooling * m);

    // Update particle arrays
    gas.vel_x[i] = vx_new;
    gas.vel_y[i] = vy_new;
    gas.vel_z[i] = vz_new;
    gas.u[i] -= u_cooling;
    gas.entropy[i] -= gas.entropy[i] * expansion_factor;

    // Maintain Dual Energy synchronization
    double spec_ke_old = ke_old / m;
    double spec_ke_new = ke_new / m;
    gas.total_energy[i] += (spec_ke_new - spec_ke_old) - u_cooling;

    return work;
}

static void apply_gas_particle_gravity_kick(GasParticleSystem& gas, double dt,
                                            double a, double H,
                                            const Config& config,
                                            bool is_kick1) {
    double step_grav_work = 0.0;
    double step_exp_work = 0.0;
    const size_t n = gas.num_particles;
    double dt_macro = 2.0 * dt;

#pragma omp parallel for reduction(+ : step_grav_work, step_exp_work) \
    schedule(static)
    for (size_t i = 0; i < n; i++) {
        if (config.individual_particle_timesteps) {
            double t_start = gas.t_end[i] - gas.dt_step[i];
            bool do_kick = false;
            if (is_kick1) {
                if (std::abs(gas.global_time - dt_macro - t_start) < 1e-10)
                    do_kick = true;
            } else {
                if (std::abs(gas.global_time - gas.t_end[i]) < 1e-10)
                    do_kick = true;
            }
            if (!do_kick) continue;
        }

        double half_dt = gas.dt_step[i] / 2.0;
        KickWork work =
            kick_single_gas_particle_gravity(gas, i, half_dt, a, H, config);

        step_grav_work += work.grav_work;
        step_exp_work += work.exp_work;
    }

    gas.accumulated_gravitational_work += step_grav_work;
    gas.accumulated_expansion_work += step_exp_work;
}

void kick_single_gas_particle_hydro(GasParticleSystem& gas, size_t i, double dt,
                                    const Config& config) {
    double gamma_minus_1 = config.gamma - 1.0;
    constexpr double min_energy = 1e-20;

    double m = gas.mass[i];
    double vx_old = gas.vel_x[i];
    double vy_old = gas.vel_y[i];
    double vz_old = gas.vel_z[i];
    double ke_old =
        0.5 * m * (vx_old * vx_old + vy_old * vy_old + vz_old * vz_old);

    gas.vel_x[i] += gas.hydro_acc_x[i] * dt;
    gas.vel_y[i] += gas.hydro_acc_y[i] * dt;
    gas.vel_z[i] += gas.hydro_acc_z[i] * dt;

    double ke_new = 0.5 * m *
                    (gas.vel_x[i] * gas.vel_x[i] + gas.vel_y[i] * gas.vel_y[i] +
                     gas.vel_z[i] * gas.vel_z[i]);
    double delta_ke = ke_new - ke_old;

    double delta_u = gas.du_dt[i] * dt;
    gas.u[i] += delta_u;

    gas.total_energy[i] += gas.de_dt[i] * dt;

    if (gas.u[i] < min_energy) gas.u[i] = min_energy;
    if (gas.total_energy[i] < min_energy) gas.total_energy[i] = min_energy;

    gas.entropy[i] =
        gamma_minus_1 * gas.u[i] / std::pow(gas.rho[i], gamma_minus_1);
    gas.pressure[i] = gamma_minus_1 * gas.rho[i] * gas.u[i];
}

static void apply_gas_particle_hydro_kick(GasParticleSystem& gas, double dt,
                                          double a, const Config& config,
                                          bool is_kick1) {
    const size_t n = gas.num_particles;
    double dt_macro = 2.0 * dt;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) {
        if (config.individual_particle_timesteps) {
            double t_start = gas.t_end[i] - gas.dt_step[i];
            bool do_kick = false;
            if (is_kick1) {
                if (std::abs(gas.global_time - dt_macro - t_start) < 1e-10)
                    do_kick = true;
            } else {
                if (std::abs(gas.global_time - gas.t_end[i]) < 1e-10)
                    do_kick = true;
            }
            if (!do_kick) continue;
        }

        double half_dt = gas.dt_step[i] / 2.0;
        kick_single_gas_particle_hydro(gas, i, half_dt, config);
    }
}

static void apply_gas_particle_drift(GasParticleSystem& gas, double dt,
                                     const Vec3& domain_size) {
    double inv_domain_x = 1.0 / domain_size.x;
    double inv_domain_y = 1.0 / domain_size.y;
    double inv_domain_z = 1.0 / domain_size.z;
    const size_t n = gas.num_particles;

    // Drift from t_current to global_time
#pragma omp parallel for simd schedule(static)
    for (size_t i = 0; i < n; i++) {
        double drift_dt = gas.global_time - gas.t_current[i];
        if (drift_dt <= 0.0) continue;  // Skip if already drifted

        double nx = gas.pos_x[i] + gas.vel_x[i] * drift_dt;
        double ny = gas.pos_y[i] + gas.vel_y[i] * drift_dt;
        double nz = gas.pos_z[i] + gas.vel_z[i] * drift_dt;

        // Periodic boundaries
        gas.pos_x[i] = nx - domain_size.x * std::floor(nx * inv_domain_x);
        gas.pos_y[i] = ny - domain_size.y * std::floor(ny * inv_domain_y);
        gas.pos_z[i] = nz - domain_size.z * std::floor(nz * inv_domain_z);

        gas.t_current[i] = gas.global_time;  // Mark as drifted
    }
}

static void calculate_max_acceleration(ParticleSystem& dm) {
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

    dm.max_accel_sq = local_max;
}

static void calculate_max_acceleration(GasParticleSystem& gas) {
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
    gas.max_accel_sq = local_max;
}

static void reset_all_accels(SimState& state, Config& config) {
    state.pm_gravity_x.setZero();
    state.pm_gravity_y.setZero();
    state.pm_gravity_z.setZero();

    std::fill(state.dm.acc_x.begin(), state.dm.acc_x.end(), 0.0);
    std::fill(state.dm.acc_y.begin(), state.dm.acc_y.end(), 0.0);
    std::fill(state.dm.acc_z.begin(), state.dm.acc_z.end(), 0.0);
    if (config.hydro_method == HydroMethod::MFM) {
        GasParticleSystem& gas = *state.mfm_gas;
        std::fill(gas.acc_x.begin(), gas.acc_x.end(), 0.0);
        std::fill(gas.acc_y.begin(), gas.acc_y.end(), 0.0);
        std::fill(gas.acc_z.begin(), gas.acc_z.end(), 0.0);
    }
}

void compute_forces(SimState& state, Config& config, Diagnostics& diag) {
    // Update softening
    GravitySolver::update_softening(state, config);

    // Build trees
    {
        ScopedTimer timer(diag, TimerRegion::Tree);
        state.dm.build_lbvh(config);
        if (config.hydro_method == HydroMethod::MFM) {
            state.mfm_gas->build_lbvh(config);
        }
    }

    // Hydro
    if (config.hydro_method == HydroMethod::MFM) {
        {
            ScopedTimer timer(diag, TimerRegion::Density);
            DensitySolver::compute_density_and_h(*state.mfm_gas, config,
                                                 state.dm);
        }
        {
            ScopedTimer hydro_timer(diag, TimerRegion::Hydro);
            state.mfm_gas->hydro_step(config, state.scale_factor);
        }
        if (config.enable_cooling) {
            ScopedTimer cooling_timer(diag, TimerRegion::Cool);
            state.mfm_gas->apply_cooling(state.scale_factor, config,
                                         state.cooling);
        }
    }

    // PM
    {
        ScopedTimer pm_timer(diag, TimerRegion::PM);
        PM::bin_and_assign_mass(state, config);

        if (config.use_PM) {
            PM::compute_acceleration(state, config);
            PM::interpolate_forces(state, config);
        } else {
            reset_all_accels(state, config);
        }
    }

    // PP
    {
        ScopedTimer pp_timer(diag, TimerRegion::PP);
        if (config.use_PP) {
            // DM-DM
            GravitySolver::compute_and_add_dm_pp_forces(
                state.dm, state.scale_factor, config, diag);

            if (config.hydro_method == HydroMethod::MFM) {
                // Gas-Gas
                GravitySolver::compute_and_add_gas_pp_forces(
                    *state.mfm_gas, state.scale_factor, config, diag);

                // Gas-DM
                GravitySolver::compute_gas_dm_cross_pp_forces(
                    state.dm, *state.mfm_gas, state.scale_factor, config, diag);
            }
        }

        calculate_max_acceleration(state.dm);
        if (config.hydro_method == HydroMethod::MFM) {
            calculate_max_acceleration(*state.mfm_gas);
        }
    }
}

static void eulerian_mesh_flux(SimState& state, Diagnostics& diag,
                               Config& config, double dt) {
    // Approximate the scale factor at the half-step (t + dt/2)
    double mid_a = state.scale_factor * (1.0 + 0.5 * state.hubble_param * dt);
    {
        ScopedTimer hydro_timer(diag, TimerRegion::Hydro);
        state.gas->hydro_step(dt, mid_a);
    }
    if (config.enable_cooling) {
        ScopedTimer cooling_timer(diag, TimerRegion::Cool);
        state.gas->apply_cooling(dt, mid_a, state.cooling);
    }
}

static void apply_kick(SimState& state, double dt, const Config& config,
                       bool is_kick1) {
    // Forward the config and kick state to the updated DM function
    apply_dm_kick(state.dm, dt, state.scale_factor, state.hubble_param, config,
                  is_kick1);

    if (config.hydro_method == HydroMethod::Eulerian) {
        apply_mesh_gas_gravity_kick(*state.gas, state.pm_gravity_x,
                                    state.pm_gravity_y, state.pm_gravity_z, dt,
                                    state.scale_factor, state.hubble_param,
                                    config);
    } else if (config.hydro_method == HydroMethod::MFM) {
        apply_gas_particle_gravity_kick(*state.mfm_gas, dt, state.scale_factor,
                                        state.hubble_param, config, is_kick1);
    }

    if (config.hydro_method == HydroMethod::MFM) {
        // Note: apply_gas_particle_hydro_kick ignores the dt passed
        // here and uses gas.dt_step[i] / 2.0 internally
        apply_gas_particle_hydro_kick(*state.mfm_gas, dt, state.scale_factor,
                                      config, is_kick1);
    }
}

static void apply_drift(SimState& state, double dt, Config& config,
                        Diagnostics& diag) {
    apply_dm_drift(state.dm, dt, config.domain_size);

    if (config.hydro_method == HydroMethod::Eulerian) {
        eulerian_mesh_flux(state, diag, config, dt);
    } else if (config.hydro_method == HydroMethod::MFM) {
        apply_gas_particle_drift(*state.mfm_gas, dt, config.domain_size);
    }
}

void KDK_step(SimState& state, TimestepInfo& ts, Config& config,
              Diagnostics& diag) {
    double dt = ts.dt_macro;

    // Clock sync and wake-ups
    state.dm.sync_and_activate(dt, config);
    if (config.hydro_method == HydroMethod::MFM) {
        state.mfm_gas->sync_and_activate(dt, state.scale_factor,
                                         state.hubble_param, config);
    }

    // KICK 1 (Half step)
    apply_kick(state, dt / 2.0, config, true);

    // DRIFT
    apply_drift(state, dt, config, diag);

    // Update cosmology and compute forces
    state.total_time += dt;
    update_cosmology(state, config);
    compute_forces(state, config, diag);

    // KICK 2 (Half step)
    apply_kick(state, dt / 2.0, config, false);

    // SYNCHRONIZATION BEFORE I/O AND DIAGNOSTICS
    // Force the primitive variables (u, pressure, entropy) to conform
    // to the kicked conservative variables (total_energy and velocity)
    if (config.hydro_method == HydroMethod::MFM) {
        state.mfm_gas->update_primitive_variables(config, state.scale_factor);
    } else if (config.hydro_method == HydroMethod::Eulerian) {
        state.gas->update_primitive_variables(state.scale_factor);
    }

    // Time bin assignment
    state.dm.update_particle_timesteps(config.fixed_dt, config);
    if (config.hydro_method == HydroMethod::MFM) {
        state.mfm_gas->update_particle_timesteps(
            config.fixed_dt, state.scale_factor, config, state.cooling);
    }
}
