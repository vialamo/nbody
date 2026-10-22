#include "integrator.h"

#include <omp.h>

#include "config.h"
#include "cooling.h"
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

namespace Integrator {

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
                if (std::abs(dm.global_time - dt_macro - t_start) < 1e-10)
                    do_kick = true;
            } else {
                if (std::abs(dm.global_time - dm.t_end[i]) < 1e-10)
                    do_kick = true;
            }
            if (!do_kick) continue;
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

KickWork kick_single_gas_particle(GasParticleSystem& gas, size_t i, double dt,
                                  double a, double H, const Config& config) {
    KickWork work;
    double gamma_minus_1 = config.gamma - 1.0;
    constexpr double min_energy = 1e-20;

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
    double ke_old =
        0.5 * m * (vx_old * vx_old + vy_old * vy_old + vz_old * vz_old);

    // Physical Acceleration (Gravity + Hydro)
    // Isolate the gravity-only kick first to compute
    // the isolated gravitational work for the diagnostics
    double vx_g = vx_old + gas.acc_x[i] * dt;
    double vy_g = vy_old + gas.acc_y[i] * dt;
    double vz_g = vz_old + gas.acc_z[i] * dt;
    double ke_g = 0.5 * m * (vx_g * vx_g + vy_g * vy_g + vz_g * vz_g);

    // Gravitational Work (dW_grav = dK_grav)
    work.grav_work = (ke_g - ke_old);

    // Add Hydro acceleration to get the total intermediate physical
    // velocity
    double vx_phys = vx_g + gas.hydro_acc_x[i] * dt;
    double vy_phys = vy_g + gas.hydro_acc_y[i] * dt;
    double vz_phys = vz_g + gas.hydro_acc_z[i] * dt;
    double ke_phys =
        0.5 * m * (vx_phys * vx_phys + vy_phys * vy_phys + vz_phys * vz_phys);

    // Cosmological Drag
    // Apply Hubble drag to the combined physical velocity
    gas.vel_x[i] = vx_phys * drag_factor;
    gas.vel_y[i] = vy_phys * drag_factor;
    gas.vel_z[i] = vz_phys * drag_factor;

    double ke_new = 0.5 * m *
                    (gas.vel_x[i] * gas.vel_x[i] + gas.vel_y[i] * gas.vel_y[i] +
                     gas.vel_z[i] * gas.vel_z[i]);
    double w_drag =
        ke_new -
        ke_phys;  // Negative value representing kinetic energy lost to drag

    // Thermodynamics (Hydro du/dt + Adiabatic Expansion Cooling)
    double u_phys = gas.u[i] + gas.du_dt[i] * dt;
    double u_cooling =
        u_phys *
        expansion_factor;  // Apply expansion to the updated physical energy

    gas.u[i] = u_phys - u_cooling;
    if (gas.u[i] < min_energy) gas.u[i] = min_energy;

    // Expansion Work Accumulation (Kinetic drag loss + Thermal expansion loss)
    work.exp_work = -w_drag + (u_cooling * m);

    // Total Energy Synchronization
    // total_energy changes due to:
    //   -de_dt (Riemann hydro flux, which intrinsically includes hydro
    //    kinetic work)
    //   -Gravity kinetic work
    //   -Cosmological losses (drag + cooling)
    gas.total_energy[i] +=
        (gas.de_dt[i] * dt) + (work.grav_work / m) + (w_drag / m) - u_cooling;

    if (gas.total_energy[i] < min_energy) gas.total_energy[i] = min_energy;

    // Update EOS Variables
    gas.entropy[i] =
        gamma_minus_1 * gas.u[i] / std::pow(gas.rho[i], gamma_minus_1);
    gas.pressure[i] = gamma_minus_1 * gas.rho[i] * gas.u[i];

    return work;
}

static void apply_gas_kick(GasParticleSystem& gas, double dt, double a,
                           double H, const Config& config, bool is_kick1) {
    double step_grav_work = 0.0;
    double step_exp_work = 0.0;
    const size_t n = gas.num_particles;

    // dt passed to this function is half of the macro step
    double dt_macro = 2.0 * dt;

#pragma omp parallel for reduction(+ : step_grav_work, step_exp_work) \
    schedule(static)
    for (size_t i = 0; i < n; i++) {
        // Determine if this particle should be kicked based on timeline
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

        // Execute the kick
        double half_dt = gas.dt_step[i] / 2.0;
        KickWork work = kick_single_gas_particle(gas, i, half_dt, a, H, config);

        // Accumulate local diagnostic work
        step_grav_work += work.grav_work;
        step_exp_work += work.exp_work;
    }

    // Update global diagnostics
    gas.accumulated_gravitational_work += step_grav_work;
    gas.accumulated_expansion_work += step_exp_work;
}

static void apply_gas_drift(GasParticleSystem& gas, double dt,
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

static void reset_all_accels(SimState& state, Config& config) {
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
#ifndef GRAVITY_OWN_MEMORY
        state.dm.build_lbvh(config);
#endif
        if (config.hydro_method == HydroMethod::MFM) {
            state.mfm_gas->build_lbvh(config);
        }
    }

    // Hydro
    if (config.hydro_method == HydroMethod::MFM) {
        {
            ScopedTimer timer(diag, TimerRegion::Density);
            DensitySolver::compute_density_and_h(*state.mfm_gas, config);
        }
        {
            ScopedTimer hydro_timer(diag, TimerRegion::Hydro);
            state.mfm_gas->hydro_step(config, state.scale_factor);
        }
        if (config.enable_cooling) {
            ScopedTimer cooling_timer(diag, TimerRegion::Cool);
            state.cooling.apply_cooling(*state.mfm_gas, config,
                                        state.scale_factor);
        }
    }

    // PM
    {
        ScopedTimer pm_timer(diag, TimerRegion::PM);
        if (config.use_PM) {
            state.pm.bin_and_assign_mass(state, config);
            state.pm.compute_acceleration(state, config);
            state.pm.interpolate_forces(state, config);
        } else {
            reset_all_accels(state, config);
        }
    }

    // PP
    if (config.use_PP) {
        ScopedTimer pp_timer(diag, TimerRegion::PP);
        GravitySolver::add_pp_accel(state.dm, state.mfm_gas.get(),
                                    state.scale_factor, config, diag);
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
        apply_mesh_gas_gravity_kick(*state.gas, state.pm.acc_x, state.pm.acc_y,
                                    state.pm.acc_z, dt, state.scale_factor,
                                    state.hubble_param, config);
    } else if (config.hydro_method == HydroMethod::MFM) {
        apply_gas_kick(*state.mfm_gas, dt, state.scale_factor,
                       state.hubble_param, config, is_kick1);
    }
}

static void apply_drift(SimState& state, double dt, Config& config,
                        Diagnostics& diag) {
    apply_dm_drift(state.dm, dt, config.domain_size);

    if (config.hydro_method == HydroMethod::Eulerian) {
        eulerian_mesh_flux(state, diag, config, dt);
    } else if (config.hydro_method == HydroMethod::MFM) {
        apply_gas_drift(*state.mfm_gas, dt, config.domain_size);
    }
}

void update_timesteps(SimState& state, const Config& config) {
    state.dm.update_particle_timesteps(config.fixed_dt, config);
    if (config.hydro_method == HydroMethod::MFM) {
        state.mfm_gas->update_particle_timesteps(
            config.fixed_dt, state.scale_factor, config, state.cooling);
    }
}

void KDK_step(SimState& state, double dt, Config& config, Diagnostics& diag) {
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
    Cosmology::update(state, config);
    compute_forces(state, config, diag);

    // KICK 2 (Half step)
    apply_kick(state, dt / 2.0, config, false);

    // Time bin assignment
    update_timesteps(state, config);
}

}  // namespace Integrator
