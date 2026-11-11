#include "integrator.h"

#include <omp.h>

#include "config.h"
#include "cooling.h"
#include "cosmology.h"
#include "density.h"
#include "diagnostics.h"
#include "gravity.h"
#include "hydro.h"
#include "particles.h"
#include "pm.h"
#include "state.h"
#include "thermodynamics.h"
#include "types.h"

namespace {
void apply_dm_kick(ParticleSystem& dm, double dt, double a, double H,
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
                if (std::abs(dm.global_time - dt_macro - t_start) <
                    TimeConstants::SYNC_TOLERANCE)
                    do_kick = true;
            } else {
                if (std::abs(dm.global_time - dm.t_end[i]) <
                    TimeConstants::SYNC_TOLERANCE)
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

    dm.budget.gravitational_work += step_grav_work;
    dm.budget.expansion_work += step_exp_work;
}

void apply_dm_drift(ParticleSystem& dm, double dt, const Vec3& domain_size) {
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

struct KickResult {
    double ke_spec_new;
    double grav_work;
    double drag_work;
};

inline KickResult apply_velocity_kicks(GasParticleSystem& gas, size_t i,
                                       double dt, double drag_factor) {
    KickResult res;
    const double m = gas.mass[i];

    // Cache old velocities
    double vx_old = gas.vel_x[i];
    double vy_old = gas.vel_y[i];
    double vz_old = gas.vel_z[i];

    double ke_old =
        0.5 * m * (vx_old * vx_old + vy_old * vy_old + vz_old * vz_old);

    // Gravity-only integration (for diagnostics)
    double vx_g = vx_old + gas.acc_x[i] * dt;
    double vy_g = vy_old + gas.acc_y[i] * dt;
    double vz_g = vz_old + gas.acc_z[i] * dt;
    double ke_g = 0.5 * m * (vx_g * vx_g + vy_g * vy_g + vz_g * vz_g);
    res.grav_work = ke_g - ke_old;

    // Total physical integration (Gravity + Hydro)
    double vx_phys = vx_g + gas.hydro_acc_x[i] * dt;
    double vy_phys = vy_g + gas.hydro_acc_y[i] * dt;
    double vz_phys = vz_g + gas.hydro_acc_z[i] * dt;
    double ke_phys =
        0.5 * m * (vx_phys * vx_phys + vy_phys * vy_phys + vz_phys * vz_phys);

    // Cosmological drag applied directly to the global arrays
    gas.vel_x[i] = vx_phys * drag_factor;
    gas.vel_y[i] = vy_phys * drag_factor;
    gas.vel_z[i] = vz_phys * drag_factor;

    double ke_new = 0.5 * m *
                    (gas.vel_x[i] * gas.vel_x[i] + gas.vel_y[i] * gas.vel_y[i] +
                     gas.vel_z[i] * gas.vel_z[i]);

    res.drag_work = ke_new - ke_phys;
    res.ke_spec_new = ke_new / m;

    return res;
}

void apply_gas_drift(GasParticleSystem& gas, double dt,
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

void apply_drift(SimState& state, double dt, const Config& config,
                 Diagnostics& diag) {
    apply_dm_drift(state.dm, dt, config.domain_size);
    apply_gas_drift(state.gas, dt, config.domain_size);
}

void update_timesteps(SimState& state, const Config& config) {
    state.dm.update_particle_timesteps(config.fixed_dt, config,
                                       state.gravity.softening_length);
    state.gas.update_particle_timesteps(config.fixed_dt, state.scale_factor,
                                        config, state.gravity.softening_length);
}
}  // namespace

Integrator::KickWork Integrator::kick_single_gas_particle(
    GasParticleSystem& gas, size_t i, double dt, double a, double H,
    const Config& config) {
    KickWork work;

    // Cosmological Factors
    double drag_factor, expansion_factor;
    Cosmology::compute_step_factors(a, H, dt, config.gamma, drag_factor,
                                    expansion_factor);

    // Kinematics (In-place velocity updates & mechanical work)
    auto kin = apply_velocity_kicks(gas, i, dt, drag_factor);
    work.grav_work = kin.grav_work;

    // Thermodynamic Update (In-place integration & arbitration)
    Thermodynamics::kick_particle_thermodynamics(
        gas, i, kin.ke_spec_new, kin.grav_work, kin.drag_work, dt,
        expansion_factor, a, config, work);

    return work;
}

void Integrator::apply_gas_kick(GasParticleSystem& gas, double dt, double a,
                                double H, const Config& config, bool is_kick1) {
    double step_grav_work = 0.0;
    double step_exp_work = 0.0;
    double step_switch_energy = 0.0;
    double step_injected_energy = 0.0;
    const size_t n = gas.num_particles;

    // dt passed to this function is half of the macro step
    double dt_macro = 2.0 * dt;

#pragma omp parallel for reduction(+ : step_grav_work, step_exp_work, \
                                       step_switch_energy,            \
                                       step_injected_energy) schedule(static)
    for (size_t i = 0; i < n; i++) {
        // Determine if this particle should be kicked based on timeline
        if (config.individual_particle_timesteps) {
            double t_start = gas.t_end[i] - gas.dt_step[i];
            bool do_kick = false;
            if (is_kick1) {
                if (std::abs(gas.global_time - dt_macro - t_start) <
                    TimeConstants::SYNC_TOLERANCE)
                    do_kick = true;
            } else {
                if (std::abs(gas.global_time - gas.t_end[i]) <
                    TimeConstants::SYNC_TOLERANCE)
                    do_kick = true;
            }
            if (!do_kick) continue;
        }

        // Execute the kick
        double half_dt = gas.dt_step[i] / 2.0;
        Integrator::KickWork work =
            kick_single_gas_particle(gas, i, half_dt, a, H, config);

        // Accumulate local diagnostic work
        step_grav_work += work.grav_work;
        step_exp_work += work.exp_work;
        step_switch_energy += work.switch_energy;
        step_injected_energy += work.injected_energy;
    }

    // Update global diagnostics
    gas.budget.gravitational_work += step_grav_work;
    gas.budget.expansion_work += step_exp_work;
    gas.budget.entropy_switch += step_switch_energy;
    gas.budget.floor_injection += step_injected_energy;
}

void Integrator::compute_forces(SimState& state, const Config& config,
                                Diagnostics& diag) {
    // Update softening
    state.gravity.update_softening(state.scale_factor, config);

    // Build trees
    {
        ScopedTimer timer(diag, TimerRegion::Tree);
#ifndef GRAVITY_GPU_MEMORY
        state.dm.build_lbvh(config);
#endif
        state.gas.build_lbvh(config);
    }

    // PM
    {
        ScopedTimer pm_timer(diag, TimerRegion::PM);
        if (config.use_PM) {
            state.pm.bin_and_assign_mass(state, config);
            state.pm.compute_acceleration(state, config);
            state.pm.interpolate_forces(state, config);
        } else {
            state.gravity.reset_accels(state, config);
        }
    }

    // PP
#ifdef GPU_CPU_PARALLEL_EXECUTION
    if (config.use_PP) {
        ScopedTimer pp_timer(diag, TimerRegion::PP);
        pp_worker.launch([&]() {
            ScopedTimer pp_timer(diag, TimerRegion::PP_GPU);
            state.gravity.add_pp_accel(state.dm, state.gas, state.scale_factor,
                                       config, diag);
        });
    }
#endif

    {
        ScopedTimer timer(diag, TimerRegion::Density);
        DensitySolver::compute_density_and_h(state.gas, config,
                                             state.gravity.softening_length);
    }
    {
        ScopedTimer hydro_timer(diag, TimerRegion::Hydro);
        HydroSolver::step(state.gas, config, state.scale_factor);
    }

    if (config.use_PP) {
        ScopedTimer pp_timer(diag, TimerRegion::PP);
#ifdef GPU_CPU_PARALLEL_EXECUTION
        pp_worker.wait();
#else
        state.gravity.add_pp_accel(state.dm, state.gas, state.scale_factor,
                                   config, diag);
#endif
        state.gravity.map_results(state.dm, state.gas, diag);
    }
}

void Integrator::apply_kick(SimState& state, double dt, const Config& config,
                            bool is_kick1) {
    // Forward the config and kick state to the updated DM function
    apply_dm_kick(state.dm, dt, state.scale_factor, state.hubble_param, config,
                  is_kick1);

    apply_gas_kick(state.gas, dt, state.scale_factor, state.hubble_param,
                   config, is_kick1);
}

void Integrator::KDK_step(SimState& state, double dt, const Config& config,
                          Diagnostics& diag) {
    // Clock sync and wake-ups
    state.dm.sync_and_activate(dt, config);
    state.gas.sync_and_activate(*this, dt, state.scale_factor,
                                state.hubble_param, config);

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

    if (config.enable_cooling) {
        ScopedTimer cooling_timer(diag, TimerRegion::Cool);
        state.cooling.apply_cooling(state.gas, config, state.scale_factor);
    }

    // Time bin assignment
    update_timesteps(state, config);
}

void Integrator::initialize(SimState& state, const Config& config) {
    state.gravity.update_softening(state.scale_factor, config);

    state.gas.build_lbvh(config);
    DensitySolver::compute_density_and_h(state.gas, config,
                                         state.gravity.softening_length);
    Thermodynamics::initialize_primitive_variables(state.gas, config,
                                                   state.scale_factor);

    state.dm.sync_and_activate(0.0, config);
    state.gas.sync_and_activate(*this, 0.0, state.scale_factor,
                                state.hubble_param, config);

    Diagnostics dummy_diag;
    compute_forces(state, config, dummy_diag);

    update_timesteps(state, config);
}
