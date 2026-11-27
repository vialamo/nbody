#include "thermodynamics.h"

#include "gas.h"
#include "integrator.h"

namespace {
inline double get_jeans_pressure_floor(double rho, double h, double a,
                                       const Config& config) {
    const double N_J = config.jeans_resolution_factor;
    const double G_eff = config.G / (a * a * a);
    return (N_J * N_J) / (M_PI * config.gamma) * G_eff * (rho * rho) * (h * h);
}
}  // namespace

namespace Thermodynamics {
void initialize_primitive_variables(GasParticleSystem& gas,
                                    const Config& config, double a) {
    size_t num_active = gas.get_num_active();
    if (num_active == 0) return;

    const double gamma_minus_1 = config.gamma - 1.0;

    // Calculate the physical internal energy floor
    const double u_floor =
        get_internal_energy_from_temp(config.temp_floor_k, a, config);

    double step_floor_heating = 0.0;

#pragma omp parallel for simd reduction(+ : step_floor_heating) schedule(static)
    for (size_t k = 0; k < num_active; ++k) {
        size_t i = gas.active_indices[k];

        // Enforce density and metallicity floors/bounds
        gas.rho[i] = std::max(gas.rho[i], config.density_floor);
        gas.metal_frac[i] = std::max(0.0, std::min(gas.metal_frac[i], 1.0));

        // Enforce thermodynamic energy floor
        double u_temp = gas.u[i];
        if (u_temp < u_floor) {
            double injected_u = u_floor - u_temp;
            step_floor_heating += injected_u * gas.mass[i];

            u_temp = u_floor;
            gas.u[i] = u_temp;
            gas.total_energy[i] += injected_u;
        }

        // Always sync entropy and pressure to the current valid state
        gas.entropy[i] =
            gamma_minus_1 * u_temp / std::pow(gas.rho[i], gamma_minus_1);
        gas.pressure[i] = gamma_minus_1 * gas.rho[i] * u_temp;
    }

    gas.budget.photoheated += step_floor_heating;
}

void predict_primitive_variables(GasParticleSystem& gas, const Config& config,
                                 double a) {
    if (gas.num_particles == 0) return;

    const double gamma_minus_1 = config.gamma - 1.0;
    const double u_floor =
        get_internal_energy_from_temp(config.temp_floor_k, a, config);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < gas.num_particles; ++i) {
        // Project from the Kick 1 half-step to the current global time
        uint64_t ti_kicked_to = gas.ti_end[i] - (gas.ti_step[i] / 2);

        // Use signed integer math, as the particle's Kick 1 target could
        // theoretically be ahead of or behind the current global time if it is
        // a sleeping particle
        int64_t ti_pred = static_cast<int64_t>(gas.ti_global_time) -
                          static_cast<int64_t>(ti_kicked_to);

        double dt_pred = ti_pred * config.timebase_to_float;

        gas.vel_x_pred[i] =
            gas.vel_x[i] + (gas.acc_x[i] + gas.hydro_acc_x[i]) * dt_pred;
        gas.vel_y_pred[i] =
            gas.vel_y[i] + (gas.acc_y[i] + gas.hydro_acc_y[i]) * dt_pred;
        gas.vel_z_pred[i] =
            gas.vel_z[i] + (gas.acc_z[i] + gas.hydro_acc_z[i]) * dt_pred;

        // Predict the new internal energy adiabatically using the new density.
        // This captures the mechanical PdV work of the drift
        double u_pred = gas.entropy[i] * std::pow(gas.rho[i], gamma_minus_1) /
                        gamma_minus_1;

        // Enforce the thermodynamic floor on the predicted state
        if (u_pred < u_floor) u_pred = u_floor;

        double p_pred = gamma_minus_1 * gas.rho[i] * u_pred;
        double p_jeans =
            get_jeans_pressure_floor(gas.rho[i], gas.min_h_cap, a, config);
        gas.pressure_pred[i] = std::max(p_pred, p_jeans);
    }
}

double get_u_floor(double a, const Config& config) {
    double u_floor =
        get_internal_energy_from_temp(config.temp_floor_k, a, config);
    return u_floor;
}

double get_pressure_floor(double a, const Config& config) {
    double u_floor = get_u_floor(a, config);
    double gamma_minus_1 = config.gamma - 1.0;
    double pressure_floor = gamma_minus_1 * config.density_floor * u_floor;
    return pressure_floor;
}

void kick_particle_thermodynamics(GasParticleSystem& gas, size_t i,
                                  double ke_spec, double grav_work,
                                  double drag_work, double dt,
                                  double expansion_factor, double a,
                                  const Config& config, bool is_kick1,
                                  Integrator::KickWork& work) {
    constexpr double min_energy = 0.0;
    const double gamma_minus_1 = config.gamma - 1.0;
    const double m = gas.mass[i];
    const double rho = gas.rho[i];

    work.entropy_energy = 0.0;
    work.def_energy = 0.0;
    work.injected_energy = 0.0;

    // Base Integration (Explicit and Conservative)
    double u_post_hydro = gas.u[i] + (is_kick1 ? 0.0 : gas.hydro_du[i]);
    double u_cooling = u_post_hydro * expansion_factor;

    double u_raw = std::max(u_post_hydro - u_cooling, min_energy);
    work.exp_work = -drag_work + (u_cooling * m);

    double e_tot_post = gas.total_energy[i] +
                        (is_kick1 ? 0.0 : gas.hydro_de[i]) + (grav_work / m) +
                        (drag_work / m) - u_cooling;
    e_tot_post = std::max(e_tot_post, min_energy);

    if (!is_kick1) {
        gas.hydro_du[i] = 0.0;
        gas.hydro_de[i] = 0.0;
    }

    double u_cons = e_tot_post - ke_spec;

    // Gravity Noise Thresholding
    double a_grav_mag =
        std::sqrt(gas.acc_x[i] * gas.acc_x[i] + gas.acc_y[i] * gas.acc_y[i] +
                  gas.acc_z[i] * gas.acc_z[i]);
    double delta_E_grav = a_grav_mag * gas.h[i];

    // Arbitration (DEF & Entropy Switch)
    double thresh_trunc =
        config.dual_energy_threshold * (ke_spec + delta_E_grav);
    double thresh_kin =
        config.entropy_switch_akin * (gas.max_rel_ke[i] + u_raw);
    double thresh_grav = config.entropy_switch_agrav * delta_E_grav;

    if (u_raw > thresh_trunc) {
        // NORMAL REGIME (Shocked/Hot): derive 'u' from the integrated total
        // energy
        gas.u[i] = std::max(u_cons, min_energy);
        gas.total_energy[i] = e_tot_post;
    } else if (u_raw < thresh_kin || u_raw < thresh_grav) {
        // ENTROPY FALLBACK (Cold/Smooth): Flow is adiabatic but
        // kinetic/gravitational noise is high. Recover 'u' from entropy
        double u_clean =
            gas.entropy[i] * std::pow(rho, gamma_minus_1) / gamma_minus_1;
        u_clean -= (u_clean * expansion_factor);

        work.entropy_energy = (u_clean - u_cons) * m;
        gas.u[i] = std::max(u_clean, min_energy);
        gas.total_energy[i] = gas.u[i] + ke_spec;
    } else {
        // DUAL ENERGY FORMALISM (Intermediate): Truncation error ruins
        // conservative subtraction. Trust explicitly integrated 'u_raw'
        work.def_energy = (u_raw - u_cons) * m;
        gas.u[i] = u_raw;
        gas.total_energy[i] = gas.u[i] + ke_spec;
    }

    // Apply Thermodynamic Floor
    double u_floor =
        get_internal_energy_from_temp(config.temp_floor_k, a, config);
    if (gas.u[i] < u_floor) {
        double injected = u_floor - gas.u[i];
        work.injected_energy = injected * m;
        gas.u[i] = u_floor;
        gas.total_energy[i] += injected;
    }

    // Finalize Entropy and Pressure to match internal energy
    gas.entropy[i] = gamma_minus_1 * gas.u[i] / std::pow(rho, gamma_minus_1);
    gas.pressure[i] = gamma_minus_1 * rho * gas.u[i];
}

void kick_wakeup_catchup(GasParticleSystem& gas, size_t i, double a,
                         double dt_actual, double expansion_factor,
                         const Config& config) {
    const double u_floor = get_u_floor(a, config);

    // Thermodynamic update using continuous rates
    gas.u[i] += gas.du_dt[i] * dt_actual;
    gas.total_energy[i] += gas.de_dt[i] * dt_actual;

    // Apply adiabatic expansion cooling
    double u_cooling = gas.u[i] * expansion_factor;
    gas.u[i] -= u_cooling;
    gas.total_energy[i] -= u_cooling;

    // Floor enforcement to prevent crashes during the catch-up
    gas.u[i] = std::max(gas.u[i], u_floor);
    gas.total_energy[i] = std::max(gas.total_energy[i], u_floor);

    // Re-sync state
    double gamma_minus_1 = config.gamma - 1.0;
    gas.entropy[i] =
        gamma_minus_1 * gas.u[i] / std::pow(gas.rho[i], gamma_minus_1);
    gas.pressure[i] = gamma_minus_1 * gas.rho[i] * gas.u[i];
}

}  // namespace Thermodynamics
