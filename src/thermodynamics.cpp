#include "thermodynamics.h"

#include "integrator.h"
#include "mfm.h"

namespace Thermodynamics {

void initialize_primitive_variables(GasParticleSystem& gas,
                                    const Config& config, double a) {
    if (gas.num_active == 0) return;

    const double gamma_minus_1 = config.gamma - 1.0;

    // Calculate the physical internal energy floor
    const double u_floor =
        get_internal_energy_from_temp(config.temp_floor_k, a, config);

    double step_floor_heating = 0.0;

    if (gas.num_active != 0) {
#pragma omp parallel for simd reduction(+ : step_floor_heating) schedule(static)
        for (size_t k = 0; k < gas.num_active; ++k) {
            size_t i = gas.active_indices[k];

            // Enforce density and metallicity floors/bounds
            gas.rho[i] = std::max(gas.rho[i], density_floor);
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
    }

    gas.accumulated_photoheating_energy += step_floor_heating;
}

void predict_primitive_variables(GasParticleSystem& gas, const Config& config,
                                 double a) {
    if (gas.num_active == 0) return;

    const double gamma_minus_1 = config.gamma - 1.0;
    const double u_floor =
        get_internal_energy_from_temp(config.temp_floor_k, a, config);

#pragma omp parallel for schedule(static)
    for (size_t k = 0; k < gas.num_active; ++k) {
        size_t i = gas.active_indices[k];

        // Project from the Kick 1 half-step to the global time (t + dt)
        double dt_pred_i =
            gas.global_time - (gas.t_end[i] - gas.dt_step[i] / 2.0);

        gas.vel_x_pred[i] =
            gas.vel_x[i] + (gas.acc_x[i] + gas.hydro_acc_x[i]) * dt_pred_i;
        gas.vel_y_pred[i] =
            gas.vel_y[i] + (gas.acc_y[i] + gas.hydro_acc_y[i]) * dt_pred_i;
        gas.vel_z_pred[i] =
            gas.vel_z[i] + (gas.acc_z[i] + gas.hydro_acc_z[i]) * dt_pred_i;

        // Predict the new internal energy adiabatically using the new density.
        // This captures the mechanical PdV work of the drift
        double u_pred_i = gas.entropy[i] * std::pow(gas.rho[i], gamma_minus_1) /
                          gamma_minus_1;

        // Enforce the thermodynamic floor on the predicted state
        if (u_pred_i < u_floor) u_pred_i = u_floor;

        gas.pressure_pred[i] = gamma_minus_1 * gas.rho[i] * u_pred_i;
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
    double pressure_floor = gamma_minus_1 * density_floor * u_floor;
    return pressure_floor;
}

void kick_particle_thermodynamics(GasParticleSystem& gas, size_t i,
                                  double ke_spec_new, double grav_work,
                                  double drag_work, double dt,
                                  double expansion_factor, double a,
                                  const Config& config,
                                  Integrator::KickWork& work) {
    constexpr double min_energy = 0.0;
    const double gamma_minus_1 = config.gamma - 1.0;
    const double m = gas.mass[i];
    const double rho = gas.rho[i];

    work.switch_energy = 0.0;
    work.injected_energy = 0.0;

    // Base Integration (Explicit and Conservative)
    double u_post_hydro = gas.u[i] + gas.du_dt[i] * dt;
    double u_cooling = u_post_hydro * expansion_factor;

    double u_raw = std::max(u_post_hydro - u_cooling, min_energy);
    work.exp_work = -drag_work + (u_cooling * m);

    double e_tot_post = gas.total_energy[i] + (gas.de_dt[i] * dt) +
                        (grav_work / m) + (drag_work / m) - u_cooling;
    e_tot_post = std::max(e_tot_post, min_energy);

    double u_cons = e_tot_post - ke_spec_new;

    // Gravity Noise Thresholding
    double a_grav_mag =
        std::sqrt(gas.acc_x[i] * gas.acc_x[i] + gas.acc_y[i] * gas.acc_y[i] +
                  gas.acc_z[i] * gas.acc_z[i]);
    double delta_E_grav = a_grav_mag * gas.h[i];

    // Arbitration (DEF & Entropy Switch)
    double thresh_trunc = config.dual_energy_threshold * ke_spec_new;
    double thresh_kin =
        config.entropy_switch_akin * (gas.max_rel_ke[i] + u_raw);
    double thresh_grav = config.entropy_switch_agrav * delta_E_grav;

    if (u_raw < thresh_kin || u_raw < thresh_grav) {
        // ENTROPY FALLBACK: Revert to pristine adiabatic state
        double u_clean =
            gas.entropy[i] * std::pow(rho, gamma_minus_1) / gamma_minus_1;
        u_clean -= (u_clean * expansion_factor);

        work.switch_energy = (u_clean - u_raw) * m;
        gas.u[i] = std::max(u_clean, min_energy);
        gas.total_energy[i] = gas.u[i] + ke_spec_new;
    } else if (u_cons < thresh_trunc) {
        // DUAL ENERGY FORMALISM: Trust explicit integration
        work.switch_energy = (u_raw - u_cons) * m;
        gas.u[i] = u_raw;
        gas.total_energy[i] = gas.u[i] + ke_spec_new;
    } else {
        // NORMAL REGIME: Trust conservative integration
        gas.u[i] = std::max(u_cons, min_energy);
        gas.total_energy[i] = e_tot_post;
    }

    // Apply Thermodynamic Floor
    const double u_floor = Thermodynamics::get_internal_energy_from_temp(
        config.temp_floor_k, a, config);
    if (gas.u[i] < u_floor) {
        double injected = u_floor - gas.u[i];
        work.injected_energy = injected * m;
        gas.u[i] = u_floor;
        gas.total_energy[i] += injected;
    }

    // Finalize Entropy and Pressure to match settled internal energy
    gas.entropy[i] = gamma_minus_1 * gas.u[i] / std::pow(rho, gamma_minus_1);
    gas.pressure[i] = gamma_minus_1 * rho * gas.u[i];
}

}  // namespace Thermodynamics
