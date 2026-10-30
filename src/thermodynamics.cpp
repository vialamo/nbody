#include "thermodynamics.h"

#include "cooling.h"
#include "mfm.h"

namespace Thermodynamics {

static inline void evaluate_particle_thermodynamics(
    const GasParticleSystem& gas, size_t i, double u_floor, double alpha_kin,
    double alpha_grav, double gamma_minus_1, double& out_u, double& out_p,
    double& out_s, double& out_e_tot, double& out_injected_u,
    double& out_switch_energy) {
    out_injected_u = 0.0;
    out_switch_energy = 0.0;

    double local_rho = std::max(gas.rho[i], density_floor);

    // Calculate specific kinetic energy
    double ke =
        0.5 * (gas.vel_x[i] * gas.vel_x[i] + gas.vel_y[i] * gas.vel_y[i] +
               gas.vel_z[i] * gas.vel_z[i]);

    // Calculate the conserved internal energy
    double u_cons = gas.total_energy[i] - ke;

    double threshold_truncation = alpha_kin * ke;
    double threshold_kin = alpha_kin * (gas.max_rel_ke[i] + gas.u[i]);
    double threshold_grav = alpha_grav * gas.delta_E_grav[i];

    double u_temp = gas.u[i];
    double s_temp = gas.entropy[i];
    double e_tot_temp = gas.total_energy[i];

// #define USE_ENTROPY_SWITCH
#ifdef USE_ENTROPY_SWITCH
    if (u_temp < threshold_kin || u_temp < threshold_grav) {
        // FALLBACK TRIGGERED: Extreme Mach number detected.
        // Discard the 'u' updated by the Riemann solver and use entropy-based
        // adiabatic evolution
        double u_new =
            s_temp * std::pow(local_rho, gamma_minus_1) / gamma_minus_1;

        out_switch_energy = (u_new - u_temp) * gas.mass[i];
        u_temp = u_new;
        e_tot_temp = u_temp + ke;
    } else
#endif
        if (u_cons < threshold_truncation) {
        // DUAL ENERGY FORMALISM: Highly supersonic/kinetic-dominated flow.
        // u_cons is corrupted by truncation error. Trust the explicitly
        // integrated u[i]
        u_temp = gas.u[i];

        out_switch_energy = (gas.u[i] - u_cons) * gas.mass[i];

        s_temp = gamma_minus_1 * u_temp / std::pow(local_rho, gamma_minus_1);
        e_tot_temp =
            u_temp + ke;  // Synchronize total energy to match the trusted u[i]

    } else {
        // NORMAL REGIME: Trust the conservative internal energy.
        // Re-sync the entropy array to match the shock-heated state
        u_temp = u_cons;
        s_temp = gamma_minus_1 * u_temp / std::pow(local_rho, gamma_minus_1);
    }

    // Enforce thermodynamic floor
    if (u_temp < u_floor) {
        // Track the added energy
        double injected_u = u_floor - u_temp;
        out_injected_u = injected_u * gas.mass[i];

        u_temp = u_floor;
        s_temp = gamma_minus_1 * u_temp / std::pow(local_rho, gamma_minus_1);
        e_tot_temp += injected_u;
    }

    // Assign final computed variables to the output references
    out_u = u_temp;
    out_p = gamma_minus_1 * local_rho * u_temp;
    out_s = s_temp;
    out_e_tot = e_tot_temp;
}

void update_primitive_variables(GasParticleSystem& gas, const Config& config,
                                double a) {
    if (gas.num_particles == 0) return;

    static int num_cycles = 0;

    const double gamma_minus_1 = config.gamma - 1.0;

    // Calculate the physical internal energy floor
    const double u_floor =
        get_internal_energy_from_temp(config.temp_floor_k, a, config);

    double alpha_kin = config.dual_energy_threshold;
    double alpha_grav = config.dual_energy_threshold;

    double step_floor_heating = 0.0;
    double step_entropy_switch = 0.0;

    if (gas.num_active != 0) {
#pragma omp parallel for simd reduction( \
        + : step_floor_heating, step_entropy_switch) schedule(static)
        for (size_t k = 0; k < gas.num_active; ++k) {
            size_t i = gas.active_indices[k];

            double local_injected_u = 0.0;
            double local_switch_energy = 0.0;

            // Update the class arrays
            evaluate_particle_thermodynamics(
                gas, i, u_floor, alpha_kin, alpha_grav, gamma_minus_1, gas.u[i],
                gas.pressure[i], gas.entropy[i], gas.total_energy[i],
                local_injected_u, local_switch_energy);

            gas.rho[i] = std::max(gas.rho[i], density_floor);
            gas.metal_frac[i] = std::max(0.0, std::min(gas.metal_frac[i], 1.0));

            step_floor_heating += local_injected_u;
            step_entropy_switch += local_switch_energy;
        }
    }
    num_cycles++;

    gas.accumulated_photoheating_energy += step_floor_heating;
    gas.accumulated_entropy_switch_energy += step_entropy_switch;
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

        double u_pred_i = gas.u[i] + gas.du_dt[i] * dt_pred_i;

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

}  // namespace Thermodynamics
