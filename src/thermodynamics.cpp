#include "thermodynamics.h"

#include "cooling.h"
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

        // Predict the new internal energy adiabatically using the fresh density
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

}  // namespace Thermodynamics
