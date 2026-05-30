#include "cooling.h"

#include <cmath>

#include "constants.h"

namespace cooling {

double compute_cooling_rate(double u_code, double rho_code, double a,
                            const Config& config) {
    double T = get_temp_from_internal_energy(u_code, a, config);

    if (T <= 10.0) return 0.0;

    // Convert code density to CGS
    double rho_cgs = (rho_code * config.UNIT_DENSITY_CGS) / (a * a * a);

    // Number density (n_e = n_i = n)
    double n = rho_cgs / (config.PRIMORDIAL_MU * constants::M_P_CGS);

    // Bremsstrahlung specific energy loss rate (erg / g / s)
    double lambda_cgs = 1.4e-27 * std::sqrt(T) * n * n;
    double du_dt_cgs = lambda_cgs / rho_cgs;

    // Convert specific energy loss rate back to code units
    return (du_dt_cgs * config.COOLING_CONVERSION_FACTOR) / (a * a);
}

double solve_cooling_implicit(double u_old, double rho_code, double a,
                              double dt, const Config& config,
                              int& iterations_taken) {
    double dynamic_u_floor =
        get_internal_energy_from_temp(config.TEMP_FLOOR_KELVIN, a, config);

    if (u_old <= dynamic_u_floor) {
        return u_old;  // Already at or below the floor, no cooling
    }

    int iter = 0;
    const double TOLERANCE = 1e-5;  // Require 0.001% convergence

    double u_guess = u_old;

    for (; iter < MAX_ITER; ++iter) {
        double lambda = compute_cooling_rate(u_guess, rho_code, a, config);

        // The root-finding function: f(u) = 0
        double f_u = u_guess - u_old + dt * lambda;

        // Numerical derivative for the Jacobian: d(Lambda)/du
        // Using central difference for better accuracy
        double eps = 1e-4 * u_guess;
        double lambda_plus =
            compute_cooling_rate(u_guess + eps, rho_code, a, config);
        double lambda_minus =
            compute_cooling_rate(u_guess - eps, rho_code, a, config);
        double dlambda_du = (lambda_plus - lambda_minus) / (2.0 * eps);

        // Derivative of the root-finding function
        double df_du = 1.0 + dt * dlambda_du;

        // Newton-Raphson step
        double du = f_u / df_du;
        u_guess -= du;

        // Enforce physical bounds (prevent negative energy or dropping below
        // the CMB/floor)
        if (u_guess < dynamic_u_floor) {
            u_guess = dynamic_u_floor;
            break;  // Hit the floor, we can safely exit
        }

        // Check for convergence
        if (std::abs(du / u_guess) < TOLERANCE) {
            break;
        }
    }

    iterations_taken = iter;

    return u_guess;
}

}  // namespace cooling