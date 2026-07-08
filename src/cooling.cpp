#include "cooling.h"

#include <cmath>

#include "constants.h"
#include "math_utils.h"

namespace cooling {

static QuadraticInterpolator A_curve;

void initialize(const Config& config) {
    // Empirical fit for subgrid clumping
    // x = resolution y = A factor
    Eigen::Vector2d p1(0.325, 0.5);
    Eigen::Vector2d p2(0.5, 1.4);
    Eigen::Vector2d p3(1.04, 10.0);

    A_curve.Set(p1, p2, p3);
}

double compute_cooling_rate(double u_code, double rho_code, double a,
                            const Config& config) {
    double T = get_temp_from_internal_energy(u_code, a, config);

    // Convert code density to CGS
    double rho_cgs = (rho_code * config.unit_density_cgs) / (a * a * a);

    // Number density (n_e = n_i = n)
    double n = rho_cgs / (config.primordial_mu * constants::M_P_CGS);

    double clumping_factor = 1.0;
    if (config.enable_subgrid_clumping) {
        double mean_rho_code =
            config.gas_total_mass /
            (config.domain_size * config.domain_size * config.domain_size);
        double overdensity = rho_code / mean_rho_code;

        // Apply clumping only in collapsing regions
        if (overdensity > 1.0) {
            double resolution = config.box_size_mpc / config.mesh_size;
            double A = config.subgrid_clumping_amplitude < 0
                           ? std::max(0.0, A_curve.evaluate(resolution))
                           : config.subgrid_clumping_amplitude;

            constexpr double B = 0.7;
            double A_dynamic = A * a;
            clumping_factor = 1.0 + A_dynamic * std::pow(overdensity, B);
        }
    }

    // Bremsstrahlung specific energy loss rate (erg / g / s)
    double lambda_cgs = 1.4e-27 * std::sqrt(T) * n * n * clumping_factor;
    double du_dt_cgs = lambda_cgs / rho_cgs;

    // Convert specific energy loss rate back to code units
    return (du_dt_cgs * config.cooling_conversion_factor) / (a * a);
}

double solve_cooling_implicit(double u_old, double rho_code, double a,
                              double dt, const Config& config,
                              int& iterations_taken) {
    // Use the highest of the physical radiation floor or the hydro
    // temp floor
    double target_floor_k =
        std::max(config.cooling_cutoff_k, config.temp_floor_k);
    double u_rad_floor =
        get_internal_energy_from_temp(target_floor_k, a, config);

    if (u_old <= u_rad_floor) {
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
        // Using central difference
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

        // Enforce physical bounds
        if (u_guess < u_rad_floor) {
            u_guess = u_rad_floor;
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
