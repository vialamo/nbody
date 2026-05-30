#pragma once
#include "config.h"

namespace cooling {

constexpr int MAX_ITER = 50;

// Convert code internal energy to Temperature (Kelvin)
inline double get_temp_from_internal_energy(double u_code, double a,
                                            const Config& config) {
    return u_code * (a * a) * config.FACTOR_U_TO_T;
}

// Convert Temperature (Kelvin) to code internal energy
inline double get_internal_energy_from_temp(double T_kelvin, double a,
                                            const Config& config) {
    return (T_kelvin * config.FACTOR_T_TO_U) / (a * a);
}

// Computes the cooling rate Lambda in code units [u_code / t_code]
// u_code is the specific internal energy.
// rho_code is the code density.
double compute_cooling_rate(double u_code, double rho_code, double a,
                            const Config& config);

// Implicitly solves for the new internal energy after cooling over timestep dt
// Uses Newton-Raphson to ensure stability for stiff cooling rates.
double solve_cooling_implicit(double u_old, double rho_code, double a,
                              double dt, const Config& config,
                              int& iterations_taken);

}  // namespace cooling