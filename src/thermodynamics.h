#pragma once
#include <vector>

#include "config.h"

class GasParticleSystem;

namespace Thermodynamics {
// Sync tracked internal energy, total energy, and compute pressure
void update_primitive_variables(GasParticleSystem& gas, const Config& config,
                                double a);

void predict_primitive_variables(GasParticleSystem& gas, const Config& config,
                                 double a);

// Computes synced snapshot arrays (u, P, S, E_tot) without modifying the
// state
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
           std::vector<double>>
compute_snapshot_primitives(const GasParticleSystem& gas, const Config& config,
                            double a);


double get_u_floor(double a, const Config& config);

// Conversion helpers
inline double get_temp_from_internal_energy(double u_code, double a,
                                            const Config& config) {
    return u_code * (a * a) * config.factor_u_to_t;
}

inline double get_internal_energy_from_temp(double T_kelvin, double a,
                                            const Config& config) {
    return (T_kelvin * config.factor_t_to_u) / (a * a);
}

}  // namespace Thermodynamics
