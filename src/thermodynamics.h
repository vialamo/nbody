#pragma once
#include <vector>

#include "config.h"
#include "integrator.h"

class GasParticleSystem;

namespace Thermodynamics {

void initialize_primitive_variables(GasParticleSystem& gas,
                                    const Config& config, double a);

void predict_primitive_variables(GasParticleSystem& gas, const Config& config,
                                 double a);

void kick_particle_thermodynamics(GasParticleSystem& gas, size_t i,
                                  double ke_spec_new, double grav_work,
                                  double drag_work, double dt,
                                  double expansion_factor, double a,
                                  const Config& config,
                                  Integrator::KickWork& work);

double get_u_floor(double a, const Config& config);
double get_pressure_floor(double a, const Config& config);

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
