#pragma once

#include <stddef.h>

struct Config;
struct Diagnostics;
struct SimState;
class GasParticleSystem;

namespace Integrator {

void compute_forces(SimState& state, Config& config, Diagnostics& diag);
void update_timesteps(SimState& state, const Config& config);

void KDK_step(SimState& state, double dt, Config& config,
              Diagnostics& diag);

struct KickWork {
    double grav_work = 0.0;
    double exp_work = 0.0;
};

KickWork kick_single_gas_particle(GasParticleSystem& gas, size_t i,
                                          double dt, double a, double H,
                                          const Config& config);
                                    
} // end namespace
