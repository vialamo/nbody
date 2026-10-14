#pragma once

#include <stddef.h>

struct Config;
struct Diagnostics;
struct SimState;
struct TimestepInfo;
class GasParticleSystem;

void compute_forces(SimState& state, Config& config, Diagnostics& diag);

void KDK_step(SimState& state, TimestepInfo& ts, Config& config,
              Diagnostics& diag);

struct KickWork {
    double grav_work = 0.0;
    double exp_work = 0.0;
};

KickWork kick_single_gas_particle_gravity(GasParticleSystem& gas, size_t i,
                                          double dt, double a, double H,
                                          const Config& config);

void kick_single_gas_particle_hydro(GasParticleSystem& gas, size_t i, double dt,
                                    const Config& config);
