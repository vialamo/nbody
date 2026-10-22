#pragma once

struct Config;
struct Diagnostics;
struct SimState;
class GasParticleSystem;
class ParticleSystem;

#define GRAVITY_OWN_MEMORY

namespace GravitySolver {

void update_softening(SimState& state, Config& config);

void add_pp_accel(ParticleSystem& dm, GasParticleSystem* gas, double a,
                  const Config& config, Diagnostics& diag);

// Timestep
double get_timestep(const ParticleSystem& dm, const Config& config);
double get_timestep(const GasParticleSystem& gas, const Config& config);

}  // namespace GravitySolver
