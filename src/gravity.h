#pragma once

struct Config;
struct Diagnostics;
struct SimState;
class GasParticleSystem;
class ParticleSystem;

namespace GravitySolver {

void update_softening(SimState& state, Config& config);

// Gas
void compute_and_add_gas_pp_forces(GasParticleSystem& gas, double a,
                               const Config& config, Diagnostics& diag);

// DM
void compute_and_add_dm_pp_forces(ParticleSystem& dm, double a,
                                  const Config& config, Diagnostics& diag);

// Gas - DM
void compute_gas_dm_cross_pp_forces(ParticleSystem& dm,
                                GasParticleSystem& gas, double a,
                                const Config& config, Diagnostics& diag);
                                

// Timestep                                
double get_timestep(const ParticleSystem& dm, const Config& config);
double get_timestep(const GasParticleSystem& gas, const Config& config);

}  // namespace GravitySolver
