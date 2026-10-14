#pragma once

#include "config.h"
#include "diagnostics.h"
#include "mfm.h"
#include "particles.h"

struct SimState;

namespace GravitySolver {

void update_softening(SimState& state, Config& config);

// Gas
void compute_and_add_gas_pp_forces(GasParticleSystem& gas, double a,
                               const Config& config, Diagnostics& diag);

double get_gravity_timestep(const GasParticleSystem& gas, const Config& config);

// DM
void compute_and_add_dm_pp_forces(ParticleSystem& dm, double a,
                                  const Config& config, Diagnostics& diag);

double get_gravity_timestep(const ParticleSystem& dm, const Config& config);

// Gas - DM
void compute_gas_dm_cross_pp_forces(ParticleSystem& dm,
                                GasParticleSystem& gas, double a,
                                const Config& config, Diagnostics& diag);

}  // namespace GravitySolver
