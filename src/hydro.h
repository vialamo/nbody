#pragma once

#include "config.h"
#include "mfm.h"  // Needed for GasParticleSystem definition

namespace HydroSolver {

// Compute the least-squares matrix gradients for all primitive variables
void compute_gradients(GasParticleSystem& gas, const Config& config);

// Solve the Riemann problem between neighbors and update hydro_acc and du_dt
void compute_hydro_forces(GasParticleSystem& gas, const Config& config,
                          double a, double dt);

}  // namespace HydroSolver
