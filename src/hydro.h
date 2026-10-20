#pragma once

struct Config;
class GasParticleSystem;

namespace HydroSolver {

// Compute the least-squares matrix gradients for all primitive variables
void compute_gradients(GasParticleSystem& gas, const Config& config);

// Solve the Riemann problem between neighbors and update hydro_acc and du_dt
void compute_forces(GasParticleSystem& gas, const Config& config,
                          double a);

}  // namespace HydroSolver
