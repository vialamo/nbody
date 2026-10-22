#pragma once

class GasParticleSystem;
struct Config;

namespace DensitySolver {

// Main entry point for computing density, smoothing lengths, and gradient
// matrices
void compute_density_and_h(GasParticleSystem& gas, const Config& config);

}  // namespace DensitySolver
