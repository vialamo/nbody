#pragma once

#include "config.h"
#include "mfm.h"  // Needed for GasParticleSystem definition
#include "particles.h"

namespace DensitySolver {

// Main entry point for computing density, smoothing lengths, and gradient
// matrices
void compute_density_and_h(GasParticleSystem& gas, const Config& config,
                           const ParticleSystem& dm);

}  // namespace DensitySolver
