#pragma once

#include "config.h"
#include "diagnostics.h"
#include "mfm.h"  // Needed for GasParticleSystem and KickWork definitions
#include "particles.h"

namespace GravitySolver {

void bin_and_assign_mass(GasParticleSystem& gas, const Config& config);

void interpolate_cic_forces(GasParticleSystem& gas, const Grid3D& ax_grid,
                            const Grid3D& ay_grid, const Grid3D& az_grid,
                            const Config& config);

void compute_and_add_pp_forces(GasParticleSystem& gas, double a,
                               const Config& config, Diagnostics& diag);

void compute_cross_pp_forces(GasParticleSystem& gas, double a,
                             ParticleSystem& dm, const Config& config,
                             Diagnostics& diag);

KickWork kick_particle_gravity(GasParticleSystem& gas, size_t i, double dt,
                               double a, double H, const Config& config);

double get_gravity_timestep(const GasParticleSystem& gas, const Config& config);

}  // namespace GravitySolver
