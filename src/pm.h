#pragma once
struct Config;
struct SimState;

namespace PM {

// Splats particle masses onto a 3D grid and caches the interpolation weights
void bin_and_assign_mass(SimState& state, const Config& config);

// Takes the density grid, solves Poisson's equation in Fourier space,
// and computes the spatial gradients to yield acceleration grids
void compute_acceleration(SimState& state, const Config& config);

// Gathers forces from 3D grids back to the particles using cached weights
void interpolate_forces(SimState& state, const Config& config);

}  // namespace PM
