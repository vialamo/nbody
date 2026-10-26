#pragma once

#include "types.h"

struct Config;
struct SimState;

class PM {
   public:
    Grid3D rho;
    Grid3D phi;

    // Comoving gravitational acceleration on the mesh
    Grid3D acc_x;
    Grid3D acc_y;
    Grid3D acc_z;

    PM(const Config& config);
    ~PM();

    // Splats particle masses onto a grid and caches the interpolation weights
    void bin_and_assign_mass(SimState& state, const Config& config);

    // Takes the density grid, solves Poisson's equation in Fourier space,
    // and computes the spatial gradients to yield acceleration grids
    void compute_acceleration(SimState& state, const Config& config);

    // Gathers forces from grids back to the particles using cached weights
    void interpolate_forces(SimState& state, const Config& config);
};
