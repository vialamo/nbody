#pragma once
#include "gas.h"
#include "particles.h"
#include "types.h"

struct SimState {
    ParticleSystem dm;
    GasGrid gas;

    Grid3D total_rho;
    Grid3D phi;

    // Global Gravitational Field on the Eulerian Mesh
    Grid3D gravity_x;
    Grid3D gravity_y;
    Grid3D gravity_z;

    double total_time = 0.0;
    double scale_factor = 1.0;
    double hubble_param = 0.0;

    SimState(const Config& config)
        : dm(config),
          gas(config),
          gravity_x(config.mesh_size),
          gravity_y(config.mesh_size),
          gravity_z(config.mesh_size),
          total_rho(config.mesh_size),
          phi(config.mesh_size) {}
};